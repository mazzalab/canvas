__author__ = "Tommaso Mazza"
__copyright__ = "Copyright 2017, The AnnotateCNV Project"
__version__ = "0.0.9"
__maintainer__ = "Tommaso Mazza"
__email__ = "bioinformatics@css-mendel.it"
__status__ = "Development"
__date__ = "30/01/2018"
__creator__ = "t.mazza"
__license__ = u"""
  Copyright (C) 20016-2017  Tommaso Mazza <t,mazza@css-mendel.it>
  Viale Regina Margherita 261, 00198 Rome, Italy

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
  02110-1301 USA
  """

""""
It annotates a tab-delimited input file with a set of BED files, with textual features in their 4th columns
"""

import argparse
import re
import time
import collections
import pandas as pd
from pandas import DataFrame
from pymongo import MongoClient
import sys
import os
from io import StringIO
import itertools
from collections import OrderedDict
from subprocess import call, Popen, PIPE
import glob
from pprint import pprint
import json


class MainApp:
    def __init__(self, args):

        self.args = args
        self.connection = MongoClient('localhost', 27017, unicode_decode_error_handler='ignore')
    
        self.db = self.connection['BED']
        self.db_gl = self.connection['genelists']
        self.db_t = self.connection['targets']
        self.d = vars(self.args)
        self.mirbase_dict = {}  #this stores the mirbase annotation that will be written as a separate json
        print(self.d)
        # input()
        

    def process(self):
    
        # sys.stdout = sys.stderr = open(self.args.out.replace('.xlsx', '_log.txt'), 'wt')
        
        for arg in self.d:
            
            if arg not in ['cnv_line', 'cnv_file', 'out', 'distance', 'all_beds', 'all_genelists', 'reference'] and \
                    (self.d[arg] or self.d['all_beds']) and not arg.startswith('__'):
                sys.stdout.write("Checking presence of {}...\n".format(arg))
                
                # Checking presence of required genelists and BEDs....
                if 'genelist' in arg and not arg.replace('_genelist', '') in self.db_gl.collection_names():
                    sys.exit("The collection {} could not be found in the genelists db. Exiting.".format(arg))
                if not 'genelist' in arg and not arg in self.db.collection_names():
                    sys.exit("The collection {} could not be found in the BED db. Exiting.".format(arg))

                if arg == 'mirbase':
                    if not ('tarbase' in self.db_t.collection_names() or
                                    'targetscan' in self.db_t.collection_names()):
                        sys.exit("Need tarbase and/or targetscan collection in targets db to annotate mirna. "
                                 "Exiting.")
                    if not self.d['mirna'] and not self.d['all_beds']:
                        sys.exit('The --mirbase option depends on the prior execution of --mirna. '
                                 'Please include this option as well.')
                        
        sys.stdout.write("Requested DBs are present. Proceeding to annotate...\n")
        
        # In case the input is a file...
        if self.args.cnv_file:
            if os.path.exists(self.args.cnv_file):
                try:
                    cnv_info = self.read_cnv_coordinates_file(self.args.cnv_file)
                except NameError:
                    return -1

            else:
                sys.exit("CNV file not found: {}".format(self.args.cnv_file))
                
        # Alternatively, with a text line input..
        elif self.args.cnv_line:
            cnv_info = self.read_cnv_coordinates_line(self.args.cnv_line)

        if self.args.reference != 'hg19':
            # Testing presence of CrossMap executable in PATH
            test = Popen(['resources/liftover_files/crossmap'], stdout=PIPE, stderr=PIPE)
            out, err = test.communicate()
    
            if not err.decode('ascii').strip().startswith('liftOver - Move'):
                sys.exit('crossmap not found in system PATH, or not executable.')
    
            print("Chosen reference is", self.args.reference)
            chainfile = os.path.join("resources/liftover_files/chains",
                                     "{}ToHg19.over.chain".format(self.args.reference))
            print(cnv_info)
    
            # For the conversion, it is necessary to save a temporary file
            # and run CrossMap on it.
            tempfile = os.path.join(os.path.dirname(self.args.out), 'tempfile.tsv')
            tempfile_converted = tempfile.replace(".tsv", "_converted.tsv")
            tempfile_unmapped = tempfile.replace(".tsv", "_unmapped.tsv")
            cnv_info.to_csv(tempfile, sep='\t',
                            columns=['CHR', 'START', 'END'], index=False, header=False)
    
            call(["resources/liftover_files/crossmap", tempfile, chainfile,
                  tempfile_converted, tempfile_unmapped])
            cnv_info_converted = pd.read_table(tempfile_converted, encoding='cp1252', sep='\t', header=None,
                                               names=["CHR", "START", "END"])
    
            # if no cnv was split during liftover, any possible extra column will be preserved.
            # otherwise, only CHR, START, END will be retained.
    
            if len(cnv_info_converted) == len(cnv_info):
                print("Conversion successful. Re-adding additional input fields if present.")
                cnv_info["CHR"] = cnv_info_converted["CHR"]
                cnv_info["START"] = cnv_info_converted["START"]
                cnv_info["END"] = cnv_info_converted["END"]
            else:
                print("One or more CNVs were split during liftover from {0} to hg19. If the input"
                      "file contained extra columns besides CHR, START, END, they will be "
                      "discarded.".format(self.args.reference))
                cnv_info = cnv_info_converted
            print("Converted:")
            print(cnv_info)
            # Cleaning
            for f in glob.glob(os.path.dirname(self.args.out) + '/tempfile*.tsv'):
                os.remove(f)
                
        self.out_dataframe = cnv_info

        # Annotation based on the options that were chosen
        count = 0
        extra_info = {}
        for arg in self.d:
            if arg not in ['cnv_line', 'cnv_file', 'out', 'distance', 'all_beds', 'all_genelists', 'reference'] and \
                    (self.d[arg] or self.d['all_beds']) and not arg.startswith('__') and not 'genelist' in arg:
                count += 1
                sys.stdout.write("Adding annotation for {}...\n".format(arg))
                f = open(os.path.dirname(self.args.out)+'/'+str(count)+'_'+arg+'.progress', 'w')
                f.close()
                if arg != 'mirbase':
                    extra_annot_info = self.add_annotation(arg)
                    extra_info.update({arg:extra_annot_info})
                else:
                    for target_db in self.db_t.collection_names():
                        if target_db != 'system.indexes':
                            self.mirbase_dict[target_db] = []  # initialized the target key (tarbase or targetscan) in the mirbase dict
                            
                            dict_inside, dict_cross, dict_distal = self.add_mirna_target(target=target_db, unique=True)
                            self.mirbase_dict[target_db].extend((dict_inside, dict_cross, dict_distal))
        
        # Genelists annotation
        sys.stdout.write("Adding genelists classifications...\n")
        if self.d['gene'] or self.d['all_beds']:
            for name in sorted(self.db_gl.collection_names()):
                if name != 'system.indexes':
                    print("checking ", name)
                    if self.d[name+'_genelist'] or self.d['all_genelists']:
                        print(name+'_genelist', self.d[name+'_genelist'])
                        print("Adding {} gene classification...\n".format(name))
                        self.add_meta_gene(name)
        else:
            sys.stderr.write("WARNING: genelists not added since --gene option was not included.\n")

        #Writing final file
        write_file(self.out_dataframe, self.mirbase_dict, self.args.out, extra_info)
        # return self.out_dataframe.reset_index().to_json(orient='records')

        #Cleaning
        for f in glob.glob(os.path.dirname(self.args.out)+'/*.progress'):
            os.remove(f)
            
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__

        
        return 0
    
    def read_cnv_coordinates_file(self, cnv_file: str) -> DataFrame:
        """
        Read and annotate the original CNV file as a DataFrame
        :param cnv_file: File path and name of the original CNV file to be annotated
        :return: A DataFrame containing the CNV to be annotate, one per line
        """
        cnv_coords = pd.read_table(cnv_file, encoding='cp1252', sep='\t')
        cnv_coords.columns = map(str.upper, cnv_coords.columns)
        if not (hasattr(cnv_coords, 'START')):
            raise NameError("The input file does not contain the required 'Start' field.")
        if not (hasattr(cnv_coords, 'END')):
            raise NameError("The input file does not contain the required 'End' field.")
        if not (hasattr(cnv_coords, 'CHR')):
            raise NameError("The input file does not contain the required 'Chr' field.")
        cnv_coords = cnv_coords[
            cnv_coords.START.notnull() & cnv_coords.END.notnull() & cnv_coords.CHR.notnull()]
        
        cnv_coords['START'] = cnv_coords['START'].astype('int')
        cnv_coords['END'] = cnv_coords['END'].astype('int')
        return cnv_coords

    def read_cnv_coordinates_line(self, cnv_line: str) -> DataFrame:
        """
        Read and annotate the original CNV file as a DataFrame
        :param cnv_file: File path and name of the original CNV file to be annotated
        :return: A DataFrame containing the CNV to be annotate, one per line
        """
        print("LINE")
        print(cnv_line)
        m = re.match(r'(?P<chr>chr[\dXYM]+):(?P<start>\d+)-(?P<end>\d+)', cnv_line)
        parsed_line = StringIO("""CHR\tSTART\tEND\n{0}\t{1}\t{2}""".format(m.group('chr'), m.group('start'),
                                                                           m.group('end')))
        cnv_coords = pd.read_csv(parsed_line, sep='\t')
        return cnv_coords
    
    def add_annotation(self, annotation_db: str):
        """
        Take a DataFrame and add in the last six columns the annotation provided in the 4th column of the BED file
        :param str annotation_db: File path and name of the annotation BED file
        """
        db = self.connection['BED'][annotation_db]
        sys.stdout.write('Annotating {}...\n'.format(annotation_db))
        distance_from_gene = self.args.distance
        
        inside_molecules = []
        """:type : list[str]"""
        inside_molecules_coords = []
        """:type : list[str]"""
        inside_molecules_count = []
        """:type : list[int]"""
        cross_molecules = []
        """:type : list[str]"""
        cross_molecules_coords = []
        """:type : list[str]"""
        cross_molecules_count = []
        """:type : list[int]"""
        distal_molecules = []
        """:type : list[str]"""
        distal_molecules_coords = []
        """:type : list[str]"""
        distal_molecules_count = []
        """:type : list[int]"""
        extra_info = {"inside": [], "cross": [], "distal": []}
        
        for row in self.out_dataframe.itertuples():
            chrom = row.CHR
            start = row.START
            end = row.END
            # print("finding inside")
            inside_molecule_data = db.find({"$and": [{"chr": chrom, "start": {"$gte": start},
                                            "end": {"$lte": end}}]})

            inside_molecule_data_list = list(inside_molecule_data)
            if inside_molecule_data_list and len(inside_molecule_data_list[0]) > 5:
                extra_info_cnv = list(inside_molecule_data_list)
                for i in range(0, len(extra_info_cnv)):
                    extra_info_cnv[i].pop('_id')
                extra_info['inside'].append({'cnv': "{0}:{1}-{2}".format(chrom, start, end), 'data': extra_info_cnv})
            elif not inside_molecule_data_list or len(inside_molecule_data_list[0]) <= 5:
                extra_info['inside'].append({'cnv': "{0}:{1}-{2}".format(chrom, start, end), 'data': []})

            inside_molecule = inside_molecule_data.distinct('info')
            # find the coordinates of these DISTINCT genes
            distinct_inside_genes_data = list(
                db.find({"info": {"$in": inside_molecule}}, {'start': 1, 'end': 1}))
            inside_molecule_coords = list(str(d['start'])+'-'+str(d['end']) for d in distinct_inside_genes_data)

            # print("finding cross")

            cross_molecule_data = db.find(
                {"$or": [
                    {"$and": [
                        {"chr": chrom, "start": {"$gte": start, "$lte": end}, "end": {"$gt": end}}
                    ]},
                    {"$and": [
                        {"chr": chrom, "start": {"$lte": start}, "end": {"$gte": start, "$lt": end}}
                    ]},
                    {"$and": [
                        {"chr": chrom}, {"start": {"$lt": start}}, {"end": {"$gt": start}},
                        {"start": {"$lt": end}}, {"end": {"$gt": end}}
                    ]}

                ]}
            )
            cross_molecule_data_list = list(cross_molecule_data)
            if cross_molecule_data_list and len(cross_molecule_data_list[0]) > 5:
                extra_info_cnv = list(cross_molecule_data_list)
                for i in range(0, len(extra_info_cnv)):
                    extra_info_cnv[i].pop('_id')
                extra_info['cross'].append({'cnv': "{0}:{1}-{2}".format(chrom, start, end), 'data': extra_info_cnv[0]})
            elif not cross_molecule_data_list:
                extra_info['cross'].append({'cnv': "{0}:{1}-{2}".format(chrom, start, end), 'data': []})

            cross_molecule = cross_molecule_data.distinct('info')
            # find the coordinates of these DISTINCT genes
            distinct_cross_genes_data = list(
                db.find({"info": {"$in": cross_molecule}}, {'start': 1, 'end': 1}))
            cross_molecule_coords = list(str(d['start'])+'-'+str(d['end']) for d in distinct_cross_genes_data)

            # print("finding distal")

            distal_molecule_data = db.find(
                {"$or": [
                    {"$and": [
                        {"chr": chrom}, {"start": {"$gt": end}}, {"start": {"$lt": distance_from_gene+end}},
                        {"end": {"$lt": distance_from_gene + end}}, {"end": {"$gt": end}}
                    ]},
                    {"$and": [
                        {"chr": chrom}, {"end": {"$lt": start}}, {"end": {"$gt": start-distance_from_gene}},
                        {"start": {"$gt": start - distance_from_gene}}, {"start": {"$lt": start}}
                    ]}
                ]
                })
            distal_molecule_data_list = list(distal_molecule_data)
            if distal_molecule_data_list and len(distal_molecule_data_list[0]) > 5:
                extra_info_cnv = list(distal_molecule_data_list)
                for i in range(0,len(extra_info_cnv)):
                    extra_info_cnv[i].pop('_id')
                extra_info['distal'].append({'cnv': "{0}:{1}-{2}".format(chrom, start, end), 'data': extra_info_cnv})
            elif not distal_molecule_data_list:
                extra_info['distal'].append({'cnv': "{0}:{1}-{2}".format(chrom, start, end), 'data': []})
            distal_molecule = distal_molecule_data.distinct('info')
            # find the coordinates of these DISTINCT genes
            distinct_distal_genes_data = list(db.find({"info":{"$in":distal_molecule}}, {'start': 1, 'end': 1}))
            distal_molecule_coords = list(str(d['start'])+'-'+str(d['end']) for d in distinct_distal_genes_data)
            
            # print("step2")
            # print(len(inside_molecule))
            if len(inside_molecule) > 0:
                inside_molecules.append(";".join(str(x) for x in inside_molecule))
                inside_molecules_count.append(len(inside_molecule))
                inside_molecules_coords.append(";".join(str(x) for x in inside_molecule_coords))
            else:
                inside_molecules.append(".")
                inside_molecules_count.append(0)
                inside_molecules_coords.append(".")

            # print(len(cross_molecule))
            if len(cross_molecule) > 0:
                cross_molecules.append(";".join(str(x) for x in cross_molecule))
                cross_molecules_count.append(len(cross_molecule))
                cross_molecules_coords.append(";".join(str(x) for x in cross_molecule_coords))

            else:
                cross_molecules.append(".")
                cross_molecules_count.append(0)
                cross_molecules_coords.append(".")
                
            # print(len(distal_molecule))
            if len(distal_molecule) > 0:
                distal_molecules.append(";".join(str(x) for x in distal_molecule))
                distal_molecules_count.append(len(distal_molecule))
                distal_molecules_coords.append(";".join(str(x) for x in distal_molecule_coords))
            else:
                distal_molecules.append(".")
                distal_molecules_count.append(0)
                distal_molecules_coords.append(".")
                
        # print("final")
        self.out_dataframe.loc[:, annotation_db + '_inside'] = inside_molecules
        self.out_dataframe.loc[:, annotation_db + '_inside_coords'] = inside_molecules_coords
        self.out_dataframe.loc[:, annotation_db + '_inside_count'] = inside_molecules_count
        # print("Inside", len(inside_molecules[0].split(';')), len(inside_molecules_coords[0].split(';')), inside_molecules_count)
        
        self.out_dataframe.loc[:, annotation_db + '_cross'] = cross_molecules
        self.out_dataframe.loc[:, annotation_db + '_cross_coords'] = cross_molecules_coords
        self.out_dataframe.loc[:, annotation_db + '_cross_count'] = cross_molecules_count
        # print("cross", len(cross_molecules[0].split(';')), len(cross_molecules_coords[0].split(';')), cross_molecules_count)

        self.out_dataframe.loc[:, annotation_db + '_distal'] = distal_molecules
        self.out_dataframe.loc[:, annotation_db + '_distal_coords'] = distal_molecules_coords
        self.out_dataframe.loc[:, annotation_db + '_distal_count'] = distal_molecules_count
        # print("distal", len(distal_molecules[0].split(';')), len(distal_molecules_coords[0].split(';')), distal_molecules_count)
        
        return extra_info
    
    def __get_genetarget(self, mirs_in_cnv: str, mirbase_dict: dict, target_dict: dict) -> list:
        """
        Take a list of miRs, formatted as from miRBase, parse it, get the corresponding mature miR symbols, and get a list
        of their target genes
        :param str mirs_in_cnv: List of miRs contained in a CNV region, formatted as from miRBase
        :param dict mirbase_dict: Dictionary containing the association MI id (stem loop) -> mirna_names (mature)
        :param dict target_dict: Dictionary associating mature mirna -> target genes
        :return: List of targeted genes by the miRs contained in the CNV
        """
        m_inside = re.findall(
            r'(?:\.|(\"ID=(?P<mi_name>MI[0-9]+);[Alias=MI[0-9]+]?;Name=(?P<mirna_name>[A-Za-z0-9\-]+)\"[,]?))',
            mirs_in_cnv)
        
        mature_mir_gene_names = []
        mature_mir_gene_names_dict = OrderedDict()

        for (other, mi_name, mirna_name) in m_inside:
            if mi_name in mirbase_dict:
                mature_mirna_names = mirbase_dict[mi_name]
                gene_targets = []
                """:type : list[str]"""
                for mature_mirna_name in mature_mirna_names:
                    if mature_mirna_name in target_dict:
                        gene_targets = gene_targets + target_dict[mature_mirna_name]
            else:
                gene_targets = []
            
            mature_mir_gene_names = mature_mir_gene_names + gene_targets
            if gene_targets:
                mature_mir_gene_names_dict[mirna_name+';'+mi_name] = gene_targets
            else:
                mature_mir_gene_names_dict[mirna_name+';'+mi_name] = ['.']
            # print(mature_mir_gene_names)
            # print(mature_mir_gene_names_dict)
            # print("now the tuples are", len(mature_mir_gene_names_dict))
            # input()
        # print(mature_mir_gene_names)
        # print("this is what is going ^^")
        # input()
        # list()
        # print(list(itertools.chain.from_iterable(list(mature_mir_gene_names_dict.values()))))
        # print("these are the keys ^^^")
        # print(mature_mir_gene_names_dict)
        # input()
        # print(mature_mir_gene_names == list(itertools.chain.from_iterable(list(mature_mir_gene_names_dict.values()))))
        # return mature_mir_gene_names_dict if len(mature_mir_gene_names) > 0 else ["."]

        return mature_mir_gene_names_dict

    def add_mirna_target(self, target: str, unique: bool):
        """
        Take a DataFrame and annotate it with targeting genes of miRs according to a target db
        :param target: Suffix string for the column names containing gene target (inside, cross and distal)
        :param bool unique: If true, only unique genes will be reported for each CNV region
        """
        
        db_mirbase = self.connection['BED']['mirbase']
        
        if unique:
            mirbase_list = db_mirbase.find({
                "$and": [
                    {"feature": {"$ne": "miRNA_primary_transcript"}}
                ]}
            ).distinct('info')
        else:
            
            mirbase_list =[el['info'] for el in db_mirbase.find({
                "$and": [
                    {"feature": {"$ne": "miRNA_primary_transcript"}}
                ]}, {"_id": 0, "info": 1}
            )]
            
        mirbase_info = {}
        for entry in mirbase_list:
            m = re.match(r'.+;Name=(?P<mirna_name>[A-Za-z0-9\-]+);Derives_from=(?P<mi>MI[0-9]+)$', entry)
            if not m.group('mi') in mirbase_info:
                mirbase_info[m.group('mi')] = [m.group('mirna_name')]
            else:
                mirbase_info[m.group('mi')].append(m.group('mirna_name'))
        
        # print(target)
        db_target = self.connection['targets'][target]
        target_info = {}
        # print(list(db_target.find()))
        for entry in list(db_target.find()):
            if entry['mirna'] not in target_info:
                target_info[entry['mirna']] = [entry['geneName']]
            else:
                target_info[entry['mirna']].append(entry['geneName'])
        #
        inside_mature_mir_gene_names = []
        """: type : list[str] """
        inside_mature_mir_gene_names_count = []
        """: type : list[int] """
        cross_mature_mir_gene_names = []
        """: type : list[str] """
        cross_mature_mir_gene_names_count = []
        """: type : list[int] """
        distal_mature_mir_gene_names = []
        """: type : list[str] """
        distal_mature_mir_gene_names_count = []
        """: type : list[int] """
        
        target_genes_dict_inside = {}
        target_genes_dict_cross = {}
        target_genes_dict_distal = {}
        
        for row in self.out_dataframe.itertuples():
            miR_inside = row.mirna_inside
            target_genes_dict_inside = self.__get_genetarget(miR_inside, mirbase_info, target_info)
            target_genes_names = list(itertools.chain.from_iterable(list(target_genes_dict_inside.values())))
            target_genes_names = list(set(target_genes_names) if unique else target_genes_names)
            if target_genes_names != ['.']:
                target_genes_names = list(filter(lambda a: a != '.', target_genes_names))
                
            inside_mature_mir_gene_names.append(";".join(target_genes_names))
            inside_mature_mir_gene_names_count.append(len(target_genes_names) if list(target_genes_names)[0] != "." else 0)

            miR_cross = row.mirna_cross
            target_genes_dict_cross = self.__get_genetarget(miR_cross, mirbase_info, target_info)
            target_genes_names = list(itertools.chain.from_iterable(list(target_genes_dict_cross.values())))
            target_genes_names = list(set(target_genes_names) if unique else target_genes_names)
            if target_genes_names != ['.']:
                target_genes_names = list(filter(lambda a: a != '.', target_genes_names))
            cross_mature_mir_gene_names.append(";".join(target_genes_names))
            cross_mature_mir_gene_names_count.append(len(target_genes_names) if list(target_genes_names)[0] != "." else 0)

            miR_distal = row.mirna_distal
            target_genes_dict_distal = self.__get_genetarget(miR_distal, mirbase_info, target_info)
            target_genes_names = list(itertools.chain.from_iterable(list(target_genes_dict_distal.values())))
            target_genes_names = list(set(target_genes_names) if unique else target_genes_names)
            if target_genes_names != ['.']:
                target_genes_names = list(filter(lambda a: a != '.', target_genes_names))
            distal_mature_mir_gene_names.append(";".join(target_genes_names))
            distal_mature_mir_gene_names_count.append(len(target_genes_names) if list(target_genes_names)[0] != "." else 0)

        self.out_dataframe.loc[:, target + '_inside'] = inside_mature_mir_gene_names
        self.out_dataframe.loc[:, target + '_inside_count'] = inside_mature_mir_gene_names_count
        self.out_dataframe.loc[:, target + '_cross'] = cross_mature_mir_gene_names
        self.out_dataframe.loc[:, target + '_cross_count'] = cross_mature_mir_gene_names_count
        self.out_dataframe.loc[:, target + '_distal'] = distal_mature_mir_gene_names
        self.out_dataframe.loc[:, target + '_distal_count'] = distal_mature_mir_gene_names_count
        
        return target_genes_dict_inside, target_genes_dict_cross, target_genes_dict_distal

    def add_meta_gene(self, genelist):
        db = self.connection['genelists'][genelist]
        genelist_l = set(db.find().distinct('gene'))
        print("Genelist")
        genes_inside = []
        genes_cross = []
        genes_distal = []
        genes_inside_count = []
        genes_cross_count = []
        genes_distal_count = []
        genes_inside_coords = []
        genes_cross_coords = []
        genes_distal_coords = []
        
        for row in self.out_dataframe.itertuples():
            querylist = set([g.split(':')[0] for g in row.gene_inside.split(';')])
            #This retrieves the gene name in the genelist. It populates a list with gene names only.
            inside_no_transcript = ';'.join(list(querylist & genelist_l))

            #This explodes each gene into all the transcripts
            inside = []
            for r in row.gene_inside.split(';'):
                    print("Q:",r)
                    if r.split(':')[0] in inside_no_transcript:
                        inside.append(r)
            inside_coords = []
            if inside:
                for i in inside:
                    inside_coords.append(row.gene_inside_coords.split(';')[row.gene_inside.split(';').index(i)])
                inside_coords = ';'.join(inside_coords)
            genes_inside.append(';'.join(inside)) if inside else genes_inside.append('.')
            genes_inside_coords.append(inside_coords) if inside_coords else genes_inside_coords.append('.')
            genes_inside_count.append(str(len(inside))) if inside else genes_inside_count.append('0')


            querylist = set([g.split(':')[0] for g in row.gene_cross.split(';')])
            #This retrieves the gene name in the genelist. It populates a list with gene names only.
            cross_no_transcript = ';'.join(list(querylist & genelist_l))

            #This explodes each gene into all the transcripts
            cross = []
            for r in row.gene_cross.split(';'):
                    print("Q:",r)
                    if r.split(':')[0] in cross_no_transcript:
                        cross.append(r)
            cross_coords = []
            if cross:
                for i in cross:
                    cross_coords.append(row.gene_cross_coords.split(';')[row.gene_cross.split(';').index(i)])
                cross_coords = ';'.join(cross_coords)
            genes_cross.append(';'.join(cross)) if cross else genes_cross.append('.')
            genes_cross_coords.append(cross_coords) if cross_coords else genes_cross_coords.append('.')
            genes_cross_count.append(str(len(cross))) if cross else genes_cross_count.append('0')



            querylist = set([g.split(':')[0] for g in row.gene_distal.split(';')])
            #This retrieves the gene name in the genelist. It populates a list with gene names only.
            distal_no_transcript = ';'.join(list(querylist & genelist_l))

            #This explodes each gene into all the transcripts
            distal = []
            for r in row.gene_distal.split(';'):
                    print("Q:",r)
                    if r.split(':')[0] in distal_no_transcript:
                        distal.append(r)
            distal_coords = []
            if distal:
                for i in distal:
                    distal_coords.append(row.gene_distal_coords.split(';')[row.gene_distal.split(';').index(i)])
                distal_coords = ';'.join(distal_coords)
            genes_distal.append(';'.join(distal)) if distal else genes_distal.append('.')
            genes_distal_coords.append(distal_coords) if distal_coords else genes_distal_coords.append('.')
            genes_distal_count.append(str(len(distal))) if distal else genes_distal_count.append('0')

            
        self.out_dataframe.loc[:, genelist + '_genelist_inside'] = genes_inside
        self.out_dataframe.loc[:, genelist + '_genelist_inside_coords'] = genes_inside_coords
        self.out_dataframe.loc[:, genelist + '_genelist_inside_count'] = genes_inside_count
        self.out_dataframe.loc[:, genelist + '_genelist_cross'] = genes_cross
        self.out_dataframe.loc[:, genelist + '_genelist_cross_coords'] = genes_cross_coords
        self.out_dataframe.loc[:, genelist + '_genelist_cross_count'] = genes_cross_count
        self.out_dataframe.loc[:, genelist + '_genelist_distal'] = genes_distal
        self.out_dataframe.loc[:, genelist + '_genelist_distal_coords'] = genes_distal_coords
        self.out_dataframe.loc[:, genelist + '_genelist_distal_count'] = genes_distal_count

    
def write_file(cnv_infolist: DataFrame, mirbase_dict: dict, out_filename: str, extra_info: dict):
    """
    Write a DataFrame to excel
    :param DataFrame cnv_infolist: Annotated pandas DataFrame to be written to xlsx file
    :param out_filename: File name of the final xlsx file
    """
    print(cnv_infolist)
    print(extra_info.keys())
    # Writing Excel file
    writer = pd.ExcelWriter(out_filename, engine='xlsxwriter')
    cnv_infolist.to_excel(writer, sheet_name='Annotated CNV - ' + time.strftime("%d-%m-%Y"), startrow=1,
                          header=False, index=False)
    
    workbook = writer.book
    worksheet = writer.sheets['Annotated CNV - ' + time.strftime("%d-%m-%Y")]
    worksheet.freeze_panes(1, 0)
    header_format = workbook.add_format({
        'bold': True,
        # 'text_wrap': True,
        'valign': 'top',
        'font_color': 'red',
        # 'fg_color': '#ffca6f',
        'border': 0})
    
    for col_num, value in enumerate(cnv_infolist.columns.values):
        worksheet.write(0, col_num, value, header_format)
    
    writer.save()
    
    # Writing JSON file
    json_data = cnv_infolist.reset_index().to_json(orient='records')
    with open(re.sub('.xlsx', '.json', out_filename), 'w') as f:
        f.write(json_data)

    # Writing CSV file
    cnv_infolist.to_csv(re.sub('.xlsx', '.csv', out_filename))
    
    #mirBase section (if --mirbase selected)
    if mirbase_dict != {}:
        with open(re.sub('.xlsx', '_mirbase.json', out_filename), 'w') as f:
            json.dump(mirbase_dict, f)

    # Writing JSON with extra info (for analyses that provide some)
    with open(re.sub('.xlsx', '_extra.json', out_filename), 'w') as f:
        json.dump(extra_info, f)


    
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--cnv-file", help="CNV file to be annotated")
    parser.add_argument("--out", required=True, help="Annotated file")
    parser.add_argument("--cnv-line", help="Single CNV to be annotated (E.g.: 'chr1:11110-11150'")
    
    parser.add_argument("--gene", action='store_true', required=False,
                        help="BED file of all RefSeq genes")
    parser.add_argument("--coding_gene", action='store_true', required=False,
                        help="BED file of all coding RefSeq genes")
    parser.add_argument("--noncoding_gene", action='store_true', required=False,
                        help="BED file of all non-coding RefSeq genes")
    parser.add_argument("--mirna", action='store_true', required=False, help="BED file of known miRNAs")
    parser.add_argument("--longNC", action='store_true', required=False, help="BED file of known long "
                                                                              "Non-coding molecules")
    parser.add_argument("--circRNA", action='store_true', required=False, help="BED file of known "
                                                                               "circular RNA molecules")
    parser.add_argument("--pseudogene", action='store_true', required=False,
                        help="BED file of known pseudogenes from GENECODE")
    parser.add_argument("--mirbase", action='store_true', required=False, help="miRBase file")
    parser.add_argument("--enhancer", action='store_true', required=False, help="miRBase file")

    parser.add_argument("--ucr", action='store_true', required=False, help="Ultra Conserved regions")
    parser.add_argument("--har", action='store_true', required=False, help="Hypervariable regions")

    parser.add_argument("--all", dest='all_beds', action='store_true', required=False,
                        help="Perform all available annotations")
    
    parser.add_argument("--all-genelists", action='store_true', required=False,
                        help="Perform all available gene classifications")
    # parser.add_argument("--ASD-genelist", dest='ASD_genelist', action='store_true', required=False,
    #                     help="Perform ASD gene classification")
    parser.add_argument("--ID-genelist", dest='ID_genelist', action='store_true', required=False,
                        help="Perform ID gene classification")
    # parser.add_argument("--IDb-genelist", dest='ID_b_genelist', action='store_true', required=False,
    #                     help="Perform ID_b gene classification")
    parser.add_argument("--dosage-sensitive-genelist", dest='dosage_sensitive_genelist', action='store_true', required=False,
                        help="Perform dosage-sensitive gene classification")
    # parser.add_argument("--epilepsy-genelist", dest='epilessia_genelist', action='store_true', required=False,
    #                     help="Perform epilepsy gene classification")
    # parser.add_argument("--malformations-genelist", dest='malformazioni_genelist', action='store_true', required=False,
    #                     help="Perform malformations gene classification")
    parser.add_argument("--mendeliome-genelist", dest='mendeliome_genelist', action='store_true', required=False,
                        help="Perform mendeliome gene classification")
    parser.add_argument("--ohnologs-genelist", dest='ohnologs_genelist', action='store_true', required=False,
                        help="Perform ohnolog genes classification")
    parser.add_argument("--imprinted-genelist", dest='imprinted_genelist', action='store_true', required=False,
                        help="Perform imprinted genes classification")
    # parser.add_argument("--pubmed-autism-genelist", dest='pubmed_autism_genelist', action='store_true', required=False,
    #                     help="Perform pubmed autism gene classification")
    # parser.add_argument("--pubmed-brain-genelist", dest='pubmed_brain_malformations_genelist',
    #                     action='store_true', required=False,
    #                     help="Perform pubmed brain malformations gene classification")
    # parser.add_argument("--pubmed-epilepsy-genelist", dest='pubmed_epilepsy_or_seizures_genelist',
    #                     action='store_true', required=False,
    #                     help="Perform pubmed epilepsy gene classification")
    # parser.add_argument("--pubmed-ID-genelist", dest='pubmed_intellectual_disability_genelist',
    #                     action='store_true', required=False,
    #                     help="Perform pubmed intellectual disability gene classification")

    parser.add_argument("-D", "--distance", type=int, default=1000000, required=False,
                        help="Distance from gene (Default 1Mb)")
    parser.add_argument("-r", "--reference", default='hg19', required=False,
                        help="Reference genome", choices=['hg19', 'hg18', 'hg38']) #hg38 coming soon
    args = parser.parse_args()
    if args.cnv_file and args.cnv_line:
        sys.exit("Input line(s) can be provided as --cnv_line OR --cnv_file. "
                 "They cannot be specified together.")
    elif not (args.cnv_file or args.cnv_line):
        sys.exit("Please provide input as either --cnv_line or --cnv_file.")
    success = MainApp(args).process()
    if success == -1:
        sys.stdout.write("There was a problem in annotating CNVs.")
    print(success)
