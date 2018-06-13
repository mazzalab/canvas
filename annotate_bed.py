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
        
    def process(self):
        
        for arg in self.d:
            
            if arg not in ['cnv_line', 'cnv_file', 'out', 'distance', 'all_beds', 'all_genelists'] and \
                    (self.d[arg] or self.d['all_beds']) and not arg.startswith('__'):
                print("Checking presence of {}...".format(arg))
                
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
        if self.args.cnv_file:
            if os.path.exists(self.args.cnv_file):
                cnv_info = self.read_cnv_coordinates_file(self.args.cnv_file)
            else:
                sys.exit("CNV file not found: {}".format(self.args.cnv_file))
        elif self.args.cnv_line:
            cnv_info = self.read_cnv_coordinates_line(self.args.cnv_line)
        
        self.out_dataframe = cnv_info

        # Annotation based on the options that were chosen
        for arg in self.d:
            if arg not in ['cnv_line', 'cnv_file', 'out', 'distance', 'all_beds', 'all_genelists'] and \
                    (self.d[arg] or self.d['all_beds']) and not arg.startswith('__') and not 'genelist' in arg:
                print("Adding annotation for", arg)
                if arg != 'mirbase':
                    self.add_annotation(arg)
                else:
                    for target_db in self.db_t.collection_names():
                        self.mirbase_dict[target_db] = []  # initialized the target key (tarbase or targetscan) in the mirbase dict
                        
                        dict_inside, dict_cross, dict_distal = self.add_mirna_target(target=target_db, unique=True)
                        self.mirbase_dict[target_db].extend((dict_inside, dict_cross, dict_distal))
        
        # Genelists annotation
        print("Adding genelists classifications...")
        if self.d['gene'] or self.d['all_genelists']:
            for name in sorted(self.db_gl.collection_names()):
                if self.d[name+'_genelist'] or self.d['all_genelists']:
                    print("Adding {} gene classification...".format(name))
                    self.add_meta_gene(name)
        else:
            sys.stderr.write("WARNING: Could not add genelists annotation since --gene option was not included.\n")
        
        #Writing final file
        write_file(self.out_dataframe, self.mirbase_dict, self.args.out)
        # return self.out_dataframe.reset_index().to_json(orient='records')
    
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
        
        m = re.match(r'(?P<chr>chr[0-9]+):(?P<start>\d+)-(?P<end>\d+)', cnv_line)
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
        print('Annotating', annotation_db)
        distance_from_gene = self.args.distance
        
        inside_molecules = []
        """:type : list[str]"""
        inside_molecules_count = []
        """:type : list[int]"""
        cross_molecules = []
        """:type : list[str]"""
        cross_molecules_count = []
        """:type : list[int]"""
        distal_molecules = []
        """:type : list[str]"""
        distal_molecules_count = []
        """:type : list[int]"""

        for row in self.out_dataframe.itertuples():
            chrom = row.CHR
            start = row.START
            end = row.END
            inside_molecule = db.find({"$and": [{"chr": chrom, "start": {"$gte": start},
                                                 "end": {"$lte": end}}]}).distinct('info')
            
            cross_molecule = db.find(
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
            ).distinct('info')

            distal_molecule = db.find(
                {"$or": [
                    {"$and": [
                        {"chr": chrom}, {"start": {"$gt": end}}, {"start": {"$lt": distance_from_gene+end}},
                        {"end": {"$lt": distance_from_gene + end}}, {"end": {"$gt": end}}
                    ]},
                    {"$and": [
                        {"chr": chrom}, {"end": {"$lt": start}}, {"end": {"$gt": start-distance_from_gene}},
                        {"start": {"$gt": start - distance_from_gene}}, {"start": {"$lt": start}}
                    ]}
                ]}
            ).distinct('info')
            
            if len(inside_molecule) > 0:
                inside_molecules.append(",".join(inside_molecule))
                inside_molecules_count.append(len(inside_molecule))
            else:
                inside_molecules.append(".")
                inside_molecules_count.append(0)

            if len(cross_molecule) > 0:
                cross_molecules.append(",".join(cross_molecule))
                cross_molecules_count.append(len(cross_molecule))

            else:
                cross_molecules.append(".")
                cross_molecules_count.append(0)

            if len(distal_molecule) > 0:
                distal_molecules.append(",".join(distal_molecule))
                distal_molecules_count.append(len(distal_molecule))
            else:
                distal_molecules.append(".")
                distal_molecules_count.append(0)

        self.out_dataframe.loc[:, annotation_db + '_inside'] = inside_molecules
        self.out_dataframe.loc[:, annotation_db + '_inside_count'] = inside_molecules_count
        self.out_dataframe.loc[:, annotation_db + '_cross'] = cross_molecules
        self.out_dataframe.loc[:, annotation_db + '_cross_count'] = cross_molecules_count
        self.out_dataframe.loc[:, annotation_db + '_distal'] = distal_molecules
        self.out_dataframe.loc[:, annotation_db + '_distal_count'] = distal_molecules_count
        
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
                mature_mir_gene_names_dict[mi_name] = gene_targets
            else:
                mature_mir_gene_names_dict[mi_name] = ['.']
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
        
        db_target = self.connection['targets'][target]
        target_info = {}
        
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
                
            inside_mature_mir_gene_names.append(",".join(target_genes_names))
            inside_mature_mir_gene_names_count.append(len(target_genes_names) if list(target_genes_names)[0] != "." else 0)

            miR_cross = row.mirna_cross
            target_genes_dict_cross = self.__get_genetarget(miR_cross, mirbase_info, target_info)
            target_genes_names = list(itertools.chain.from_iterable(list(target_genes_dict_cross.values())))
            target_genes_names = list(set(target_genes_names) if unique else target_genes_names)
            if target_genes_names != ['.']:
                target_genes_names = list(filter(lambda a: a != '.', target_genes_names))
            cross_mature_mir_gene_names.append(",".join(target_genes_names))
            cross_mature_mir_gene_names_count.append(len(target_genes_names) if list(target_genes_names)[0] != "." else 0)

            miR_distal = row.mirna_distal
            target_genes_dict_distal = self.__get_genetarget(miR_distal, mirbase_info, target_info)
            target_genes_names = list(itertools.chain.from_iterable(list(target_genes_dict_distal.values())))
            target_genes_names = list(set(target_genes_names) if unique else target_genes_names)
            if target_genes_names != ['.']:
                target_genes_names = list(filter(lambda a: a != '.', target_genes_names))
            distal_mature_mir_gene_names.append(",".join(target_genes_names))
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
        
        genes_inside = []
        genes_cross = []
        genes_distal = []
        genes_inside_count = []
        genes_cross_count = []
        genes_distal_count = []
        
        for row in self.out_dataframe.itertuples():
            inside = ','.join(list(set(row.gene_inside.split(',')) & genelist_l))
            genes_inside.append(inside) if inside else genes_inside.append('.')
            genes_inside_count.append(str(len(inside.split(',')))) if inside else genes_inside_count.append('0')

            cross = ','.join(list(set(row.gene_cross.split(',')) & genelist_l))
            genes_cross.append(cross) if cross else genes_cross.append('.')
            genes_cross_count.append(str(len(cross.split(',')))) if cross else genes_cross_count.append('0')

            distal = ','.join(list(set(row.gene_distal.split(',')) & genelist_l))
            genes_distal.append(distal) if distal else genes_distal.append('.')
            genes_distal_count.append(str(len(distal.split(',')))) if distal else genes_distal_count.append('0')
            
        self.out_dataframe.loc[:, genelist + '_inside'] = genes_inside
        self.out_dataframe.loc[:, genelist + '_inside_count'] = genes_inside_count
        self.out_dataframe.loc[:, genelist + '_cross'] = genes_cross
        self.out_dataframe.loc[:, genelist + '_cross_count'] = genes_cross_count
        self.out_dataframe.loc[:, genelist + '_distal'] = genes_distal
        self.out_dataframe.loc[:, genelist + '_distal_count'] = genes_distal_count

    
def write_file(cnv_infolist: DataFrame, mirbase_dict: dict, out_filename: str):
    """
    Write a DataFrame to excel
    :param DataFrame cnv_infolist: Annotated pandas DataFrame to be written to xlsx file
    :param out_filename: File name of the final xlsx file
    """
    
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
    json = cnv_infolist.reset_index().to_json(orient='records')
    with open(re.sub('.xlsx', '.json', out_filename), 'w') as f:
        f.write(json)

    # Writing CSV file
    cnv_infolist.to_csv(re.sub('.xlsx', '.csv', out_filename))
    
    #mirBase section (if --mirbase selected)
    if mirbase_dict != {}:
        with open(re.sub('.xlsx', '_mirbase.json', out_filename), 'w') as f:
            import json
            json.dump(mirbase_dict, f)
    
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
    parser.add_argument("--all", dest='all_beds', action='store_true', required=False,
                        help="Perform all available annotations")
    
    parser.add_argument("--all-genelists", action='store_true', required=False,
                        help="Perform all available gene classifications")
    parser.add_argument("--ASD-genelist", dest='ASD_genelist', action='store_true', required=False,
                        help="Perform ASD gene classification")
    parser.add_argument("--IDa-genelist", dest='ID_a_genelist', action='store_true', required=False,
                        help="Perform ID_a gene classification")
    parser.add_argument("--IDb-genelist", dest='ID_b_genelist', action='store_true', required=False,
                        help="Perform ID_b gene classification")
    parser.add_argument("--dosage-sensitive-genelist", dest='dosage_sensitive_genelist', action='store_true', required=False,
                        help="Perform dosage-sensitive gene classification")
    parser.add_argument("--epilepsy-genelist", dest='epilessia_genelist', action='store_true', required=False,
                        help="Perform epilepsy gene classification")
    parser.add_argument("--malformations-genelist", dest='malformazioni_genelist', action='store_true', required=False,
                        help="Perform malformations gene classification")
    parser.add_argument("--mendeliome-genelist", dest='mendeliome_genelist', action='store_true', required=False,
                        help="Perform mendeliome gene classification")
    parser.add_argument("--onologs-genelist", dest='onologhi_genelist', action='store_true', required=False,
                        help="Perform onologs gene classification")
    parser.add_argument("--pubmed-autism-genelist", dest='pubmed_autism_genelist', action='store_true', required=False,
                        help="Perform pubmed autism gene classification")
    parser.add_argument("--pubmed-brain-genelist", dest='pubmed_brain_malformations_genelist',
                        action='store_true', required=False,
                        help="Perform pubmed brain malformations gene classification")
    parser.add_argument("--pubmed-epilepsy-genelist", dest='pubmed_epilepsy_or_seizures_genelist',
                        action='store_true', required=False,
                        help="Perform pubmed epilepsy gene classification")
    parser.add_argument("--pubmed-ID-genelist", dest='pubmed_intellectual_disability_genelist',
                        action='store_true', required=False,
                        help="Perform pubmed intellectual disability gene classification")

    parser.add_argument("-D", "--distance", type=int, default=1000000, required=False,
                        help="Distance from gene (Default 1Mb)")
        
    args = parser.parse_args()
    if args.cnv_file and args.cnv_line:
        sys.exit("Input line(s) can be provided as --cnv_line OR --cnv_file. "
                 "They cannot be specified together.")
    elif not (args.cnv_file or args.cnv_line):
        sys.exit("Please provide input as either --cnv_line or --cnv_file.")
    MainApp(args).process()
