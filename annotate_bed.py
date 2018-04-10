""""
It annotates a tab-delimited input file with a set of BED files, with textual features in their 4th columns
"""

import argparse
import re
import time
import collections
import pandas as pd
from pandas import DataFrame

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



def read_cnv_coordinates(cnv_file: str) -> DataFrame:
    """
    Read and annotate the original CNV file as a DataFrame
    :param cnv_file: File path and name of the original CNV file to be annotated
    :return: A DataFrame containing the CNV to be annotate, one per line
    """

    cnv_coords = pd.read_table(cnv_file, encoding='cp1252')
    # cnv_coords = cnv_coords[['CHR', 'START', 'END']]
    # cnv_coords = cnv_coords[(pd.isnull(cnv_coords.START)) & (pd.isnull(cnv_coords.END)) & (pd.isnull(cnv_coords.CHR))]
    cnv_coords = cnv_coords[cnv_coords.START.notnull() & cnv_coords.END.notnull() & cnv_coords.CHR.notnull()]

    cnv_coords['START'] = cnv_coords['START'].astype('int')
    cnv_coords['END'] = cnv_coords['END'].astype('int')
    # cnv_coords = cnv_coords.reset_index(drop=True)

    return cnv_coords


def add_annotation(cnv_tobe_annotated: DataFrame, annotation_bedfile: str, column_name_suffix: str) -> DataFrame:
    """
    Take a DataFrame and add in the last six columns the annotation provided in the 4th column of the BED file
    :param str column_name_suffix: Name of the six columns to be added
    :param str annotation_bedfile: File path and name of the annotation BED file
    :param DataFrame cnv_tobe_annotated: DataFrame to be annotated
    :return: A DataFrame containing the input DataFrame with the new annotations as last 6 columns
    """

    distance_from_gene = 1000000

    annotation_info = []
    """:type : list[list[str]] """
    with open(annotation_bedfile)as f:
        for line in f:
            L = line.strip().split("\t")
            if L[0] != "chr":
                annotation_info.append(L)

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
    cnv_annotated = cnv_tobe_annotated
    """:type : DataFrame """
    for index, row in cnv_tobe_annotated.iterrows():
        chrom = row['CHR']
        start = row['START']
        end = row['END']

        inside_molecule = []
        """:type : list[str]"""

        cross_molecule = []
        """:type : list[str]"""

        distal_molecule = []
        """:type : list[str]"""

        for annotation_coord in annotation_info:
            """:type : list[str]"""
            chrom_bed = annotation_coord[0]

            try:
                start_bed = int(annotation_coord[1])
                end_bed = int(annotation_coord[2])
            except ValueError:
                print("Oops! The BED file does not contain valid genomic coordinates. Let's check!")
                break

            if chrom == chrom_bed and start <= start_bed and end >= end_bed:
                # let's check this fourth field, which must contain the feature to be added (e.g., miRNA Symbol)
                inside_molecule.append(annotation_coord[3])

            elif (chrom == chrom_bed and start <= start_bed <= end and end < end_bed) or \
                    (chrom == chrom_bed and start_bed <= start <= end_bed and end > end_bed) or (
                        chrom == chrom_bed and start_bed < start < end_bed and start_bed < end < end_bed):
                cross_molecule.append(annotation_coord[3])

            elif (chrom == chrom_bed and end < start_bed and start_bed - end < distance_from_gene) or \
                    (chrom == chrom_bed and start > end_bed and start - end_bed < distance_from_gene):
                distal_molecule.append(annotation_coord[3])

        if len(inside_molecule) > 0:
            inside_molecule = set(inside_molecule)
            inside_molecules.append(",".join(inside_molecule))
            inside_molecules_count.append(len(inside_molecule))
        else:
            inside_molecules.append(".")
            inside_molecules_count.append(0)

        if len(cross_molecule) > 0:
            cross_molecule = set(cross_molecule)
            cross_molecules.append(",".join(cross_molecule))
            cross_molecules_count.append(len(cross_molecule))
        else:
            cross_molecules.append(".")
            cross_molecules_count.append(0)

        if len(distal_molecule) > 0:
            distal_molecule = set(distal_molecule)
            distal_molecules.append(",".join(distal_molecule))
            distal_molecules_count.append(len(distal_molecule))
        else:
            distal_molecules.append(".")
            distal_molecules_count.append(0)

    cnv_annotated.loc[:, column_name_suffix + '_inside'] = inside_molecules
    cnv_annotated.loc[:, column_name_suffix + '_inside_count'] = inside_molecules_count
    cnv_annotated.loc[:, column_name_suffix + '_cross'] = cross_molecules
    cnv_annotated.loc[:, column_name_suffix + '_cross_count'] = cross_molecules_count
    cnv_annotated.loc[:, column_name_suffix + '_distal'] = distal_molecules
    cnv_annotated.loc[:, column_name_suffix + '_distal_count'] = distal_molecules_count

    return cnv_annotated


def __get_genetarget(mirs_in_cnv: str, mirbase_dict: dict, target_dict: dict) -> list:
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

    return mature_mir_gene_names if len(mature_mir_gene_names) > 0 else ["."]


def add_mirna_target(cnv_tobe_annotated: DataFrame, mirbase_file: str, db_target_file: str,
                     mir_gene_index: collections.namedtuple, unique: bool,
                     miR_colname_suffix: str,
                     target_colname_suffix: str) -> DataFrame:
    """
    Take a DataFrame and annotate it with targeting genes of miRs according to a db_target_file
    :param target_colname_suffix: Suffix string for the column names containing gene target (inside, cross and distal)
    :param namedtuple mir_gene_index: 0-based column indices of miR and Gene targets
    :param DataFrame cnv_tobe_annotated: DataFrame to be annotated
    :param str mirbase_file: Original flat file of miRBase
    :param str db_target_file: Flat file of a database reporting "Gene Symbol" and targeting "miRNA" plus an header
    :param bool unique: If true, only unique genes will be reported for each CNV region
    :param str miR_colname_suffix: Suffix string for the column name containing annotated miRs in the CNVs
    :return DataFrame: Annotated DataFrame
    """

    mirbase_info = {}
    """:type : dict[str, list[str]] """
    with open(mirbase_file)as mb:
        for line in mb:
            if not line.startswith("#"):
                mbL = line.rstrip().split("\t")
                if mbL[2] != "miRNA_primary_transcript":
                    m = re.match(r'.+;Name=(?P<mirna_name>[A-Za-z0-9\-]+);Derives_from=(?P<mi>MI[0-9]+)$', mbL[8])
                    if not m.group('mi') in mirbase_info:
                        mirbase_info[m.group('mi')] = [m.group('mirna_name')]
                    else:
                        mirbase_info[m.group('mi')].append(m.group('mirna_name'))

    target_info = {}
    """:type : dict[str, list[str]] """
    with open(db_target_file)as t:
        t.readline()  # Skip header line
        mir_index = mir_gene_index.miR
        gene_index = mir_gene_index.gene

        for line in t:
            tL = line.rstrip().split("\t")

            if not tL[mir_index] in target_info:
                target_info[tL[mir_index]] = [tL[gene_index]]
            else:
                target_info[tL[mir_index]].append(tL[gene_index])

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
    for index, row in cnv_tobe_annotated.iterrows():
        miR_inside = row[miR_colname_suffix + '_inside']
        target_genes = __get_genetarget(miR_inside, mirbase_info, target_info)
        target_genes = set(target_genes) if unique else target_genes
        inside_mature_mir_gene_names.append(",".join(target_genes))
        inside_mature_mir_gene_names_count.append(len(target_genes) if list(target_genes)[0] != "." else 0)

        miR_cross = row[miR_colname_suffix + '_cross']
        target_genes = __get_genetarget(miR_cross, mirbase_info, target_info)
        target_genes = set(target_genes) if unique else target_genes
        cross_mature_mir_gene_names.append(",".join(target_genes))
        cross_mature_mir_gene_names_count.append(len(target_genes) if list(target_genes)[0] != "." else 0)

        miR_distal = row[miR_colname_suffix + '_distal']
        target_genes = __get_genetarget(miR_distal, mirbase_info, target_info)
        target_genes = set(target_genes) if unique else target_genes
        distal_mature_mir_gene_names.append(",".join(target_genes))
        distal_mature_mir_gene_names_count.append(len(target_genes) if list(target_genes)[0] != "." else 0)

    cnv_annotated = cnv_tobe_annotated
    cnv_annotated.loc[:, target_colname_suffix + '_inside'] = inside_mature_mir_gene_names
    cnv_annotated.loc[:, target_colname_suffix + '_inside_count'] = inside_mature_mir_gene_names_count
    cnv_annotated.loc[:, target_colname_suffix + '_cross'] = cross_mature_mir_gene_names
    cnv_annotated.loc[:, target_colname_suffix + '_cross_count'] = cross_mature_mir_gene_names_count
    cnv_annotated.loc[:, target_colname_suffix + '_distal'] = distal_mature_mir_gene_names
    cnv_annotated.loc[:, target_colname_suffix + '_distal_count'] = distal_mature_mir_gene_names_count

    return cnv_annotated


def write_file(cnv_infolist: DataFrame, out_filename: str):
    """
    Write a DataFrame to excel
    :param DataFrame cnv_infolist: Annotated pandas DataFrame to be written to xlsx file
    :param out_filename: File name of the final xlsx file
    """

    # cnv_infolist.to_excel(out_filename, sheet_name='Annotated CNV', header=True)
    writer = pd.ExcelWriter(out_filename, engine='xlsxwriter')
    cnv_infolist.to_excel(writer, sheet_name='Annotated CNV - '+time.strftime("%d-%m-%Y"), startrow=1, header=False, index=False)

    workbook = writer.book
    worksheet = writer.sheets['Annotated CNV - '+time.strftime("%d-%m-%Y")]
    # worksheet.set_column('A:D', 10)
    # worksheet.set_column('E:ZZ', 15)
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--cnv", required=True, help="CNV file to be annotated")
    parser.add_argument("--out", required=True, help="Annotated file")

    parser.add_argument("--gene", required=False, help="BED file of all RefSeq genes")
    parser.add_argument("--coding_gene", required=False, help="BED file of all coding RefSeq genes")
    parser.add_argument("--noncoding_gene", required=False, help="BED file of all non-coding RefSeq genes")
    parser.add_argument("--mirna", required=False, help="BED file of known miRNAs")
    parser.add_argument("--longNC", required=False, help="BED file of known long Non-coding molecules")
    parser.add_argument("--circRNA", required=False, help="BED file of known circular RNA molecules")
    parser.add_argument("--pseudogene", required=False, help="BED file of known pseudogenes from GENECODE")
    parser.add_argument("--mirbase", required=False, help="miRBase file")
    parser.add_argument("--targetscan", required=False, help="TargetScan file")
    parser.add_argument("--tarbase", required=False, help="Tarbase file")

    args = parser.parse_args()
    cnv_bedfile = args.cnv
    gene_bedfile = args.gene
    coding_gene_bedfile = args.coding_gene
    noncoding_gene_bedfile = args.noncoding_gene
    mirna_bedfile = args.mirna
    longNC_bedfile = args.longNC
    circRNA_bedfile = args.circRNA
    preudogene_bedfile = args.pseudogene
    mirbase_file = args.mirbase
    targetscan_file = args.targetscan
    tarbase_file = args.tarbase
    out_file = args.out

    cnv_info = read_cnv_coordinates(cnv_bedfile)

    out_dataframe = cnv_info
    """: type : DataFrame """
    if gene_bedfile:
        out_dataframe = add_annotation(out_dataframe, gene_bedfile, "gene")
    if coding_gene_bedfile:
        out_dataframe = add_annotation(out_dataframe, coding_gene_bedfile, "coding_gene")
    if noncoding_gene_bedfile:
        out_dataframe = add_annotation(out_dataframe, noncoding_gene_bedfile, "noncoding_gene")
    if preudogene_bedfile:
        out_dataframe = add_annotation(out_dataframe, preudogene_bedfile, "pseudogene")
    if mirna_bedfile:
        out_dataframe = add_annotation(out_dataframe, mirna_bedfile, "miR")
    if mirbase_file and targetscan_file or tarbase_file:
        mir_gene_ind = collections.namedtuple('mir_gene_ind', 'miR gene')
        out_dataframe = add_mirna_target(out_dataframe, mirbase_file, targetscan_file, mir_gene_ind(miR=1, gene=0),
                                         True, "miR", "Targetscan")
        out_dataframe = add_mirna_target(out_dataframe, mirbase_file, tarbase_file, mir_gene_ind(miR=2, gene=1),
                                         True, "miR", "Tarbase")
    if longNC_bedfile:
        out_dataframe = add_annotation(out_dataframe, longNC_bedfile, "lnc")
    if circRNA_bedfile:
        out_dataframe = add_annotation(out_dataframe, circRNA_bedfile, "circ")

    write_file(out_dataframe, out_file)
