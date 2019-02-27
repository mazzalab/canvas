""""
File description
"""

import argparse
from collections import Counter

import pandas as pd

__author__ = "Tommaso Mazza"
__copyright__ = "Copyright 2017, The pairwiseCNV Project"
__version__ = "0.0.1"
__email__ = "bioinformatics@css-mendel.it"
__status__ = "Development"
__date__ = "17/10/2017"
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


def percent_overlap(row1, row2, padding: int) -> (float, float):
    """
    Calculate the overlap between CNV in row1 and in row2. Eventually add a padding region to both sizes of the CNVs.
    :param row1: Row containing the first CNV to be intersected
    :param row2: Row containing the second CNV to be intersected
    :param padding: Padding region to be added to both sizes of the CNVs
    :return: A tuple of percentages of overlap between both CNVs, normalized on the length of each CNV
    """
    chr1 = row1[1]
    start1 = row1[2] - padding
    end1 = row1[3] + padding

    chr2 = row2[1]
    start2 = row2[2] - padding
    end2 = row2[3] + padding

    if chr1 != chr2:
        return 0, 0
    else:
        row1_range = range(start1, end1)
        row2_range = range(start2, end2)
        row1_range_set = set(row1_range)
        row2_range_set = set(row2_range)

        inters = row1_range_set.intersection(row2_range_set)

        if not inters:
            return 0, 0
        else:
            extreme_inters = max(inters) - min(inters) + 1
            return round((extreme_inters / (end1 - start1)), 2), round((extreme_inters / (end2 - start2)), 2)


def pairwise_match(cnv_file: str, padding: int) -> object:
    """

    :param cnv_file:
    :param padding:
    :return:
    """
    raw_table = pd.read_table(cnv_file, usecols=['PAZIENTE', 'CHR', 'START', 'END'], encoding="latin1")
    raw_table = raw_table[raw_table['START'] != "."]
    raw_table = raw_table[raw_table['END'] != "."]
    raw_table = raw_table[raw_table['CHR'] != "."]
    raw_table = raw_table.reset_index(drop=True)

    # Casting START and END columns to INT32
    raw_table['START'] = raw_table['START'].astype(int)
    raw_table['END'] = raw_table['END'].astype(int)

    # Making indies unique appending to replicate patient names an incremental number
    patients_names = raw_table['PAZIENTE'].tolist()
    patients_unique_names = [s + str(suffix) if num > 1 else s for s, num in Counter(patients_names).items() for suffix in
                       range(1, num + 1)]
    out_matrix = pd.DataFrame(columns=patients_unique_names, index=patients_unique_names)

    for outer_index, outer_row in raw_table.iterrows():
        for inner_index, inner_row in raw_table.iloc[outer_index:].iterrows():
            if outer_index == inner_index:
                out_matrix.iloc[outer_index, inner_index] = 1.0
            else:
                paired_percent_ovl = percent_overlap(outer_row, inner_row, padding)
                out_matrix.iloc[outer_index, inner_index] = paired_percent_ovl[0]
                out_matrix.iloc[inner_index, outer_index] = paired_percent_ovl[1]

    return out_matrix


if __name__ == "__main__":
    import time

    parser = argparse.ArgumentParser()
    parser.add_argument("--cnv", required=True, help="CNV file to be annotated")
    parser.add_argument("--out", required=True, help="Annotated file")
    parser.add_argument("--pad", required=True, help="Padding extension", type=int)

    args = parser.parse_args()
    cnv_file = args.cnv
    out_file = args.out
    padding_size = args.pad

    start = time.perf_counter()
    matrix = pairwise_match(cnv_file, padding_size)
    writer = pd.ExcelWriter(out_file)
    matrix.to_excel(writer, 'Sheet1', index=True)
    writer.save()
    end = time.perf_counter()
    print("--- Elapsed time: {:.2f} seconds ---".format(end - start))
