import argparse
import itertools
import sys
import pandas as pd
from CNVOverlap.cnv_operations import CNVOperations
from CNVOverlap.CNV import CNV
from xlsxwriter.utility import xl_rowcol_to_cell
import xlsxwriter
from pymongo import MongoClient

connection = MongoClient('localhost', 27017, unicode_decode_error_handler='ignore')
db = connection['Incas_overlaps']

import os



def load_db(collection: str):
    cnv_list = []
    for record in db[collection].find():
        # print(record)
        cnv_list.append(CNV(record['chr'], record['start'], record['end'], record['type']))
    # for index, row in data.iterrows():
    #     try:
    #         yield CNV(row['chr'], row['start'], row['end'], row['type'])
    #     except ValueError:
    #         print("CNV {}:{}-{} [{}] skipped".format(row["chr"], row['start'], row['end'], row['type']))
    return cnv_list

def load_cnv(cnv_file: str):
    data = pd.read_table(cnv_file, sep='\t', dtype={'CHR':str})
    data.columns = map(str.upper, data.columns)
    notype = False
    if not (hasattr(data, 'START')):
        raise NameError("The input file does not contain the required 'Start' field.")
    if not (hasattr(data, 'END')):
        raise NameError("The input file does not contain the required 'End' field.")
    if not (hasattr(data, 'CHR')):
        raise NameError("The input file does not contain the required 'Chr' field.")
    if not (hasattr(data, 'TYPE')):
        sys.stdout.write("The 'type' field could not be found in the cnv file. It will be omitted in the output.\n")
        notype = True
    data = data[data.START.notnull() & data.END.notnull() & data.CHR.notnull()]
    data['ID'] = data['CHR'].map(str)+':'+data['START'].map(str)+'-'+data['END'].map(str)
    if notype:
        data = data[['ID', 'CHR', 'START', 'END']]
    else:
        data = data[['ID', 'CHR', 'START', 'END', 'TYPE']]
    data['START'] = data['START'].astype('int')
    data['END'] = data['END'].astype('int')

    seen = set()
    cnv_list = []
    for index, row in data.iterrows():
        row_id = str(row['CHR'])+' '+ str(row['START'])+ ' ' + str(row['END'])
        if row_id not in seen:
            if notype:
                cnv_list.append(CNV(row['CHR'], row['START'], row['END'], ''))
            else:
                cnv_list.append(CNV(row['CHR'], row['START'], row['END'], row['TYPE']))
            seen.add(row_id)
    return cnv_list
    

class OverlapApp:
    
    def __init__(self, args):
        self.args = args
        self.search_db = False
        print(self.args)
        self.d = vars(self.args)

        
    def process(self):
        
        if not os.path.exists(self.args.input_1):
            
            sys.stderr.write("Input file 1 not found. Exiting\n")
            return -1
    
        if not os.path.exists(self.args.input_2):
            if self.args.input_2 in db.collection_names():
                print("Using collection {} as target".format(self.args.input_2))
                self.search_db = True
            else:
                sys.stderr.write("Input file 2 not found. Exiting\n")
                return -1
 
        pathout = os.path.abspath(os.path.join(self.args.output_prefix, os.pardir))
        print("PATHOUT", pathout)
        if not os.path.exists(pathout):
            sys.stderr.write("Trying to write results into a non existing directory ({})\n".format(pathout))
            return -1
        
        sys.stdout.write("Input 1: {0}\nInput 2: {1}\nOutput: {2}\n".format(self.args.input_1, self.args.input_2,
                                                                    pathout))
        cnv_list1 = list(set(load_cnv(self.args.input_1)))
        names_1 = []

        if self.search_db == False:
            cnv_list2 = list(set(load_cnv(self.args.input_2)))
        else:
            cnv_list2 = list(set(load_db(self.args.input_2)))
        names_2 = []
        for c1 in cnv_list1:
            names_1.append(str(c1).split()[0])
        for c2 in cnv_list2:
            names_2.append(str(c2).split()[0])
        
        if self.args.mode == 'reciprocal':

            if self.args.combine_mode == 'combination':
                print("COMBINATION!")

                cnv_gen = itertools.combinations(set(cnv_list1 + cnv_list2), 2)
                names = names_1 = names_2 = list(set(list(names_1 + names_2)))
                out_matrix = pd.DataFrame(columns=names, index=names)
                tot_comparisons = max(len(names) * len(names)-len(names),
                                  len(list(set(cnv_list1 + cnv_list2))) * len(list(set(cnv_list1 + cnv_list2)))-len(list(set(cnv_list1 + cnv_list2))))
            
            elif self.args.combine_mode == 'product':
                print("PRODUCT!")
                cnv_gen = itertools.product(cnv_list1, cnv_list2)
                out_dict = dict.fromkeys(names_1, {})
                tot_comparisons = len(names_1) * len(names_2)
            
            count = 0

            with open(self.args.output_prefix+'_list.csv', 'w') as out_file:
                out_file.write("QUERY\tTARGET\tOVERLAP %\n".format(self.args.min_overlap))
                for res in CNVOperations.reciprocal_overlap(cnv_gen, self.args.padding, self.args.combine_mode):
                    
                    if res[2] >= self.args.min_overlap:
                        out_file.write(str(res[0]).replace(' []','') + "\t" + str(res[1]).replace(' []','') + "\t" + str(round(float(res[2]),2)) + "\n")
                    if res[3] >= self.args.min_overlap:
                        out_file.write(str(res[1]).replace(' []','') + "\t" + str(res[0]).replace(' []','') + "\t" + str(round(float(res[3]),2)) + "\n")
                    if self.args.combine_mode == 'combination':
                        out_matrix.loc[str(res[0]).split()[0], str(res[1]).split()[0]] = round(float(res[2]), 2)
                        out_matrix.loc[str(res[1]).split()[0], str(res[0]).split()[0]] = round(float(res[3]), 2)
                        count += 1
                    elif self.args.combine_mode == 'product':
                        
                        if out_dict[str(res[0]).split()[0]] == {}:
                            out_dict[str(res[0]).split()[0]] = {str(res[1]).split()[0]: round(float(res[2]), 2)}
                        else:
                            out_dict[str(res[0]).split()[0]].update({str(res[1]).split()[0]:round(float(res[2]),2)})
                    count += 1
                    print('{:.2%}%'.format(count/tot_comparisons), end="\r", flush=True)

            print('\n')
            print("percentage", count, tot_comparisons)
            
            if self.args.combine_mode == 'combination':
                writer = pd.ExcelWriter(self.args.output_prefix+'_matrix.xlsx', engine='xlsxwriter')
                out_matrix.to_excel(writer, 'Sheet1', index=True)
                workbook = writer.book
                worksheet = writer.sheets['Sheet1']
                headercells_format = workbook.add_format({'bold': True, 'font_color': 'green', 'align':'center'})
                greycell_format = workbook.add_format({'bg_color': '#dddddd'})
                worksheet.conditional_format('B2:{0}'.format(xl_rowcol_to_cell(len(names_1), len(names_2))), {'type': '2_color_scale', 'min_value': 0, 'max_value': 1, 'min_color': '#FFFFFF', 'max_color': '#00b60c', 'min_type':'num', 'max_type':'num'})
                
                for i in range(1, len(names_1)+1):
                    worksheet.write("{0}".format(xl_rowcol_to_cell(i, i)), '', greycell_format)
                    print("{0}".format(xl_rowcol_to_cell(i, i)))
                worksheet.freeze_panes(1, 1)
    
                # Row of counts
                cell_format = workbook.add_format({'bold': True, 'font_color': 'green', 'align':'right'})
    
                for col_num in range(0, len(names_2)):
                    worksheet.write_formula(len(names_1)+1, col_num+1,
                                            '= COUNTIF(%s:%s,"<>%d")-1' % (xl_rowcol_to_cell(1, col_num+1), xl_rowcol_to_cell(len(names_1), col_num+1), 0), cell_format)
                # Column of counts
                for row_num in range(0, len(names_1)):
                    worksheet.write_formula(row_num+1, len(names_2)+1,
                                            '= COUNTIF(%s:%s,"<>%d")-1' % (xl_rowcol_to_cell(row_num+1, 1), xl_rowcol_to_cell(row_num+1, len(names_2)), 0), cell_format)
    
                worksheet.write('A{}'.format(len(names_2)+2), "Overlap counts", headercells_format)
                worksheet.write('{}'.format(xl_rowcol_to_cell(0, len(names_1)+1)), "Overlap counts", headercells_format)
                worksheet.set_column(0, len(names_2)+1, 30)
                writer.save()
                
            elif self.args.combine_mode == 'product':
                # print(out_dict)
                col = 0
                workbook = xlsxwriter.Workbook(self.args.output_prefix+'_matrix.xlsx')
                worksheet = workbook.add_worksheet()
                headercells_format = workbook.add_format(
                    {'bold': True, 'align': 'center'})
                worksheet.conditional_format('B2:{0}'.format(xl_rowcol_to_cell(len(names_2), len(names_1))),
                                             {'type': '2_color_scale', 'min_value': 0, 'max_value': 1,
                                              'min_color': '#FFFFFF', 'max_color': '#00b60c',
                                              'min_type': 'num', 'max_type': 'num'})
                
                for key in sorted(out_dict.keys()):
                    col += 1
                    worksheet.write(0, col, key, headercells_format)
                    row = 0
                    for item in sorted(out_dict[key]):
                        worksheet.write(row+1, 0, item, headercells_format)
                        worksheet.write(row + 1, col, out_dict[key][item])
                        row += 1

                worksheet.freeze_panes(1, 1)
                
                # Row of counts
                overlap_cell_format = workbook.add_format({'bold': True, 'font_color': 'green', 'align': 'right'})
                for col_num in range(0, len(names_1)):
                    worksheet.write_formula(len(names_2) + 1, col_num + 1,
                                            '= COUNTIF(%s:%s,"<>%d")' % (xl_rowcol_to_cell(1, col_num + 1),
                                                                           xl_rowcol_to_cell(len(names_2),
                                                                                             col_num + 1), 0),
                                            overlap_cell_format)
                # Column of counts
                for row_num in range(0, len(names_2)):
                    worksheet.write_formula(row_num + 1, len(names_1) + 1,
                                            '= COUNTIF(%s:%s,"<>%d")' % (xl_rowcol_to_cell(row_num + 1, 1),
                                                                           xl_rowcol_to_cell(row_num + 1,
                                                                                             len(names_1)),
                                                                           0), overlap_cell_format)

                worksheet.write('A{}'.format(len(names_2)+2), "Overlap counts", overlap_cell_format)
                worksheet.write('{}'.format(xl_rowcol_to_cell(0, len(names_1) + 1)), "Overlap counts",
                                overlap_cell_format)

                worksheet.set_column(0, len(names_1)+1, 30)
                
                workbook.close()

            
        elif self.args.mode == 'spanning':
            names_1 = list(set(names_1))
            names_2 = list(set(names_2))
            cnv_gen = itertools.product(set(cnv_list1), set(cnv_list2))

            tot_comparisons = len(names_1) * len(names_2)
            count = 0
            with open(self.args.output_prefix+'_list.csv', 'w') as out_file:
                out_file.write("QUERY\tTARGET\tOVERLAP %\n".format(self.args.min_overlap, self.args.padding))
                for res in CNVOperations.spanning_overlap(cnv_gen, self.args.padding):
                    if res[2] >= self.args.min_overlap and res[3] <= self.args.span:
                        out_file.write(str(res[0]).replace(' []','') + "\t" + str(res[1]).replace(' []','') + "\t" + str(round(float(res[2]), 2)) + "\n")
                    count += 1
                    print('{:.2%}%'.format(count/tot_comparisons), end="\r", flush=True)
            print('\n')
        
        return 0

        
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", help="Overlapping mode", choices=['reciprocal', 'spanning'])
    parser.add_argument("--combine-mode",
                        help="The itertools algorithm to be used to generate comparisons. Combination is "
                             "useful for self comparisons and to compare two small datasets reciprocally; "
                             "Product is indicated for comparisons against large DBs.",
                        choices=['combination', 'product'], default='combination')
    parser.add_argument("--min-overlap", type=float, default=50, help="Min. percentage of overlap in order "
                                                                      "to consider two CNVs as overlapping")
    parser.add_argument("--padding", type=int, default=0, help="Tolerance for overlap")
    parser.add_argument("--span", type=int, default=100000, help="Maximum length that does not overlap in order "
                                                            "to consider two CNVs as overlapping")
    parser.add_argument("--input_1", help="First CNV file")
    parser.add_argument("--input_2", help="Second CNV file")
    parser.add_argument("--output_prefix", help="Output prefix for the resulting files.")
    args = parser.parse_args()
    
    
    success = OverlapApp(args).process()
    if success == -1:
        sys.stdout.write("There was a problem in annotating CNVs.\n")
        
    print(success)
    # cnv1 = CNV("chr1", 10, 20, "loss")
    # cnv2 = CNV("chr1", 8, 15, "loss")
    # cnv3 = CNV("chr1", 1, 30, "gain")
    # cnv_list = [cnv1, cnv2, cnv3]
    # min_ovl = 0.1
    # print("IMPORT CNV1")
    # cnv_list = load_cnv("db/lista CNV.csv")
    # print("CNV 1 IMPORTED")
    # dgv1 = CNV("chr1", 3, 7, "loss")
    # dgv2 = CNV("chr2", 100, 250, "gain")
    # dgv3 = CNV("chr1", 5, 15, "loss")
    # dgv4 = CNV("chr4", 1000000, 2000000, "loss")
    # dgv5 = CNV("chr4", 75, 410, "loss")
    # dgv6 = CNV("chr4", 403, 430, "gain")
    # dgv_list = [dgv1, dgv2, dgv3, dgv4, dgv5, dgv6]
    # print("IMPORT CNV2")
    # dgv_list = load_cnv("db/DGV_GRCh37_hg19_variants_2016-05-15.txt") # load_dgv("F:\Dropbox\Applicazioni\CNVOverlap\db\DGV_GRCh37_hg19_variants_2016-05-15.txt")
    # print("CNV2 IMPORTED")

    # # Let's intersect all combination of CNV in the same set
    # cnv_gen = itertools.combinations(cnv_list, 2)
    # print("Reciprocal overlap within CNV: {}".format(list(C
    # NVOperations.reciprocal_overlap(cnv_gen, 0.5))))
    #
    # # Let's intersect all combination of CNV in the same set, considering a spanning range of 100k bp
    # cnv_gen = itertools.combinations(cnv_list, 2)
    # print("Spanning overlap within CNV (span = 100k): {}".format(list(CNVOperations.spanning_overlap(cnv_gen, 0.5, 100000))))

    # Let's intersect any element of one set with any elements of another set
    # cnv_dgv_gen = itertools.product(cnv_list, dgv_list)
    # print(cnv_list)
    # print(dgv_list)
    # # print("Reciprocal overlap within CNV x DGV: {}".format(list(CNVOperations.reciprocal_overlap(cnv_dgv_gen, 0.5))))
    # with open("db/cnv_dgv_match.txt", 'w') as the_file:
    #     the_file.write("CNV\tDGV\tPERC\tMIN_OVERLAP_{}%\n".format(min_ovl*100))
    #
    #     for res in CNVOperations.reciprocal_overlap(cnv_dgv_gen):
    #         if res[2] >= min_ovl:
    #             the_file.write(str(res[0]) + "\t" + str(res[1]) + "\t" + str(res[2]) + "\t" + str(res[2] >= min_ovl) + "\n")
    #         if res[3] >= min_ovl:
    #             the_file.write(str(res[1]) + "\t" + str(res[0]) + "\t" + str(res[3]) + "\t" + str(res[3] >= min_ovl) + "\n")

    sys.exit()

    # # Let's melt the CNV list
    # cnv_gen = itertools.combinations(cnv_list, 2)
    # print("Melting CNV list: {}".format(list(CNVOperations.melt(cnv_gen))))
    #
    # # Let's melt the DGV list
    # dgv_gen = itertools.combinations(dgv_list, 2)
    # print("Melting DGV list: {}".format(list(CNVOperations.melt(dgv_gen))))
