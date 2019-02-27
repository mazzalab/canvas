import itertools
import sys


class CNVOperations:

    @staticmethod
    def reciprocal_overlap(cnv_pairs_gen, padding=0, combine='combination'):
        # print(cnv_pairs_gen, padding, combine)
        for cnv_pair in cnv_pairs_gen:
            # print(cnv_pair)
            cnv1 = cnv_pair[0]
            cnv2 = cnv_pair[1]
            cnv2_cnv1 = 0
            
            if padding != 0 and padding != '0':
                sys.stdout.write("Extending CNVs with a padding of {}".format(padding))
                cnv1._CNV__start = cnv1._CNV__start - padding
                cnv1._CNV__end = cnv1._CNV__end + padding
                cnv2._CNV__start = cnv2._CNV__start - padding
                cnv2._CNV__end = cnv2._CNV__end + padding

            if cnv2.length == 0:
                sys.stdout.write("CNV of length 1, pass: "+str(cnv2)+'\n')
                cnv2_cnv1 = 0
            elif cnv1.length == 0:
                cnv1_cnv2 = 0
                sys.stdout.write("CNV of length 1, pass: "+str(cnv1)+'\n')
            else:
                cnv1_cnv2 = 100 * cnv1.intersects_with(cnv2) / cnv1.length
                if combine == 'combination':
                    cnv2_cnv1 = 100 * cnv2.intersects_with(cnv1) / cnv2.length
                # print("RESULT", cnv1, cnv2, cnv1_cnv2, cnv2_cnv1)
            yield (cnv1, cnv2, cnv1_cnv2, cnv2_cnv1)

    @staticmethod
    def spanning_overlap(cnv_pairs_gen, min_ovl: float, padding=0):
        for cnv_pair in cnv_pairs_gen:
            cnv1 = cnv_pair[0]
            cnv2 = cnv_pair[1]
            
            if padding != 0:
                sys.stdout.write("Extending CNVs with a padding of {}".format(padding))
                cnv1._CNV__start = cnv1._CNV__start - padding
                cnv1._CNV__end = cnv1._CNV__end + padding
                cnv2._CNV__start = cnv2._CNV__start - padding
                cnv2._CNV__end = cnv2._CNV__end + padding
            
            if cnv1.length == 0:
                cnv1_cnv2 = 0
                cnv1_notin_cnv2 = 0
                sys.stdout.write("CNV of length 1, pass: "+str(cnv1)+'\n')
            elif cnv1 == cnv2:
                cnv1_cnv2 = 0
                cnv1_notin_cnv2 = 0
            else:
                cnv1_cnv2 = cnv1.intersects_with(cnv2)
                cnv1_notin_cnv2 = cnv1.length - cnv1_cnv2
                cnv1_cnv2 = 100*cnv1_cnv2 / cnv1.length

            # if cnv1_cnv2 >= min_ovl and cnv1_notin_cnv2 <= span:
            yield (cnv1, cnv2, cnv1_cnv2, cnv1_notin_cnv2)

    @staticmethod
    def melt(cnv_pairs_gen):
        melted = set()
        tobe_removed = set()
        not_melted = set()
        melt_down = False

        for cnv_pair in cnv_pairs_gen:
            cnv1 = cnv_pair[0]
            cnv2 = cnv_pair[1]

            if cnv1.intersects_with(cnv2) and cnv1.type == cnv2.type:
                melted.add(cnv1.melts_with(cnv2))

                tobe_removed.add(cnv1)
                tobe_removed.add(cnv2)
                melt_down = True
            else:
                not_melted.add(cnv1)
                not_melted.add(cnv2)

        # melted U (not_melted - tobe_removed)
        not_melted.difference_update(tobe_removed)
        melted = melted.union(not_melted)

        if not melt_down or len(melted) == 1:
            return melted
        else:
            return CNVOperations.melt(itertools.combinations(melted, 2))
