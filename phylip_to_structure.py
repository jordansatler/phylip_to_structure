"""
Convert phylip SNP file to Structure input file.
SNPs are assumed to be unlinked.

author: J. Satler
date: 15 Dec 2020
version: 1

usage: python phylip_to_structure.py infile.snps

"""
import os
import sys

#dictionary of nucleotide and ambiguity codes for Structure
codes = {"A":["1","1"], "C":["2","2"], "G":["3","3"], "T":["4","4"],
         "Y":["2", "4"], "R":["1", "3"], "W":["1", "4"],
         "S":["2", "3"], "K":["3","4"], "M":["1", "2"]}

def get_snp_data(file):
    """read SNP data file"""
    with open(file, "r") as data:
        return [i.strip() for i in data if not i.split()[1][0].isdigit()]

def convert_data_to_numbers(data):
    """convert SNPs to values"""
    res = {ind.split()[0]:[] for ind in data}
    for line in data:
        name = line.split()[0]
        seq = line.split()[1]
        for bp in seq:
            if bp in codes:
                res[name].append(codes[bp])
            else:
                res[name].append(["-9","-9"])
    return res

def get_missing_data_per_ind(matrix):
    """print amount of missing data per individual"""
    for ind, seq in sorted(matrix.items()):
        bp = [i for i in seq if i[0] != "-9"]
        print "{0}: {1:.02f}% present".format(ind,
                                             (len(bp)/float(len(seq)) * 100))

def remove_ind_missing_data(matrix):
    """remove individuals missing over 90% of data"""
    newMat = {}
    for ind, seq in sorted(matrix.items()):
        bp = [i for i in seq if i[0] != "-9"]
        sam = len(bp)/float(len(seq))
        if sam >= 0.1:
            newMat[ind] = seq
        else:
            print "{0} sampled for {1:.02}% of snps - REMOVED".format(ind,
                                                                      sam * 100)
    return newMat

def write_to_file(fname, matrix, num_ind, num_snp):
    """write data to new file"""
    f = "str_{0}_{1}ind_{2}snp.str".format(fname,
                                           str(num_ind),
                                           str(num_snp))
    with open(f, "w") as out:
        for k, v in sorted(matrix.items()):
            #get individual alleles
            h1 = [a[0] for a in v]
            h2 = [a[1] for a in v]
            st = "{0}_a\t{1}\n{2}_b\t{3}\n".format(k, '\t'.join(h1),
                                                   k, '\t'.join(h2))
            out.write(st)

def main():
    if len(sys.argv) != 2:
        print "python phylip_to_structure.py infile.snps"
        sys.exit()

    d = get_snp_data(sys.argv[1])
    dc = convert_data_to_numbers(d)
    dc_comp = remove_ind_missing_data(dc)
    fname = os.path.basename(sys.argv[1]).split(".")[0]
    write_to_file(fname, dc_comp, len(dc_comp), len(dc[d[0].split()[0]]))

    #print percent of data present for sample
    get_missing_data_per_ind(dc)

if __name__ == '__main__':
    main()
