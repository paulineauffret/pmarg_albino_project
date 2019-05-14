# -*- coding: utf-8 -*-

#Import libraries
import sys
import argparse
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO
import difflib 

##############################################################################
#WARNING : run with Python2
#How to use : see "cmd_get_prot_seq_custom.sh"
#2018 Pauline Auffret paulineauffret88@gmail.com
#-----------------------------------------------------------------------------

##############################################################################
parser = argparse.ArgumentParser(description = 'Roughly annotates SNPs and predict if they are synonymous substitutions')
parser.add_argument('snp_file', type = str, help = 'SNP info text file (subset from vcf file, ie 4 columns : CHROM    POS    REF    ALT)')
parser.add_argument('orf_file', type = str, help = 'orf fasta file (longest_orfs.cds obtained from Transdecoder.LongOrfs)')
parser.add_argument('-o', '--output', dest = 'output_file', default = "annotate_snp.out", help = 'output file name')
args = parser.parse_args()
#-----------------------------------------------------------------------------

##############################################################################
#Get SNP infos and store into dictionnary
def store_snp(line, snp_infos):
    #not header nor empty lines
    if not line.startswith("CHROM") and not len(line) == 0:
        line = line.replace("\n", "").split("\t")
        #only takes unique nucleotide variations    
        if len(str(line[2])) == 1 and len(str(line[3])) == 1:
            snp_infos[line[0]] = [line[1], line[2], line[3]]
    return snp_infos

#-----------------------------------------------------------------------------
#Predict if SNPS are synonymous substitutions
def annotate_snp(outfile, snp_infos, orf_infos):
    for elem in orf_infos: 
        syn = "SYN"
        orf_seq = str(elem.seq)
        desc = str(elem.description).split(" ")
        id = str(elem.id.split('.')[0])
        if id in snp_infos:
            start = int(desc[3].split(":")[1].split("-", 1)[0])
            end = int(desc[3].split(":")[1].split("-", 1)[1].split("(")[0])
            strand = desc[3].split(":")[1].split("-", 1)[1].split("(")[1].split(")")[0]
            if end<start:
                tmp = start
                start = end
                end = tmp
            pos = int(snp_infos[id][0])
            ref = snp_infos[id][1]
            alt = snp_infos[id][2]
            if int(pos) < int(start) or int(pos) > int(end):
                syn = "NOT IN ORF"
                orf_seq = "NA"
                prot_orf = "NA"
                seq_mut = "NA"
                prot_seq_mut = "NA"
                aa_ref = "NA"
                aa_mut = "NA" 
            else:
                orf_pos = int(pos)-int(start)+1
                orf_pos_prot = orf_pos/3
                orf_seq = str(orf_seq)
                #REF and ALT are always given on forward strand so need to rev_complement if the ORF is given on the - strand
                if strand == "-" :
                    orf_seq = Seq(orf_seq)
                    orf_seq = orf_seq.reverse_complement()
                orf_pos_nt = orf_seq[orf_pos-1]
                if orf_pos_nt !=  ref:
                    print("Error in determining the ref allele in ORF sequence.")
                seq_mut = list(orf_seq)
                seq_mut[orf_pos-1] = str(alt)
                seq_mut = ''.join(seq_mut)
                orf_seq = Seq(str(orf_seq), generic_dna)
                seq_mut = Seq(str(seq_mut), generic_dna)
                #Back to rev comp with ALT in correct position
                if strand == "-":
                    seq_mut = seq_mut.reverse_complement()
                    orf_seq = orf_seq.reverse_complement()
                prot_orf = orf_seq.translate()
                prot_seq_mut = seq_mut.translate()
                aa_ref = prot_orf[orf_pos_prot-1]
                aa_mut = prot_seq_mut[orf_pos_prot-1]
                if aa_ref !=  aa_mut:
                    syn = "NSYN"
            the_big_list = [str(id), str(start), str(end), str(syn), str(strand), str(orf_seq), str(prot_orf), str(seq_mut), str(prot_seq_mut), str(aa_ref), str(aa_mut)]    
            for desc in snp_infos:
                the_big_list.append(str(desc))
                the_big_list.append(str(snp_infos[desc][0]))
                the_big_list.append(str(snp_infos[desc][1]))
                the_big_list.append(str(snp_infos[desc][2]))
            the_big_list = map(lambda x: x + '\t', the_big_list)
            outfile.writelines(the_big_list)            
            outfile.write("\n")
    return
#-----------------------------------------------------------------------------


##############################################################################
def main():
    try:
        with open(args.snp_file, "rU") as snp_file, open(args.orf_file, "rU") as orf_file, open(args.output_file, "a") as out_file:
            snp_infos = dict()
            for line in snp_file:
                snp_infos = store_snp(line, snp_infos)
            orf_infos = SeqIO.parse(orf_file, "fasta")
            annotate_snp(out_file, snp_infos, orf_infos)
    except IOError, e:
        print "File opening failed: %s" % e.strerror

if __name__ == "__main__":
    main()

##############################################################################
##############################################################################




