'''
Author: Brian G. Peterson (bgpete@umich.edu)
Purpose:
This script determines if a protein coding sequence contains any mutations compared to a reference protein sequence. 
Will count up single and multiple mutations and export data into table format of amino acid position vs amino acids. 
Multiple mutations are also exported as tab-delimited text file of all linked mutations (same read) WTaa###Mutaa
Arguments:
    1) str name to fasta file containing sequences that will be compared
    2) str name to fasta file of reference sequence (expecting only one sequence)
    3) str name for base file name. Will append suffixes for files being written
Output:
    output 5 kinds of files
    single codon changes mutational counts "_Single.tsv"
    multi codon changes mutational counts "_Multi.tsv"
    all codon changes (singles and multiple) mutational count table "_All.tsv"
    multi codon changes as tsv with each column representing amino acid changes "_MultiMutationList.tsv"
    tsv of mutation frequency "_MutationFrequency.tsv"
Assumptions:
    Assume that the fasta file has explicit gaps written out using "-" for length of protein.
    Meant to work downstream of "BAM_Translation_v2.5.py" 
'''

import pandas as pd
import sys
from Bio import SeqIO


def make_protein(length):
    '''
    Takes in length of protein and returns list of dictionaries that will be used to count the frequency of each amino acid on a per position basis.
    Dictionary also tracks what is the WT amino acid "Ref_AA"; and tracks # of reads that are WT, single mutation, multiple mutations, and total read count per position.
    '''
    ## contains the 20 standard amino acids; stop codon "*"; any amino acid (result of poor read quality giving "N") "X"; count gaps "-"" ;total # of reads "total"
    aa_dict= {"A": 0 ,"G": 0 , "I": 0 , "L": 0 , "P": 0 , "V": 0 , "F": 0 ,
        "W": 0 , "Y": 0 , "D": 0 , "E": 0 , "R": 0 ,"H": 0 , "K": 0 ,
        "S": 0 ,"T": 0 , "C": 0 , "M": 0 ,"N": 0 ,"Q": 0 , "*":0, "X":0, "-":0,
        "Ref_AA": '',"WT_read_count":0,"1-mutation_read_count":0,
        "Multi-mutation_read_count":0, "Total_read_count" : 0}
    out=[]
    x = 0
    while x< length:
        out.append(aa_dict.copy()) 
        x+=1
    return out

def call_mutations(in_file, ref_file,out_file,include_WT='Y'):
    '''
    in_file is  str location to fasta containing file with sequences going to determine location of mutations
    ref_file is str location to fasta format reference amino acid sequence (only 1 protein)
    out_file is a base name used for outputting data
        output 5 kinds of files
        single codon changes mutational counts
        multi codon changes mutational counts
        all codon changes (singles and multiple) mutational count table
        multi codon changes as tsv with each column representing amino acid changes
        tsv of mutation frequency
    include_WT == optional parameter if the WT amino acid counts will be outputted; "Y" is Yes, "N" is No
    '''
    ## include WT amino acid counts default is True
    iWT = True
    if include_WT.upper() == 'N':
        iWT = False

    wt_seq = '' ## wt protein sequence
    ## get wt_seq from reference fasta (second fasta)
    with open(ref_file,'r') as file:
        wt_recs = list(SeqIO.parse(file,'fasta'))
        if len(wt_recs) > 1 or len(wt_recs) ==0 :
            print('Multiple Sequences detected in Reference Fasta or No sequences detected \n \
            Exiting program: reference fasta should only contain 1 record')
            exit()
        else:
            wt_seq = wt_recs[0].seq.upper()

    single_mut = make_protein(len(wt_seq)) ## codon table for single mutations
    multi_mut = make_protein(len(wt_seq)) ## codon table for multiple mutations
    all_mut = make_protein(len(wt_seq)) ## codon table for all mutations singles and multiple
    mutation_freq = {} ## will be an dictionary of # of mutations and counts. Key is # (int) of mutations per contig; value is counts (int) (# of reads with those mutation)

    multi_mutation_list_out_file = open(str(out_file+"_MultiMutationList.tsv"),'w')


    x = 0
    while x < len(wt_seq):
        single_mut[x]['Ref_AA'] = wt_seq[x]
        multi_mut[x]['Ref_AA'] = wt_seq[x]
        all_mut[x]['Ref_AA'] = wt_seq[x]
        x+=1

    ## loop through and call mutations from all fasta records
    with open(in_file,'r') as file:
        for record in SeqIO.parse(file,'fasta'):

            r_seq = record.seq.upper()
            ## check to make sure reads match length if they don't print reads that don't match
            if len(r_seq) != len(wt_seq):
                print("Warning: Length of reference and experimental read seq don't match for %s" %(r_seq.id))
                continue


            ## determine # of mutations; loop through fasta sequences and determine # of mutations exclude gaps from analysis ('-')
            mutation_num = 0 ## number of mutations in a read
            mutations = [] ## list of lists that will contain read mutant position [WT codon, mutation position (start codon = 0), mutation]
            wt_loc =[] ## any position that matches wt will be appended with # (location of WT); note start codon = 0
            gap_loc =[] ## any position with gaps

            x=0 ## counter
            while x < len(wt_seq):
                wt_aa = wt_seq[x]
                r_aa = r_seq[x]
                if r_aa == "-":
                    ##ignore these as mutations so don't increase mutation_num, but will output where there are gaps in tables
                    gap_loc.append(x)

                elif wt_aa == r_aa:
                    wt_loc.append(x)

                elif wt_aa != r_aa:
                    mutation_num+=1
                    mutations.append([wt_aa,x,r_aa])

                x+=1
            ## update # of mutations in mutation frequency
            if mutation_num in mutation_freq:
                mutation_freq[mutation_num]+=1
            else:
                mutation_freq[mutation_num] =1


            ## write out mutations to tables for fasta 
            if mutation_num == 0:
                ## WT sequence increase WT reads and total # of reads; two options to either to include WT read counts in output or just in column of WT_read_count

                if iWT:
                    ## include the WT reads in codon count tables
                    for x in wt_loc:

                        single_mut[x]['WT_read_count']+=1
                        single_mut[x]['Total_read_count'] +=1
                        ##use wt_seq to increase WTcount
                        single_mut[x][wt_seq[x]]+=1


                        multi_mut[x]['WT_read_count']+=1
                        multi_mut[x]['Total_read_count'] +=1
                        ##use wt_seq to increase WTcount
                        multi_mut[x][wt_seq[x]]+=1

                        all_mut[x]['WT_read_count']+=1
                        all_mut[x]['Total_read_count'] +=1
                        ##use wt_seq to increase WTcount
                        all_mut[x][wt_seq[x]]+=1

                else:
                    for x in wt_loc:
                        ## x is location of WT amino acid. Won't update the WT amino acid in codon count table, just record # of WT amino acids reads.
                        single_mut[x]['WT_read_count']+=1
                        single_mut[x]['Total_read_count'] +=1

                        multi_mut[x]['WT_read_count']+=1
                        multi_mut[x]['Total_read_count'] +=1

                        all_mut[x]['WT_read_count']+=1
                        all_mut[x]['Total_read_count'] +=1

                ## fill in gap location

                for x in gap_loc:
                    ## don't count as a read count as gaps are lack of coverage
                    single_mut[x]["-"] +=1
                    multi_mut[x]["-"] +=1
                    all_mut[x]["-"] +=1

            elif mutation_num == 1:
                ## single mutation
                ## add up read count first; while it is odd seeing wt_loc; these reads still contain WT sequences at all positions except 1, but count read as single mutation_num
                for x in wt_loc:
                    single_mut[x]['1-mutation_read_count']+=1
                    single_mut[x]['Total_read_count'] +=1

                    multi_mut[x]['1-mutation_read_count']+=1
                    multi_mut[x]['Total_read_count'] +=1

                    all_mut[x]['1-mutation_read_count']+=1
                    all_mut[x]['Total_read_count'] +=1

                for x in mutations:  ## list of list that will contain read mutant position [wt aa, aa loc so x, mutant aa]
                    ## x is list [wt aa, amino acid position (start codon = 0), mutated amino acids]
                    ##x[0] = WT amino acid
                    ##x[1] = amino acid position
                    ##x[2] = mutated amino acid
                    ## there should only be 1 list in mutations given only 1 mutations

                    ## increase count of amino mutations in single and all mutation tables
                    single_mut[x[1]][x[2]]+=1
                    all_mut[x[1]][x[2]]+=1

                    ## increase read counts
                    single_mut[x[1]]['1-mutation_read_count']+=1
                    single_mut[x[1]]['Total_read_count'] +=1

                    multi_mut[x[1]]['1-mutation_read_count']+=1
                    multi_mut[x[1]]['Total_read_count'] +=1

                    all_mut[x[1]]['1-mutation_read_count']+=1
                    all_mut[x[1]]['Total_read_count'] +=1

                ## fill in gaps location
                for x in gap_loc:
                    single_mut[x]["-"] +=1
                    multi_mut[x]["-"] +=1
                    all_mut[x]["-"] +=1

            elif mutation_num >1:
                ## multiple mutations
                for x in wt_loc:
                    single_mut[x]['Multi-mutation_read_count']+=1
                    single_mut[x]['Total_read_count'] +=1

                    multi_mut[x]['Multi-mutation_read_count']+=1
                    multi_mut[x]['Total_read_count'] +=1

                    all_mut[x]['Multi-mutation_read_count']+=1
                    all_mut[x]['Total_read_count'] +=1


                ## add multi codon changes to table and print out the multi mutation to tsv so know which mutations are linked
                line_out = str(record.id+'\t') ##record ID outputs with multi mutations
                for x in mutations:
                    ## x is list [wt aa, amino acid position (start codon = 0), mutated amino acids]
                    ##x[0] = WT amino acid
                    ##x[1] = amino acid position
                    ##x[2] = mutated amino acid

                    ## update the multi mutation table and all mutation table; don't update single table
                    multi_mut[x[1]][x[2]]+=1
                    all_mut[x[1]][x[2]]+=1

                    ## increase read counts
                    single_mut[x[1]]['Multi-mutation_read_count']+=1
                    single_mut[x[1]]['Total_read_count'] +=1

                    multi_mut[x[1]]['Multi-mutation_read_count']+=1
                    multi_mut[x[1]]['Total_read_count'] +=1

                    all_mut[x[1]]['Multi-mutation_read_count']+=1
                    all_mut[x[1]]['Total_read_count'] +=1

                    ##line will be tsv out WTaa###Mutaa ### will be corrected so start codon is #1
                    line_out+= str(x[0]+str(int(x[1]+1))+x[2]+"\t")


                ## write out mutations to out_file; change final tab character to newline
                multi_mutation_list_out_file.write(str(line_out[:-1]+"\n"))

                ## fill in gaps location
                for x in gap_loc:
                    ## don't count as a read count though as we don't know if gap comes from indel (should be excluded)
                    single_mut[x]["-"] +=1
                    multi_mut[x]["-"] +=1
                    all_mut[x]["-"] +=1


    multi_mutation_list_out_file.close()
    ## write out codon tables; 4 files total to write; for the amino acid counts change index so first row = 1 so start Met is actually 1 (not zero) (see the rename axis)
    df_out1 = pd.DataFrame(single_mut).set_index(pd.Index(list(range(1,len(single_mut)+1)))).rename_axis('Amino Acid Position-Start=1',axis=0)
    df_out1.to_csv(str(out_file+"_Single.tsv"),sep="\t")

    df_out_mult = pd.DataFrame(multi_mut).set_index(pd.Index(list(range(1,len(multi_mut)+1)))).rename_axis('Amino Acid Position-Start=1',axis=0)
    df_out_mult.to_csv(str(out_file+"_Multi.tsv"),sep="\t")

    df_out_all= pd.DataFrame(all_mut).set_index(pd.Index(list(range(1,len(all_mut)+1)))).rename_axis('Amino Acid Position-Start=1',axis=0)
    df_out_all.to_csv(str(out_file+"_All.tsv"),sep="\t")

    df_out_freq = pd.DataFrame.from_dict(mutation_freq, orient='index',columns=['# of counts']).rename_axis('# of mutations',axis=0)
    df_out_freq.to_csv(str(out_file+"_MutationFrequency.tsv"),sep= "\t")

def main():
    '''
    Arguments:
    1) str name to fasta file containing sequences that will be compared
    2) str name to fasta file of reference sequence (expecting only one sequence)
    3) str name for base file name. Will append suffixes for files being written
        ## single mutation codon table
        ## multi mutation codon table
        ## multi mutation list
        ## total mutations codon table
        ## mutation frequency
    4) "Y" or "N" Optional 4th argument to include WT codon counts ("Y) in table (default is True if argument is not given)
    '''

    if len(sys.argv) == 4:
        call_mutations(sys.argv[1],sys.argv[2],sys.argv[3])
    elif len(sys.argv) == 5:
        call_mutations(sys.argv[1],sys.argv[2],sys.argv[3],include_WT=sys.argv[4])


if __name__ == '__main__':
    main()