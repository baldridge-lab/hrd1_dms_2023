'''
Author: Brian G. Peterson (bgpete@umich.edu)
Purpose: 
The purpose of this script is to translate DNA sequences that have been aligned to a reference DNA fragment from a BAM (binary version of Sequence Alignment Map (SAM)) sorted by name. 
The script is designed to only translate one open reading frame (containing no introns) at a time and is useful for a deep mutational scanning (DMS). 
Assumptions: 
Only 1 open reading frame/ 1 translation product for the reference. There are no  soft-clipping of the aligned sequences. 
Will not translate any deletion/ insertion (indels) even if indel produces an in frame mutation.
Read is only a single contig (merged paired reads; or only single read)
Arguments:
    0) = python file name
    1) str file name input of bam alignment sorted on name ".bam"
    2) str file name output of translation products (fasta format) ".fa"
    3) str file name of reads that fail to translate (short, fall in promoter, term, indels) ".tsv"
    4) int nucleotide location of "A" in "ATG" of 1Met in reference sequence note first nucleotide is 0.
    5) int nucleotide location of third (wobble) basepair in stop codon in reference sequence note first nucleotide is 0.
    6) int minimal length of sequence to translate (default is 150)
Output:
    1) Fasta file ".fa" containing translated sequences with sequence name referring to sequence name
    2) tab-delimited file ".tsv containing any reads that did not translate and reason for not translating
'''

import sys
import pysam
import Bio.Seq as BS
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

def read_translate(seq, ref,seq_name, t_start, t_end,min_len=150):
    '''
    Takes 6 arguments (with 1 optional) and returns "keyword" about where read aligns and then a BioPython SeqRecord object
        The keyword can be 'prom' if read aligns in promoter; 'term' if aligns in terminator ; 'short' if not long enough (less than min_len)
        SeqRecord object will have protein sequence as long as not prom,term,short and SeqRecord id = seq_name
    Arguments:
    sequence= "str" of forward sequence that aligned to reference sequence
    ref= tuple of first and last alignment location of sequence passed (start  alignment loc, last alignmnet loc)
        note if the aligner sends in softcips (local alignment used) this translation function won't work.
    seq_name = "str" sequence name (query_name from fastq header name)
    t_start= "int" location of A in ATG of Met1 used for alignment. Note first nucleotide in reference is zero
    t_end = "int" location of the wobble nucleotide in stop codon. Note first nucleotide in reference sequence is zero
        nnnATGtttTAGnnn : t_start would be 3 ; t_end = 11
    min_len is length of smallest read to translate after trimming leading/ trailing nucleotides. default is 150 nucleotides
    '''
    ## 6 situations that could be occurring
    ## seq in promoter only; return 'prom'
    ## seq in terminator only; return 'term'
    ## seq in promoter into ORF
    ## seq in ORF ends in ORF
    ## seq in ORF ends in terminator
    ## seq in promoter ends in terminator (very long read/ short protein)

    ##will be the nucleotide sequence variable, but trimmed so front/ back of sequence are in triplet coding frame
    trimmed = ''

    ## position of first amino acid once trimmed is translated; start codon (Met) = 0
    aa_start = int()

    ## need +1 as on t_end on wobble position and want a full stop codon
    len_protein = int(t_end+1-t_start)/3 

    if ref[1] <= int(t_start+1):
        ## should get rid of anything in the promoter or partial start codon
        return 'prom', SeqRecord(None,id=seq_name)

    elif ref[0] >= int(t_end -1):
        ## should get rid of anything in the termiantor/ partial stop codon
        return 'term', SeqRecord(None,id=seq_name)

    elif ref[0] <= t_start and ref[1] <= t_end:
    ## promoter into ORF trim up until the start codon
    ## then trim back end so triplet
        trim_front = t_start -ref[0]
        trim_back = int(len(seq)-trim_front)%3

        ## if the seq is inframe %3 returns zero and [:-trim_back] ([:-0])  returns empty string
        if trim_back ==0:
            trimmed = seq[trim_front:]
        else:
            trimmed = seq[trim_front:-trim_back]

        aa_start = int(0)

    elif ref[0] >= t_start and ref[1] <= t_end:
        ## in ORF for start and end nucleotides also if falls on A of ATG and wobble of stop codon
        ## can be  inframe/ off by 1, off by 2 see below

        ##determines # of nucleotides that need to be trimmed
        ##012345678 ;; so t_start would be 0
        ##ATGtcaATC ;; below are nucleotide being
        ##---tcaACT ;; (3-0)%3 = 0 (inframe)
        ##----caACT ;; (4-0)%3 = 1 (trim two nucleotides to get in frame)
        ##-----aACT ;; (5-0)%3 = 2 (trim one nucleotide to get in frame)

        front_offset = int(ref[0]-t_start)%3
        if front_offset==0:
            ##in coding frame
            trim_back=int(len(seq)%3)
            ## need to catch the negative zero trim_back
            if trim_back ==0:
                trimmed=seq
            else:
                trimmed=seq[:-trim_back]
            aa_start = int(int(ref[0]-t_start)/3)

        elif front_offset==1:
            ##offset of 1 so there are two nucleotides to trim from front
            trim_back=int(len(seq)-2)%3

            if trim_back == 0:
                trimmed=seq[2:]
            else:
                trimmed=seq[2:-trim_back]
            aa_start = int(int(ref[0]+2-t_start)/3)
        elif front_offset==2:
            ## trim 1 nucleotide
            trim_back=int(len(seq)-1)%3

            if trim_back ==0:
                trimmed=seq[1:]
            else:
                trimmed=seq[1:-trim_back]
            aa_start = int(int(ref[0]+1-t_start)/3)

    elif ref[0] >= t_start and ref[1] >= t_end:
        ## starts in ORF and goes into terminator
        ## front offset code is same as above
        ## trim the backend based on t_end
        front_offset = int(ref[0]-t_start)%3
        trim_back=int(ref[1]-t_end) ## number of nucleotides aligned in terminator
        if front_offset==0:
            ##in coding frame
            if trim_back ==0:
                trimmed=seq
            else:
                trimmed=seq[:-trim_back]
            aa_start = int(int(ref[0]-t_start)/3)
        elif front_offset==1:
            ## trim 2 nucleotide on front
            if trim_back ==0:
                trimmed=seq
            else:
                trimmed=seq[2:-trim_back]
            aa_start = int(int(ref[0]+2-t_start)/3)
        elif front_offset==2:
            ## trim 1 nucleotide on front
            if trim_back ==0:
                trimmed=seq
            else:
                trimmed=seq[1:-trim_back]
            aa_start = int(int(ref[0]+1-t_start)/3)
    elif ref[0] < t_start and ref[1] > t_end:
        ## read in promoter and ends in terminator

        trim_front = int(t_start -ref[0])
        trim_back = int(ref[1]-t_end)

        if trim_back ==0:
            ## this shouldn't be called as the ref[1] needs to be greater than wobble base of stop codon
            trimmed= seq[trim_front:]
        else:
            trimmed= seq[trim_front:-trim_back]
        aa_start = int(0)




    if len(trimmed) <min_len:
        return 'short', SeqRecord(None,id=seq_name)
    else:

        trans_seq = BS.translate(trimmed)
        ##full protein coding with explicit gaps using "-"
        full_protein = str(aa_start*'-'+
                            trans_seq +
                            int(len_protein-aa_start-len(trans_seq))*'-'
                            )
        ## True indicates read is good (not short, prom, term etc)
        return True, SeqRecord(BS.Seq(full_protein),id=seq_name)


def process_bam(bam_in,out_name,bad_read_out, trans_start,trans_end,min_read_len):
    '''
    function that will loop through bam alignment file using pysam and output translated sequences 
    Arguments:
    bam_in = str input of bam alignment file (needs to be sorted on name)
    out_name = str of fasta file name where translated sequences are going to be outputted
    bad_red_out = str of text file that will contain names of reads that failed to translate (indels, falls in prom/term, too short)
    trans_start= int location of A in ATG of Met1 amino acid of protein used for alignment. Note first nucleotide in reference is zero
    trans_end = int location of the wobble (third) nucleotide in stop codon. note first nucleotide in reference sequence is zero
    min_read_len = int minimal length of read to translate (post trimming leading/ trailing nucleotides to get inframe) default is 150
    Output:
    Fasta file ".fa" with name "out_name" containing translated sequences with sequence name referring to sequence name
    tab-delimited file ".tsv with name "bad_read_out" containing any reads that did not translate and reason for not translating

    '''
    ## print the inputs that were sent in for documentation in standard out
    print(bam_in,out_name,bad_read_out, trans_start,trans_end,min_read_len)
    ##bam_file_operator is object that allows iteration through bam alignment file
    bam_file_operator = pysam.AlignmentFile(bam_in,'rb')

    ##seq_outs list that will contain all the translated sequences in seqrecord format
    first_write_out = True ## for checking if first seq_out write to file
    num_write_out = 100000 ## number of bio seq records in seq_outs to write to file
    seq_outs = []

    ##track # of reads. will print to standard out after processing bam
    total_reads =0 ##total objects in bam file
    total_not_align =0 ## total objects without cigar statement (so not aligned still a read object though)
    total_indel = 0 ## total objects with indel present (based on cigar statement)
    total_softclip = 0 ##can't handle softclips ignore from analysis
    total_indel_softclip =0 ## reads with indel/ softclips
    total_prom_term = 0 ## reads that aligned in promoter/ terminator
    total_short = 0 ## reads that are too short for translation

    file_bad_reads = open(bad_read_out,'w')

    ## loop through object of bam_file "x" will be individual alignment object
    for x in bam_file_operator:
        total_reads +=1 ## increase read count with next object

        ##check to make sure read is good
        ##these are predefined and aren't editable from cmd line would have to edit the script to change what excluding
        ##if read is good proceed to writing to output
        cig = x.cigartuples
        if cig != None:
            ## read has a cigar statement so aligned
            read_indel = False
            read_softclip = False
            ##find if indel/ softclip present; loop through the tuples.
            for y in cig:
                ## see pysam documentation for how cigartuples set up and what #s mean
                ## 1 = insertion; 2 = deletion; 4 = soft clipping
                if y[0] == 1 or y[0] == 2:
                    read_indel = True
                elif y[0] == 4:
                    read_softclip == True

            seq_name = x.query_name
            ## once gone through cig statement ready to translated if passes indel/ softclip filters
            if read_indel or read_softclip:
                if read_indel and read_softclip:
                    total_indel_softclip +=1
                    file_bad_reads.write(str(seq_name+'\t'+'indel_softclip\n'))
                elif read_indel:
                    total_indel +=1
                    file_bad_reads.write(str(seq_name+'\t'+'indel\n'))
                elif read_softclip:
                    total_softclip +=1
                    file_bad_reads.write(str(seq_name+'\t'+'softclip\n'))
            else:
                ##send to read_translate function  as not indel/ softclipped (conditions read_translate can't handle)
                ##in case aligner maps to reverse strand pysam will send out reverse complement to give "original read"
                if x.is_reverse:
                    ## read is reverse (stored as reverse complement in bam, but pysam gives reverse complement back)
                    seq_for = BS.Seq(x.get_forward_sequence()).reverse_complement()
                else:
                    seq_for = x.get_forward_sequence() ##str of nucletodies that aligned 


                ## these aren't flipped if sequence is reverse
                seq_ref = x.get_reference_positions() ## list of positions the nucleotides aligned

                ##translate will return a biopython seqrecord called seq_trans
                seq_good_bad, seq_trans = read_translate(seq_for,(seq_ref[0],seq_ref[-1]),seq_name, trans_start,trans_end,min_len = min_read_len)

                if seq_good_bad== True:
                    seq_outs.append(seq_trans)
                elif seq_good_bad == 'prom':
                    file_bad_reads.write(str(seq_name+'\t'+'prom\n'))
                    total_prom_term+=1
                elif seq_good_bad == 'term':
                    file_bad_reads.write(str(seq_name+'\t'+'term\n'))
                    total_prom_term+=1
                elif seq_good_bad == 'short':
                    file_bad_reads.write(str(seq_name+'\t'+'short_read\n'))
                    total_short+=1
        else:
            ## no cigar statement so not aligned
            total_not_align +=1

        ## To prevent holding all seqrecords in memory write out every 100,000 records)
        if len(seq_outs) == num_write_out:
            if first_write_out:
                ## first time writing out data
                ## open and 'write' not append to outfile in case rerunning analysis in same folder don't want data appended
                with open(out_name,'w') as out_file:
                    SeqIO.write(seq_outs,out_file,'fasta')
                first_write_out = False
            else:
                ## second time (or more) writing to file so append to out_file
                with open(out_name,'a') as out_file:
                    SeqIO.write(seq_outs,out_file,'fasta')
            ## clear seq_outs
            seq_outs=[]

    ## write out data seqrecords to fasta; if less than 100,000 records this will be the first write; if wrote out before will append finals reads to fine.
    if first_write_out:
        with open(out_name,'w') as out_file:
            SeqIO.write(seq_outs,out_file,'fasta')
    else:
        with open(out_name,'a') as out_file:
            SeqIO.write(seq_outs,out_file,'fasta')

    ## close files
    file_bad_reads.close()
    ## print out read stats
    print('bam_in', 'total_reads','total_not_align','total_indel','total_softclip','total_indel_softclip','total_prom_term','total_short')
    print(bam_in, total_reads,total_not_align,total_indel,total_softclip,total_indel_softclip,total_prom_term,total_short)



def main():
    '''
    arguments come from sys.argv 1-5; optional argument 6 is the minimum length that DNA contig needs to be (default is 150)
    0) = python file name
    1) str file name input of bam alignment sorted on name ".bam"
    2) str file name output (fasta format) ".fa"
    3) str file name of reads that fail to translate (short, fall in promoter, term, indels) ".tsv"
    4) int nucleotide location of "A" in ATG of 1Met in reference sequence note first nucleotide is 0
    5) int nucleotide location of wobble base in stop codon in reference sequence note first nucleotide is 0
    6) int minimal length of sequence to translate (default is 150)
    '''
    if len(sys.argv) == 7:
        process_bam(str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[3]),int(sys.argv[4]),int(sys.argv[5]),int(sys.argv[6]))
    elif len(sys.argv) ==6:
        ## didn't specify minimal translation legnth in default of 150
        process_bam(str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[3]),int(sys.argv[4]),int(sys.argv[5]),150)
    else:
        print('Error: Not enough arguments')



if __name__ == '__main__':
	main()