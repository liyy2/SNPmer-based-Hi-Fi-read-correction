from utils import printf
from Bio.Seq import Seq
def correction(target_dict_tsv, retList_tsv, ret_dict_txt, fasta_reads):
    for error in target_dict_tsv.keys():
        if len(target_dict_tsv[error] == 1):
            printf('Correcting ' + error + ' with naive strategy')
            target = target_dict_tsv[error][0]
            naive_correction(error, target, snpmer2read_dict, fasta_reads)
        elif len(target_dict_tsv[error] == 2):
            printf('Correcting ' + error + ' with neighbor finding strategy')
            correction_with_neighbor(target, retList_tsv, ret_dict_txt, fasta_reads)


def naive_correction(error_identifier, target, snpmer2read_dict, fasta_reads, k = 181):
    '''
    error is the low_frequency sequence 
'''
    # ============= Grab the read name, and directionality and position of the to be corrected base ============= #
    read_name = snpmer2read_dict[error_identifier][0]['name']
    read_directionality = snpmer2read_dict[error_identifier][0]['strand']
    # corrected position =  current position + floor(k/2)
    position_in_read = snpmer2read_dict[error_identifier][0]['position'] + int(k/2)
    # collect the base of which our snpmer should be corrected to
    target_base = target['CorrectSeq']
    # reverse complement the strand if necessary
    read_seq = fasta_reads[read_name]
    if read_directionality == '-':
        read_seq = fasta_reads[read_name].complement()
    # new read = left + position of the error + right
    new_fasta_reads = read_seq[:position_in_read] + Seq(target_base) + read_seq[position_in_read + 1:]
    # reverse back the sequence if directionality is negative
    if read_directionality == '-':
        new_fasta_reads = new_fasta_reads.complement()
    # assign corrected read to its name
    fasta_reads[read_name] = new_fasta_reads
    
def correction_with_neighbor
    
