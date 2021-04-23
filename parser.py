from snpmer import SnpmerParser
from collections import defaultdict
from Bio import SeqIO

def parse_dump(filename):
    '''
    Parse the file end with .dump, which stores the snper and their frequencies in all of the reads
    Input: *.dump
    Output: {key: value} of which the key is the sequence of the snpmer, and value is their frequencies
'''
    ret_dictionary = {}
    with open(filename, 'r') as f:
        for line in f:
            output = line.split('\t')
            key = output[0]
            value = int(output[1][:-1])
            ret_dictionary[key] = value
    return ret_dictionary


def parse_tsv(filename, lower_threshold = 1, upper_threshold = 3):
    '''
    Parse the tsv file, which stores the snpmer pairs with one base difference in the middle
    and return a list of to-be-corrected target snpmers (low_fre <= lower_threshold and 
    high_fre >= upper_threshold)
    Input: *.tsv
    Output: One list stores all the tsv snpmer entries in an indexed order.
    One dictionary of which the key is to be corrected sequence, and the 
    value is a dict (index of the snpmer pair, the base pair it should 
    be corrected to in the middle)

'''
    retList = []
    targetDict = defaultdict(list)
    with open(filename, 'r') as f:
        for index, line in enumerate(f):
            output = line.split('\t')
            ref_freq = int(output[0])
            alt_freq = int(output[1])
            ref_seq = output[2]
            alt_seq = output[3][:-1]
            retList.append((ref_seq, ref_freq, alt_seq, alt_freq))
            if ref_freq <= lower_threshold:
                if alt_freq >= upper_threshold:
                    identifier = str(index) + 'R' 
                    target = {'index': index, 'ref_type': 'R', 'CorrectSeq': alt_seq[int(len(alt_seq)/2)]}
                    targetDict[identifier].append(target)
            if alt_freq <= lower_threshold:
                if ref_freq >= upper_threshold:
                    identifier = str(index) + 'A'
                    target = {'index': index, 'ref_type': 'A', 'CorrectSeq': ref_seq[int(len(ref_seq)/2)]}
                    targetDict[identifier].append(target)
    return retList, targetDict







def parse_txt(filename, retList, threshold = 3):
    '''
    Input: txt file, of which the first column is the read name, 
    and the following column conforms to the following format
    {the snpmer id}[RA]:[-]{the position within that read}
    Output: a mapping from snpmer-pair id to a tuple 
    ({read_id, directionality, position})
'''
    read2snpmer_dict = defaultdict(list)
    snpmer2read_dict = defaultdict(list)
    read2signal_dict = defaultdict(list)
    with open(filename, 'r') as f:
        for line in f:
            output = line.split('\t')
            read_name = output[0]
            output = output[1].split(' ')
            output[-1] = output[-1][:-1] 
            for snpmer in output:
                parsed_snpmer = SnpmerParser(snpmer)  
                read2snpmer_dict[read_name].append(parsed_snpmer)
                snp_index = parsed_snpmer.snpmer_id
                freq_1, freq_2 = retList[snp_index][1::2]
                if ((freq_1 > threshold) and (freq_2 > threshold)):
                    read2signal_dict.append(parsed_snpmer)
                snp_id = parsed_snpmer.generate_indentifier()
                snpmer2read_dict[snp_id].append({'name': read_name, 'strand': parsed_snpmer.directionality, 
                                                'position', parsed_snpmer.position_in_read})        
    return read2snpmer_dict, snpmer2read_dict

def parse_fasta(filename):
    dict_of_reads = {}
    for seq in SeqIO.parse(filename, "fasta"):
        dict_of_reads[seq.id] = seq
    return dict_of_reads

            