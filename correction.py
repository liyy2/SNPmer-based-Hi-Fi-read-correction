
def correction(targetDict_tsv, retList_tsv, ret_dict_txt, fasta_reads):
    for target in targetDict_tsv.keys():
        if len(targetDict_tsv[target] == 1):
            corrected_base = targetDict_tsv[target][0][1]
            naive_correction(target, snpmer2read_dict, fasta_reads)
        elif len(targetDict_tsv[target] == 2):
            correction_with_neighbor(target, retList_tsv, ret_dict_txt, fasta_reads)

def naive_correction(target, ret_dict_txt, fasta_reads):
    corrected_base = targetDict_tsv[target][0][1]
    pos = find_target_location(target, snpmer2read_dict, fasta_reads)
    correct(pos, corrected_base)

    