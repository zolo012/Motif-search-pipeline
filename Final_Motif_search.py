import requests
import copy
import math
import re
import os
from iupred3_lib import * # for iupred scoring
import subprocess
import pandas as pd
from bs4 import BeautifulSoup


########################################################################################################################
# 						Visualisation
########################################################################################################################

# ----------------------------------------------------------------------------------------------------------------------

# display matrix but not + signs before numbers -> looks matrix badly if there are any signs
def nice_mtx_show(mtx):
    for i, j in mtx.items():
        print(i, "   ".join("{:.4f}".format(value) for value in j))
    return -1

# ----------------------------------------------------------------------------------------------------------------------

# display such matrices that have + and - numbers too
def nicer_mtx_show(mtx):
    for i, j in mtx.items():
        print(i, "    ".join("{:+.4f}".format(value) for value in j))
    return 1

# ----------------------------------------------------------------------------------------------------------------------

# Convert matrix into tsv to display in excel
def mtx_to_tsv(mtx, file_path):
    file = open(file_path, 'w')
    cont_list = []
    cont_string = ''
    # loop through rows
    for aa in mtx:
        cont_string += aa + '\t'
        # loop through columms
        for cnt in mtx[aa]:
            cont_string += str(cnt) + '\t'
        cont_string += '\n'
    # append string and write it into outputfile
    cont_list.append(cont_string)
    file.writelines(cont_list)
    print('Convertation was successful!')
    return 1

# ----------------------------------------------------------------------------------------------------------------------

########################################################################################################################
# 						FASTA file processing						
########################################################################################################################

# read FASTA file from computer and extract sequences and headers in 2 different lists in dictionary
def multi_FASTA_seq_extractor(fasta_path):
    # open fasta file and close it
    with open(fasta_path, "r") as f:
        # make dictionary which consist of headers and sequences
        fasta_dict = {"header" : [],
                      "sequence" : []}
        # read lines
        file_cont = f.readlines()
    # set counter to 0 and an empty string to store actual sequence later
    i = 0
    seq = ''
    # go across all lines
    while i <= len(file_cont)-1:
        # if line is header
        if file_cont[i].startswith('>'):
            # save header
            fasta_dict['header'].append(file_cont[i])
            # if not first line then there are already stored entire sequence from previous consecutive cycles, thus must be saved into dictionary
            if i != 0:
                fasta_dict['sequence'].append(seq)
                seq = ''
        #if line is not header then add sequence part to the seq string for complete sequence
        else:
            seq += file_cont[i].strip()
        i += 1
    # last sequence must be added in this way
    fasta_dict['sequence'].append(seq)
    return fasta_dict

# ----------------------------------------------------------------------------------------------------------------------

# Extracts all information about proteins in multi fasta file and place them in dictionary within keys = header tags: {'ID' : id, 'OS' : os, ....}
def impr_multi_FASTA_seq_extractor(fasta_path):
    # create variables
    fasta_dict = {}  # store the proteins' information: {protein ID : {name: ..., OS: ..., OX: ..., GN: ..., PE: ..., SV: ...}}
    ID_list = [] # store ID of proteins to be able to place them as keys into dictionary

    # Open multifasta file and close it after using
    with open(fasta_path, "r") as f:
        # read lines
        f_lines = f.readlines()
        for line in f_lines:    # loop over all lines
            if line.startswith(">"):  # if it is header (starts with ">")
                # extract ID
                ID_start_index = line.index("|") + 1
                ID_end_index = line.index("|", ID_start_index)
                ID = line[ID_start_index:ID_end_index]
                ID_list.append(ID)
                fasta_dict.update({ID : {}})

                # extract name
                name_start_index = ID_end_index + 1
                name_end_index = line.index("OS=", name_start_index)
                name = line[name_start_index:name_end_index]
                fasta_dict[ID].update({"Name" : name})


                # extract OS
                OS_start_index = line.find("OS=")
                # if OS exists then extract it
                if OS_start_index != -1:
                    OS_start_index += 3
                    OS_center_index = line.index(" ", OS_start_index)
                    OS_end_index = line.index(" ", OS_center_index + 1)
                    OS = line[OS_start_index:OS_center_index] + line[OS_center_index : OS_end_index]
                # if OS not exist then assign 'Unknown'
                else:
                    OS = "Unknown"
                # store OS value
                fasta_dict[ID].update({"OS" : OS})

                # extract OX
                OX_start_index = line.find("OX=")
                # if OX exist then extract it
                if OX_start_index != -1:
                    OX_start_index += 3
                    OX_end_index = line.index(" ", OX_start_index)
                    OX = int(line[OX_start_index:OX_end_index])
                # if OX not exist then assign 'Unknown'
                else:
                    OX = "Unknown"
                # store OX value
                fasta_dict[ID].update({"OX" : OX})

                # extract GN
                GN_start_index = line.find("GN=")
                # if GN exist then extract it
                if GN_start_index != -1:
                    GN_start_index += 3
                    GN_end_index = line.index(" ", GN_start_index)
                    GN = line[GN_start_index:GN_end_index]
                # if GN not exist then assign 'Unknown'
                else:
                    GN = "Unknown"
                # store GN value
                fasta_dict[ID].update({"GN" : GN})

                # extract PE
                PE_start_index = line.find("PE=")
                # if PE exists then extract it
                if PE_start_index != -1:
                    PE_start_index += 3
                    PE_end_index = line.index(" ", PE_start_index)
                    PE = int(line[PE_start_index:PE_end_index])
                # if PE not exist then assign 'Unknown'
                else:
                    PE = "Unknown"
                # store PE value
                fasta_dict[ID].update({"PE" : PE})

                # extract SV
                SV_start_index = line.find("SV=")
                # if SV exists then extract it
                if SV_start_index != -1:
                    SV_start_index += 3
                    SV = int(line[SV_start_index:].strip())
                # if SV not exist then assign 'Unknown'
                else:
                    SV = "Unknown"
                # store SV value
                fasta_dict[ID].update({"SV" : SV})
                seq = ""

            else:   # if it is not header then add sequence lines until next header and store it into nested dictionary with 'sequence' key
                seq += line.strip()
                fasta_dict[ID].update({"sequence" : seq})
    return fasta_dict

# ----------------------------------------------------------------------------------------------------------------------

# Create FASTA file for selected proteins separately by using a multi-fasta file and a sequence library created by impr_multi_FASTA_seq_extractor (see later)
# inputs: protein ids to store in fasta file, search file which contains the sequences of proteins (multifasta),file_path is the path of output directory,
# sequence dictionary created by impr_multi_FASTA_seq_extractor
def fasta_creator(id, search_file, file_path, seq_dict):
    # open file in which there are the required proteins' headers
    with open(search_file, 'r') as f:
        # loop through lines of search_file
        for line in f:
            # if given line contains one of the wanted proteins' ids then save the header as a string (fasta_string)
            if id in line:
                fasta_string = '' + line
                # search its sequence in sequence library
                seq = seq_dict[id]['sequence']
                # add corresponding sequence to the fasta_string
                fasta_string += seq
                #print(fasta_string)
                # create output filename then its path
                file_name = '/' + id + '.fasta'
                file_path += file_name
                # open output file and write the fasta_string into it
                with open(file_path, 'w') as nf:
                    nf.write(fasta_string)
    print('FASTA file was successfully created.')
    return 1

# ----------------------------------------------------------------------------------------------------------------------

# Concatenate coherent headers and sequences from multi-fasta file and return with a list
def head_and_neck_merge(fasta_path):
    new_fasta_cont_list = []
    # open and read multifasta file
    fasta_file = open(fasta_path, 'r')
    fasta_cont_list = fasta_file.readlines()
    # loop through line of fasta file
    for id in range(len(fasta_cont_list)):
        # if list index is odd then it's sequence and go to the next line
        if id % 2 != 0:
            continue
        # if list index is even then it's header
        else:
            # join the current (header) and previous (sequence) lines and add to the list
            merged_head_and_neck = str(fasta_cont_list[id]) + str(fasta_cont_list[id + 1])
            new_fasta_cont_list.append(merged_head_and_neck)
    fasta_file.close()
    return new_fasta_cont_list

# ----------------------------------------------------------------------------------------------------------------------

########################################################################################################################
# 						Motif search with PSSM						
########################################################################################################################

# ----------------------------------------------------------------------------------------------------------------------

# Create count matrix from given same length motifs
def count_mtx_maker(instance_motifs_list):
    # initialize count matrix with pseudo-count (1 in this case)
    count_mtx = {aa : [1]*len(instance_motifs_list[0]) for aa in ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']}
    #print(motifs_list)
    # loop through all motifs in motif list
    for motif in instance_motifs_list:
        # loop through each residue in given motif and increase the corresponding matrix entry by one
        for idx, aa in enumerate(motif):
            count_mtx[aa][idx] += 1
    return count_mtx

# ----------------------------------------------------------------------------------------------------------------------

# Convert count matrix into probability matrix
def count_to_prob_mtx_converter(cnt_mtx, instances_number, amino_acid = True, normalized = True, norm_dict = False):
    unnorm_prob_mtx = copy.deepcopy(cnt_mtx) # Try with either .deepcopy() or list() constructor
    norm_prob_mtx = copy.deepcopy(cnt_mtx)
    # if passed instances are amino-acids
    if amino_acid == True:
        # if want normalized version
        if normalized == True:
            # normalized version
            for pos in range(len(cnt_mtx["A"])):  # go over column (referred to A because it's in both DNA and protein/polypeptide as well)
                for aa in cnt_mtx: # go over rows
                    norm_prob = cnt_mtx[aa][pos] / (20 + instances_number) / norm_dict[aa][0] #calculate probability of each aa.: number of current aa. / all instances number / frequency of actual aa. in eukaryotes
                    #print("norm_prob", norm_prob)
                    norm_prob_mtx[aa][pos] = norm_prob # change the original count value at current position to probability value
            return norm_prob_mtx

        # if want not normalized version
        elif normalized == False:
            for pos in range(len(cnt_mtx["A"])): # go over column (referred to A because it's in both DNA and protein/polypeptide as well)
                for aa in cnt_mtx: # go over rows
                    unnorm_prob = cnt_mtx[aa][pos] / (20 + instances_number) # calculate probability of each aa.: number of current aa. / all instances number
                    unnorm_prob_mtx[aa][pos] = unnorm_prob
            return unnorm_prob_mtx
        # not proper normalized value (neither 'True' nor 'False')
        else:
            print("Error. Invalid normalized argument!")
            return -1

    # if passed instances are DNA
    elif amino_acid == False:
        # if want normalized version
        if normalized == True:
            for pos in range(len(cnt_mtx["A"])): # go over column (referred to A because it's in both DNA and protein/polypeptide as well)
                for aa in cnt_mtx:
                    norm_prob = cnt_mtx[aa][pos] / (4 + instances_number) / norm_dict[aa][0] #calculate probability of each nucleotides.: number of current nucl. / all instances number / normalization values
                    norm_prob_mtx[aa][pos] = norm_prob
            return norm_prob_mtx

        # if want not normalized version
        elif normalized == False:
            for pos in range(len(cnt_mtx["A"])): # go over column (referred to A because it's in both DNA and protein/polypeptide as well)
                for aa in cnt_mtx:
                    unnorm_prob = cnt_mtx[aa][pos] / (4 + instances_number) # calculate probability of each nucleotides.: number of current nucl. / all instances number
                    unnorm_prob_mtx[aa][pos] = unnorm_prob
            return unnorm_prob_mtx
        # not proper normalized value (neither 'True' nor 'False')
        else:
            print("Error. Invalid normalized argument!")
            return -1
    # not proper amino_acid value (neither 'True' nor 'False')
    else:
        print("Error. Invalid amino_acid argument!")
        return -1

# ----------------------------------------------------------------------------------------------------------------------

# It converts the values of any kind of matrix into logarithm form
def mtx_log_converter(mtx):
    # go over all keys (aa) in got matrix and take the logarithm of values and put results in a list which will replace the values of matrix
    for aa in mtx.keys():
        log_values_list = [math.log10(value) for value in mtx[aa]]
        mtx[aa] = log_values_list
    return mtx

# ----------------------------------------------------------------------------------------------------------------------

# calculate the score of passed kmer based on (PSSM) matrix
def score(kmer, mtx):
    return sum([mtx[aa][pos] for pos, aa in enumerate(kmer)])

# ----------------------------------------------------------------------------------------------------------------------

# Count the number of motifs in proteome and stores information of good motifs in a dictionary with their scores and positions in separated lists:
# {'Motifs' : [motif1, motif2...], 'Motifs_scores' : [score1, score2...], 'Motifs_position' : [position1, position2...]}
def solid_PSSM_motif_searcher(proteome, kmer_len, mtx):
    # create variables for results storing
    good_motifs_list, good_motifs_score_list, good_motifs_pos_list, = [], [], []
    pos = 0

    while pos <= len(proteome) - kmer_len: # loop over kmers in proteome
        test_kmer = proteome[pos : pos + kmer_len]
        if test_kmer.find("U") != -1 or test_kmer.find("X") != -1:  # skip the U (selenocytein) and X (any aminoacid)
            pos += 1 # increase position by one
            continue
        else:
            kmer_score = score(test_kmer, mtx) # calculate score of kmer by using own score function
            #print('Test kmer and its score:', test_kmer, kmer_score)
            # keep only positive scored kmers
            if kmer_score > 0:
                #print(test_kmer, kmer_score)
                good_motifs_score_list.append(kmer_score) # stores scores of accepted motifs as list
                good_motifs_list.append(test_kmer) # stores accepted motifs as list
                good_motifs_pos_list.append(str(pos+1) + ":" + str(pos + kmer_len + 1)) # Must add 1 for both start and end index to make interpreation easier for humans, since Python uses 0-indexing
                                                                                    # but for indexing later in python subtract 1 from the results in dictionary!
        pos += 1 # increase position by one

    # Make dictionary of obtained data
    good_motifs_dict = {"Motifs" : good_motifs_list,
                        "Motifs_scores" : good_motifs_score_list,
                        "Motifs_position" : good_motifs_pos_list}
    # Display the number of passed kmers -> possible motifs
    print(len(good_motifs_score_list))
    return good_motifs_dict

# ----------------------------------------------------------------------------------------------------------------------

# Improved version of solid_PSSM_motif_searcher that stores information in nested dictionary: {id : {'motifs' : {app_numb : {'position' : position, 'score' : score}}}}
# inputs: proteome dictionary that contains the protein sequences (id : sequence), motif length, matrix to score, threshold to distinguish good and bad hits, 'All' = True means all kmers even bad ones will be store, equal_too = False mean just motif with greater score than threshold will be accepted, keep non = False means delete empty proteins from dictionary so when no any motif of it
def multifasta_PSSM_motif_searcher(proteome_dict, motif_len, mtx, threshold = 0, All = True, equal_too = False, keep_none = True):
    # create variables
    motif_dict = {ID : {"motifs" : {}} for ID in proteome_dict}  # store the result
    found_motifs, all_kmers = 0, 0   # for counting all possible kmers and the found motifs
    found_motifs_set = set()
    protein_id_set = set()
    # loop through proteins
    for ID in proteome_dict.keys():
        #print('ID:', ID)
        # extract the sequence from given dict
        current_seq = proteome_dict[ID]['sequence']
        #print('current_seq:', current_seq)
        # go through the kmers in current sequence
        for pos in range(len(current_seq) - motif_len + 1):
            start_index = pos  # store the start index of actual kmer         |
            #print('start_index:', start_index)                             # | for store the kmer position in protein in the result dict at the end of function
            end_index = pos + motif_len  # store the end index of actual kmer |
            #print('end_index:', end_index)
            poss_motif = current_seq[start_index : end_index]    # generate kmer
            #print('poss_motif', poss_motif)
            # check if current kmer contains U (selenocysteine) or X (any amino acid)
            if poss_motif.find("U") != -1 or poss_motif.find("X") != -1:
                # if true then skip that kmer and continue with the next kmer
                continue
            else:
                # increase the number of all kmers by one in each iteration if it doesn't have U and X
                all_kmers += 1
                # in false calculate the score of kmer based on got matrix with score function
                motif_score = score(poss_motif, mtx)
                #print('motif_score', motif_score)
                start_index += 1
                motif_pos = str(start_index) + "-" + str(end_index)    # store motif position
                if equal_too == False:
                    if motif_score > threshold:
                        # increase number of found motifs by one in each iteration
                        found_motifs += 1
                        found_motifs_set.add(poss_motif)
                        motif_pos = str(start_index) + "-" + str(end_index)
                        protein_id_set.add(ID)

                        # check the actual motif already present among motifs in result dict (needed for proper same motif occurrence number generation)
                        if bool(poss_motif in motif_dict[ID]["motifs"].keys()) == True:
                            # if true increase the occurrence of given motif by one then store gained information
                            #print('True, poss_motif', ID, poss_motif)
                            cntr = len(list(motif_dict[ID]["motifs"][poss_motif])) + 1
                            #print('cntr:', cntr)
                            motif_dict[ID]["motifs"][poss_motif].update({str(cntr): {"position": motif_pos,
                                                                                     "score": motif_score}})
                        else:
                            # if false generate one occurrence to that motif and store gained information
                            #print("Current motif doesn't present in dictionary yet.")
                            cntr = 1
                            motif_pos = str(start_index) + "-" + str(end_index)
                            motif_dict[ID]["motifs"].update({poss_motif : {str(cntr) : {"position" : motif_pos,
                                                                                        "score" : motif_score}}})
                elif equal_too == True:
                    if motif_score >= threshold:
                        # increase number of found motifs by one in each iteration
                        found_motifs += 1
                        found_motifs_set.add(poss_motif)
                        motif_pos = str(start_index) + "-" + str(end_index)
                        protein_id_set.add(ID)

                        # check the actual motif already present among motifs in result dict (needed for proper same motif occurrence number generation)
                        if bool(poss_motif in motif_dict[ID]["motifs"].keys()) == True:
                            # if true increase the occurrence of given motif by one then store gained information
                            #print('True, poss_motif', ID, poss_motif)
                            cntr = len(list(motif_dict[ID]["motifs"][poss_motif])) + 1
                            #print('cntr:', cntr)
                            motif_dict[ID]["motifs"][poss_motif].update({str(cntr): {"position": motif_pos,
                                                                                     "score": motif_score}})
                        else:
                            # if false generate one occurrence to that motif and store gained information
                            #print("Current motif doesn't present in dictionary yet.")
                            cntr = 1
                            motif_pos = str(start_index) + "-" + str(end_index)
                            motif_dict[ID]["motifs"].update({poss_motif : {str(cntr) : {"position" : motif_pos,
                                                                                        "score" : motif_score}}})

                    #print('motif_dict:', motif_dict)

    # put not found value inside such proteins keys which don't have any accepted kmers as motifs
        if bool(motif_dict[ID]["motifs"]) == False:
            if keep_none == True:
                #print('Not found motif in actual protein.')
                motif_dict[ID].update({"motifs" : "Not found"})
            elif keep_none == False:
                del motif_dict[ID]

    # display the number of all kmers and found motifs
    print("\nAll kmers that don't have X or U: {:,}".format(all_kmers), "\nFound motifs with repetition (all good motif even more in same protein): {:,}".format(found_motifs), "\nFound motifs without repetition (number of kind of good motifs): {:,}".format(len(found_motifs_set)), '\nNumber of proteins that consist at least one good motif: {:,}'.format(len(protein_id_set))) # fix: found motifs shows the number of motifs with repetition
    return motif_dict

# ----------------------------------------------------------------------------------------------------------------------

# Converts given dictionary from multifasta_PSSM_motif_searcher into such dictionary within likely to check quickly which proteins a specific motif presents in? {poss_motif : [protein_id1, protein_id2, protein_id3...]...}
def ID_to_motif_mk_dict_converter(motif_dict):  # mk = main key
    new_motif_dict = {}
    for ID in motif_dict:  # loop over the protein IDs
        #print(ID)
        for motif in motif_dict[ID]['motifs']: # loop over the motifs that are present in current protein
            if motif not in new_motif_dict:
                new_motif_dict[motif] = []
            new_motif_dict[motif].append(ID)
    #print(motif_dict)
    return new_motif_dict

# ----------------------------------------------------------------------------------------------------------------------

# extract and places score of motif candidates from dictionary created by multifasta_PSSM_motif_searcher into a dictionary: {ID1 : [scores], ID2 : [scores]}
# + places a keys ('All_scores') with scores of all motif candidates
def dict_score_extractor(dict):
    # create dictionary for score storing in later: {ID1 : [scores], ID2 : [scores]}
    score_dict = {ID : [] for ID in dict}

    # extract and put all scores in list corresponding to ID
    for ID in dict: # loop over all proteins
        for motif in dict[ID]['motifs']:  # loop over all motifs of actual protein
            for app_numb in dict[ID]['motifs'][motif]:  # go through all occurrences of actual motif of actual protein
                score = dict[ID]['motifs'][motif][app_numb]['score']
                score_dict[ID].append(score)

    # put sum of all values regardless of corresponding IDs and motifs as a value into a new key
    All_scores = []
    for ID in score_dict:
        All_scores.extend(score_dict[ID])
    score_dict['All_scores'] = All_scores

    return score_dict

# ----------------------------------------------------------------------------------------------------------------------

# Score all kmers in multi fasta file by using PSSM and store only the scores
# inputs: sequence library, motif length, matrix to score, threshold to distinguish good and bad motifs, 'All = True' means scores of all motifs regardless of their scores are stored
def PSSM_motif_score_and_store(proteome_dict, motif_len, mtx, threshold = 0, All = True):
    scores = []
    # loop through all sequences in sequence library
    for ID in proteome_dict:
        current_seq = proteome_dict[ID]['sequence']
        # loop through all possible kmers
        for pos in range(len(current_seq) - motif_len + 1):
            poss_motif = current_seq[pos : pos + motif_len]
            # if given possible motif contains unknown or special character then skip it
            if poss_motif.find('U') != -1 or poss_motif.find('X') != -1:
                continue
            # calculate score of possible motif based on matrix
            poss_motif_score = score(poss_motif, mtx)
            # if want keep scores of all motif candidates then save it
            if All == True:
                scores.append(poss_motif_score)
            # if not want keep scores of all motif candidates then just motif candidates with greater score then certain threshold value are accepted and stored
            else:
                if poss_motif_score > threshold:
                    scores.append(poss_motif_score)
    return scores

# ----------------------------------------------------------------------------------------------------------------------

# calculate the score of passed kmer based on (PSSM) matrix, but with gaps; inputs: kmer, matrix, gap score
def score_wth_gaps(kmer, mtx, gap):
    score_list = []
    print(kmer)
    # go through actual kmer
    for pos, aa in enumerate(kmer):
        # if given residue is unknown or gap then give gap score
        if aa == '-' or aa == 'X':
            score_list.append(gap)
        # else give the matrix score to it
        else:
            score_list.append(mtx[aa][pos])
    return sum(score_list) # return with sum of values of score list

# ----------------------------------------------------------------------------------------------------------------------

########################################################################################################################
# 						Motif search with Regular expression						
########################################################################################################################

# ----------------------------------------------------------------------------------------------------------------------

# Search motifs by using regex in proteins then create files with data and return with a nested dictionary: {id : {motif : {appearance_number : {'start' : position}}}, protein set, motif set with distinct motifs and motif list with multiplicated motifs
# inputs: sequence library, regex motif, motif length, output directory
def motif_regex_search(sequences_dicty, motif, motif_len, result_dir):
    motif_dicty = {}
    # loop through protins
    for id in sequences_dicty:
        #print(id)
        # use given protein's sequence
        sequence = sequences_dicty[id]['sequence']
        # loop through positions in sequence
        for pos in range(len(sequence) - motif_len + 1):
            # create actual putative motif
            put_motif = sequence[pos:pos+motif_len]
            # check whether is it correct based on regex
            x = re.search(motif, put_motif)
            # if it is correct possible motif then save that motif
            if x:
                #print(put_motif)
                # check whether given protein is already exist in motif_dicty
                if bool(id not in motif_dicty.keys()) == True:
                    # if given protein id not in motif_dicty, then add its id
                    motif_dicty.update({id : {}})
                # Check whether given putative motif already exists in motif_dicty
                if put_motif not in motif_dicty[id].keys():
                    # if given motif not present in motif_dicty of given protein, then add it !!!!!NOTE: MUST ADD 1 TO START POSITION FOR HUMANS!!!!!!
                    motif_dicty[id].update({put_motif : {1 : {'start' : pos}}})
                # if given putative motif not present in motif_dicty
                elif motif_dicty[id].keys():
                    #print(motif_dicty[id][put_motif].keys())
                    cntr = max(motif_dicty[id][put_motif].keys()) + 1
                    motif_dicty[id][put_motif].update({cntr : {'start' : pos}})
        #print(motif_dicty)

    # Count number of motifs and motif-containing proteins
    protein_set, motif_set, motif_all_list = set(), set(), list()
    for id in motif_dicty:
        protein_set.add(id)
        for motif in motif_dicty[id]:
            for cntr in motif_dicty[id][motif]:
                #print(motif)
                motif_set.add(motif)
                motif_all_list.append(motif)

    # Generate result files for all proteins that contain at least one motif
    for id in motif_dicty:
        with open(result_dir + id + '.res', 'w') as f:
            for motif in motif_dicty[id]:
                for cntr in motif_dicty[id][motif]:
                    motif_start = str(motif_dicty[id][motif][cntr]['start'])
                    f.write(motif_start + ' ')
                    f.write(motif)
                    f.write('\n')
    print('\nNumber of proteins: {}\nNumber of distinct motifs: {}\nNumber of all proteins: {}\n'.format(len(protein_set), len(motif_set), len(motif_all_list)))
    return motif_dicty, protein_set, motif_set, motif_all_list

# ----------------------------------------------------------------------------------------------------------------------

# Open and read motif_regex_search output files and create a dictionary for them: {start : [motif], start : [motif]}
def motif_search_result_read_dict_create(file_path):
    protein_motif_dicty = {}
    # open file
    with open(file_path, 'r') as f:
        # read its lines
        for line in f.readlines():
            print(line)
            # extract motif and its position
            motif = line.split()[1]
            start = int(line.split()[0])
            # place motif and its position into output dictionary
            protein_motif_dicty.update({str(start) : motif})
    #print(protein_motif_dicty)
    return protein_motif_dicty

# ----------------------------------------------------------------------------------------------------------------------

########################################################################################################################
# 						Motif candidates analysis
#						based on (1) disorder content (IuPred)
########################################################################################################################

# ----------------------------------------------------------------------------------------------------------------------

# estimate disorder content of given protein with the aid of IuPred software
def iupred_score(seq):
    return iupred(seq)

# ----------------------------------------------------------------------------------------------------------------------

# Calculate iupred score of all motifs in given dictionary generated by multifasta_PSSM_motif_searcher then keep all of them or remove those that have less iupred score than the threshold
# and return with a dictionary: {id : {motif : {app_numb : {'start' : start_index, 'end' : end_index, 'iupred_score' : iupred_score, 'PSSM_score' : score}}}
# inputs: dictionary created by multifasta_PSSM_motif_searcher, sequence library created by impr_multi_FASTA_seq_extractor, keep_all=True if not want remove bad motifs,
# threshold to distinguish good and bad motifs, equal_too=True if want allow equation as well (score >= threshold)
def iupred_filter_use(dicty, human_proteom_dict, keep_all=True, threshold=0.5, equal_too=True):
    # create and place protein ids into output dictionary
    extended_good_motif_dict = {id: {} for id in dicty}
    motifs_numb = 0
    good_motifs_numb = 0
    protein_id_set, motif_types_set = set(), set()
    for id in dicty:  # loop over proteins
        for motif in dicty[id]['motifs']:  # loop through their motifs
            extended_good_motif_dict[id].update({motif: {}})
            for app_numb in dicty[id]['motifs'][motif]:  # loop through their occurrence numbers
                position = dicty[id]['motifs'][motif][app_numb]['position']  # loop through each motif position
                # split and convert the string type position into integer start and end index
                splitter_index = position.index('-')
                start_index = int(position[:splitter_index])
                end_index = int(position[splitter_index + 1:])
                # calculate iupred score by using iupred module from iupred3 package
                protein_seq = human_proteom_dict[id]['sequence']
                protein_iupred_score = iupred_score(protein_seq)
                motif_iupred_score = protein_iupred_score[0][start_index - 1:end_index]
                # take average of motif iupred score
                avg_motif_iupred_value = avg(motif_iupred_score)
                # PSSM score
                PSSM_score = dicty[id]['motifs'][motif][app_numb]['score']
                # put information into result dictionary
                # print(id)
                # print(motif_iupred_score)
                # print(avg_motif_iupred_value)
                motifs_numb += 1
                # if want to keep bad motifs as well, then just add iupred score besides the other information
                if keep_all == True:
                    extended_good_motif_dict[id][motif].update({app_numb : {'start': start_index,
                                                                            'end': end_index,
                                                                            'iupred_score': avg_motif_iupred_value,
                                                                            'PSSM_score': PSSM_score}})
                # if not want to keep bad motifs
                elif keep_all == False:
                    # if allow equation
                    if equal_too == True:
                        # if iupred score of motif is greater or equal to threshold add it to output dictionary
                        if avg_motif_iupred_value >= threshold:
                            protein_id_set.add(id)
                            motif_types_set.add(motif)
                            good_motifs_numb += 1
                            extended_good_motif_dict[id][motif].update({app_numb: {'start': start_index,
                                                                                   'end': end_index,
                                                                                   'iupred_score': avg_motif_iupred_value,
                                                                                   'PSSM_score': PSSM_score}})
                    # if not allow equation
                    elif equal_too == False:
                        # if iupred score of motif is greater than threshold add it to output dictionary
                        if avg_motif_iupred_value > threshold:
                            protein_id_set.add(id)
                            motif_types_set.add(motif)
                            good_motifs_numb += 1
                            extended_good_motif_dict[id][motif].update({app_numb: {'start': start_index,
                                                                                   'end': end_index,
                                                                                   'iupred_score': avg_motif_iupred_value,
                                                                                   'PSSM_score': PSSM_score}})
            # if there is not occurrence of motif then delete it
            if bool(extended_good_motif_dict[id][motif]) == False:
                del extended_good_motif_dict[id][motif]
        # if there is not remaining motif in a protein then remove it
        if bool(extended_good_motif_dict[id]) == False:
            del extended_good_motif_dict[id]
    # if removal was allowed then must calculate remained number of protein, motifs, and occurrences of them and print them out
    if keep_all == False:
        print('\nNumber of remained proteins: {}'.format(len(protein_id_set)))
        print('\nNumber of kind of good motifs: {}'.format(len(motif_types_set)))
        print('\nNumber of all good motifs even in the same protein: {}'.format(good_motifs_numb))


    print('\nNumber of all motifs even in same protein: {}'.format(motifs_numb))
    return extended_good_motif_dict

# ----------------------------------------------------------------------------------------------------------------------

# Create a dictionary {protein_id : [iupred_score]} within there're iupred score of all proteins from a dictionary with this structure {id : {'sequence' : sequence}, id2 : {'sequence' : sequence}}
# or generated by impr_multi_FASTA_seq_extractor and also creates separated files for all proteins with their iupred scores
def iupred_dicty_generate(dicty, outfile_path):
    protein_iupred_dicty = {}
    # loop through proteins
    for protein_id in dicty:
        #print(protein_id)
        # try to calculate iupred score of given protein sequence, if can't then print out the general reason of the fail and give 'Unknown' value instead of score
        try:
            protein_iupred_dicty.update({protein_id : iupred_score(dicty[protein_id]['sequence'])[0]})  # [0] is needed because iupred_score function returns with 2 element list: iupred_position_score and '', but I only need the scores
            #print(iupred_score(dicty[protein_id]['sequence']))
        except:
            print(protein_id)
            print('Iupred can not predict its disorder. It is too short: {}'.format(len(dicty[protein_id]['sequence'])))
            protein_iupred_dicty.update({protein_id : 'Unknown'})
        # Generate output files for each protein
        with open(outfile_path + protein_id + '.iupred', 'w') as f:
            f.writelines([str(score) + ' ' for score in protein_iupred_dicty[protein_id]])
    return protein_iupred_dicty

# ----------------------------------------------------------------------------------------------------------------------

# Create a tsv file from dictionary after IuPred filter carried out by iupred_filter
# inputs: dictionary generated by iupred_filter, outputfile path
def extgm_dict_to_tsv_convert(dicty, file_path):
    # create and open outputfile for writing
    new_file = open(file_path, 'w')
    # Create header: 'Protein id\tmotif\tappearance number\tposition\tPSSM score\tiupred score' as the first row
    cont_list = ['Protein id\tmotif\tappearance number\tposition\tPSSM score\tiupred score\n']
    # loop through proteins
    for id in dicty:
        # loop through motifs
        for motif in dicty[id]:
            # loop through occurrence of motifs
            for app_numb in dicty[id][motif]:
                # extract motif's score and position and iupred score
                motif_start_position = str(dicty[id][motif][app_numb]['start'])
                motif_PSSM_score = str(dicty[id][motif][app_numb]['PSSM_score'])
                motif_iupred_score = str(dicty[id][motif][app_numb]['iupred_score'])
                # assemble and give row representing current motif's information
                cont_string = str(id) + '\t' + str(motif) + '\t' + str(app_numb) + '\t' + motif_start_position + '\t' + motif_PSSM_score + '\t' + motif_iupred_score + '\n'
                cont_list.append(cont_string)
    # write list into the output file
    new_file.writelines(cont_list)
    new_file.close()
    print('Convertation was successful!')
    return 1

# ----------------------------------------------------------------------------------------------------------------------

########################################################################################################################
# 						Motif candidates analysis
#						with (2) Pfam to identify mainly domains
########################################################################################################################

# ----------------------------------------------------------------------------------------------------------------------

# Create pfam result files into given directory from input fasta files
# inputs: pfam program's path, input fasta files' directory path, outputfiles path
def pfam_use(pfam_path, input_dir_path, output_dir_path):
    # Convert the pfam_path into pfam_dir path
    pfam_dir_path = pfam_path.replace('pfam_scan.pl', '')
    print(pfam_dir_path)
    # change actual directory into directory within the pfam program is in
    os.chdir(pfam_dir_path)
    # loop through fasta files
    for file_name in os.listdir(input_dir_path):
        # assemble input and outputfile path
        input_file_path = input_dir_path + file_name
        output_file_name = file_name.replace('.fasta', '_pfam_res.txt')
        output_file_path = output_dir_path + output_file_name
        print('new')
        #print(output_file_name)
        # run pfam
        command_list = ['perl', '-I', '.', pfam_path, '-fasta', input_file_path, '-dir', pfam_dir_path]
        print(command_list)
        # save pfam result
        c1 = subprocess.Popen(command_list, stdout=subprocess.PIPE)
        out, err = c1.communicate()
        c1.stdout.close()
        # if there was an error then print it out
        if err != None:
            print('There was a problem during running of ' + file_name + '!')
            break
            print('Pfam is stuck.')
            return -1
        out = out.decode('utf-8')
        print(out)
        # create and open outputfile for writing and write pfam result into it
        file = open(output_file_path, 'w')
        file.writelines(out)
        file.close()
    print('Pfam ran successfully for all given files.')
    return 1


# ----------------------------------------------------------------------------------------------------------------------

# Removes those motifs from dict that are overlapped with not good 'words' (arbitrary) by using pfam results after applied Iupred filter (!)
# and return with a dictionary: {id : {motif : {app_numb : {'start' : position, 'end' : position, 'PSSM_score' : score, iupred_score : score}}}}
# inputs: dictionary generated by iupred_filter_use, list containing forbidden pfam categories, pfam result files containing directory path
def pfam_filter_bad_words_kill(dicty, bad_words, pfam_files_path):
    killed_motifs_numb = 0
    new_dicty = {}
    cntr = 0
    # create list of all pfam files
    pfam_files_list = os.listdir(pfam_files_path)
    print(pfam_files_list)

    # Compute the number of proteins, motif types and motifs in input directory
    motif_types_org_set = set()
    protein_numb_org, motif_numb_org = 0, 0
    for id in dicty:
        #print('Run')
        protein_numb_org += 1
        for motif in dicty[id]:
            motif_types_org_set.add(motif)
            for app_numb in dicty[id][motif]:
                motif_numb_org += 1

    print('Finish')

    # loop through ids in dictionary
    for id in dicty:
        #print('RUN2')
        new_dicty.update({id : {}})
        # add _pfam_res.txt to id since pfam_files_list contains id's pfam results with that extension
        id_file_name = id + '_pfam_res.txt'
        # find the position of current id in pfam_files_list
        id_file_pos = pfam_files_list.index(id_file_name)
        # index the file that contains current id's pfam result
        pfam_file = pfam_files_list[id_file_pos]
        pfam_file_path = pfam_files_path + pfam_file
        # open the proper pfam output file
        with open(pfam_file_path, 'r') as f:
            # read lines of current file
            pfam_file_cont_list = f.readlines()
            # create a list containing only the important line which don't start with '#' and replace the spaces with commas (',')
            imp_lines_list = [line.replace(' ', ',').split(',') for line in pfam_file_cont_list if
                              not line.startswith('#') and not line.startswith('\n')]
            # if the file of pfam result of current protein is empty, then add all motifs of it to new dictionary
            if bool(imp_lines_list) == False:
                #print('new file {}'.format(id), imp_lines_list)
                for motif in dicty[id]:
                    #print(motif)
                    new_dicty[id].update({motif : {}})
                    for app_numb in dicty[id][motif]:
                        cntr += 1
                        #print(app_numb)
                        new_dicty[id][motif].update({app_numb : {'start' : dicty[id][motif][app_numb]['start'],
                                                                 'end' : dicty[id][motif][app_numb]['end'],
                                                                 'PSSM_score' : dicty[id][motif][app_numb]['PSSM_score'],
                                                                 'iupred_score' : dicty[id][motif][app_numb]['iupred_score']}})


            # If pfam output file is not empty
            elif bool(imp_lines_list) == True:
                # removes the spaces as items from list
                filtrd_imp_lines_list = [list(filter(lambda x : x != '', line)) for line in imp_lines_list]
                print('Original list of {}'.format(id), '\n', imp_lines_list)
                print('The list of {}'.format(id), '\n', filtrd_imp_lines_list)
                # find those lines that contain bad word
                bad_lines_pos_list = []
                # loop through the lines of filtered important lines
                for idx, line in enumerate(filtrd_imp_lines_list):
                    for word in line:
                        if word in bad_words:
                            bad_lines_pos_list.append(idx)
                print('Bad lines position: ', bad_lines_pos_list)
                # if bad lines pos list is not empty
                if bool(bad_lines_pos_list) == True:
                    # find and extract position of motifs of actual protein
                    for motif in dicty[id]:
                        new_dicty[id].update({motif: {}})
                        print(motif)
                        new_app_numb = 0
                        for app_numb in dicty[id][motif]:
                            print(app_numb)
                            start_pos = dicty[id][motif][app_numb]['start']
                            end_pos = dicty[id][motif][app_numb]['end'] + 1
                            motif_range_set = set(range(start_pos, end_pos))
                            # create new app numb since killing can make mess
                            print(motif_range_set)
                            # Create a list containing 0 or 1 based on that if motif is overlapped with domain (1) or not (0)
                            overlap_list = []
                            # Create range for bad words as well
                            for idx in bad_lines_pos_list:
                                print('These are them:\n', filtrd_imp_lines_list[idx])
                                bad_word_range_set = set(range(int(filtrd_imp_lines_list[idx][3]), int(filtrd_imp_lines_list[idx][4]) + 1))
                                print(bad_word_range_set)
                                # If motif range and bad word range are overlapped add 1 to overlap_list
                                if bool(motif_range_set.intersection(bad_word_range_set)) == True:
                                    overlap_list.append(1)
                                # If motif range and bad word range are not overlapped add 0 to overlap list
                                elif bool(motif_range_set.intersection(bad_word_range_set)) == False:
                                    overlap_list.append(0)
                            print('Overlap list: \n', overlap_list)
                            print(sum(overlap_list))
                            # if there're overlap with any bad word
                            if sum(overlap_list) != 0:
                                print(id, motif, app_numb, 'was killed')
                                killed_motifs_numb += 1
                            # if there aren't overlapping
                            elif sum(overlap_list) == 0:
                                new_app_numb += 1
                                new_dicty[id][motif].update({new_app_numb: {'start': start_pos,
                                                                            'end': end_pos,
                                                                            'PSSM_score': dicty[id][motif][app_numb]['PSSM_score'],
                                                                            'iupred_score': dicty[id][motif][app_numb]['iupred_score']}})
                        # if motif not occur in a protein then remove it
                        if bool(new_dicty[id][motif]) == False:
                            del new_dicty[id][motif]

                # if bad lines pos list is empty then save all motifs of current protein id
                elif bool(bad_lines_pos_list) == False:
                    for motif in dicty[id]:
                        new_dicty[id].update({motif : {}})
                        for app_numb in dicty[id][motif]:
                            new_dicty[id][motif].update({app_numb: {'start': dicty[id][motif][app_numb]['start'],
                                                                    'end': dicty[id][motif][app_numb]['end'],
                                                                    'PSSM_score': dicty[id][motif][app_numb]['PSSM_score'],
                                                                    'iupred_score': dicty[id][motif][app_numb]['iupred_score']}})
        # if not motifs in a protein at all remove protein
        if bool(new_dicty[id]) == False:
            del new_dicty[id]


    # Compute the number of motifs, motif types and proteins
    motif_types_new_set = set()
    protein_numb_new, motif_numb_new = 0, 0
    for id in new_dicty:
        protein_numb_new += 1
        for motif in new_dicty[id]:
            motif_types_new_set.add(motif)
            for app_numb in new_dicty[id][motif]:
                motif_numb_new += 1


    print('\nNumber of proteins of input: {}\nNumber of motif types of input: {}\nNumber of motifs of input: {}'.format(protein_numb_org, len(motif_types_org_set), motif_numb_org))
    print('\nRemoved motifs: {}'.format(killed_motifs_numb))
    print('\nRemained number of proteins: {}\nRemained number of motif types: {}\nRemained number of motifs: {}\n'.format(protein_numb_new, len(motif_types_new_set), motif_numb_new))
    print(motif_numb_org-motif_numb_new)

    return new_dicty



# ----------------------------------------------------------------------------------------------------------------------



# ----------------------------------------------------------------------------------------------------------------------



# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------

# Convert ox (taxonomic id) to phylogenetic level then to order number
# inputs: prot = protein list in which there are their headers and sequences generated by head_and_neck_merge,
# sort_dict = dictionary determining the taxonomic level of given protein sequence: {ID : taxonomic level, ID2 : taxonomic level 2}
def ox_to_phyl_lev_to_ord_numb_convert(prot, sort_dict):
    # extract ox from
    ox_start_id = prot.find('OX=') + 3
    ox_end_id = prot.find(' ', ox_start_id)
    ox = prot[ox_start_id:ox_end_id]
    # if it's human, then must be first
    if ox == '9606':
        ord_numb = 0
    # if it's not human, then must convert its ox into phyl_level
    else:
        phyl_lev = sort_dict[ox]
        if phyl_lev == 'Mammalia':
            ord_numb = 1
        elif phyl_lev == 'Vertebrata':
            ord_numb = 2
        elif phyl_lev == 'Eumetazoa':
            ord_numb = 3
        elif phyl_lev == 'Opisthokonta':
            ord_numb = 4
        elif phyl_lev == 'Plant':
            ord_numb = 5
        elif phyl_lev == 'Eukaryota':
            ord_numb = 6
    return ord_numb

# ----------------------------------------------------------------------------------------------------------------------

# Create a new multi-fasta file within the order of fasta files based on phylogenetic levels from given list by using ox_to_phyl_lev_to_ord_numb_convert
# inputs: fasta_list generated from multifasta by head_and_neck_merge, sort_dict = dictionary determining the taxonomic level of given protein sequence: {ID : taxonomic level, ID2 : taxonomic level 2}
# fasta_path = path of output sorted multifasta file
def multi_fasta_cont_sort_file_creator(fasta_list, sort_dict, fasta_path):
    # Sort them by using ox_to_phyl_lev_to_ord_numb_convert
    fasta_list.sort(key = lambda x : ox_to_phyl_lev_to_ord_numb_convert(x, sort_dict))
    # Create new file and write sorted list content into that
    new_file = open(fasta_path, 'w')
    for fasta in fasta_list:
        new_file.writelines(fasta)
    new_file.close()
    print('\nSorted multi-fasta files creation was successful!')
    return 1

# ----------------------------------------------------------------------------------------------------------------------

# Cut off gaps in sequences based on got reference
def indels_free_create(seq, ref_seq):
    # determine length of sequence and reference sequence
    seq_len = len(seq)
    ref_seq_len = len(ref_seq)
    ref_seq2, seq2 = '', ''
    print(seq_len)
    print(ref_seq_len)
    # If reference sequence is longer or same as the sequence use the indices of shorter one (sequence) to loop through
    if ref_seq_len >= seq_len:
        for id, aa in enumerate(seq):
            #print(aa, ref_seq[id])
            # if given position is not gap in reference sequence then add residue to the sequence2 amd current residue of ref. seq. into ref. seq. 2 (skip gaps)
            if ref_seq[id] != '-':
                seq2 += aa
                ref_seq2 += ref_seq[id]
        # add remanining residues except gaps from reference sequence into reference sequence 2 since sequence is shorter
        id += 1
        while id != len(ref_seq):
            if ref_seq[id] != '-':
                ref_seq2 += ref_seq[id]
            id += 1
    # If reference sequence is shorter than the sequence use the indices of shorter one (reference sequence) to loop through
    else:
        for id, aa in enumerate(ref_seq):
            #print(aa, seq[id])
            # if given position is not gap in reference sequence then add residue to the ref.seq.2 and current residue of seq. into seq. 2 (skip gaps)
            if aa != '-':
                seq2 += seq[id]
                ref_seq2 += aa
        id += 1
        # add the remaining amino acids to seq2
        while id != len(seq):
            seq2 += seq[id]
            id += 1
    print('\nOriginal sequence:', seq)
    print('\nOriginal reference sequence:', ref_seq)
    print('\nModified sequence:', seq2)
    print('\nModified reference sequence:', ref_seq2)
    return seq2, ref_seq2

# ----------------------------------------------------------------------------------------------------------------------

# Removes gaps based on reference sequence (human) in a bunch of fasta file in a directory and creates the ungapped versions of them into a select directory
def gap_kill(sorted_mafftalign_dir, output_gapped_sorted_mafft_dir):
    # loop through the input files
    for sorted_alignm_file_name in os.listdir(sorted_mafftalign_dir):
        #print('\nSorted alignment file\'s name:', sorted_alignm_file_name)
        # store the sorted, ungapped file content
        new_sorted_alingm_file_cont = []
        sorted_alignm_file_path = sorted_mafftalign_dir + sorted_alignm_file_name
        # open the current file
        sorted_alignm_file = open(sorted_alignm_file_path, 'r')
        # read lines of current file
        sorted_alignm_file_cont = sorted_alignm_file.readlines()
        #print('\n', merged_align_fasta_list)
        # loop through the lines
        for id in range(len(sorted_alignm_file_cont)):
            # find the human sequence
            if '9606' in sorted_alignm_file_cont[id]:
                # save its sequence (id+1 means that)
                human_header = sorted_alignm_file_cont[id]
                human_seq = sorted_alignm_file_cont[id+1]
                #print('\nHuman header:', human_header)
                new_sorted_alingm_file_cont.append(human_header)
                #print('\nHuman sequence:', human_seq)
                # Remove the gaps from human sequence
                human_seq2, human_ref_seq2 = indels_free_create(human_seq, human_seq)
                new_sorted_alingm_file_cont.append(human_ref_seq2)

    # compare all sequences to reference human sequence and cut off amino acid based on the latter and also remove gaps from it then write the results into file
        for id in range(len(sorted_alignm_file_cont)):
            # extract all sequences except human one and the same ones
            if id % 2 != 0 and not '9606' in sorted_alignm_file_cont[id-1]:
                seq = sorted_alignm_file_cont[id]
                #print('\nSequence:', seq)
                # remove gaps from current sequence based on reference sequence (human)
                seq2, ref_seq2 = indels_free_create(seq, human_seq)
                new_sorted_alingm_file_cont.append(seq2)
            elif id % 2 == 0 and not '9606' in sorted_alignm_file_cont[id]:
                header = sorted_alignm_file_cont[id]
                #print('\nHeader:', header)
                new_sorted_alingm_file_cont.append(header)
        #print('\nNew sorted alignmnent file content:', new_sorted_alingm_file_cont)

    # Create the file and put modified content into that
        output_gapped_sorted_mafft_align_path = output_gapped_sorted_mafft_dir + sorted_alignm_file_name
        new_multifasta_file = open(output_gapped_sorted_mafft_align_path, 'w')
        for line in new_sorted_alingm_file_cont:
            new_multifasta_file.writelines(line)
        new_multifasta_file.close()
        sorted_alignm_file.close()
        print('Gap was removed succesfully and new file was created: {}'.format(new_multifasta_file))
    print('Gap kill ran perfectly for all mafft files!')
    return 1

# ----------------------------------------------------------------------------------------------------------------------

# It does as the same as previous algorithm but for that case when Drosi reference protein is used
def gap_kill_drosi(sorted_mafftalign_dir, output_gapped_sorted_mafft_dir):

    # loop through the input files
    for sorted_alignm_file_name in os.listdir(sorted_mafftalign_dir):
        #print('\nSorted alignment file\'s name:', sorted_alignm_file_name)
        # store the sorted, ungapped file content
        new_sorted_alingm_file_cont = []
        sorted_alignm_file_path = sorted_mafftalign_dir + sorted_alignm_file_name
        # open the current file
        sorted_alignm_file = open(sorted_alignm_file_path, 'r')
        # read lines of current file
        sorted_alignm_file_cont = sorted_alignm_file.readlines()
        #print('\n', merged_align_fasta_list)
        # loop through the lines
        for id in range(len(sorted_alignm_file_cont)):
            # find the human sequence
            if '7227' in sorted_alignm_file_cont[id]:
                # save its sequence (id+1 means that)
                human_header = sorted_alignm_file_cont[id]
                human_seq = sorted_alignm_file_cont[id+1]
                #print('\nHuman header:', human_header)
                new_sorted_alingm_file_cont.append(human_header)
                #print('\nHuman sequence:', human_seq)
                # Remove the gaps from human sequence
                human_seq2, human_ref_seq2 = indels_free_create(human_seq, human_seq)
                new_sorted_alingm_file_cont.append(human_ref_seq2)

    # compare all sequences to reference human sequence and cut off amino acid based on the latter and also remove gaps from it then write the results into file

        for id in range(len(sorted_alignm_file_cont)):
            # extract all sequences except human one and the same ones
            if id % 2 != 0 and not '7227' in sorted_alignm_file_cont[id-1]:
                seq = sorted_alignm_file_cont[id]
                #print('\nSequence:', seq)
                # remove gaps from current sequence based on reference sequence (human)
                seq2, ref_seq2 = indels_free_create(seq, human_seq)
                new_sorted_alingm_file_cont.append(seq2)
            elif id % 2 == 0 and not '7227' in sorted_alignm_file_cont[id]:
                header = sorted_alignm_file_cont[id]
                #print('\nHeader:', header)
                new_sorted_alingm_file_cont.append(header)
        #print('\nNew sorted alignmnent file content:', new_sorted_alingm_file_cont)

    # Create the file and put modified content into that
        output_gapped_sorted_mafft_align_path = output_gapped_sorted_mafft_dir + sorted_alignm_file_name
        new_multifasta_file = open(output_gapped_sorted_mafft_align_path, 'w')
        for line in new_sorted_alingm_file_cont:
            new_multifasta_file.writelines(line)
        new_multifasta_file.close()
        sorted_alignm_file.close()
        print('Gap was removed succesfully and new file was created: {}'.format(new_multifasta_file))
    print('Gap kill ran perfectly for all mafft files!')
    return 1

# ----------------------------------------------------------------------------------------------------------------------

# Score the motif at the same place in other species' proteins where it's in human protein by using PSSM and return with a dictionary: {id : {original_motif : {app_numb : {OX : {motif : PSSM_score, 'start' : start}}}}}
# inputs: gapped and ordered multifasta files containing directory, filtered motif candidate containing dictionary {ID : {}}, matrix to score, gap score
def species_motifs_score(gapped_ordered_mafft_dir, surv_motifs_dir, norm_prob_log_mtx, gap = 0 ):
    motifs_in_species_dict = {}
    # loop through the gapped and ordered mafft alignments
    for gp_ord_mafftalign_file_name in os.listdir(gapped_ordered_mafft_dir):
        # assemble current file and open the it
        gp_ord_mafftalign_file_path = gapped_ordered_mafft_dir + gp_ord_mafftalign_file_name
        gp_ord_mafftalign_file = open(gp_ord_mafftalign_file_path, 'r')
        # read its lines
        gp_ord_mafftalign_cont = gp_ord_mafftalign_file.readlines()
        # create the proper protein id key for extended_good_motifs_signalp_phb_filtrd_dict indexing from filename
        id = gp_ord_mafftalign_file_name[:gp_ord_mafftalign_file_name.index('.')]
        #print('ID:', id)
        # place id into output dictionary
        motifs_in_species_dict.update({id: {}})
        # loop through current protein's all motifs
        for motif in surv_motifs_dir[id]:
            motifs_in_species_dict[id].update({motif : {}})
            print(motif)
            for app_numb in surv_motifs_dir[id][motif]:
                print(app_numb)
                motifs_in_species_dict[id][motif].update({app_numb : {}})
                # extract start and end positions of motif
                start_id = int(surv_motifs_dir[id][motif][app_numb]['start']) - 1   # -1 because positions were stored in human way but Python uses 0-indexing
                print(start_id, start_id + 9)
                # loop through the lines in the current multi-fasta file
                for line_id in range(len(gp_ord_mafftalign_cont)):
                    # if this true, then it's the sequence
                    if line_id % 2 != 0:
                        seq = gp_ord_mafftalign_cont[line_id]
                        # select the motif in the positions of good disordered motif
                        found_motif = seq[start_id:start_id+9]
                        print('Found motif: ', found_motif)
                        # score it with gaps
                        PSSM_score = score_wth_gaps(found_motif, norm_prob_log_mtx, gap)
                        motifs_in_species_dict[id][motif][app_numb][OX].update({found_motif: PSSM_score,
                                                                                'start' : start_id + 1})
                        # print('\nSearched motif:', motif)
                        # print('\nFound motif:   ', found_motif)
                        # print('\nPSSM score:    ', PSSM_score)
                    # if actual line is header
                    else:
                        # extract headers as well
                        header = gp_ord_mafftalign_cont[line_id]
                        # print('Header:', header)
                        # extract ox from header to known the motif is in which species
                        OX_start_id = header.index('OX=') + 3
                        OX_end_id = header.index(' ', OX_start_id)
                        OX = header[OX_start_id:OX_end_id]
                        # print('\n', OX)
                        # Save all information into dictionary
                        motifs_in_species_dict[id][motif][app_numb].update({OX: {}})

    return motifs_in_species_dict

# ----------------------------------------------------------------------------------------------------------------------

# Create a dictionary within the motifs' occurence in other species (OXs) are placed in phylogenetic groups: {id : {original_motif : {app_numb : {group : [ox1, ox2]}}}}}}
# inputs: dictionary created by species_motifs_score algorithm, dictionary to know given OX belongs to which of the taxonomic groups {OX1 : group1, OX2 : group2, OX3 : group1....}
def phyl_group_dict_create(motifs_in_species_dict, qfo_dict):
    phyl_group_motifs_dict = {}
    # loop through motifs in species dict until reach scores
    for id in motifs_in_species_dict:
        print('\nID:', id)
        phyl_group_motifs_dict.update({id : {}})
        for orig_motif in motifs_in_species_dict[id]:
            print('\nOriginal motif:', orig_motif)
            phyl_group_motifs_dict[id].update({orig_motif: {}})
            for app_numb in motifs_in_species_dict[id][orig_motif]:
                print('\napp_numb:', app_numb)
                phyl_group_motifs_dict[id][orig_motif].update({app_numb : {}})
                print('See:', motifs_in_species_dict[id][orig_motif][app_numb].keys())
                for ox in motifs_in_species_dict[id][orig_motif][app_numb]:
                    print('\nOX:', ox)
                    # find the corresponding phylogenetic group based on qfo_dict
                    phyl_group = qfo_dict[ox]
                    print(phyl_group)
                    # append given ox to corresponding taxonomical group
                    if phyl_group not in phyl_group_motifs_dict[id][orig_motif][app_numb].keys():
                        phyl_group_motifs_dict[id][orig_motif][app_numb].update({phyl_group : []})
                    phyl_group_motifs_dict[id][orig_motif][app_numb][phyl_group].append(ox)
    return phyl_group_motifs_dict

# ----------------------------------------------------------------------------------------------------------------------

# Calculate evolutional conservation = number of accepted motifs / number of species from given phylg. group with found motifs in given protein
# then return with a dictionary: {id : {motif : {app_numb : {'start' : position, group1 : [conservation_score, fraction], group2 : [conservation_score, fraction]...}}}}
# input: dictionary generated by phyl_group_dict_create algorithm, dictionary created by species_motifs_score algorithm, threshold to distinguish good and bad motifs
def evol_cons_calculate(phyl_group_motifs_dict, motifs_in_species_dict, threshold = 0):
    evol_cons_dict = {}
    # loop through proteins
    for id in phyl_group_motifs_dict:
        print('\n', id)
        evol_cons_dict.update({id : {}})
        # loop through original motifs in current protein
        for orig_motif in phyl_group_motifs_dict[id]:
            print('\n', orig_motif)
            evol_cons_dict[id].update({orig_motif : {}})
            # loop through all occurrences
            for app_numb in phyl_group_motifs_dict[id][orig_motif]:
                print('\n', app_numb)
                evol_cons_dict[id][orig_motif].update({app_numb : {}})
                # loop through phylogenetic groups
                for phyl_group in phyl_group_motifs_dict[id][orig_motif][app_numb]:
                    print('\n', phyl_group)
                    cntr = 0
                    score_list = []
                    # loop through all species
                    for ox in phyl_group_motifs_dict[id][orig_motif][app_numb][phyl_group]:
                        print('\n', ox)
                        # increase species counter at the case of all species except human one since that was the reference
                        if ox != '9606':
                            cntr += 1
                        # find score of motifs
                        motif_and_start_list = list(motifs_in_species_dict[id][orig_motif][app_numb][ox].keys())
                        print('\nMotif:', motif_and_start_list)
                        score = motifs_in_species_dict[id][orig_motif][app_numb][ox][motif_and_start_list[0]]
                        evol_cons_dict[id][orig_motif][app_numb].update({'start': motifs_in_species_dict[id][orig_motif][app_numb][ox]['start']})
                        print('\nScore:', score)
                        # if motif score is above threshold and not come from human then append its score into the list
                        if score >= threshold and ox != '9606':
                            print('Yes')
                            score_list.append(score)
                        print('\n', score_list)
                    # if no score then no orthologs
                    if not score_list:
                        evol_cons_dict[id][orig_motif][app_numb].update({phyl_group: [0, '0/{}'.format(cntr)]})
                    # if there are score then there are orthologs
                    elif score_list:
                        evol_cons_percent = len(score_list)/cntr
                        print('\nEvol cons percent:', evol_cons_percent)
                        evol_cons_dict[id][orig_motif][app_numb].update({phyl_group : [evol_cons_percent, '{}/{}'.format(len(score_list), cntr)]})

    return evol_cons_dict

# ----------------------------------------------------------------------------------------------------------------------

# Create a tsv file from given evol_con_dict with all phylgroup display
# inputs: dictionary created by evol_cons_calculate algorithm, output filepath, list contains groups of OXs, proportion say inclued or not the fraction (accepted motifs/all motifs)
def evol_dict_to_tsv_convert(dicty, file_path, phyl_group_list, proportion = True):
    # create and open output file for writing
    new_file = open(file_path, 'w')
    cont_string = ''
    # Create the header as firt item of list: Protein id\tmotif\tappearance number\tposition\tphyl_group1\tphyl_group2...
    cont_list = ['Protein id\tmotif\tappearance number\tposition']
    for phyl_group in phyl_group_list:
        cont_list[0] += '\t' + phyl_group + '\t'
    cont_list[0] += '\n'
    print(cont_list)

    # Create the rows after header
    # loop through proteins
    for id in dicty:
        # loop through motifs
        for motif in dicty[id]:
            # loop through appearance number of motif
            for app_numb in dicty[id][motif]:
                cont_string = str(id) + '\t' + str(motif) + '\t' + str(app_numb) + '\t' + str(dicty[id][motif][app_numb]['start'])
                # # loop through phylogenetic groups
                for phyl_group in phyl_group_list:
                    # proportion == False then the fraction are not copied from given dictionary
                    if proportion == False:
                        # extract all evol cons. value if not exist then it's 0
                        if phyl_group not in list(dicty[id][motif][app_numb]):
                            cont_string += '\t0.0'
                        elif phyl_group in list(dicty[id][motif][app_numb]):
                            cont_string += '\t' + str(dicty[id][motif][app_numb][phyl_group][0])
                    # proportion == True then the fraction are copied from given dictionary
                    elif proportion == True:
                        if phyl_group not in list(dicty[id][motif][app_numb]):
                            cont_string += '\t0.0' + '\t0/0'
                        elif phyl_group in list(dicty[id][motif][app_numb]):
                            cont_string += '\t' + str(dicty[id][motif][app_numb][phyl_group][0]) + '\t' + str(dicty[id][motif][app_numb][phyl_group][1])
                cont_string += '\n'
                print(cont_string)
                # add current row to the content of output file
                cont_list.append(cont_string)
    print(cont_list)
    # write content into output file
    new_file.writelines(cont_list)
    new_file.close()
    print('Convertation was successful!')
    return 1

# ----------------------------------------------------------------------------------------------------------------------

# create a dictionary in which the species are seen in their groups: {group1 : [OX, OX2, OX3], group2 : [OX4, OX5]...}
# inputs: qfo_dictionary {OX1 : group1, OX2 : group2, OX3 : group1.....}
def qfo_group_ox_dict(qfo_dict):
    result_qfo_dict = {phyl_group : [] for phyl_group in set(qfo_dict.values())}
    print('\n', result_qfo_dict)
    for ox in qfo_dict:
        phyl_group = qfo_dict[ox]
        result_qfo_dict[phyl_group].append(ox)
    print('\n', result_qfo_dict)
    return result_qfo_dict

# ----------------------------------------------------------------------------------------------------------------------



# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------



# ----------------------------------------------------------------------------------------------------------------------

# Create tsv file from dictionary generated by species_motifs_score algorithm
# inputs: dictionary generated by species_motifs_score algorithm, outputfile_path, dictionary containing sequence library created by impr_multi_FASTA_seq_extractor
def motifs_in_species_dict_to_tsv_convert(dicty, file_path, human_human_dicty):
    # create and open outputfile for writing
    file = open(file_path, 'w')
    cont_list =  []
    # loop through proteins
    for id in dicty:
        cntr = 0
        # extract gene name that encodes given protein
        gene_name = human_human_dicty[id]['GN']
        # loop through motifs
        for orig_motif in dicty[id]:
            # loop through appearance of motif
            for app_numb in dicty[id][orig_motif]:
                # # loop through species that have motif
                for ox in dicty[id][orig_motif][app_numb]:
                    # extract motifs and start position of them from a given species
                    motifs_set = set(motif for motif in dicty[id][orig_motif][app_numb][ox].keys() if motif != 'start')
                    start_pos = dicty[id][orig_motif][app_numb][ox]['start']
                    # loop through motifs
                    for motif in motifs_set:
                        # try to extract PSSM score of given motif if it exist if not then skip this step
                        try:
                            PSSM_score = dicty[id][orig_motif][app_numb][ox][motif]
                        except:
                            continue
                        # if given motif comes from human
                        if ox == '9606':
                            # create header and information of motif to the file (header is only created in this case, because human sequence must be the first)
                            cont_string = 'Protein id\tgene name\toriginal motif\tappearance\tPosition\tPSSM_score\n'
                            cont_list.append(cont_string)
                            cont_string = '' + str(id) + '\t' + str(gene_name) + '\t' + str(orig_motif) + '\t' + str(app_numb) + '\t' + str(start_pos) + '\t' + str(PSSM_score) + '\n'
                            cont_list.append(cont_string)
                        # if not human originated motif then just assemble the information of motif and place into the output file
                        elif ox != '9606':
                            cont_string = '\t\t' + str(motif) + '\t' + str(ox) + '\t\t' + str(PSSM_score) + '\n'
                            cont_list.append(cont_string)
    # write information into outputfile
    file.writelines(cont_list)
    print('Convertation was successful!')
    return 1

# ----------------------------------------------------------------------------------------------------------------------



# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------



# ----------------------------------------------------------------------------------------------------------------------

# It runs signalp for each protein based on a dictionary if given protein fasta not exist in signalp_output_dir then it makes its fasta file into fasta directory and run signalp for it
# inputs: dictionary within at least the protein ids must be in, path of directory that contains fasta files, outputdirectory in which the results of SignalP are placed
def signalp_run(dicty, fasta_dir, signalp_output_dir):
    # Check filenames of signalp_output_dir
    filenames = os.listdir(signalp_output_dir)
    # if signalp_output_dir not empty
    if filenames:
        print('There are files in signalp_output_dir!')
        # loop through protein IDs
        for id in dicty:
            # if given protein id not present in signalp_output_dir
            if id not in filenames:
                # Make fasta file name and path into fasta_dir
                id_file_name = id + '.fasta'
                fasta_file_path = fasta_dir + id_file_name
                # Make a directory (protein_id) for output directory of signalp6, since signalp has 6 outputs for a protein which not present its name
                out_file_dir = signalp_output_dir + id
                subprocess.run(['mkdir', out_file_dir])
                # Run signalp6
                command_list = ['signalp6', '--fastafile', fasta_file_path, '--output_dir', out_file_dir, '--organism', 'eukarya']
                subprocess.run(command_list)
            # if given protein id is already present in signalp_output_dir
            else:
                print('File is already existed for {} protein!'.format(id))

    # if signalp_output_dir is empty
    elif not filenames:
        # loop through protein IDs
        print('Not files in signalp_output_dir!')
        for id in dicty:
            # Make a directory (protein_id) for output directory of signalp6, since signalp has 6 outputs for a protein which not present its name
            out_file_dir = signalp_output_dir + id
            subprocess.run(['mkdir', out_file_dir])
            # Make fasta file path
            fasta_file_path = fasta_dir + id + '.fasta'
            # Run signalp6
            command_list = ['signalp6', '--fastafile', fasta_file_path, '--output_dir', out_file_dir, '--organism', 'eukarya']
            subprocess.run(command_list)
            print('SignalP ran successfully for {} protein.'.format(id))
    print('Signalp6 ran successfully for all given files.')
    return 1

# ----------------------------------------------------------------------------------------------------------------------

# Removes motifs that located extracellular or transmembrane part of protein
# Returns with 2 directories: 1, contains survived motifs; 2, contains died motifs
# inputs: dictionary with at least the protein ids, path to phobius program, directory containing fasta files, phobius output directory
def phobius_run(dicty, phobius_path, fasta_dir, phb_output_dir):
    # loop through proteins
    for id in dicty:
        # Make fasta file name and path
        id_file_name = id + '.fasta'
        id_file_path = fasta_dir + id_file_name
        # Run phobius and save its result into file
        command_list = ['perl', phobius_path, id_file_path]
        c1 = subprocess.Popen(command_list, stdout=subprocess.PIPE)
        out, err = c1.communicate()
        c1.stdout.close()
        if err != None:
            print('There was a problem during running this file:', id_file_name)
            print('Phobius is stuck!')
            break
            return -1
        out = out.decode('utf-8')
        phb_output_path = phb_output_dir + id + '.phb_output'
        file = open(phb_output_path, 'w')
        file.writelines(out)
        file.close()
    print('Phobius ran successfully for all given files.')
    return 1

# ----------------------------------------------------------------------------------------------------------------------

# Removes those motifs from given dictionary which overlap with any unallowed phobius category
# Returns with 2 directories: 1, contains survived motifs; 2, contains died motifs;
# both directories: {id : {motif : {app_numb : {'start' : position, 'end' : position, 'PSSM_score' : score, iupred_score : score}}}}
# inputs: dictionary generated by iupred_filter_use, directory path that contain phobius results, list contains forbidden phobius categories
def phobius_kill(dicty, phb_output_dir, bad_words):
    killed_motifs_numb = 0
    # copy input dicty
    dicty_copy = copy.deepcopy(dicty)
    # Make dictionary to store removed motifs from original one
    deleted_dicty = {}

    # Compute the number of all motifs, motif types and proteins
    protein_numb_old, motif_types_old_set, motif_numb_old = 0, set(), 0
    for id in dicty_copy:
        protein_numb_old += 1
        for motif in dicty_copy[id]:
            motif_types_old_set.add(motif)
            for app_numb in dicty_copy[id][motif]:
                motif_numb_old += 1

    # Main part: removing motifs that are in extracellular or transmembrane region of proteins
    # loop through proteins
    for id in dicty:
        print('ID:', id)
        deleted_dicty.update({id : {}})
        # open the corresponding phobius output
        file_name = id + '.phb_output'
        file_path = phb_output_dir + file_name
        bad_word_range_set_list = []
        with open(file_path, 'r') as f:
            for line in f.readlines():
                print('Line:', line)
                line_words_list = line.split()
                # Assemble NON and CYTOPLASMIC words since they are split away during line.split() method
                if len(line_words_list) == 6:
                    if line_words_list[4] == 'NON':
                        corr_word = line_words_list[4] + ' ' + line_words_list[5]
                        del line_words_list[4:6]
                        line_words_list.append(corr_word)

                # if line contains one bad word
                for word in line_words_list:
                    print('Word:', word)
                    if word in bad_words:
                        # extract positions of bad words
                        start_pos = int(line_words_list[2])
                        end_pos = int(line_words_list[3])
                        bad_word_range = set(range(start_pos, end_pos + 1))
                        bad_word_range_set_list.append(bad_word_range)
                        print('Bad word range:', bad_word_range, '\nBad word range set list:', bad_word_range_set_list)
        # loop through motifs
        for motif in dicty[id]:
            print('Motif:', motif)
            deleted_dicty[id].update({motif : {}})
            # loop through its occurrence
            for app_numb in dicty[id][motif]:
                print('App_numb:', app_numb)
                # extract motif's position
                motif_start = int(dicty[id][motif][app_numb]['start'])
                motif_end = int(dicty[id][motif][app_numb]['end'])
                motif_range = set(range(motif_start, motif_end + 1))
                # check overlap between bad words and current motif
                interception_list = [motif_range.intersection(rng) for rng in bad_word_range_set_list]
                print('Motif range:', motif_range)
                # if motif overlaps with any bad word range
                if bool(interception_list) == True:
                    #print('Yes!')
                    deleted_dicty[id][motif].update({app_numb : {'start' : motif_start,
                                                                 'end' : motif_end,
                                                                 'iupred_score' : dicty[id][motif][app_numb]['iupred_score'],
                                                                 'PSSM_score' : dicty[id][motif][app_numb]['PSSM_score']}})
                    killed_motifs_numb += 1
                    # remove the bad motif from copy dictionary
                    del dicty_copy[id][motif][app_numb]
                    print('Removed motif:', id, motif, app_numb)
        # if not remaining occurrences of motif then remove it
        if bool(dicty_copy[id][motif]) == False:
            del dicty_copy[id][motif]
        if bool(deleted_dicty[id][motif]) == False:
            del deleted_dicty[id][motif]
    # if not remaining motif in a protein then remove it
    if bool(dicty_copy[id]) == False:
        del dicty_copy[id]
    if bool(deleted_dicty[id]) == False:
        del deleted_dicty[id]

    # Compute the number of all proteins, motifs and motif types
    protein_numb_new, motif_types_new_set, motif_numb_new = 0, set(), 0
    for id in dicty_copy:
        protein_numb_new += 1
        for motif in dicty_copy[id]:
            motif_types_new_set.add(motif)
            for app_numb in dicty_copy[id][motif]:
                motif_numb_new += 1

    print('\nNumber of proteins of input: {}\nNumber of motif types of input: {}\nNumber of motifs of input: {}'.format(protein_numb_old, len(motif_types_old_set), motif_numb_old))
    print('\nRemoved motifs: {}'.format(killed_motifs_numb))
    print('\nRemained number of proteins: {}\nRemained number of motif types: {}\nRemained number of motifs: {}\n'.format(protein_numb_new, len(motif_types_new_set), motif_numb_new))
    print(motif_numb_old - motif_numb_new)

    return dicty_copy, deleted_dicty

# ----------------------------------------------------------------------------------------------------------------------

# Remove those motifs that have signalpeptide and are in non cytoplasmatic or transmembrane region of protein
# Returns with 2 directories: 1, contains survived motifs; 2, contains died motifs;
# both directories: {id : {motif : {app_numb : {'start' : position, 'end' : position, 'PSSM_score' : score, iupred_score : score}}}}
def signalp_phobius_kill(dicty, signalp_dir, phobius_dir, bad_words_list):
    dicty_copy = copy.deepcopy(dicty)
    deleted_dicty = {}
    killed_motifs_numb = 0
    # Compute the number of all proteins, motifs and motif types
    protein_numb_old, motif_types_old_set, motif_numb_old = 0, set(), 0
    for id in dicty_copy:
        protein_numb_old += 1
        for motif in dicty_copy[id]:
            motif_types_old_set.add(motif)
            for app_numb in dicty_copy[id][motif]:
                motif_numb_old += 1

    # Main part
    # loop through proteins
    for id in dicty:
        # print(id)
        deleted_dicty.update({id: {}})
        # find its signalp file
        signalp_file_path = signalp_dir + id + '/output.gff3'
        # open and extract information from it
        with open(signalp_file_path, 'r') as signalp_file:
            good_lines_list = [line for line in signalp_file.readlines() if not line.startswith('##')]
        # If there's signal peptide then check the phobius output
        if good_lines_list:
            # print('True', good_lines_list)
            # find given protein's phobius output
            phobius_file_path = phobius_dir + id + '.phb_output'
            # open and read its phobius result file
            with open(phobius_file_path, 'r') as phobius_file:
                bad_word_range_set_list = []
                for line in phobius_file.readlines():
                    # print('Line:', line)
                    line_words_list = line.split()
                    # Correct the NON CYTOPLASMIC. if it's in the file since line.split() splits NON CYTOPLASMIC. into two distinct list items
                    if len(line_words_list) == 6:
                        if line_words_list[4] == 'NON':
                            corr_word = line_words_list[4] + ' ' + line_words_list[5]
                            del line_words_list[4:6]
                            line_words_list.append(corr_word)

                    # if line contains one bad word
                    for word in line_words_list:
                        # print('Word:', word)
                        if word in bad_words_list:
                            # print('Yes bad word!!', word)
                            # extract positions
                            start_pos = int(line_words_list[2])
                            end_pos = int(line_words_list[3])
                            bad_word_range = set(range(start_pos, end_pos + 1))
                            bad_word_range_set_list.append(bad_word_range)
                    # print('\nBad word range set list:', bad_word_range_set_list)
            # loop through motifs
            for motif in dicty[id]:
                # print('Motif:', motif)
                deleted_dicty[id].update({motif: {}})
                for app_numb in dicty[id][motif]:
                    # print('App_numb:', app_numb)
                    # extract position of actual motif
                    motif_start = int(dicty[id][motif][app_numb]['start'])
                    motif_end = int(dicty[id][motif][app_numb]['end'])
                    motif_range = set(range(motif_start, motif_end + 1))
                    # check overlap between motif's and eny bad word's range
                    interception_list = [motif_range.intersection(rng) for rng in bad_word_range_set_list]
                    interception_list_result = [1 for inters in interception_list if inters]
                    # print('Motif range:', motif_range)
                    # print('Interception list:', interception_list)
                    # if motif overlaps with any bad word range remove it and place into deleted_dicty
                    if interception_list_result:
                        # print('Yes!')
                        # print('Interception list result:', interception_list_result)

                        deleted_dicty[id][motif].update({app_numb: {'start': motif_start,
                                                                    'end': motif_end,
                                                                    'iupred_score': dicty[id][motif][app_numb][
                                                                        'iupred_score'],
                                                                    'PSSM_score': dicty[id][motif][app_numb][
                                                                        'PSSM_score']}})
                        killed_motifs_numb += 1
                        del dicty_copy[id][motif][app_numb]
                        # print('Removed motif:', id, motif, app_numb)
                # Check for emtpy motif keys in both dictionary
                if not deleted_dicty[id][motif]: del deleted_dicty[id][motif]
                if not dicty_copy[id][motif]: del dicty_copy[id][motif]
        # Check for emtpy protein keys in both dictionary
        if not deleted_dicty[id]: del deleted_dicty[id]
        if not dicty_copy[id]: del dicty_copy[id]

    # Compute the number of all proteins, motifs and motif types
    protein_numb_new, motif_types_new_set, motif_numb_new = 0, set(), 0
    for id in dicty_copy:
        protein_numb_new += 1
        for motif in dicty_copy[id]:
            motif_types_new_set.add(motif)
            for app_numb in dicty_copy[id][motif]:
                motif_numb_new += 1

    print('\nNumber of proteins of input: {}\nNumber of motif types of input: {}\nNumber of motifs of input: {}'.format(
        protein_numb_old, len(motif_types_old_set), motif_numb_old))
    print('\nRemoved motifs: {}'.format(killed_motifs_numb))
    print(
        '\nRemained number of proteins: {}\nRemained number of motif types: {}\nRemained number of motifs: {}\n'.format(
            protein_numb_new, len(motif_types_new_set), motif_numb_new))
    print(motif_numb_old - motif_numb_new)

    return dicty_copy, deleted_dicty

# ----------------------------------------------------------------------------------------------------------------------

# Run gopher for non-existing files in a directory within there are already ones the gopher have already ran for by using dictionary
# inputs: dictionary that contains at least proteins' ids, path to fasta files containing directory, path to gopher program, path to directory of gopher results, path to database for gopher, path to directory of already existing gopher results
def gopher_run(dicty, fasta_file_dir, gopher_path, gopher_output_dir, database, copied_maftalign_dir):
    # read filenames in directory copied_maftalign_dir
    filenames = os.listdir(copied_maftalign_dir)
    # loop through proteins in dictionary
    for id in dicty:
        # assemble corresponding gopher result filename
        mafft_filename = id + '.orthaln.fas'
        print('File_name: ', mafft_filename)
        # if result for given protein not existing then run gopher for it
        if mafft_filename not in filenames:
            fasta_file_path = fasta_file_dir + id + '.fasta'
            subprocess.run(['python', gopher_path, 'seqin={}'.format(fasta_file_path), 'orthdb={}'.format(database),'gopherdir={}'.format(gopher_output_dir), 'orthblast', 'orthfas', 'orthalign'])
        # if result for given protein has already been existing then not run gopher for it, just print the text (see below in code)
        elif mafft_filename in filenames:
            print('The file is already existing for this protein: {}'.format(mafft_filename))
        print('Gopher successfully ran on {}'.format(id))
    print('Gopher ran successfully all selected files!')
    return 1

# ----------------------------------------------------------------------------------------------------------------------

# Run gopher for selected all files in given directory
# inputs: path to fasta files containing directory, path to gopher program, path to directory of gopher results, path to database for gopher
def gopher_dir_run(fasta_file_dir, gopher_path, gopher_output_dir, database):
    # List file names of fasta file containing directory
    file_names = os.listdir(fasta_file_dir)
    # loop through file names and generate file path for all of them
    for file_name in file_names:
        print('File name: ', file_name)
        file_path = fasta_file_dir + file_name
        # run gopher for all proteins
        subprocess.run(['python', gopher_path, 'seqin={}'.format(file_path), 'orthdb={}'.format(database), 'gopherdir={}'.format(gopher_output_dir), 'orthblast', 'orthfas', 'orthalign'])
    print('Gopher ran successfully all selected files!')
    return 1

# ----------------------------------------------------------------------------------------------------------------------

# Counts number of species which its proteome is presented in given fasta file
def species_count(fasta_file_path):
    # open given fasta file, read it, then extract each distinct OX values referring to different species and print result out
    with open(fasta_file_path, 'r') as f:
            sp_codes_set = set(line[line.index('OX=')+3:line.index(' ', line.index('OX='))] for line in f.readlines() if line.startswith('>'))
            print('Species: {}'.format(sp_codes_set))
            print('Number of species: {}'.format(len(sp_codes_set)))
    return sp_codes_set

# ----------------------------------------------------------------------------------------------------------------------

# Retrieve scientific names of taxon IDs from given list and create a tsv file (data retrieved from: NCBI taxonomy browser: https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi)
def species_decode(species_list, output_file_path):
    # loop through all species in species list and retrieve their scientific names which are stored as the values of corresponding 'tax_id'
    species_dicty = {tax_id : BeautifulSoup(requests.get('https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=ownid&lvl=3&lin=f&keep=1&srchmode=1&unlock'.replace('ownid', tax_id)).text, 'html.parser').find('strong').get_text() for tax_id in species_list}
    print('Number of species: ', len(list(species_dicty.keys())))
    # generate tsv file from dictionary
    pd.DataFrame.from_dict(species_dicty, orient='index', columns=['Scientific name']).to_csv(output_file_path, sep='\t', index_label='Taxon ID')
    print('Tsv file was created successfully!')
    return 1

# ----------------------------------------------------------------------------------------------------------------------

# Create tsv file in which all ortholog proteins' peptides where motif is saved are
# inputs: dictionary {protein id : {start : motif...}...} created by motif_search_result_read_dict_create, sequence library generated by impr_multi_FASTA_seq_extractor, path to directory of gopher results, path to directory of ungapped gopher results, path to output tsv le
def species_motif_tsv_generate(dicty, sequence_library_dicty, gopher_result_dir, gopher_gapped_result_dir, outputfile_path):
    file_cont_list = []
    file_cont_string = ''
    # create header
    header_part = 'protein id\tgene id\tstart\tspecies\tmotif\n'
    # loop through proteins
    for protein_id in dicty:
        cntr = 0
        file_cont_string += header_part
        print('Protein id: ', protein_id)
        # extract gene id of given protein
        gene_id = sequence_library_dicty[protein_id]['GN']
        # assemble protein id and gene id in a tsv format
        file_cont_string += protein_id + '\t' + gene_id + '\t'
        print(dicty[protein_id])
        # loop through start positions of motifs from given proteins
        for start in dicty[protein_id]:
            print('\nMotif start: ', start)
            # add start position into list that contains content of output tsv file
            file_cont_string += start + '\t'
            # extract length and the motif of itself from dictionary
            motif_len = len(dicty[protein_id][start][0])
            motif = dicty[protein_id][start][0]
            print('\nMotif: ', motif)
            print(motif_len)
            # open corresponding gopher result file and read it
            with open(gopher_result_dir + protein_id + '.orthaln.fas', 'r') as f:
                gopher_outputfile_cont = f.readlines()
                # if motif is found in first sequence in gopher result file then extract its start and end positions
                if gopher_outputfile_cont[1].find(motif) != -1:
                    motif_start_in_seq = gopher_outputfile_cont[1].find(motif)
                    motif_end_in_seq =  motif_start_in_seq + motif_len
                    # loop through rest of the lines in gopher result file
                    for line in gopher_outputfile_cont:
                        # if line starts with '>' then it is header and extract species name from it
                        if line.startswith('>'):
                            sp_name_start = line.index('OS=') + 3
                            sp_name_end = line.index(' ', line.index(' ', sp_name_start)+1)
                            sp_name = line[sp_name_start:sp_name_end]
                            # if the motif from first sequence then it is different than others (not need tabulator before it)
                            if cntr == 0:
                                file_cont_string += sp_name + '\t'
                            # if the motif not from first sequence then 3 tabulators are needed before it to visualise all motifs below each other
                            elif cntr != 0:
                                file_cont_string +=  '\t' + '\t' + '\t' + sp_name + '\t'
                        # if line not starts with '>' then extract peptide from sequence where motif is in reference
                        else:
                            # if the motif from first sequence was in this cycle then must increase counter by one otherwise not
                            if cntr == 0:
                                file_cont_string += line[motif_start_in_seq:motif_start_in_seq+motif_len] + '\n'
                                cntr += 1
                            elif cntr != 0:
                                file_cont_string += line[motif_start_in_seq:motif_start_in_seq+motif_len] + '\n'
                # if not motif in the first sequence in gopher result file, then it should contain gaps
                elif gopher_outputfile_cont[1].find(motif) == -1:
                    # open ungapped version of sequences represented in another file and read it
                    with open(gopher_gapped_result_dir + protein_id + '.orthaln.fas', 'r') as f:
                        gopher_gapped_cont = f.readlines()
                        # extract start position and end position of given motif
                        motif_start_in_seq = gopher_gapped_cont[1].find(motif)
                        motif_end_in_seq = motif_start_in_seq + motif_len
                        # loop through lines
                        for line in gopher_gapped_cont:
                            # if line starts with '>' then its a header, extract OS value
                            if line.startswith('>'):
                                sp_name_start = line.index('OS=') + 3
                                sp_name_end = line.index(' ', line.index(' ', sp_name_start) + 1)
                                sp_name = line[sp_name_start:sp_name_end]
                                file_cont_string += '\t' + '\t' + '\t' + sp_name + '\t'
                            # if line not start with '>' then its sequence must extract motif by using start position and length of reference motif
                            else:
                                file_cont_string += line[motif_start_in_seq:motif_start_in_seq + motif_len] + '\t' + '\n'
                print(file_cont_string)
    # create outputfile and write its content into it
    with open(outputfile_path, 'w') as outf:
        outf.writelines(file_cont_string)

# ----------------------------------------------------------------------------------------------------------------------



# ----------------------------------------------------------------------------------------------------------------------

# Integrate all information of motifs into files (iupred, pfam, phobius) depending on which of them selected and return with dictionary: {start : [motif, iupred, pfam, phobius]} New dicty: {start : [motif, iupred, pfam, phobius, uniprot subcell, uniprot topology]}
def integrate_all(dicty, iupred_file_path=False, pfam_file_path=False, phobius_file_path=False):
    # Integrate Iupred score if path was given
    if iupred_file_path:
        lines_list = []
        # Open and read iupred_dictionary
        with open(iupred_file_path, 'r') as f:
        # It contains only one line and I can go through all characters in the string by just one line convert strings into floats
            for line in f.readlines()[0].split():
		# convert strings into floats
                lines_list.append(float(line))
    # Loop through dictionary to find actual motif and extract its position which is used to extract Iupred scores belonging to that motif position
        for start in dicty:
            motif = dicty[start]
            #print(lines_list[int(start):int(start)+len(motif)])
	    # take average of corresponding position-specific Iupred scores
            motif_iupred_score = np.mean(lines_list[int(start):int(start)+len(motif)])
            # place average iupred score as an list item
            dicty[start] = [motif, motif_iupred_score]

    # Integrate Pfam results if path was given
    if pfam_file_path:
        # Open pfam results of given protein
        with open(pfam_file_path, 'r') as f:
            # read lines of current file
            pfam_file_cont_list = f.readlines()
            # create a list containing only the important lines which don't start with '#' or new line character (\n), and replace the spaces with commas (',')
            imp_lines_list = [line.replace(' ', ',').split(',') for line in pfam_file_cont_list if
                              not line.startswith('#') and not line.startswith('\n')]
            # loop through motif of given protein
            for start in dicty:
		# extract motif itself
                motif = dicty[start][0]
                # if the file of pfam result of current protein is empty, then add all motifs of it to new dictionary
                if bool(imp_lines_list) == False:
                    # add '-' to dictionary list
                    dicty[start].append('-')
		# if the file of pfam result is not empty
                elif bool(imp_lines_list) == True:
                    # removes the spaces as items from list
                    filtrd_imp_lines_list = [list(filter(lambda x: x != '', line)) for line in imp_lines_list]
                    #print('Original list, '\n', imp_lines_list)
                    #print('The filtered list', '\n', filtrd_imp_lines_list)
                    # Create motif range set. Set is good choice since of its interception function
                    # Add 1 to start because it represented 0-based indexing, but pfam uses 1-based indexing
                    motif_range_set = set(range(int(start)+1, int(start)+1+len(motif)))
                    #print(motif_range_set)
                    # Create overlap_list variable for later usage (see later)
                    overlap_list = []
                    # loop through all important lines one by one
                    for line_id, line in enumerate(filtrd_imp_lines_list):
                        #print(line)
                        # extract start and end of given pfam object
                        pfam_object_start = int(line[1])
                        pfam_object_end = int(line[2])
                        # create range of given pfam object
                        pfam_object_range_set = set(range(pfam_object_start, pfam_object_end+1))
                        #print(pfam_object_range_set)
                        # Check for overlaps between motif_range_set and pfam_object_range_set
                        overlaps = motif_range_set.intersection(pfam_object_range_set)
                        #print('Overlaps: ', overlaps)
                        # If motif_range_set and pfam_object_range_set are overlapped add 1 to overlap_list
                        if bool(motif_range_set.intersection(pfam_object_range_set)) == True:
                            overlap_list.append(1)
                        # If motif_range_set and pfam_object_range_set  are not overlapped add 0 to overlap list
                        elif bool(motif_range_set.intersection(pfam_object_range_set)) == False:
                            overlap_list.append(0)
                    # if sum of overlap_list equals to 0 then add '-' sign to dicty representing given motif not overlap with any pfam object
                    #print(overlap_list)
                    if sum(overlap_list) == 0:
                        dicty[start].append('-')
                    # if there is overlap between motif_range_set and pfam_object_range_set, then extract the pfam object's name from given line and add it to dicty
                    elif sum(overlap_list) != 0:
                        pfam_object_line_id_array = np.where(np.array(overlap_list) == 1)[0]
                        #print(pfam_object_line_id_array)
                        for line_id in pfam_object_line_id_array:
                            #print(filtrd_imp_lines_list[line_id])
                            dicty[start].append(filtrd_imp_lines_list[line_id][7])

    # Integrate Phobius results if path was given
    # EXTRACT DATA FROM PHOBIUS OUTPUT
    if phobius_file_path:
        # open Phobius result
        with open(phobius_file_path, 'r') as f:
            # read lines
            lines_list = f.readlines()
            # loop through lines
            for line_id, line in enumerate(lines_list):
                # if line not start with '//'
                if not line.startswith('//'):
                    #print('Line:', line)
                    # split line where spaces are and return with a list
                    line_words_list = line.split()
                    # if given list contains 6 words and 5th word is NON then it must correct it by joining with next word resulting NON CYTOPLASMATIC otherwise NON and CITOPLASMIC would be two distinct list items
                    if len(line_words_list) == 6:
                        if line_words_list[4] == 'NON':
                            corr_word = line_words_list[4] + ' ' + line_words_list[5]
                            del line_words_list[4:6]
                            line_words_list.append(corr_word)
                            lines_list[line_id] = line_words_list
                    # replace each line of lines_list with split version words list for easier search and indexing later
                    lines_list[line_id] = line_words_list
                # Remove '//' as last line
                elif line.startswith('//'):
                    del lines_list[line_id]
            #print(lines_list)

            # Run this section if want to extract Phobius output information and uniprot subcellular location and topology if that is also available
            # check for signal AND transmem in protein
            decision_help_list = [line[1] for line in lines_list if line[1] == 'SIGNAL' or line[1] == 'TRANSMEM']
            # loop through motifs of given protein
            for start in dicty:
                # create motif_range_set. Set is good choice since of its interception function
                motif = dicty[start][0]
                # EXTRACT PHOBIUS OUTPUT OF MOTIF
                # Add 1 to start because it represented 0-based indexing, but phobius uses 1-based indexing
                motif_range_set = set(range(int(start) + 1, int(start) + 1 + len(motif)))
                # print('See line: ', lines_list)
                line_id = 1
                # extract start and end position of given protein part from first line which is not header
                protein_part_range_set = set(range(int(lines_list[line_id][2]), int(lines_list[line_id][3]) + 1))
                # print(protein_part_range_set)
                # Go line by line until there is overlap between motif_range_set and protein_part_range_set
                while not motif_range_set.intersection(protein_part_range_set):
                    line_id += 1
                    protein_part_range_set = set(range(int(lines_list[line_id][2]), int(lines_list[line_id][3]) + 1))
                # 1, CASE if there is signal in protein then add SIGNAL into dicty and then add extracted topology of motif from Phobius output as a string to the SIGNAL string
                if 'SIGNAL' in decision_help_list:
                    dicty[start].append('SIGNAL ')
		    # if motif is in signal part, then there is not 4th element of given line, thereby must extract 1st element of it
                    try:
                        dicty[start][3] += lines_list[line_id][4]
                    except:
                        dicty[start][3] += lines_list[line_id][1] 
                if 'SIGNAL' not in Phobius output then save only the Phobius category within the motif is without 'SIGNAL' category into the dictionary
		elif 'SIGNAL' not in decision_help_list:
                    dicty[start].append(lines_list[line_id][4])

                # RETRIEVE UNIPROT DATA OF GIVEN PROTEIN
                BASE = "http://www.uniprot.org"
                KB_ENDPOINT = '/uniprot/'
                # TOOL_ENDPOINT = '/uploadlists/'
                payload = {'query': protein_id,
                           'format': 'tab',
                           # 'columns': 'id,entry_name,reviewed,protein_names,organism,ec,keywords'
                           'columns': 'id,entry_name,feature(TRANSMEMBRANE),feature(TOPOLOGICAL DOMAIN),comment(SUBCELLULAR LOCATION)'
                           # comment the following line to exclude isoforms
                           # 'include': 'yes',
                           }
                protein_uniprot_topology = requests.get(BASE + KB_ENDPOINT, params=payload).content
                protein_uniprot_topology = protein_uniprot_topology.decode('ascii')
                # Extract important line because it can sometimes retrieve multiple proteins' topology (reason not known)
                important_protein_topology = protein_uniprot_topology[protein_uniprot_topology.find(protein_id):protein_uniprot_topology.find('\n', protein_uniprot_topology.find(protein_id))]
                #print('\nOriginal important protein topology: ', important_protein_topology)
                # split it at separaters: '\t'
                important_protein_topology = important_protein_topology.split('\t')
                print('\nSplit important protein topology: ', important_protein_topology)
                
		# EXTRACT SUBCELLULAR LOCATION OF GIVEN PROTEIN FROM UNIPROT DATA
                # if there is subcellular location
                if important_protein_topology[4]:
                    # extract subcellular location from uniprot output without 'NOTE='
                    subcell_loc = important_protein_topology[4][:important_protein_topology[4].find('Note=')]
                    #print('\nSubcell loc: ', subcell_loc)
                    # loop through words )list items) in subcell_loc list and remove 'SUBCELLULAR LOCATION: '
                    subcell_predata = [re.sub('SUBCELLULAR LOCATION: ', '', subcell_word) for subcell_word in subcell_loc.split(';')]
                    #print('\nSubcell predata: ', subcell_predata)
		    # loop through words in subcell_predata list which not contain 'SUBCELLULAR LOCATION: '
                    for word_id, word in enumerate(subcell_predata):
			# if there is period in given word, then there are multiple locations after each other --> must be split them into single list items + remove the original word to avoid duplications
                        if word.find('.') != -1:
                            #print('word: ', word)
                            subcell_predata += word.split('.')
                            del subcell_predata[word_id]
                    #print('\nSubcell predata again: ', subcell_predata)
		    # assign properly perpared word list to subcell_data list
                    subcell_data = subcell_predata
                    #print('\nSubcellular location data: ', subcell_data)
		    # loop through words in subcell_data list to remove more meaningless extra information which are found with the presence of asterisk
                    subcell_data_list = [re.sub('{.*.}', '', subcell_word).strip() for subcell_word in subcell_data]
                    #print('\nBetter subcellular location data: ', subcell_data_list)
		    # remove created spaces from the list
                    subcell_data_list = [subcell_word for subcell_word in subcell_data_list if subcell_word != ' ']
		    # merge and convert list items into a single string
                    subcell_data_string = ', '.join(subcell_data_list).strip()
		    # if there is ',' at the end of string then remove it
                    if subcell_data_string[-1] == ',':
                        subcell_data_string = subcell_data_string[:-1]
                    print('\nSubcellular location string: ', subcell_data_string)
                    # Save subcellular location of given motif in certain protein
                    dicty[start].append(subcell_data_string)
                # if no uniprot subcellular location then save '-' into dicty
                elif not important_protein_topology[4]:
                    dicty[start].append('-')

                # EXTRACT UNIPROT TOPOLOGY OF GIVEN PROTEIN
                # if 3rd element is not empty string, then there is topology, thus extract it
                if important_protein_topology[2]:
                    # Create a variable to be able to decide was it successfull to get the topology of given motif in protein by using uniprot retrieved data (see later in use within if condition out of for loop)
                    success = 0
                    # select topological data
                    for protein_topology in important_protein_topology[2:4]:
                        # extract positions from protein_uniprot_topology string: find positions of all '..' since protein parts' start and end values separated by '..' ex.: 133..153
                        matches = re.finditer('\.\.', protein_topology)
                        matches_pos = [match.start() for match in matches]
                        #print(matches_pos)
                        # loop through all matches and generate protein_part_range_set until no overlap with motif_range_set or till the last match
                        for period_pos in matches_pos:
                            # if the uniprot range is wrong like ?..2001 then skip it since it is not known where the given category starts or ends
                            print(protein_topology[period_pos])
                            if protein_topology[protein_topology.rfind(' ', 0, period_pos) + 1:period_pos] == '?':
                                continue
                            #print(matches_pos)
                            #print(protein_uniprot_topology.rfind(' ', 0, period_pos))
                            protein_part_start = int(protein_topology[protein_topology.rfind(' ', 0, period_pos) + 1:period_pos])
                            # print('Protein part start: ', protein_part_start)
                            protein_part_end = int(protein_topology[period_pos + 2:protein_topology.find(';', period_pos)]) + 1
                            # print('Protein part end: ', protein_part_end)
                            protein_part_range_set = set(range(protein_part_start, protein_part_end))
                            # print(protein_part_range_set)
                            # if there is intersection between motif_range_set and protein_part_range_set then extract protein part topology
                            if motif_range_set.intersection(protein_part_range_set):
                                # find position of /note=... because protein part topology is after
                                protein_part_topology_start = protein_topology.find('/note=', period_pos) + 7
                                protein_part_topology_end = protein_topology.find('";', protein_part_topology_start)
                                protein_part_topology = protein_topology[protein_part_topology_start:protein_part_topology_end]
                                # print(protein_part_topology)
                                # save protein_part_topology in dicty for given motif and end for loop because there are not intersections
                                print(protein_part_topology)
                                dicty[start].append(protein_part_topology.upper())
                                # Set success variable to 1 representing success = there was data about given motif's topology in given uniprot retrieved data
                                success = 1
                                break
                    # if there was not intersection between motif_range_set and protein_part_range_set then extract subcellular location of given protein
                    if success == 0:
                        # add '-' into dicty as representing no uniprot topology
                        dicty[start].append('-')
                # if no topology for given protein in uniprot
                elif not important_protein_topology[2]:
                    # add '-' into dicty as representing no uniprot topology
                    dicty[start].append('-')

                # Create response (accepted or rejected) based on these data
                good_topology_list = ['CYTOPLASMIC', 'NUCLEAR', 'SIGNAL CYTOPLASMIC']
                good_subcell_loc_list = ['Nucleus', 'Cytoplasm', 'Cytoplasmic side']
                # 1. Try to use uniprot topology
                if dicty[start][5] != '-':
                    print('\nRun 1')
                    # if it is 'Cytoplasmic' then accept it else reject it
                    dicty[start].append('ACCEPTED') if dicty[start][5] in good_topology_list else dicty[start].append('REJECTED')
                # 2. if transmembrane region in Phobius output, then consider Phobius output and believe it
                elif 'TRANSMEM' in decision_help_list:
                    print('\nRun 2')
                    dicty[start].append('ACCEPTED') if dicty[start][3].strip('.') in good_topology_list else dicty[start].append('REJECTED')
                # 3. try to use uniprot subcellular location
                elif dicty[start][4] != '-':
                    print('\nRun 3')
                    #use subcell good word list and accept it if any word is presented in subcell data
                    res = sum([1 if dicty[start][4].find(word) != -1 else 0 for word in good_subcell_loc_list])
                    if dicty[start][4].find('Nucleus membrane') != -1:
                        res -= 1
                    dicty[start].append('ACCEPTED') if res != 0 else dicty[start].append('REJECTED')
                # 4. try to use Phobius output
                elif dicty[start][3] != '-':
                    print('\nRun 4')
                    # if there is SIGNAL and no transmembrane in Phobius output, then reject it
                    if 'SIGNAL' in decision_help_list:
                        dicty[start].append('REJECTED')
                    # if not SIGNAL and no transmembrane then accept it
                    else:
                        dicty[start].append('ACCEPTED')
    print(dicty)
    return dicty

# ----------------------------------------------------------------------------------------------------------------------

# Retrieve SlimPrints information of given protein to observe island-like conservation and create a file for it
def slimprints_retrieve(acc, output_file_path):
    # retrieve data of given protein
    url = "http://bioware.ucd.ie/rest/slimprints?uniprotid={}".format(acc)
    r = requests.get(url)
    # encode it to change byte code ('u')
    slimprint_output = r.content.decode('ascii')
    print('\n', acc)
    print(slimprint_output)
    # create outputfile from retrieved data
    outf = open(output_file_path, 'w')
    outf.writelines(slimprint_output)

# ----------------------------------------------------------------------------------------------------------------------

# Integrate gopher result in such way calculate Shannon value and regex values (see code below) and return with dictionary: {start : [motif, iupred, pfam, phobius, 'mammalia', shannon_val, regex_val, regex_fraction, 'vertebrata', shannon_val, regex_val, regex_fraction, 'Eumetazoa, shannon_val, regex_val, regex_fraction, 'Opisthokonta', shannon_val, regex_val, regex_fraction, 'Plants', shannon_val, regex_val, regex_fraction, 'Eukaryota', shannon_val, regex_val, regex_fraction]}
# Moreover integrate slimprints result to see islands of conservation
def gopher_slimprint_res_integrate(dicty, slimprint_result_path, gopher_ordered_result_path, ordering_base_path, sorted_gapped_mafftalign_output_path, accepted_gaps):
    # open slimprint result
    with open(slimprint_result_path, 'r') as f:
	# read its lines
        slimprint_lines = f.readlines()
        print('\nRetrieved slimprints data: ', slimprint_lines)
        # loop through motifs
        for start in dicty:
	    # extract motif itself
            motif = dicty[start][0]
	    # create range of motif
            motif_range_set = set(range(int(start)+1, int(start)+len(motif)+1))
            # if there was not data of certain protein in Slimprint server then give 'Unknown', '-' and 'REJECTED' into dictionary
            if len(slimprint_lines) == 1 and slimprint_lines[0].find('<h1>Not Found</h1>') != -1:
                print('\nSlimprints data was not retrieved from server.')
                dicty[start].append('Unknown')
                dicty[start].append('-')
                dicty[start].append('REJECTED')
	    # if there was data about given protein in Slimprint server
            else:
                slimprint_position_list = []
                # split lines and loop through lines from second line since header is not needed
                for line in slimprint_lines[1:]:
                    # split lines to get words and extract start and end position of conservation island detected by SlimPrints and append to slimprint_position_list list
                    split_line = line.split('\t')
                    start_pos = int(split_line[5].split(':')[0])
                    end_pos = int(split_line[5].split(':')[1])+1
                    print('Start and end position: ', start_pos, end_pos)
                    slimprint_position_list.append(set(range(start_pos, end_pos+1)))
                print('\nSlimprint position list: ', slimprint_position_list)
                print('\nMotif range set: ', motif_range_set)
                success = 0
                good_slimprint_list = []
                good_slimprint_score_list = []
		# loop through conservation islands from list
                for line_id, slimprint_range_set in enumerate(slimprint_position_list):
                    # extract score to save it later
                    slimprint_score = float(slimprint_lines[line_id+1].split('\t')[1])
		    # if there is intersection between motif and slimprint islands' range sets and their score is above or equal to a certain threshold
                    if motif_range_set.intersection(slimprint_range_set) and slimprint_score <= slimprint_threshold:
                        # Save slimprint score into a list as a string
                        good_slimprint_score_list.append(str(slimprint_score))
                        print(float(slimprint_lines[line_id + 1].split('\t')[1]))
                        print('\nOverlapt slimprint range set:', slimprint_position_list[line_id])
                        # calculate the fraction of overlap between motif and slimprint island
                        overlap_fraction_str = str(len(motif_range_set.intersection(slimprint_range_set))) + '/' + str(len(slimprint_range_set))
                        # extract SlimPrints' motif regex
                        slimprint_motif = slimprint_lines[line_id+1].split('\t')[2].strip()
                        good_slimprint_list.append(slimprint_motif + ' ; ' + overlap_fraction_str)
			# increment by one the success variable to sign the successfull overlap between motif and one of the slimprint islands
                        success += 1
                # if there was not overlap then add 'Not found', '-' and 'REJECTED' to corresponding key in dictionary
		if success == 0:
                    dicty[start].append('Not found')
                    dicty[start].append('-')
                    dicty[start].append('REJECTED')
                # if there was overlap then add conservation islands' regexes and their scores and 'ACCEPTED' to the dictionary
		elif success != 0:
                    dicty[start].append(' | '.join(good_slimprint_list))
                    dicty[start].append(' | '.join(good_slimprint_score_list))
                    dicty[start].append('ACCEPTED') if sum([1 if float(good_score) <= slimprint_threshold else 0 for good_score in good_slimprint_score_list]) != 0 else dicty[start].append('REJECTED')

    # open ordering_base_path to create a dictionary which fasta files of multi fast mafftalign result are sorted based on: Mammalia -> Vertebrata -> Eumetazoa -> Opisthokonta -> Plant -> Eukaryota
    QFO_wth_phylgen = open(ordering_base_path, 'r')
    # Convert the taxonomy id and phylum into dictionary: {id: phylum}
    lines = QFO_wth_phylgen.readlines()
    QFO_phylgen_dict = {}
    for line in lines:
        id = line.split('\t')[0]
        phylum = line.split('\t')[2]
        QFO_phylgen_dict.update({id: phylum})
    #print(QFO_phylgen_dict)

    # open corresponding protein's ordered gopher result
    with open(gopher_ordered_result_path, 'r') as f:
        file_cont = f.readlines()
        #print(file_cont)
    # loop through motifs
    for start in dicty:
        # create evol_consv_dicty to save the motifs later where human motif's position is in ortholog proteins; goal: {evolutionary_level : [motif1, motif2...]...}
        evol_consv_dicty = {evol_level : [] for evol_level in set(QFO_phylgen_dict.values())}
	# extract motif itself
        motif = dicty[start][0]
        #print('Current motif: ', motif)
        # extract motif position from human sequence which is the first sequence, but I have to refer to 2nd item in list because the first one is the header in the file
        motif_in_seq_start = file_cont[1].find(motif)
        # If motif exist in it without gaps then print it out and use the original sorted output instead of gapped ones
        if motif_in_seq_start != -1:
            #print('Found motif in first sequence: ', file_cont[1][motif_in_seq_start:motif_in_seq_start+len(motif)])
            # loop through each odd positioned lines and extract peptide sequence from same position: start from 3rd line id since human is not counted in evolution conservativity calculation
            for line_id in range(3, len(file_cont), 2):
                # extract ox (organism ID) from header (previous line) to be able to search in previously created QFO_phylgen_dict which says the type of organism
                ox = file_cont[line_id-1][file_cont[line_id-1].find('OX=')+3:file_cont[line_id-1].find(' ', file_cont[line_id-1].find('OX=')+3)]
                # search the evolutionary level of given organism's proteins
                current_evol_level = QFO_phylgen_dict[ox]
                # add motif to list of evol_consv_dicty
                evol_consv_dicty[current_evol_level].append(file_cont[line_id][motif_in_seq_start:motif_in_seq_start+len(motif)])
                #print('New line: ', file_cont[line_id][motif_in_seq_start:motif_in_seq_start+len(motif)])
            print(evol_consv_dicty)
        # if motif contains gap, thus not possible to find it, then removes gaps and search motif position in human sequence then as did it previously, add peptides of ortholog proteins in same position to evol_consv_dicty
        elif motif_in_seq_start == -1:
            print('There are gaps in reference motif!')
            # open gapped version of mafftalign result
            with open(sorted_gapped_mafftalign_output_path, 'r') as f:
                file_gaped_cont = f.readlines()
                # search motif in human sequence of ungapped mafftalign result
                motif_in_seq_start = file_gaped_cont[1].find(motif)
                print('Found motif in first sequence: ', file_gaped_cont[1][motif_in_seq_start:motif_in_seq_start + len(motif)])
                # loop through each odd positioned lines and extract peptide from same position: start from 3rd line id since human is not counted in evolution conservativity calculation
                for line_id in range(3, len(file_cont), 2):
                    # extract ox (organism ID) from header (previous line) to be able to search in previously created QFO_phylgen_dict which says the type of organism
                    ox = file_gaped_cont[line_id - 1][file_gaped_cont[line_id - 1].find('OX=') + 3:file_gaped_cont[line_id - 1].find(' ', file_gaped_cont[line_id - 1].find('OX=') + 3)]
                    # search the evolutionary level of given organism's proteins
                    current_evol_level = QFO_phylgen_dict[ox]
                    # add motif to list of evol_consv_dicty
                    evol_consv_dicty[current_evol_level].append(file_gaped_cont[line_id][motif_in_seq_start:motif_in_seq_start + len(motif)])
                    #print('New line: ', file_cont[line_id][motif_in_seq_start:motif_in_seq_start + len(motif)])
                print(evol_consv_dicty)


        # Calculate evolutionary conservativity by using Shannon index and regex separately; Shannon entropy: -sum(n1/N*ln(n1/N))
        # loop through evolutionary levels
        for evol_level in ['Mammalia', 'Vertebrata', 'Eumetazoa', 'Opisthokonta', 'Plant', 'Eukaryota']:
            print(evol_level)
            # if no orthologs for given evol_level then give 0 values as converted average Shannon conservativity, 0 and 0/0 rate as binary evol_conservativity
            if not evol_consv_dicty[evol_level]:
                #print('Yes')
                dicty[start].append(evol_level)
                dicty[start].append('0/0')
                dicty[start].append('-')
                dicty[start].append('-')
                dicty[start].append('-')
            # if there are orthologs for given evol_level then calculate Shannon index with those sequences that not have gaps
            elif evol_consv_dicty[evol_level]:
                # Add evolutionary level's name to dicty
                dicty[start].append(evol_level)
                # Check for gaps and 'X' amino acids (see section below till the end of elif part)
		removed_seq = 0  # to calculate the number of rejected motif candidates
		# count number of all motif candidates that are in same position than human one through species in actual taxonomical level
                all_seq = len(evol_consv_dicty[evol_level])
                print('\nAll ortholog sequences on evol level: ', all_seq)
                current_evol_level_list = copy.deepcopy(evol_consv_dicty[evol_level])
                # loop through sequences in given evolutionary level
		for seq in current_evol_level_list:
                    print(seq)
                    print('Gap number in {}: {}'.format(seq, seq.count('-')))
		    # if number of gaps is greater than the allowed number of gaps
                    if seq.count('-') > accepted_gaps:
			# remove given motif candidate from evol_consv_dicty
                        evol_consv_dicty[evol_level].remove(seq)
			# increase number of removed sequence variable by one
                        removed_seq += 1
		    # if there is any unknown or special amino acids in given motif candidate
                    elif seq.count('X') > 0:
			# remove given motif candidate from evol_consv_dicty
                        evol_consv_dicty[evol_level].remove(seq)
			# increase number of removed sequence variable by one
                        removed_seq += 1
                # after the removal of motif candidates that contain any unknown or special amino acid or more gaps than allowed 
		# add a fraction into dictionary to represent the number of removed candidates according to all: removed candidates number / all candidates number
		dicty[start].append(str(removed_seq) + '/' + str(all_seq))
                print('Gap-deleted evol level list: ', evol_consv_dicty[evol_level])
		#calculate the number of remained sequences
                survived_seq = len(evol_consv_dicty[evol_level])
                print('\nNumber of survived seq on evol level: ', survived_seq)
                # if more than half of sequences were remained and there are at least 2 survived sequences then calculate Shannon otherwise Shannon is useless
                if survived_seq/all_seq >= 0.5 and survived_seq >= 2:
                    print('Accepted')
                    # create shannon_pos_list to save shannon values of each position since it is calculated for all position in given motif candidate
                    shannon_pos_list = []
                    # go through each position of motif
                    for motif_pos in range(len(motif)):
                        # create amino acid dictionary for counting frequency of amino acid in each position
                        aa_dicty = {aa: 0 for aa in ['A', 'R', 'N', 'D', 'C', 'Q', 'G', 'E', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']}
                        print(motif_pos)
                        # loop through motifs of given evolutionary level and calculate frequency of amino acids
                        for actual_motif in evol_consv_dicty[evol_level]:
                            aa_dicty[actual_motif[motif_pos]] += 1
                        print(aa_dicty)
                        # Calculate Shannon index for each keys in aa_dicty, then take absolute values of sum of shannon values, but it is subtracted from 1 to reverse logic of its index such way: larger values mean higher conservativity. Originally smaller one means higher conservativity. (see below)
                        print('Separated shannon prevalues: ', [float(aa_dicty[aa])/ float(sum(aa_dicty.values())) * (np.log2(float(aa_dicty[aa])/ float(sum(aa_dicty.values())))) for aa in aa_dicty if aa_dicty[aa] > 0])
                        shannon_val = abs(sum([float(aa_dicty[aa])/ float(sum(aa_dicty.values())) * (np.log2(float(aa_dicty[aa])/ float(sum(aa_dicty.values())))) for aa in aa_dicty if aa_dicty[aa] > 0]))
                        print('\nShannon value: ', shannon_val)
                        # normalization by division of log2(min(N, K)) -> N = number of residues in given column, K = number of  symbol types
                        shannon_val_norm = shannon_val / (np.log2(min(survived_seq, len(aa_dicty.keys()))))
                        print('\nNormalized Shannon value: ', shannon_val_norm)
                        # Reverse Shannon entropy logic: higher score means more conservativity
                        shannon_val_norm_conv = 1 - shannon_val_norm
                        print('\nConverted normalized shannon value: ', shannon_val_norm_conv)
	        # add score into list shannon_pos_list
                        shannon_pos_list.append(shannon_val_norm_conv)
                    # take average of shannon values in each position
                    print('Average of Shannon: ', np.mean(shannon_pos_list))
                    dicty[start].append(np.mean(shannon_pos_list))
		# if less than half of sequences were remained and there are not at least 2 survived sequences then add '-' into dictionary
                elif survived_seq/all_seq < 0.5 or survived_seq < 2:
                    dicty[start].append('-')

                # Calculate evolutionary conservativity in binary way via regex comparison for given motif; only consider previously accepted and kept motifs (not that has either special or unknown amino acid or more gaps than allowed)
                # If regex matches to given motif then 1 otherwise 0 is added to list
                regex_motif_comparison_list = [1 if re.match(used_regex, motif) else 0 for motif in current_evol_level_list]
                # Make fraction where numerator is number of matched motifs while denominator presents number of all motifs in given evolutionary level
                regex_evol_conserv_value = sum(regex_motif_comparison_list)/len(regex_motif_comparison_list)
                regex_evol_conserv_fraction = str(sum(regex_motif_comparison_list)) + '/' + str(len(regex_motif_comparison_list))
                # Add above fraction and percentage of it into dictionary
                dicty[start].append(regex_evol_conserv_value)
                dicty[start].append(regex_evol_conserv_fraction)

    print(dicty)
    return dicty

# ----------------------------------------------------------------------------------------------------------------------

# Create a file for given protein within there is all integrated data: motif start, motif, iupred, pfam, phobius, evolutionary level, corresponding 'shannon value', corresponding regex matching percentage and fraction
# inputs: dictionary {'start' : [motif, iupred_score,]}
def info_integr_output_generate(dicty, output_file_path):
    # create output file content
    output_file_cont_list = []
    # loop through dicty to extract all information
    for start in dicty:
        # add motif start position to output list
        output_file_cont_list.append(str(int(start) +1) + '\t')   # add 1 to start position because start positions were saves in a python-preferred way (0 based indexing)
        # loop through all data of given motif and add them to output list with tabulator separation
        for data in dicty[start]:
            output_file_cont_list.append(str(data) + '\t')
        # add a new line character at the end of a given motif's information containing line
        output_file_cont_list.append('\n')
    print(output_file_cont_list)
    # open output file and write content of output list into it
    with open(output_file_path, 'w') as f:
        f.writelines(output_file_cont_list)
    print('File was successfully created!')
    return 1

# ----------------------------------------------------------------------------------------------------------------------

##########################################################################################################
How to run filtering programs
##########################################################################################################

# Read given protein's motif_regex_search result file and create a dictionary: {start : [motif], start : [motif]}
#protein_motif_dicty = motif_search_result_read_dict_create(motif_search_res_dir + protein_id + '.res')

# Integrate of all Pfam, Iupred, Phobius data of motifs in given protein and create a dictionary: {start : [motif, iupred, pfam, phobius]}
#motif_dicty = integrate_all(protein_motif_dicty, iupred_analysis_dir + protein_id + '.iupred', pfam_analysis_dir + protein_id + '_pfam_res.txt', phobius_analysis_dir + protein_id + '.phb_output')

# Integrate gopher result in such way calculate Shannon value and regex values (see code above) and return with dictionary: {start : [motif, iupred, pfam, phobius, 'mammalia', shannon_val, regex_val, regex_fraction, 'vertebrata', shannon_val, regex_val, regex_fraction, 'Eumetazoa, shannon_val, regex_val, regex_fraction, 'Opisthokonta', shannon_val, regex_val, regex_fraction, 'Plants', shannon_val, regex_val, regex_fraction, 'Eukaryota', shannon_val, regex_val, regex_fraction]}
#motif_dicty = gopher_slimprint_res_integrate(motif_dicty, slimprints_data + protein_id + '.slimprints', gopher_analysis_dir + 'Sorted_mafft_align/' + protein_id + '.orthaln.fas', eukaryota_taxonomy_list, gopher_analysis_dir + 'Gapped_sorted_mafftalign/' + protein_id + '.orthaln.fas', 0)

# Use the filter to add 'ACCEPTED' or 'REJECTED' as a new line
#motif_filter(motif_dicty, iupred_cutoff, pfam_limit_list)

# Create a file for given protein within there is all integrated data: motif start, motif, iupred, pfam, phobius, evolutionary level, corresponding 'shannon value', corresponding regex matching percentage and fraction
#info_integr_output_generate(motif_dicty, integrated_data_output_dir + protein_id + '.integr.txt')






# Used this to get protein list of motif containing proteins for GNU parallel run this script in terminal
#ls *.res | sed 's/.res//' > protein_id_list

# Used to save terminal outputs in a file when I'm running this script in gnu parallel
# script /home/zolo012/KELCH/Motif_search_v1/Protein_results/Integrated_data/Log_file


##########################################################################################################
Protein-protein interaction (PPI) section
##########################################################################################################

# Open retreived intact PPI data containing file and create a dictionary: {interacting_protein_id : {detection_method : [MI_score, MI_score]}}
# input: file path to PPI data containing file; protein which interaction partners are looked for
# output dictionary: {protein_id : {method : [score, score...]...}...}
def (input_file_path, searched_protein):
    # Create output dictionary to save
    dicty = {}
    # open text file
    with open(input_file_path, 'r') as f:
        # loop through each line
        for line in f.readlines():
            # find those lines within got searched_protein and another protein are interacted
            words_list = line.split('|')
            #print(words_list)
            # extract partner1 and 2 in current line
            partner1 = words_list[0][words_list[0].find(':')+1:words_list[0].find('\t')]
            partner2 = words_list[0][words_list[0].rfind('uniprotkb:')+10:words_list[0].rfind('\t')]
            #print(partner1, partner2)
            # if partner1 is not the same as the other partner
            if partner1 != partner2:
                # select partner which is not main interacting protein (not given searched protein)
                partner = partner1 if partner1 != searched_protein else partner2
                # extract method
                method_id = [word_id for word_id, words in enumerate(words_list) if '"(' in words][0]
                method = words_list[method_id][words_list[method_id].find('"(')+2:words_list[method_id].find(')', words_list[method_id].find('"('))]
                # extract intact MI score
                intact_miscore = [words[words.find('intact-miscore:')+15:words.find('\t', words.find('intact-miscore:'))] for words in words_list if 'intact-miscore:' in words][0]
                #print(partner)
                #print(method)
                #print(intact_miscore)
	# if partner not have key in dictionary then create one
                if partner not in dicty.keys():
                    dicty.update({partner : {}})
	# if method is not among the values of given key in dictionary then add it
                if method not in dicty[partner].keys():
                    dicty[partner].update({method : []})
	# add PPI score to dictionary
                dicty[partner][method].append(intact_miscore)
        #print(dicty)
    return dicty

# ----------------------------------------------------------------------------------------------------------------------

# Open PPI data from Biogrid which was only downloaded by hand from website (https://thebiogrid.org/) and create a dictionary based on that: {interacting_protein_id : {method : occurrence}}}
# input: path to Biogrid PPI data containing file, protein which interaction partners are looked for
# output dictionary: {partner : {method : occurrence_number}}
def biogrid_ppi_dicty_generate(input_file_path, searched_protein):
    dicty = {}
    # open file
    with open(input_file_path, 'r') as f:
        # read file, but start from 2nd line. 1st one is the header
        for line in f.readlines()[1:]:
            #print(line)
            # extract both taxon id to only run for human proteins
            taxid1 = line[line.find('taxid:')+6:line.find('\t', line.find('taxid:'))]
            taxid2 = line[line.rfind('taxid:')+6:line.find('\t', line.rfind('taxid:'))]
            #print(taxid1, taxid2)
            # if both proteins are human ones
            if taxid1 == '9606' and taxid2 == '9606':
                # extract partner1 and partner2
                partner1 = line[line.find('uniprot/swiss-prot:')+19:line.find('|', line.find('uniprot/swiss-prot:'))]
                partner2 = line[line.rfind('uniprot/swiss-prot:')+19:line.find('|', line.rfind('uniprot/swiss-prot:'))]
                #print(partner1, partner2)
                # if partners are not same
                if partner1 != partner2:
                    # extract not searched partner
                    partner = partner1 if partner1 != searched_protein else partner2
                    #print(partner)
                    # extract method
                    method = line[line.find('"(')+2:line.find(')', line.find('"('))]
                    #print(method)
                    # save everything into dicty
                    if partner not in dicty.keys():
                        dicty[partner] = {}
                    if method not in dicty[partner].keys():
                        dicty[partner][method] = 1
                    elif method in dicty[partner].keys():
                        dicty[partner][method] += 1
        #print(dicty)
    return dicty

# ----------------------------------------------------------------------------------------------------------------------

# Retrieves PPI data of given protein from Intact database: https://www.ebi.ac.uk/intact/home into an outputdir | accesion_dicty ex.: {KLHL12 : ['Q53G59']}
def retrieve_data_intact(accession_dicty, outputdir):
    # go across a dictionary
    for name, accs in accession_dicty.items():
        #print(name, end="\t")
        # go through protein accession ids
        for acc in accs:
            #set and assemble URL to retrieve the data of actual protein
            url = "http://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/interactor/{0}?format=tab27".format(acc)
            r = requests.get(url)
            # write result into a file which place into given output directory
            with open(outputdir + acc, 'w') as f:
                f.write(r.content.decode('utf-8'))
    return 1

# ----------------------------------------------------------------------------------------------------------------------

# Read integrated output files into dictionary and complete it with gained PPI data
# input dictionary: {protein_id : {start : [motif, iupred, pfam, phobius, evolutionary level, corresponding 'shannon value', corresponding regex matching percentage and fraction]}
# output dictionary: {protein_id : {start : [motif, iupred, pfam, phobius, evolutionary level, corresponding 'shannon value', corresponding regex matching percentage and fraction, PPI_method, PPI score]}} 
def integrate_data_dicty_generate(integrated_proteins_data_dir):
    # Create dictionary to save data
    dicty = {protein_id : {'start' : []}}
    # loop through proteins in given directory
    for file_name in os.listdir(integrated_proteins_data_dir):
        if file_name.find('.integr.txt') != -1:
            # create protein id from filename
            protein_id = file_name[:file_name.find('.')]
            # add protein id into dicty
            dicty.update({protein_id : {}})
            # open given file
            with open(integrated_proteins_data_dir + file_name, 'r') as f:
                # Read file
                for line in f.readlines():
                    line_list = line.strip().split('\t')
                    # extract start and save it as 2nd keys
                    start = line_list[0]
                    dicty[protein_id].update({start : []})
                    # Fill dictionary up
                    dicty[protein_id][start] = line_list[1:]
    #print(dicty)
    return dicty

# ----------------------------------------------------------------------------------------------------------------------


# Use retrieved PPI data and add method and PPI those proteins that contain any putative KELCH motif(s)
# inputs: dictionary generated by integrate_data_dicty_generate algorithm, sequence library/dictionary created by impr_multi_FASTA_seq_extractor, dictionary created by intact_ppi_dicty_create, dictionary generated by biogrid_ppi_dicty_generate, path to where output file generate, list with the names of columns (header)
def intact_biogrid_ppi_add(protein_dicty, sequence_lib_dicty, intact_data_dicty, biogrid_data_dicty, output_file_path, header_names):
    # loop through protein_dicty
    for protein_id in protein_dicty:
        for start in protein_dicty[protein_id]:
            # check given protein id in the intact_data_dicty
            if protein_id in intact_data_dicty.keys():
                # add method and MI score to list of given protein in protein_dicty as lists
                method_list, mi_score_list = [], []
                for method, mi_score in intact_data_dicty[protein_id].items():
                    method_list.append(method)
                    mi_score_list += mi_score
                protein_dicty[protein_id][start].append(method_list)
                protein_dicty[protein_id][start].append(mi_score_list)
            # if given protein id not in intact_data_dicty check biogrid_data_dicty
            elif protein_id in biogrid_data_dicty.keys():
                for method, occurence in biogrid_data_dicty[protein_id].items():
                    method_list.append(method)
                    mi_score_list.append(str(occurence))
                protein_dicty[protein_id][start].append(method_list)
                protein_dicty[protein_id][start].append(mi_score_list)
            # if given protein id not in intact_data_dicty and not in biogrid_data_dicty then give '-' '-' as no method and no score
            else:
                protein_dicty[protein_id][start].append('-')
                protein_dicty[protein_id][start].append('-')
    #print(protein_dicty)

    # Create final result file with all result and a filtered version
    output_file_cont = []
    # add header with new line character at its end to the list
    output_file_cont.append('\t'.join(header_names_list) + '\n')
    # loop through proteins in dictionary protein_dicty
    for protein_id in protein_dicty:
	# extract corresponding gene name of protein
        gene_id = sequence_lib_dicty[protein_id]['GN']
	# loop through motifs in dictionary protein_dicty
        for start in protein_dicty[protein_id]:
	    # create output_file_string and add protein name, gene name and motif's start position with tabulator separation
            output_file_string = protein_id + '\t' + gene_id + '\t' + start + '\t'
	    # if there are multiple PPI data for a given motif containing protein then add all of them to the output_file_string with tabulator separation
            if type(protein_dicty[protein_id][start][-1]) is list:
               # add all information except PPI data of given protein since in this case multiple PPI data are stored in lists
	       for word in protein_dicty[protein_id][start][:-2]:
                   output_file_string += word + '\t'
	       # add PPI method data with space separation
               for word in protein_dicty[protein_id][start][-2]:
                   #print(protein_dicty[protein_id][start][:-2])
                   output_file_string += word + ' '
	       # add a tabulator then
               output_file_string += '\t'
	       # add PPI scores with space separation
               for word in protein_dicty[protein_id][start][-1]:
                   output_file_string += word + ' '
            # if there are not multiple PPI data for a given protein then add it to the output_file_string with tabulator separation
	    elif not type(protein_dicty[protein_id][start][-1]) is list:
                for word in protein_dicty[protein_id][start][:-1]:
                    output_file_string += word + '\t'
                output_file_string += protein_dicty[protein_id][start][-1]
            # add new line character
	    output_file_string += '\n'
            output_file_cont.append(output_file_string)
    print(output_file_cont)
    # open output file and write information into it
    with open(output_file_path, 'w') as f:
        f.writelines(output_file_cont)

    # Save protein dictionary with Pickle
    print('Final file was successfully created!')
    outf = open(output_file_path + '_dictionary', 'wb')
    pickle.dump(protein_dicty, outf)
    outf.close()

    return 1

# ----------------------------------------------------------------------------------------------------------------------




