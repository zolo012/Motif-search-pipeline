from Final_motif_search import *
import sys
import numpy as np
import pickle
import requests
import copy
import os
import re

# Place variables
motif_search_res_dir = '/home/zolo012/KELCH/Motif_search_v1/Protein_results/'
iupred_analysis_dir = '/home/zolo012/Human_proteom/Iupred_analysis/'
pfam_analysis_dir = '/home/zolo012/Human_proteom/Pfam_analysis/Pfam_output/'
phobius_analysis_dir = '/home/zolo012/Human_proteom/Phobius_analysis/'
gopher_analysis_dir = '/home/zolo012/Human_proteom/Gopher_results/'
integrated_data_output_dir = '/home/zolo012/KELCH/Motif_search_v1/Protein_results/Integrated_data/'
eukaryota_taxonomy_list = '/home/zolo012/QfO_data/QFO_eukaryota_Taxonomy.list'
slimprints_data = '/home/zolo012/KELCH/Motif_search_v1/Slimprints/'



# Read given protein's motif_regex_search result file and create a dictionary: {start : [motif], start : [motif]}
protein_motif_dicty = motif_search_result_read_dict_create(motif_search_res_dir + protein_id + '.res')

# Integrate of all Pfam, Iupred, Phobius data of motifs in given protein and create a dictionary: {start : [motif, iupred, pfam, phobius]}
motif_dicty = integrate_all(protein_motif_dicty, iupred_analysis_dir + protein_id + '.iupred', pfam_analysis_dir + protein_id + '_pfam_res.txt', phobius_analysis_dir + protein_id + '.phb_output')

# Integrate gopher result in such way calculate Shannon value and regex values (see code above) and return with dictionary: {start : [motif, iupred, pfam, phobius, 'mammalia', shannon_val, regex_val, regex_fraction, 'vertebrata', shannon_val, regex_val, regex_fraction, 'Eumetazoa, shannon_val, regex_val, regex_fraction, 'Opisthokonta', shannon_val, regex_val, regex_fraction, 'Plants', shannon_val, regex_val, regex_fraction, 'Eukaryota', shannon_val, regex_val, regex_fraction]}
motif_dicty = gopher_slimprint_res_integrate(motif_dicty, slimprints_data + protein_id + '.slimprints', gopher_analysis_dir + 'Sorted_mafft_align/' + protein_id + '.orthaln.fas', eukaryota_taxonomy_list, gopher_analysis_dir + 'Gapped_sorted_mafftalign/' + protein_id + '.orthaln.fas', 0)

# Use the filter to add 'ACCEPTED' or 'REJECTED' as a new line
motif_filter(motif_dicty, iupred_cutoff, pfam_limit_list)

# Create a file for given protein within there is all integrated data: motif start, motif, iupred, pfam, phobius, evolutionary level, corresponding 'shannon value', corresponding regex matching percentage and fraction
info_integr_output_generate(motif_dicty, integrated_data_output_dir + protein_id + '.integr.txt')
