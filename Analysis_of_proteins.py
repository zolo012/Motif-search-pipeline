from Motif_search import *
import requests
import pickle
from bs4 import BeautifulSoup
import pandas as pd
import os

# Place variables
human_proteom_place = '/home/zolo012/Human_proteom/human_proteom_no_isoforms_dictionary'
human_proteom_raw_place = '/home/zolo012/Human_proteom/human_proteom_no_isoforms.fasta'
iupred_result_dir = '/home/zolo012/Human_proteom/Iupred_analysis/'
pfam_input = '/home/zolo012/Human_proteom/Pfam_analysis/Input_files/'
pfam_path = '/home/zolo012/Pfam/PfamScan/pfam_scan.pl'
pfam_output = '/home/zolo012/Human_proteom/Pfam_analysis/Pfam_output/'
signalp_output_dir = '/home/zolo012/Human_proteom/SignalP_analysis/'
phobius_output_dir = '/home/zolo012/Human_proteom/Phobius_analysis/'
QfO_proteom_48_filepath = '/home/zolo012/QfO_data/Eukaryota/eukaryota_48_proteomes.fasta'
gopher_path = '/home/zolo012/Slimsuite/SLiMSuite-master/tools/gopher.py'
gopher_output_dir = '/home/zolo012/Python/Gopher_results'
QFO_eukaryota_Taxonomy_list = '/home/zolo012/QfO_data/QFO_eukaryota_Taxonomy.list'
mafftalign_base_path = '/home/zolo012/Human_proteom/Gopher_results/Supplement/eukaryota_48_proteomes/HUMAN/mafftALN/'
sorted_mafftalign_base_path = '/home/zolo012/Human_proteom/Gopher_results/Supplement/Sorted_mafft_align/'
sorted_gapped_mafft_dir = '/home/zolo012/Human_proteom/Gopher_results/Supplement/Gapped_sorted_mafft_align/'


#################################################################
# A, Run all analysis (Pfam, Iupred, Phobius, SignalP) for each human protein
#################################################################

# 1. Open human proteom without isoforms file and create a variable to be able to use and save it too
human_proteom_dicty = impr_multi_FASTA_seq_extractor(human_proteom_raw_place)

# Save human proteom
outfile = open(human_proteom_place, 'wb')
pickle.dump(human_proteom_dicty, outfile)
outfile.close()


# Open human proteom
infile = open(human_proteom_place, 'rb')
human_proteom_dicty = pickle.load(infile)

print('\nThere are {} human proteins in human proteom!'.format(len(human_proteom_dicty.keys())))

# 2. Run iupred for all proteins in human proteom -- there's a warning:   raise ValueError("If mode is 'interp', window_length must be less "
# ValueError: If mode is 'interp', window_length must be less than or equal to the size of x.
human_proteom_iupred_dicty = iupred_dicty_generate(human_proteom_dicty, iupred_result_dir)

# Save human_proteom_iupred_dicty
outfile = open('/home/zolo012/Human_proteom/Iupred_analysis/human_proteom_iupred_dictionary', 'wb')
pickle.dump(human_proteom_iupred_dicty, outfile)


# Open human_proteom_iupred_dicty
infile = open('/home/zolo012/Human_proteom/Iupred_analysis/human_proteom_iupred_dictionary', 'rb')
human_proteom_iupred_dicty = pickle.load(infile)


# 3. Run Pfam analysis for all proteins in human proteom
# 3.a: create fasta files for all proteins 
for id in human_proteom_dicty:
  fasta_creator(id, human_proteom_raw_place, pfam_input, human_proteom_dicty)


# 3.b: do pfam analysis
pfam_use(pfam_path, pfam_input, pfam_output)

# 4. Run SingalP for all proteins in human proteom
signalp_run(human_proteom_dicty, pfam_input, signalp_output_dir)

# 5. Run phobius for all proteins in human proteom
phobius_run(human_proteom_dicty, phobius_path, pfam_input, )

#################################################################
# B, Evolution of proteins
#################################################################

# 1, Prepare QFO database for Gopher
# 1.a.: Create QfO dictionary
QfO_proteom_48_dict = impr_multi_FASTA_seq_extractor(QfO_proteom_48_filepath)
print(QfO_proteom_48_dict, '\nNumber of proteins {:,}'.format(len(QfO_proteom_48_dict)))

# 1.b: See the organisms
QfO_organisms_set = set([inf['OS'] for inf in QfO_proteom_48_dict.values()])
#print('\nOrganisms in QfO data', QfO_organisms_set)

# 2, run Gopher for all human proteins
gopher_dir_run(pfam_input, gopher_path, gopher_output_dir, QfO_proteom_48_filepath)

# 3. Order ortholog proteins from human to lower phylogenetic levels in multifasta files
# Open the table from QFO_eukaryota_taxonomy list
QFO_wth_phylgen = open(QFO_eukaryota_Taxonomy_list, 'r')

# Convert the taxonomy id and phylum into dictionary: {id : phylum}
lines = QFO_wth_phylgen.readlines()
QFO_phylgen_dict = {}
for line in lines:
    id = line.split('\t')[0]
    phylum = line.split('\t')[2]
    QFO_phylgen_dict.update({id : phylum})
QFO_wth_phylgen.close()    
print('\nQFO phylgen dict:', QFO_phylgen_dict)

# Loop through alignment files
for alignm_file_name in os.listdir(mafftalign_base_path):
    alignm_file_path = mafftalign_base_path + alignm_file_name
    # Merge coherent headers and sequences
    merged_alignm_file_cont_list = head_and_neck_merge(alignm_file_path)
    # Sort them and create new files
    output_path = sorted_mafftalign_base_path + alignm_file_name
    multi_fasta_cont_sort_file_creator(merged_alignm_file_cont_list, QFO_phylgen_dict, output_path)

# 4. Generate ungapped version of orthologs of all proteins; ungapped version means removal of gaps in ortholod proteins based on reference protein (itc.: human ones)
gap_kill(sorted_mafftalign_base_path, sorted_gapped_mafft_dir)

#################################################################
# Run motif search in whole human proteom by using regex
#################################################################

motif_dicty, protein_set, motif_set, motif_all_list = motif_regex_search(human_proteom_dicty, 'PG.P[P,E]', 5, '/home/zolo012/KELCH/Motif_search_v1/Protein_results/')

# Save motif things
outfile = open('/home/zolo012/KELCH/Motif_search_v1/motif_dictionary', 'wb')
pickle.dump(motif_dicty, outfile)

outfile = open('/home/zolo012/KELCH/Motif_search_v1/protein_set', 'wb')
pickle.dump(protein_set, outfile)

outfile = open('/home/zolo012/KELCH/Motif_search_v1/motif_set', 'wb')
pickle.dump(motif_set, outfile)

outfile = open('/home/zolo012/KELCH/Motif_search_v1/motif_all_list', 'wb')
pickle.dump(motif_all_list, outfile)








