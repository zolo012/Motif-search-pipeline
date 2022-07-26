from Motif_search import *
import requests
import pickle
from bs4 import BeautifulSoup
import pandas as pd
import os

# Place variables
human_proteom_place = '/home/zolo012/Human_proteom/human_proteom_no_isoforms_dictionary'
human_proteom_raw_place = '/home/zolo012/Human_proteom/human_proteom_no_isoforms.fasta'



# A, Run all analysis (Pfam, Iupred, Phobius, SignalP) for each human protein

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
human_proteom_iupred_dicty = iupred_dicty_generate(human_proteom_dicty, '/home/zolo012/Human_proteom/Iupred_analysis/')

# Save human_proteom_iupred_dicty
outfile = open('/home/zolo012/Human_proteom/Iupred_analysis/human_proteom_iupred_dictionary', 'wb')
pickle.dump(human_proteom_iupred_dicty, outfile)


# Open human_proteom_iupred_dicty
infile = open('/home/zolo012/Human_proteom/Iupred_analysis/human_proteom_iupred_dictionary', 'rb')
human_proteom_iupred_dicty = pickle.load(infile)


# 3. Run Pfam analysis for all proteins in human proteom
# 3.a: create fasta files for all proteins 
pfam_input = '/home/zolo012/Human_proteom/Pfam_analysis/Input_files/'
for id in human_proteom_dicty:
  fasta_creator(id, human_proteom_raw_place, pfam_input, human_proteom_dicty)


# 3.b: do pfam analysis
pfam_path = '/home/zolo012/Pfam/PfamScan/pfam_scan.pl'
pfam_output = '/home/zolo012/Human_proteom/Pfam_analysis/Pfam_output/'
pfam_use(pfam_path, pfam_input, pfam_output)

# 4. Run SingalP for all proteins in human proteom
signalp_output_dir = '/home/zolo012/Human_proteom/SignalP_analysis/'
signalp_run(human_proteom_dicty, pfam_input, signalp_output_dir)

# 5. Run phobius for all proteins in human proteom
phobius_output_dir = '/home/zolo012/Human_proteom/Phobius_analysis/'
phobius_run(human_proteom_dicty, phobius_path, pfam_input, )
