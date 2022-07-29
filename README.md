# Motif-search-pipeline

In this repository there are functions that can find putative motifs in proteins by using regular expressions (RegEx) or position-specific scoring matrices (PSSM) generated from already knowns motifs in interaction partners of a given protein. Then these candidates are filtered based on IuPred, Pfam, location, evolution and protein-protein interaction (PPI) data. 

1. Final_motif_search.py file contains all functions that I have written during motif search project.
2. Analysis_of_proteins.py file makes all analysis for proteins. Analysis includes Pfam, IUPred, Phobius, SingalP, evolutional analysis by Gopher. There is also the function that looks for moitf candidates in proteins based on regular expression (regex).
3. Motif_filtering.py file read the motif candidates with their corresponding information and filter out them based on set thresholds for IUPred, Pfam, evolutional filters.
4. Slimprints_retrieve.py create distinct SLiMPrints output files for given protein. (this is another evolutional filter to make a high confidence list)
5. PPI_data_retrieve.py returns with protein-protein interaction (PPI) data from Intact and Biogrid databases, then it applies two additional filters on motif candidates that survived aforementioned criteria (see point 3.) to remove unlikely hits. Two additional filters are based on SLiMPrints results and PPI data.
