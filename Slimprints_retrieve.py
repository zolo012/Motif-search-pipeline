from Final_motif_search import *
import requests
import sys

slimprints_data = '/home/zolo012/KELCH/Motif_search_v1/Slimprints/'

protein_id = sys.argv[1]

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

slimprints_retrieve(protein_id, slimprints_data+protein_id + '.slimprints')
