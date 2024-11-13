import oloutils as ut
import re

import sys
import argparse

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="")

    parser.add_argument('-i', '--input',
                        dest= "infile",
                        action= "store",
                        default= "./",
                        help="")
    
    parser.add_argument('-o', '--ouput',
                        dest= "outfile",
                        action= "store",
                        default= "./",
                        help="")
    
    parser.add_argument('-rp', '--remove_pathways',
                        dest= "remove_pathways",
                        action= "store",
                        default= "./",
                        help="")
    
    options = parser.parse_args()

    if len(sys.argv) <= 3:
        sys.stderr.write("No input, output or pathways to remove provided. Please, try again. \n")
        exit()

    pathways = {}

    f = open(options.infile)
    sys.stderr.write("Loaded file: %s \n" % options.infile)

    pathw_remove = ut.load_json(options.remove_pathways, compressed=True)
    sys.stderr.write("Loaded file: %s \n" % options.remove_pathways)

    sys.stderr.write("Parsing files... \n")

    for line in f:
        line = line.strip()
        line = line.split("\t")

        organism = line[5]    
        chebid = line[0]
        pathid = line[1]

        species = ['Homo sapiens', 'Mus musculus']
        if organism in species:            
            if organism not in pathways.keys():
                pathways[organism] = {}
            
            if chebid not in pathways[organism].keys():
                pathways[organism][chebid] = []
            
            if pathid not in pathw_remove[organism]["difference"]:
                pathways[organism][chebid].append(pathid)
    
    for org in pathways.keys():
        sys.stderr.write("%d molecules present in %s pathways \n" %
                         (len(pathways[org]), org))
        
    ut.save_json(pathways, options.outfile, compressed=True)
    sys.stderr.write("File saved as %s \n" % options.outfile)
