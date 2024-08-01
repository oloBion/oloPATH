import utils as ut
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
    
    options = parser.parse_args()

    if len(sys.argv) <= 2:
        sys.stderr.write("No input or output provided. Please, try again. \n")
        exit()

    pathways = {}

    f = open(options.infile)
    sys.stderr.write("Loaded file: %s \n" % options.infile)

    sys.stderr.write("Parsing file... \n")

    for line in f:
        line = line.strip()
        line = line.split("\t")

        organism = line[2]    
        pathid = line[0]
        pathnm = line[1]
        
        if organism not in pathways.keys():
            pathways[organism] = {}
        
        pathways[organism][pathid] = {'diplay_name': pathnm}
    
    for org in pathways.keys():
        sys.stderr.write("%d pathways found in %s \n" %
                         (len(pathways[org]), org))
        
    ut.save_json(pathways, options.outfile, compressed=False)
    sys.stderr.write("File saved as %s \n" % options.outfile)
