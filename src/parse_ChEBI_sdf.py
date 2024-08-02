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
    
    molecules = {}

    f = open(options.infile)
    sys.stderr.write("Loaded file: %s \n" % options.infile)

    sys.stderr.write("Parsing file... \n")

    for num, line in enumerate(f):
        if line.startswith("> <ChEBI ID>"):
            chebid = f.readline(num+1).strip()
            chebid = re.findall(r'\d+', chebid)[0]
            molecules[chebid] = {}
        if line.startswith("> <ChEBI Name>"):
            name = f.readline(num+1).strip()
            molecules[chebid]['name'] = name
        if line.startswith("> <InChIKey>"):
            inchk = f.readline(num+1).strip().split('-', 1)[0]
            molecules[chebid]['inchikey'] = inchk

    sys.stderr.write("%d molecules found \n" % len(molecules))
    
    ut.save_json(molecules, options.outfile, compressed=True)
    sys.stderr.write("File saved as %s \n" % options.outfile)
