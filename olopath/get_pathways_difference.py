import oloutils as ut
import re

import sys
import argparse

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="")

    parser.add_argument('-pep', '--pathways_endpoints',
                        dest= "pathways_endpoints",
                        action= "store",
                        default= "./",
                        help="")
    
    parser.add_argument('-pal', '--pathways_allevels',
                        dest= "pathways_allevels",
                        action= "store",
                        default= "./",
                        help="")
    
    parser.add_argument('-o', '--ouput',
                        dest= "outfile",
                        action= "store",
                        default= "./",
                        help="")
    
    options = parser.parse_args()

    if len(sys.argv) <= 3:
        sys.stderr.write("No input or output provided. Please, try again. \n")
        exit()

    pathways = {}
    species = ['Homo sapiens', 'Mus musculus']

    # End points
    endp_file = open(options.pathways_endpoints)
    
    sys.stderr.write("Loaded file: %s \n" % options.pathways_endpoints)
    sys.stderr.write("Parsing file... \n")

    for line in endp_file:
        line = line.strip()
        line = line.split("\t")

        organism = line[5]
        pathid = line[1]

        if organism in species:            
            if organism not in pathways.keys():
                pathways[organism] = {}
            
            if "endpoint" not in pathways[organism].keys():
                pathways[organism]["endpoint"] = set()

            pathways[organism]["endpoint"].add(pathid)
    
    # All levels
    allvl_file = open(options.pathways_allevels)

    sys.stderr.write("Loaded file: %s \n" % options.pathways_allevels)
    sys.stderr.write("Parsing file... \n")

    for line in allvl_file:
        line = line.strip()
        line = line.split("\t")

        organism = line[5]
        pathid = line[1]

        if organism in species:            
            if organism not in pathways.keys():
                pathways[organism] = {}
            
            if "allevels" not in pathways[organism].keys():
                pathways[organism]["allevels"] = set()
            
            pathways[organism]["allevels"].add(pathid)
    
    # Difference
    sys.stderr.write('''Getting the pathways difference between
                     end points and all levels... \n''')
    
    for org, patw in pathways.items():
        endp_ptw = patw["endpoint"]
        allv_ptw = patw["allevels"]
        diff = allv_ptw.difference(endp_ptw)
        
        pathways[org]["endpoint"] = list(endp_ptw)
        pathways[org]["allevels"] = list(allv_ptw)
        pathways[org]["difference"] = list(diff)
    
    ut.save_json(pathways, options.outfile, compressed=True)
    sys.stderr.write("File saved as %s \n" % options.outfile)

    

