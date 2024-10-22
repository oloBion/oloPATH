import oloutils as ut
from bisect import insort

import sys
import argparse

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="")

    parser.add_argument('-pcm', '--plantcyc_molecules',
                        dest= "plantcyc_molecules",
                        action= "store",
                        required= True,
                        help="")
    
    parser.add_argument('-pcp', '--plantcyc_pathways',
                        dest= "plantcyc_pathways",
                        action= "store",
                        required= True,
                        help="")
    
    parser.add_argument('-o', '--ouput',
                        dest= "outfile",
                        action= "store",
                        required= True,
                        help="")
    
    options = parser.parse_args()

    pcmolecules = ut.load_json(options.plantcyc_molecules, compressed=False)
    sys.stderr.write("Loaded file: %s \n" % options.plantcyc_molecules)

    plantcyc_pathways = ut.load_json(options.plantcyc_pathways, compressed=False)
    sys.stderr.write("Loaded file: %s \n" % options.plantcyc_pathways)
    
    sys.stderr.write("Total molecules: %d \n" % len(pcmolecules))

    sys.stderr.write("Total pathways: %d \n" % len(plantcyc_pathways))

    sys.stderr.write("Creating plantcyc database... \n")

    organism_database = {'molecules': {},
                         'pathways': {}}
    
    for mol in pcmolecules.keys():
        if 'inchikey' in pcmolecules[mol].keys():
            molnm = pcmolecules[mol]['names'][0]
            molinchi = pcmolecules[mol]['inchikey']
            organism_database["molecules"][mol] = {'name': molnm,
                                                   'inchikey': molinchi,
                                                   'pathways': []}
            
    sys.stderr.write("Total molecules in db: %d \n" % len(organism_database['molecules']))

    for path in plantcyc_pathways:
        pathid = path['id']
        pathnm = path['name']
        organism_database["pathways"][pathid] = {'name': pathnm}
        pathmol = path["reaction_metaboltites"]
        for mol in pathmol:
            try:
                organism_database["molecules"][mol]['pathways'].append(pathid)
            except KeyError:
                continue
    
    sys.stderr.write("Database created \n")

    ut.save_json(organism_database, options.outfile, compressed=True)
    sys.stderr.write("File saved as %s \n" % options.outfile)
