import utils as ut
from bisect import insort

import sys
import argparse

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="")

    parser.add_argument('-chm', '--chebi_molecules',
                        dest= "chebi_molecules",
                        action= "store",
                        required= True,
                        help="")
    
    parser.add_argument('-chp', '--chebi2reactome',
                        dest= "chebi2reactome",
                        action= "store",
                        required= True,
                        help="")
    
    parser.add_argument('-rp', '--reactome_pathways',
                        dest= "reactome_pathways",
                        action= "store",
                        required= True,
                        help="")
    
    parser.add_argument('-s', '--species',
                        dest= "species",
                        action= "store",
                        required= True,
                        choices=['Homo sapiens', 'Mus musculus'],
                        help="")
    
    parser.add_argument('-o', '--ouput',
                        dest= "outfile",
                        action= "store",
                        required= True,
                        help="")
    
    options = parser.parse_args()

    species = options.species

    chmolecules = ut.load_json(options.chebi_molecules, compressed=True)
    sys.stderr.write("Loaded file: %s \n" % options.chebi_molecules)

    ch2reactome = ut.load_json(options.chebi2reactome, compressed=True)
    ch2reactome = ch2reactome[species]
    sys.stderr.write("Loaded file: %s \n" % options.chebi2reactome)

    reactome_pathways = ut.load_json(options.reactome_pathways, compressed=True)
    reactome_pathways = reactome_pathways[species]
    sys.stderr.write("Loaded file: %s \n" % options.reactome_pathways)

    sys.stderr.write("ChEBI molecules: %d \n" % len(chmolecules))
    sys.stderr.write("%s molecules: %d \n" % (species, len(ch2reactome)))
    
    diff = set(ch2reactome.keys()).difference(set(chmolecules.keys()))
    sys.stderr.write("Difference between %s and total molecules: %d \n" % (species, len(diff)))

    sys.stderr.write("%s pathways: %d \n" % (species, len(reactome_pathways)))

    sys.stderr.write("Creating %s database... \n" % species)

    organism_database = {'molecules': {},
                         'inchikey': {},
                         'pathways': reactome_pathways}
    
    for mol in ch2reactome.keys():
        try:
            organism_database['molecules'][mol] = chmolecules[mol]
            organism_database['molecules'][mol]['pathways'] = ch2reactome[mol]
            inch = chmolecules[mol]['inchikey']
            if inch not in organism_database['inchikey'].keys():
                organism_database['inchikey'][inch] = []
            insort(organism_database['inchikey'][inch], mol)
        except KeyError:
            continue
    
    sys.stderr.write("Database created \n")

    ut.save_json(organism_database, options.outfile, compressed=True)
    sys.stderr.write("File saved as %s \n" % options.outfile)
