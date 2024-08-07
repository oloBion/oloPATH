import os
import olopath.oloutils as ut
from olopath.variables import MOLID, INCHIKEY, MOLNM
from collections import defaultdict


class Database(object):
    
    def __init__(self, species):

        self.species = species
        directory = os.path.join(os.path.dirname(__file__), 'data')
        json_file = os.path.abspath(os.path.join(directory,
                                                 '%s.json.zip' % species))
        self.database = ut.load_json(json_file, compressed=True)

    def load(self):
        return self.database


class DataSource(object):

    def __init__(self, intensity_df, annotation_df, study_design, species):
        """
        Creates a data source for oloPATH analysis
        :param intesities_df: a dataframe of peak intensities, where
        index = alignment id and columns = sample name
        :param annotation_df: a dataframe where index = alignment id
        and column = InChiKey
        :param study_design: a dictionary specifying the study design, where:
            - Keys are 'case' and 'control'
            - Values are a dictionary containing the name of the group ('name')
            and a list of samples ('samples')
        :param species: the species name (Homo sapiens or Mus musculus)
        """

        self.intensity_df = intensity_df
        self.annotation_df = annotation_df
        self.study_design = study_design
        self.species = species.capitalize()

        self.database = Database(species).load()

        self.molecules = self.database['molecules']
        self.pathways = self.database['pathways']

        self.annotation_molid_df = self.get_molid_from_inchikey()

        self.mols_in_pathways = self.molecules_in_pathways()

        self.pathways_in_data = self.pathways_in_dataset()

        self.species_molid_count = len(self.database['molecules'])
        self.dataset_molid_count = len(self.annotation_molid_df)


    def get_molid_from_inchikey(self):
        self.inchikey = self.database['inchikey']

        annotation_molid_df = self.annotation_df.copy()

        for alignid, row in annotation_molid_df.iterrows():
            inchk = row[INCHIKEY]
            molnm = row[MOLNM]
            try:
                molid = self.inchikey[inchk][0]
            except KeyError:
                continue

            annotation_molid_df.loc[alignid, MOLID] = molid
            annotation_molid_df.loc[alignid, MOLNM] = molnm
        annotation_molid_df = annotation_molid_df[[MOLID, MOLNM]].dropna()
        return annotation_molid_df


    def molecules_in_pathways(self):
        molecules_in_pathways = defaultdict(set)
        for molid in self.molecules.keys():
            pathways = self.molecules[molid]['pathways']
            for path in pathways:
                molecules_in_pathways[path].add(molid)
        return molecules_in_pathways


    def pathways_in_dataset(self):
        pathways_in_dataset = {}
        for alignid, row in self.annotation_molid_df.iterrows():
            molid = row[MOLID]
            pathways = self.molecules[molid]['pathways']
            for path in pathways:
                if path not in pathways_in_dataset.keys():
                    name = self.pathways[path]['name']
                    pathways_in_dataset[path] = {'name': name,
                                                 'alignid': []}
                pathways_in_dataset[path]['alignid'].append(alignid)
        return pathways_in_dataset


