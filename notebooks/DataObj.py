import os
import src.utils as ut
from collections import defaultdict


INCHIKEY = 'inchikey'
CHEBID = 'chemi_id'


class Database(object):
    
    def __init__(self, species):

        self.species = species
        json_file = os.path.abspath(
            '../data/metabolic_pathways/%s.json.zip' % species)
        
        self.database = ut.load_json(json_file, compressed=True)

    def load(self):
        return self.database


class DataSource(object):

    def __init__(self, intensities_df, annotation_df, study_design, species):
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

        self.intensities_df = intensities_df
        self.annotation_df = annotation_df
        self.study_design = study_design
        self.species = species.capitalize()

        self.database = Database(species).load()

        self.molecules = self.database['molecules']
        self.pathways = self.database['pathways']

        self.annotation_chebid_df = self.get_chebid_from_inchikey()

        self.mols_in_pathways = self.molecules_in_pathways()

        self.paths_in_data, self.alignid_in_paths = self.pathways_in_dataset()

        self.species_chebid_count = len(self.database['molecules'])
        self.dataset_chebid_count = len(self.annotation_df)


    def get_chebid_from_inchikey(self):
        self.inchikey = self.database['inchikey']

        annotation_chebid_df = self.annotation_df.copy()

        for alignid, row in annotation_chebid_df.iterrows():
            inchk = row['inchikey']
            try:
                chebid = self.inchikey[inchk][0]
            except KeyError:
                continue

            annotation_chebid_df.loc[alignid, CHEBID] = chebid
        annotation_chebid_df = annotation_chebid_df[CHEBID].dropna()
        return annotation_chebid_df


    def molecules_in_pathways(self):
        molecules_in_pathways = defaultdict(set)
        for chebid in self.molecules.keys():
            pathways = self.molecules[chebid]['pathways']
            for path in pathways:
                molecules_in_pathways[path].add(chebid)
        return molecules_in_pathways


    def pathways_in_dataset(self):
        alignid_in_pathways = defaultdict(list)
        for alignid, chebid in self.annotation_chebid_df.items():
            pathways = self.molecules[chebid]['pathways']
            for path in pathways:
                alignid_in_pathways[path].append(alignid)
        pathways_in_dataset = set(alignid_in_pathways.key())
        return pathways_in_dataset, alignid_in_pathways

