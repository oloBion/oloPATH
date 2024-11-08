import os
import olopath.oloutils as ut
from olopath.variables import MOLID, INCHIKEY, MOLNM, ORIG_INCHIK, DB_MOLNM, \
    ALIGNID, PVALUE, FC2
from collections import defaultdict
import olopath.preprocessing as pcss
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind


class Database(object):
    
    def __init__(self, species):

        self.species = species
        directory = os.path.join(os.path.dirname(__file__), 'data', 'metabolic_pathways')
        json_file = os.path.abspath(os.path.join(directory,
                                                 '%s.json.zip' % species))
        self.database = ut.load_json(json_file, compressed=True)

    def load(self):
        return self.database


class DataSource(object):

    def __init__(self, intensity_df, annotation_df, study_design, species,
                 pvalue=0.05, logscale=False, mode='1/10'):
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

        self.study_design = study_design
        intensity_df = self.filter_intesities_by_groups(intensity_df,
                                                        study_design)
        if self.check_missing_values(intensity_df):
            self.intensity_df = self.preprocess_data(intensity_df, mode)
        else:
            self.intensity_df = intensity_df
        self.annotation_df = self.filter_annotation_with_preprocessed_data(annotation_df)
        self.species = species.capitalize()

        self.statistics_df = self.compute_pvalues(logscale)

        self.database = Database(species).load()

        self.molecules = self.database['molecules']
        self.pathways = self.database['pathways']
        self.inchikey = self.create_inchikey_db()

        self.annotation_molid_df = self.get_molid_from_inchikey()

        self.mols_in_pathways = self.molecules_in_pathways()

        self.pathways_in_data = self.pathways_in_dataset(pvalue)

        self.species_molid_count = len(self.molecules)
        self.dataset_molid_count = len(self.annotation_molid_df)


    def check_missing_values(self, df):
        is0 = df.isin([0]).sum().sum()
        isnan = df.isna().sum().sum()
        return is0 != 0 or isnan != 0
    
    def filter_intesities_by_groups(self, df, study_design):
        case_samples = study_design['case']['samples']
        control_samples = study_design['control']['samples']
        samples = case_samples + control_samples

        df = df.loc[:, samples]

        return df
    

    def preprocess_data(self, df, mode):
        study_design = self.study_design

        df = pcss.ZeroAndNegativeReplace().process(df)
        df = pcss.MissingValueImputation(study_design).process(df)
        df = pcss.RowAverageImputation(study_design).process(df, mode)

        return df
    

    def filter_annotation_with_preprocessed_data(self, df):
        alignids = self.intensity_df.index

        df = df.loc[alignids]

        return df
    

    def compute_pvalues(self, logscale):
        study_design = self.study_design
        annotation_df = self.annotation_df
        intensity_df = self.intensity_df
        

        case_samples = study_design['case']['samples']
        case_df = intensity_df.loc[:, case_samples]

        control_samples = study_design['control']['samples']
        control_df = intensity_df.loc[:, control_samples]

        fc = case_df.mean(axis=1) / control_df.mean(axis=1)
        fc2 = np.log2(fc)

        if logscale:
            case = pcss.LogNormalisation().process(case_df)
            control = pcss.LogNormalisation().process(control_df)
        else:
            case = case_df
            control = control_df

        pvalues = ttest_ind(case, control, equal_var=False, axis=1)

        stats_df = pd.DataFrame(data={FC2: fc2, PVALUE: pvalues.pvalue},
                                index=intensity_df.index)
        
        stats_df = pd.merge(annotation_df, stats_df, on=ALIGNID)

        return stats_df


    def get_molid_from_inchikey(self):
        annotation_molid_df = self.annotation_df.copy()

        for alignid, row in annotation_molid_df.iterrows():
            orig_inchik = row[ORIG_INCHIK]
            inchk = row[INCHIKEY]
            molid_first = list(self.inchikey[inchk]['molid'])
            if len(molid_first) == 1:
                molid = molid_first[0]
                dbnme = self.molecules[molid]["name"]
                annotation_molid_df.loc[alignid, MOLID] = molid
                annotation_molid_df.loc[alignid, DB_MOLNM] = dbnme
            elif len(molid_first) > 1:
                twop_inchik = "-".join(orig_inchik.split("-", 2)[:2])
                molid_second = []
                for id in molid_first:
                    db_inchik = self.molecules[id]["inchikey"]
                    db_twop_inchik = "-".join(db_inchik.split("-", 2)[:2])
                    if db_twop_inchik == twop_inchik:
                        molid_second.append(id)
                if len(molid_second) == 1:
                    molid = molid_second[0]
                    dbnme = self.molecules[molid]["name"]
                    annotation_molid_df.loc[alignid, MOLID] = molid
                    annotation_molid_df.loc[alignid, DB_MOLNM] = dbnme
                elif len(molid_second) > 1:
                        molid_third = []
                        for id in molid_second:
                            db_inchik = self.molecules[id]["inchikey"]
                            if db_inchik == orig_inchik:
                                molid_third.append(id)
                        if len(molid_third) == 1:
                            molid = molid_third[0]
                            dbnme = self.molecules[molid]["name"]
                            annotation_molid_df.loc[alignid, MOLID] = molid
                            annotation_molid_df.loc[alignid, DB_MOLNM] = dbnme
                        else:
                            pathnum = 0
                            for id in molid_second:
                                molpathws = len(self.molecules[id]["pathways"])
                                if molpathws > pathnum:
                                    pathnum = molpathws
                                    molid = id
                                    dbnme = self.molecules[molid]["name"]
                            annotation_molid_df.loc[alignid, MOLID] = molid
                            annotation_molid_df.loc[alignid, DB_MOLNM] = dbnme
                else:
                    pathnum = 0
                    for id in molid_first:
                        molpathws = len(self.molecules[id]["pathways"])
                        if molpathws > pathnum:
                            pathnum = molpathws
                            molid = id
                            dbnme = self.molecules[molid]["name"]
                    annotation_molid_df.loc[alignid, MOLID] = molid
                    annotation_molid_df.loc[alignid, DB_MOLNM] = dbnme

        annotation_molid_df = annotation_molid_df[[MOLID, MOLNM,
                                                   DB_MOLNM]].dropna()
        return annotation_molid_df
    

    def create_inchikey_db(self):
        pathways_in_inchikeys = defaultdict(lambda: defaultdict(set),
                                                    defaultdict(set))
        for molid in self.molecules.keys():
            try:
                inchikey = self.molecules[molid]['inchikey'].split('-', 1)[0]
                pathways = self.molecules[molid]['pathways']
                pathways_in_inchikeys[inchikey]['pathways'].update(pathways)
                pathways_in_inchikeys[inchikey]['molid'].add(molid)
            except KeyError:
                continue
        return pathways_in_inchikeys


    def molecules_in_pathways(self):
        molecules_in_pathways = defaultdict(set)
        for molid in self.molecules.keys():
            pathways = self.molecules[molid]['pathways']
            for path in pathways:
                molecules_in_pathways[path].add(molid)
        return molecules_in_pathways


    def pathways_in_dataset(self, pvalue):
        pathways_in_dataset = {}
        for alignid, row in self.statistics_df.iterrows():
            inchik = row[INCHIKEY]
            alignid_pvalue = row[PVALUE]
            pathways = self.inchikey[inchik]['pathways']
            for path in pathways:
                if path not in pathways_in_dataset.keys():
                    name = self.pathways[path]['name']
                    pathways_in_dataset[path] = {'name': name,
                                                 'alignid': [],
                                                 'inchikey': set(),
                                                 'sign_inchikey': set()}
                pathways_in_dataset[path]['inchikey'].add(inchik)
                if alignid_pvalue <= pvalue:
                    pathways_in_dataset[path]['alignid'].append(alignid)
                    pathways_in_dataset[path]['sign_inchikey'].add(inchik)

        return pathways_in_dataset

