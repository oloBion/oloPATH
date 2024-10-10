import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
from scipy.stats import hypergeom
import olopath.preprocessing as pcss
from olopath.variables import PATHID, PATHNM, PVALUE, PATH_COM, HITS, PATH_COV, \
    PATH_SIG, ALIGNID


class PATHAnalysis(object):

    def __init__(self, data):
        """
        Creates PATHAnalysis object
        :param data: a DataSource object
        """

        self.data = data
        self.filtered_pathways = self.data.pathways_in_data.copy()


    def get_results(self, filter_by_hits=1):
        intensity_df = self.data.intensity_df
        processed_intensity_df = self.preprocess_data(intensity_df)
        activity_df = self.calculate_pathway_activity(processed_intensity_df)
        pvalues_df = self.compute_pvalues(activity_df)
        coverage_df = self.compute_pathways_coverage(pvalues_df)
        pathway_df = self.compute_pathway_significance(coverage_df)

        pathway_df = self.filter_by_min_hits(pathway_df, filter_by_hits)

        metabolites_df = self.get_metabolites_df()

        return pathway_df, metabolites_df
        

    def preprocess_data(self, df):
        study_design = self.data.study_design

        df = pcss.ZeroAndNegativeReplace().process(df)
        df = pcss.MissingValueImputation(study_design).process(df)
        df = pcss.RowAverageImputation(study_design).process(df)
        df = pcss.LogNormalisation().process(df)
        df = pcss.ZScoreNormalisation().process(df)

        return df

    def calculate_pathway_activity(self, intensity_df):
        """
        Calculates pathway activity given a dataframe of standardized
        intensities
        :param intensity_df: a standardized dataframe of peak intensites
        :return: a dataframe with pathway names (rows) and the SVD activity
        levels for the samples (columns)
        """
        activity_df = []
        for pathid, values in self.data.pathways_in_data.items():
            pathnm = values['name']
            alignid = values['alignid']
            data = intensity_df.loc[alignid]
            try:
                w, d, c = np.linalg.svd(np.array(data))
            except np.linalg.LinAlgError:
                w, d, c = np.linalg.svd(np.array(data))
            transformed_data = [pathid, pathnm]
            transformed_data.extend(c[0])
            activity_df.append(transformed_data)
        column_names = [PATHID, PATHNM]
        column_names.extend(intensity_df.columns)
        activity_df = pd.DataFrame(activity_df,
                                   columns=column_names).set_index(PATHID)
        return activity_df

    def compute_pvalues(self, activity_df):
        """
        Obtains p-values dataframe
        :param activity_df: activity dataframe
        :return: a df containing pathway id, pathway names and p-values
        """
        study_design = self.data.study_design

        case_samples = study_design['case']['samples']
        case_df = activity_df.loc[:, case_samples]

        control_samples = study_design['control']['samples']
        control_df = activity_df.loc[:, control_samples]

        pvalues = ttest_ind(case_df, control_df, axis=1)

        pvalues_df = activity_df[PATHNM].to_frame()
        pvalues_df[PVALUE] = pvalues.pvalue

        return pvalues_df

    def all_pathways_molecules_counts(self, pathway_ids):
        counts = []
        for pathid in pathway_ids:
            mols = len(self.data.mols_in_pathways[pathid])
            counts.append(mols)
        return counts
    
    def identified_pathways_molecules_counts(self, pathway_ids):
        counts = []
        for pathid in pathway_ids:
            mols = len(self.data.pathways_in_data[pathid]['alignid'])
            counts.append(mols)
        return counts

    def compute_pathways_coverage(self, df):
        pathways_ids = df.index

        all_molecules_per_pathway =\
            self.all_pathways_molecules_counts(pathways_ids)
        df[PATH_COM] = all_molecules_per_pathway

        identified_molecules_per_pathway =\
            self.identified_pathways_molecules_counts(pathways_ids)
        df[HITS] = identified_molecules_per_pathway

        df[PATH_COV] = ((df[HITS]/df[PATH_COM]) * 100).round(2)

        return df
    

    def compute_pathway_significance(self, df):
        pathway_significance = []

        pathway_ids = df.index

        for pathid in pathway_ids:
            M = self.data.species_molid_count
            n = df.loc[pathid, PATH_COM]

            N = self.data.dataset_molid_count
            k = df.loc[pathid, HITS]
            
            sf = hypergeom.sf(k - 1, M, n, N)
            pathway_significance.append(sf)
        
        df[PATH_SIG] = pathway_significance

        return df


    def filter_by_min_hits(self, df, filter_by_hits):
        paths_to_remove = df.loc[df[HITS] <= filter_by_hits, :].index.to_list()

        df = df.drop(index=paths_to_remove)

        for pathid in paths_to_remove:
            del self.filtered_pathways[pathid]

        return df


    def get_metabolites_df(self):
        metabolites_df = pd.DataFrame(self.filtered_pathways).T
        metabolites_df.index.name = PATHID
        metabolites_df.rename(columns={'name': PATHNM,
                                       'alignid': ALIGNID}, inplace=True)
        metabolites_df.reset_index(inplace=True)
        metabolites_df = metabolites_df.explode(ALIGNID)

        annotation_molid_df = self.data.annotation_molid_df

        metabolites_df = pd.merge(metabolites_df, annotation_molid_df,
                                  on=ALIGNID)
        metabolites_df.rename(columns={'name': PATHNM}, inplace=True)

        return metabolites_df
