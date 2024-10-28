import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
from scipy.stats import hypergeom
import olopath.preprocessing as pcss
from olopath.variables import PATHID, PATHNM, PVALUE, PATH_COM, HITS, PATH_COV, \
    PATH_SIG, ALIGNID, INCHIKEY, MOLID, MOLNM


class PATHAnalysis(object):

    def __init__(self, data):
        """
        Creates PATHAnalysis object
        :param data: a DataSource object
        """

        self.data = data
        self.filtered_pathways = self.data.pathways_in_data.copy()


    def get_results(self, filter_by_hits=1):
        pathway_df = self.create_pathway_dataframe()
        pathway_df = self.compute_pathways_coverage(pathway_df)
        pathway_df = self.compute_pathway_significance(pathway_df)

        pathway_df = self.filter_by_min_hits(pathway_df, filter_by_hits)

        if len(pathway_df) != 0:
            metabolites_df = self.get_metabolites_df()
        else:
            metabolites_df = pd.DataFrame(columns=[PATHID, PATHNM, ALIGNID,
                                                   INCHIKEY, MOLNM, MOLID])

        return pathway_df, metabolites_df
        

    def preprocess_data(self, df):
        df = pcss.LogNormalisation().process(df)
        df = pcss.ZScoreNormalisation().process(df)

        return df


    def create_pathway_dataframe(self):
        df = []
        for pathid, values in self.data.pathways_in_data.items():
            pathnm = values['name']
            hits = len(values['alignid'])
            pathcpds = len(self.data.mols_in_pathways[pathid])
            data = [pathid, pathnm, hits, pathcpds]
            
            df.append(data)
        column_names = [PATHID, PATHNM, HITS, PATH_COM]
        df = pd.DataFrame(df, columns=column_names).set_index(PATHID)
        return df


    def compute_pathways_coverage(self, df):
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
        
        df[PVALUE] = pathway_significance

        return df


    def filter_by_min_hits(self, df, filter_by_hits):
        paths_to_remove = df.loc[df[HITS] <= filter_by_hits, :].index.to_list()

        df = df.drop(index=paths_to_remove)

        for pathid in paths_to_remove:
            del self.filtered_pathways[pathid]

        return df


    def get_molid_from_inchikey(self, df):
        for indx, row in df.iterrows():
            pathid = row[PATHID]
            inchk = row[INCHIKEY]
            inch_chebids = self.data.inchikey[inchk]['molid']
            if len(inch_chebids) > 1:
                filtered_chebid = []
                for chebid in inch_chebids:
                    if chebid in self.data.mols_in_pathways[pathid]:
                        filtered_chebid.append(chebid)
                df.loc[(df[PATHID] == pathid) &
                       (df[INCHIKEY] == inchk), MOLID] = filtered_chebid[0]
            else:
                df.loc[(df[PATHID] == pathid) &
                       (df[INCHIKEY] == inchk), MOLID] = list(inch_chebids)[0]
        return df
    

    def get_metabolites_df(self):
        metabolites_df = pd.DataFrame(self.filtered_pathways).T
        metabolites_df.index.name = PATHID
        metabolites_df.rename(columns={'name': PATHNM,
                                       'alignid': ALIGNID}, inplace=True)
        metabolites_df.reset_index(inplace=True)
        metabolites_df = metabolites_df.explode(ALIGNID)

        annotation_df = self.data.annotation_df

        metabolites_df = pd.merge(metabolites_df, annotation_df,
                                  on=ALIGNID)
        metabolites_df.rename(columns={'name': PATHNM}, inplace=True)

        metabolites_df = self.get_molid_from_inchikey(metabolites_df)

        return metabolites_df
