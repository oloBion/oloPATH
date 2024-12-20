import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
from scipy.stats import hypergeom, fisher_exact
from statsmodels.stats.multitest import multipletests
import olopath.preprocessing as pcss
from olopath.variables import PATHID, PATHNM, PVALUE, PATH_COM, HITS, PATH_COV, \
    FDR, ALIGNID, INCHIKEY, ORIG_INCHIK, MOLID, MOLNM, DB_MOLNM, TOTAL_HITS, \
    FC2, HOLM


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
        pathway_df = self.get_pathway_fold_change(pathway_df)
        pathway_df = self.compute_pathway_significance(pathway_df)
        pathway_df = self.correct_pvalue(pathway_df)

        pathway_df = self.filter_by_min_hits(pathway_df, filter_by_hits)

        if len(pathway_df) != 0:
            metabolites_df = self.get_metabolites_df()
        else:
            metabolites_df = pd.DataFrame(columns=[PATHID, PATHNM, ALIGNID,
                                                   INCHIKEY, MOLNM, MOLID])

        return pathway_df, metabolites_df


    def create_pathway_dataframe(self):
        df = []
        for pathid, values in self.data.pathways_in_data.items():
            pathnm = values['name']
            significant_hits = len(values['sign_inchikey'])
            total_hits = len(values['inchikey'])
            pathcpds = len(self.data.mols_in_pathways[pathid])
            data = [pathid, pathnm, significant_hits, total_hits, pathcpds]
            
            df.append(data)
        column_names = [PATHID, PATHNM, HITS, TOTAL_HITS, PATH_COM]
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
    

    def correct_pvalue(self, df, alpha=0.05):
        fdr_padj = multipletests(df[PVALUE], alpha, "fdr_bh")
        df[FDR] = fdr_padj[1]

        holm_padj = multipletests(df[PVALUE], alpha, "holm")
        df[HOLM] = holm_padj[1]

        return df


    def get_pathway_fold_change(self, df):
        log_fold_change = []

        pathway_ids = df.index

        for pathid in pathway_ids:
            alignid = self.data.pathways_in_data[pathid]["alignid"]
            fc2 = self.data.statistics_df.loc[alignid, FC2].mean()

            log_fold_change.append(fc2)
        
        df["%s mean" % FC2] = log_fold_change

        return df


    def filter_by_min_hits(self, df, filter_by_hits):
        paths_to_remove = df.loc[df[HITS] <= filter_by_hits, :].index.to_list()

        df = df.drop(index=paths_to_remove)

        for pathid in paths_to_remove:
            del self.filtered_pathways[pathid]

        return df


    def get_molid_from_inchikey(self, df):
        for _, row in df.iterrows():
            pathid = row[PATHID]
            inchk = row[INCHIKEY]
            inch_molids = self.data.inchikey[inchk]['molid']
            if len(inch_molids) > 1:
                filtered_molid = []
                for id in inch_molids:
                    if id in self.data.mols_in_pathways[pathid]:
                        filtered_molid.append(id)
                molid = filtered_molid[0]
                molnm = self.data.molecules[molid]['name']
                df.loc[(df[PATHID] == pathid) &
                       (df[INCHIKEY] == inchk), (MOLID, DB_MOLNM)] = (molid,
                                                                      molnm)
            else:
                molid = list(inch_molids)[0]
                molnm = self.data.molecules[molid]['name']
                df.loc[(df[PATHID] == pathid) &
                       (df[INCHIKEY] == inchk), (MOLID, DB_MOLNM)] = (molid,
                                                                      molnm)
        return df
    

    def get_metabolites_df(self):
        metabolites_df = pd.DataFrame(self.filtered_pathways).T
        metabolites_df.index.name = PATHID
        metabolites_df = metabolites_df.drop(columns=['inchikey', "sign_inchikey"])
        metabolites_df.rename(columns={'name': PATHNM,
                                       'alignid': ALIGNID}, inplace=True)
        metabolites_df.reset_index(inplace=True)
        metabolites_df = metabolites_df.explode(ALIGNID)

        annotation_df = self.data.annotation_df

        metabolites_df = pd.merge(metabolites_df, annotation_df,
                                  on=ALIGNID)
        metabolites_df.rename(columns={'name': PATHNM}, inplace=True)

        metabolites_df = self.get_molid_from_inchikey(metabolites_df)

        metabolites_df[INCHIKEY] = metabolites_df[ORIG_INCHIK]
        metabolites_df = metabolites_df.drop(columns=[ORIG_INCHIK])

        statistics_df = self.data.statistics_df[[FC2, PVALUE]]

        metabolites_df = pd.merge(metabolites_df, statistics_df,
                                  on=ALIGNID)

        return metabolites_df
