import warnings

import numpy as np
import pandas as pd
from sklearn import preprocessing
from sklearn.impute import KNNImputer


class Preprocessing(object):
    def process(self):
        raise NotImplementedError()


class ZeroAndNegativeReplace(Preprocessing):
    def process(self, df):
        df[df <= 0.0] = np.nan
        return df


class MissingValueImputation(Preprocessing):
    def __init__(self, study_design):
        self.study_design = study_design

    def process(self, df):
        df = df.dropna(axis = 0, how = 'all')
        
        for grp1, grp2 in [('case', 'control'), ('control', 'case')]:
            grp1_spl = self.study_design[grp1]['samples']
            grp1_df = df.loc[:, grp1_spl].T
            
            isna = grp1_df.isna().sum()
            nan_prop = isna.loc[isna != 0]/len(grp1_df)
            alignid = nan_prop.loc[nan_prop == 1].index.values
            
            if len(alignid) != 0:
                grp2_spl = self.study_design[grp2]['samples']
                replacement_value = df.loc[alignid, grp2_spl].min(axis=1) * 0.1

                df.loc[alignid, grp1_spl] = replacement_value

        return df


class RowAverageImputation(Preprocessing):
    def __init__(self, study_design):
        self.study_design = study_design

    def process(self, df, mode='1/10'):
        for grp in self.study_design.values():
            samples = grp['samples']
            # if group values not all zeros, replace the zeros with mean of group
            group_df = df.loc[:, samples]
            if mode == '1/10':
                replacement_value = 0.1 * group_df.min(axis=1)
            elif mode == '1/5':
                replacement_value = 0.2 * group_df.min(axis=1)
            elif mode == '1':
                replacement_value = 1
            df.loc[:, samples] = group_df.mask(group_df.isna(),
                                               replacement_value, axis=0)

        return df


class KNNImputation(Preprocessing):
    def __init__(self, study_design, K=5):
        self.study_design = study_design
        self.K = K

    def process(self, df):
        for grp in self.study_design.values():
            samples = grp['samples']
            group_df = df.loc[:, samples]
            # if group values not all zeros, perform KNN imputation
            if np.sum(group_df.isnull().values) > 0:
                imputer = KNNImputer(n_neighbors=self.K)
                imputed_df = pd.DataFrame(imputer.fit_transform(group_df), index=group_df.index,
                                          columns=group_df.columns)
                df.loc[:, samples] = imputed_df
        return df


class LogNormalisation(Preprocessing):
    def process(self, df):
        return np.log(df)


class ZScoreNormalisation(Preprocessing):
    def process(self, df):
        with warnings.catch_warnings():
            warnings.filterwarnings('error')
            try:
                scaler = preprocessing.StandardScaler()
                scaled_df = scaler.fit_transform(df.transpose())  # transpose into the right shape for StandardScaler
                df = pd.DataFrame(scaled_df.transpose(), columns=df.columns,
                                  index=df.index)  # return the original shape
                return df
            except UserWarning as e:
                raise (e)
