import pandas as pd
import numpy as np

class PATHAnalysis(object):

    def __init__(self, data):
        """
        Creates PATHAnalysis object
        :param data: a DataSource object
        """

        self.data = data

    def calculate_pathway_activity(self, intensity_df):
        """
        Calculates pathway activity given a dataframe of standardized
        intensities
        :param intensity_df: a standardized dataframe of peak intensites
        :return: a dataframe with pathway names (rows) and the SVD activity
        levels for the samples (columns)
        """
        activity_df = pd.DataFrame()
        for pathw, alignid in self.data.alignid_in_paths.items():
            data = intensity_df.loc[alignid]
            try:
                w, d, c = np.linalg.svd(np.array(data))
            except np.linalg.LinAlgError:
                w, d, c = np.linalg.svd(np.array(data))
            
            # ... ; cambio en el dictionary alignid_in_paths to:
            # pathways_in_dataset = {pathid: {"name": "pw_name", 'alignid': "row ids"}}