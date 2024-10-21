import olopath.oloutils as ut
from olopath.DataObj import DataSource
from olopath.PATHAnalysis import PATHAnalysis
import pandas as pd


def run():
    data = pd.read_excel("../tests/hsa_example_1.xlsx", "AlignmentTable")
    study_design = pd.read_excel("../tests/hsa_example_1.xlsx", "StudyDesign")

    annotation_df, intensity_df, study_design = ut.load_data(data,
                                                             study_design,
                                                             case='Group1',
                                                             control='Group2')
    
    ds = DataSource(intensity_df, annotation_df, study_design, "Homo sapiens")

    PATH = PATHAnalysis(ds)
    pathway_df, metabolites_df = PATH.get_results(filter_by_hits=1)

    return pathway_df, metabolites_df


def metabolites_in_path(df, pathid, metids):
    metids_df = df.loc[df["Pathway ID"] == pathid, "Metabolite ID"].to_list()
    metids_df = sorted(metids_df)
    metids_tt = sorted(metids)
    assert metids_df == metids_tt


def test():
    pathway_df, metabolites_df = run()
    metabolites_in_path(metabolites_df, "R-HSA-71240",
                        ["18344", "16946", "57427", "58315",
                         "57912", "58095", "57844", "58095",
                         "57762", "57972", "32513"])
    metabolites_in_path(metabolites_df, "R-HSA-418594",
                        ["16296", "58432", "20067",
                         "60119", "17992", "17992",
                         "58885", "456216"])
