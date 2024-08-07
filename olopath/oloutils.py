import gzip
import json
import pandas as pd
from src.variables import INCHIKEY, MOLNM, ALIGNID, GROUP, SAMPLE


def load_json(json_file, compressed=False):
    if compressed:
        with gzip.GzipFile(json_file, 'r') as f:
            data = json.loads(f.read().decode('utf-8'))
    else:
        with open(json_file, 'r') as f:
            data = json.load(f)
    return data

def save_json(data, json_file, compressed=False):
    if compressed:
        with gzip.GzipFile(json_file, 'w') as f:
            f.write(json.dumps(data).encode('utf-8'))
    else:
        with open(json_file, 'w') as f:
            json.dump(data, f)

def load_data(data, study_design, case, control):
    if not isinstance(data, pd.DataFrame):
        data = pd.read_csv(data)

    annotation_df = data[[ALIGNID, INCHIKEY, MOLNM]]
    annotation_df = annotation_df.set_index(ALIGNID)
    annotation_df[INCHIKEY] =\
        annotation_df[INCHIKEY].apply(lambda x: x.split("-")[0])
    intensity_df = data.drop([INCHIKEY, MOLNM], axis=1)
    intensity_df = intensity_df.set_index(ALIGNID)

    if not isinstance(study_design, pd.DataFrame):
        study_design = pd.read_csv(study_design)
    control_samples = study_design.loc[study_design[GROUP] == control,
                                       SAMPLE].values
    case_samples = study_design.loc[study_design[GROUP] == case,
                                    SAMPLE].values
    study_design = {'control': {'name': control,
                                'samples': list(control_samples)},
                    'case': {'name': case,
                             'samples': list(case_samples)}
                   }
    return annotation_df, intensity_df, study_design
    