![oloPATH](images/olobion-logo.png)

# oloPATH

Package for metabolomic pathways analysis.

## Installation

oloPATH requires Python 3.10 and it can be installed from GitHub repository via:

```sh
pip install olopath@git+https://github.com/oloBion/oloPATH.git
```

## Input data

oloPATH requires the following input data in CSV format:

1. `data`: is a matrix containing metabolite information, where rows represent metabolites and columns correspond to metabolite IDs (`alignid`), metabolite names (`metabolite_name`), metabolite InChiKeys (`inchikey`) and individual samples with their associated intensities.

<center>

|alignid|metabolite_name|inchikey|C1 |C2 |C3 |K1 |K2 |K3 |
|:-----:|:-------------:|:------:|:-:|:-:|:-:|:-:|:-:|:-:|
|1|Phenol sulfate|CTYRPMDGLDAWRQ-UHFFFAOYSA-N|37833|19019|2648536|127311|2368521|19525|
|2|Nicotinamide|DFPAKSUCGFBDDF-UHFFFAOYSA-N|122401|54811|418613|95612|101095|153269

</center>

<br>

2. `study_design`: is a table of two columns that associates study groups (`group`) with individual samples (`sample`).

<center>

|group|sample|
|:---:|:----:|
|Control|C1|
|Control|C2|
|Control|C3|
|KO|K1|
|KO|K2|
|KO|K3|

</center>

## Species databases structure

Each database consists of two dictionaries: `molecules` and `pathways`.

- `molecules` dictionary has ChEBI or PlantCyc molecules identifiers as keys and dictionaries containing the molecule name, short InChiKey and associated Reactome or PlantCyc pathways identfiers as values.

```python
"molecules": {"10055": {"name": "Xamoterol",
                        "inchikey": "DXPOSRCHIDYWHW",
                        "pathways": ["R-HSA-162582",
                                     "R-HSA-372790",
                                     "R-HSA-373076",
                                     "R-HSA-375280",
                                     "R-HSA-390696",
                                     "R-HSA-500792"]
                       }
             }
```

<br>

- `pathways` dictionary has Reactome or PlantCyc pathways identifiers as keys and dictionaries containing the pathway name as values.

```python
"pathways": {"R-HSA-162582": {"name": "Signal Transduction"},
             "R-HSA-372790": {"name": "Signaling by GPCR"},
             "R-HSA-373076": {"name": "Class A/1 (Rhodopsin-like receptors)"},
             "R-HSA-375280": {"name": "Amine ligand-binding receptors"},
             "R-HSA-390696": {"name": "Adrenoceptors"},
             "R-HSA-500792": {"name": "GPCR ligand binding"}
            }
```

## Running oloPATH

Import oloPATH as a Python library and perform your pathway analysis following the instructions:

```python
import olopath.oloutils as ut
from olopath.DataObj import DataSource
from olopath.PATHAnalysis import PATHAnalysis
```

<br>

1. **Load your data**

Ensure your input files are in the proper format as specified in the [Input Data section](#input-data) and run the `load_data` function, specifying the `case` and `control` groups for analysis. The function will return `annotation_df`, `intensity_df` and `std_design` objects.

```python
annotation_df, intensity_df, study_design =\
    ut.load_data(data = "data.csv", study_design = "study_design.csv",
                 case='KO', control='Control')
```

<br>

2. **Initialize `DataSource` object**

Create the `DataSource` object with the loaded data in the previous step, specifiying the `species`, `pvalue` and `foldchange2`, `logscale` and `mode` parameters.

```python
ds = DataSource(intensity_df=intesity_df,
                annotation_df=annotation_df,
                study_design=study_design,
                species='Homo sapiens',
                pvalue=0.05,
                foldchange2=[-0.5, 0.5],
                logscale=False,
                mode='1/10')
```

Considerations for the object parameters:

- `species`: indicates the species for which the analysis is conducted. oloPATH supports humans (`Homo sapiens`), mouse (`Mus musculus`) and plants (`PlantCyc`).
- `pvalue`: specifies the threshold for statistical significance. It is set to `0.05` by default.
- `foldchange2`: specifies the cutoff values for log2(fold change). It is set to `[-0.5, 0.5]` by default.
- `logscale`: indicates if data will be log-transformed for p-values computation. It set to `False` by default.
- `mode`: indicates the approach for handling missing values, with three options: `1/10` (replace missing values with 1/10 of minimum positive value of each analyte), `1/5` (replace with 1/5 of minimum positive value of each analyte) and `1` (replace with 1).

<br>

3. **Run `PATHAnalysis`**

Run `PATHAnalysis` with the `DataSource` object created in the previous step to perform pathway analysis. The class will return `pathway_df` and `metabolites_df` objects.

```python
PATH = PATHAnalysis(ds)
pathway_df, metabolites_df = PATH.get_results(filter_by_hits=1)
```

Pathways can be filtered by hits using the `filter_by_hits` parameter, which is set to `1` by default.
