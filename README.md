![oloPATH](images/olobion-logo.png)

# oloPATH

Package for metabolomic pathways analysis.

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

```
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

```
"pathways": {"R-HSA-162582": {"name": "Signal Transduction"},
             "R-HSA-372790": {"name": "Signaling by GPCR"},
             "R-HSA-373076": {"name": "Class A/1 (Rhodopsin-like receptors)"},
             "R-HSA-375280": {"name": "Amine ligand-binding receptors"},
             "R-HSA-390696": {"name": "Adrenoceptors"},
             "R-HSA-500792": {"name": "GPCR ligand binding"}
            }
```
