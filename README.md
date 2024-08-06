![oloPATH](images/olobion-logo.png)

# oloPATH

Package for metabolomic pathways analysis.

## Input data

oloPATH requires the following input data:

1. `intensity_df`: is a matrix of metabolites intensities, where rows represent metabolites and columns represent individual samples. The first column must be the metabolite id (`alignid`) and needs to be set as index.

2. `annotation_df`: is a table of two columns that associates metabolite id (`alignid`) with metabolite InChiKey (`inchikey`). The `alignid` column has to be set as index.

3. `study_design`: is a dictionary containing information about groups as follows:

```
{'control': {'name': 'Control',
             'samples: ['C1', 'C2', 'C3', 'C4']
            },
 'case': {'name': 'KO',
          'samples: ['K1', 'K2', 'K3', 'K4']
         },
}
```


## Species databases structure

Each database consists of three dictionaries: `molecules`, `inchikey` and `pathways`.

- `molecules` dictionary has ChEMI molecules identifiers as keys and dictionaries containing the molecule name, short InChiKey and associated Reactome pathways identfiers as values.

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

- `inchikey` dictionary has short InChiKey as keys and lists containing the ChEMI identifiers associated as values.

```
"inchikey": {"DXPOSRCHIDYWHW": ["10055"]}
```

<br>

- `pathways` dictionary has Reactome pathways identifiers as keys and dictionaries containing the pathway name as values.

```
"pathways": {"R-HSA-162582": {"name": "Signal Transduction"},
             "R-HSA-372790": {"name": "Signaling by GPCR"},
             "R-HSA-373076": {"name": "Class A/1 (Rhodopsin-like receptors)"},
             "R-HSA-375280": {"name": "Amine ligand-binding receptors"},
             "R-HSA-390696": {"name": "Adrenoceptors"},
             "R-HSA-500792": {"name": "GPCR ligand binding"}
            }
```
