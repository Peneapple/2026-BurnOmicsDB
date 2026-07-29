# BurnOmicsDB

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21648675.svg)](https://doi.org/10.5281/zenodo.21648675)

**BurnOmicsDB** is a curated, standardized, and gene-centered database of human burn injury transcriptomics across tissues, post-injury time points, wound stages, and clinical contexts.

- **Live website:** https://peneapple.github.io/2026-BurnOmicsDB/
- **Archived release:** https://doi.org/10.5281/zenodo.21648675
- **GitHub repository:** https://github.com/Peneapple/2026-BurnOmicsDB

## Overview

BurnOmicsDB enables clinicians and researchers to search an official gene symbol, alias, NCBI Gene ID, Ensembl Gene ID, HGNC ID, or full gene name and examine standardized differential-expression evidence across human burn transcriptomic datasets.

Version 1.0.0 integrates six public NCBI Gene Expression Omnibus datasets:

| GEO accession | Main biological context | Data type |
|---|---|---|
| GSE139028 | Burn eschar versus normal skin | RNA-seq |
| GSE178411 | Uninjured skin, burn wounds, and hypertrophic scars | RNA-seq |
| GSE8056 | Burn wound margins across post-injury intervals | Affymetrix microarray |
| GSE19743 | Early and intermediate systemic blood responses | Affymetrix microarray |
| GSE37069 | Longitudinal systemic blood responses after severe burn injury | Affymetrix microarray |
| GSE77791 | Burn shock and hydrocortisone intervention | Affymetrix microarray |

The initial release contains:

- 6 GEO datasets;
- 22 predefined differential-expression contrasts;
- 448,735 gene–contrast result records;
- standardized gene identifiers, metadata, effect sizes, P values, FDR values, and significance labels;
- a static gene-search website;
- reproducible R, Python, and Jupyter analysis code.

## Website

The public website is deployed from the [`docs/`](docs/) directory using GitHub Pages:

**https://peneapple.github.io/2026-BurnOmicsDB/**

The website provides:

- alias-aware gene search;
- gene-centered comparison across tissues and time points;
- log2 fold change, fold change, P value, FDR, direction, and significance status;
- GEO-specific study summaries and analysis figures;
- links to the archived data and code release on Zenodo.

## Repository structure

```text
2026-BurnOmicsDB/
├── README.md
├── LICENSE
├── CITATION.cff
├── Code/
│   ├── GEO-specific R analysis scripts
│   ├── gene-reference processing notebook
│   └── website-input preparation code
└── docs/
    ├── index.html
    ├── 404.html
    ├── .nojekyll
    ├── assets/
    ├── data/
    └── figures/
```

## Code

The [`Code/`](Code/) directory contains:

- independent analysis scripts for all six GEO datasets;
- NCBI human gene-reference and alias-index processing;
- website-input preparation and validation;
- project inventory utilities.

Each GEO dataset was processed independently with a platform-appropriate workflow. RNA-seq and microarray expression matrices from different studies were not directly merged or subjected to a single cross-study batch correction.

## Data availability

The versioned BurnOmicsDB data and code release is archived on Zenodo:

**https://doi.org/10.5281/zenodo.21648675**

Raw GEO source files are not redistributed in this repository or the Zenodo release. Original data can be obtained from NCBI GEO using the accessions listed above.

## Interpretation notes

- Positive `log2FC` consistently indicates higher expression in the designated case group than in the designated control group.
- Significant differential expression is defined using both `FDR < 0.05` and `|log2FC| ≥ 1`.
- Bulk skin and blood transcriptomic differences may reflect both intracellular regulation and changes in cellular composition or tissue structure.
- GSE19743 and GSE37069 contain substantially overlapping cohorts and should not be treated as independent replication datasets.
- BurnOmicsDB is intended for research and hypothesis generation, not as a diagnostic or clinical decision-making tool.

## Citation

Please cite the archived dataset:

> Geng, P. X. (2026). *BurnOmicsDB: a curated and standardized database of human burn injury transcriptomics across tissues, time points, and clinical contexts* (Version 1.0.0) [Dataset]. Zenodo. https://doi.org/10.5281/zenodo.21648675

Machine-readable citation metadata are provided in [`CITATION.cff`](CITATION.cff).

## Authors

**Peter X. Geng**

- Peking University
- University of Colorado Boulder


## Contact

For questions about the database or repository:

**Peter X. Geng**  
Peking University  
Email: peter_geng@stu.pku.edu.cn

## License

- Software source code is released under the **MIT License**.
- Data tables, metadata, figures, and documentation are released under the **Creative Commons Attribution 4.0 International License (CC BY 4.0)**.

See [`LICENSE`](LICENSE) for details.
