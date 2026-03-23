# FRDATranscriptomicAtlas <img src="/images/dottori_lab_pentagon.svg" alt="Logo" align="right" width="180">

**FRDA Transcriptomic Atlas** is an interactive R/Shiny application for exploring and comparing transcriptomic datasets related to **Friedreich’s ataxia (FRDA)**.  This app was produced as an aid to the publication XXXXXXXXXXXXXXXXXXXXX.

The Atlas integrates multiple publicly available RNA-seq studies. Features include:

- Principal Component Analysis (PCA) - Inspect global sample structure and separation by condition, tissue/cell type, or other metadata.
- Differential Expression (DEG/DEI) Tables - Gene- and isoform-level results with log2 fold-change and adjusted p-values, filterable by direction and significance.
- Volcano Plots - Fold-change vs. significance visualisation with options to highlight genes and adjust thresholds.
- Venn Diagrams - Overlap of differentially expressed genes or isoforms across datasets.
- Expression Heatmaps - Compare transcript abundance patterns across conditions and datasets.
- Gene TPM Plots (Gene Plots) - Visualise expression of selected genes (TPM) by condition, replicate, and dataset for direct comparison.
- Cross-Dataset Forest Plots - Compare per-dataset effect sizes for selected genes to assess concordance and heterogeneity across studies.
- Download and Export - Tables and plots can be exported for reporting and reuse.

---

## Recommended installation (local)

The FRDA Transcriptomic Atlas runs through R. If you are new to R, please use the detailed installation instructions [here](https://marniemaddock.github.io/FRDATranscriptomicAtlas/). 

> ⚠️ **Initial installation may take up to ~15 minutes**  
> This is expected and occurs because the package installs large cached datasets and dependencies required for interactive analyses. These datasets are cached locally for offline reuse, so this long download is only performed once. Subsequent launches are fast and do not require re-downloading data and will load instantly. The app is ~ 600 MB.

```r
install.packages('FRDATranscriptomicAtlas', repos = c('https://marniemaddock.r-universe.dev', 'https://cloud.r-project.org'))
FRDATranscriptomicAtlas::install_deps()
```

### To run app

```
library(FRDATranscriptomicAtlas)
run_app()
```

## New to R? No problem.

The FRDA Transcriptomic Atlas is designed for users with no programming experience. All analyses and visualisations are accessible through an intuitive point-and-click interface. Trust me - you do not need to know how to use R!

[![r-universe version](https://marniemaddock.r-universe.dev/FRDATranscriptomicAtlas/badges/version)](https://marniemaddock.r-universe.dev/FRDATranscriptomicAtlas)
[![r-universe status](https://marniemaddock.r-universe.dev/FRDATranscriptomicAtlas/badges/checks)](https://marniemaddock.r-universe.dev/FRDATranscriptomicAtlas)

## Citation and Acknowledgements

We thank all researchers who generated and made publicly available the original datasets used in this meta-analysis. If you use this Atlas or any underlying data/analyses in your work, please cite the associated publication:

Maddock, M. et al. XXXXXXXXXXXXXXXXXX. [Journal details forthcoming].

The manuscript is currently on **BioRxiv**!

Friedreich ataxia transcriptomic dysregulation and identification of cell type-specific biomarkers: A systematic review and meta-analysis
Marnie L Maddock, Sara Miellet, Anjila Dongol, Amy J Hulme, Chloe K Kennedy, Louise A Corben, Rocio K Finol-Urdaneta, Alberto Nettel-Aguirre, Chiara Dionsi, Martin B Delatycki, Joel M Gottesfeld, Massimo Pandolfo, Elisabetta Soragni, Sanjay I Bidichandani, Jarmon G Lees, Shiang Y Lim, Jill S Napierala, Marek Napierala, Mirella Dottori
bioRxiv 2026.03.18.712785; doi: [https://doi.org/10.64898/2026.03.18.712785]

Users must cite the publication in any derivative analyses, figures, or reports generated from this resource. Users must also cite the original studies from which the data were derived. Please refer to the 'Datasets' section within the app for details on each dataset and its original publication.

## Future Updates and Contributions

This resource may be updated as new FRDA RNA-seq datasets become available. Researchers who wish to contribute data or suggest additions are encouraged to [contact the authors via email](mailto:mlm715@uowmail.edu.au,mmaddock@uow.edu.au,mdottori@uow.edu.au?subject=FRDA%20Transcriptomic%20Atlas%20Contribution).

## Feedback and Support
If you encounter any issues or have suggestions, feel free to:

- Open an issue on this repository: [Bug Reports/Feature Requests](https://github.com/MarnieMaddock/FRDATranscriptomicAtlas/issues)
- [Email Us](mlm715@uowmail.edu.au)

  
## License
FRDATranscriptomicAtlas is licensed under the MIT License. See [LICENSE](https://github.com/MarnieMaddock/FRDATranscriptomicAtlas/blob/main/LICENSE) for details.

---- 

![Footer](/images/footer.svg)
