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


## Feedback and Support
If you encounter any issues or have suggestions, feel free to:

- Open an issue on this repository: [Bug Reports/Feature Requests](https://github.com/MarnieMaddock/FRDATranscriptomicAtlas/issues)
- [Email Us](mlm715@uowmail.edu.au)

  
## License
FRDATranscriptomicAtlas is licensed under the MIT License. See [LICENSE](https://github.com/MarnieMaddock/FRDATranscriptomicAtlas/blob/main/LICENSE) for details.

---- 

![Footer](/images/footer.svg)
