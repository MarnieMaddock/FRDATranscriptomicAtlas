# FRDATranscriptomicAtlas <img src="/images/dottori_lab_pentagon.svg" alt="Logo" align="right" width="180">

**FRDA Transcriptomic Atlas** is an interactive R/Shiny application for exploring and comparing transcriptomic datasets related to **Friedreich’s ataxia (FRDA)**.  This app was produced as an aid to the publication XXXXXXXXXXXXXXXXXXXXX.

The Atlas integrates multiple publicly available RNA-seq studies and provides harmonised visualisation and analysis of:

- Differential gene expression (DEG)
- Gene set enrichment analysis (GSEA)
- Cross-dataset overlap, concordance, and biomarker discovery

---

## Recommended installation (local)

The FRDA Transcriptomic Atlas runs through R. If you are new to R, please use the detailed installation instructions [here](https://marniemaddock.github.io/FRDATranscriptomicAtlas/).

> ⚠️ **Initial installation may take up to ~30 minutes**  
> This is expected and occurs because the package installs large cached datasets and dependencies required for interactive analyses. These datasets are cached locally for offline reuse, so this long download is only performed once. Subsequent launches are fast and do not require re-downloading data and will load instantly.

```r
install.packages('FRDATranscriptomicAtlas', repos = c('https://marniemaddock.r-universe.dev', 'https://cloud.r-project.org'))
```

## To run app

```
library(FRDATranscriptomicAtlas)
run_app()
```

[![r-universe version](https://marniemaddock.r-universe.dev/FRDATranscriptomicAtlas/badges/version)](https://marniemaddock.r-universe.dev/FRDATranscriptomicAtlas)
[![r-universe status](https://marniemaddock.r-universe.dev/FRDATranscriptomicAtlas/badges/checks)](https://marniemaddock.r-universe.dev/FRDATranscriptomicAtlas)


## Feedback and Support
If you encounter any issues or have suggestions, feel free to:

- Open an issue on this repository
- [Email Us](mlm715@uowmail.edu.au)

  
## License
FRDATranscriptomicAtlas is licensed under the MIT License. See [LICENSE](https://github.com/MarnieMaddock/FRDATranscriptomicAtlas/blob/main/LICENSE) for details.

---- 

![Footer](/images/footer.svg)
