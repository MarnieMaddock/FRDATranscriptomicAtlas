# FRDATranscriptomicAtlas <img src="/images/dottori_lab_pentagon.svg" alt="Logo" align="right" width="180">

**FRDA Transcriptomic Atlas** is an interactive R/Shiny application for exploring and comparing transcriptomic datasets related to **Friedreich’s ataxia (FRDA)**.  This app was produced as an aid to the publication XXXXXXXXXXXXXXXXXXXXX.

The Atlas integrates multiple publicly available RNA-seq studies and provides harmonised visualisation and analysis of:

- Differential gene expression (DEG)
- Differential transcript usage (DTU)
- Isoform switching and predicted functional consequences
- Gene set enrichment analysis (GSEA)
- Cross-dataset overlap, concordance, and biomarker discovery

---

## Recommended installation (local)

> ⚠️ **Initial installation may take up to ~30 minutes**  
> This is expected and occurs because the package installs large cached datasets and dependencies required for interactive analyses.

```r
install.packages('FRDATranscriptomicAtlas', repos = c('https://marniemaddock.r-universe.dev', 'https://cloud.r-project.org'))
```

These datasets are cached locally for offline reuse, so this long download is only performed once. Subsequent launches are fast and do not require re-downloading data.

# To run app

```
library(FRDATranscriptomicAtlas)
run_app()
```


[![r-universe version](https://r-lib.r-universe.dev/FRDATranscriptomicAtlas/badges/version)](https://r-lib.r-universe.dev/FRDATranscriptomicAtlas)

[![r-universe status](https://r-lib.r-universe.dev/FRDATranscriptomicAtlas/badges/checks)](https://r-lib.r-universe.dev/FRDATranscriptomicAtlas)



## Feedback and Support
If you encounter any issues or have suggestions, feel free to:

- Open an issue on this repository
- [Email Us](mlm715@uowmail.edu.au)

  
## License
FRDATranscriptomicAtlas is licensed under the MIT License. See [LICENSE](https://github.com/MarnieMaddock/FRDATranscriptomicAtlas/blob/main/LICENSE) for details.

---- 

![Footer](/images/footer.svg)
