---
output:
  html_document: default
  pdf_document: default
---
# Replication Package

## Forecasting U.S. State-Level Volatilities with Local and Global EPUs: A Mixed-Frequency Shrinkage Combination Approach

Candila V., Cepni O., Gallo G.M. & R. Gupta (2026),  
*International Journal of Forecasting.*

---

## 📅 Assembly Date and Contact

**Package assembled:** 27 August 2026  
**Contact:** Vincenzo Candila — vcandila@unisa.it  
*(Please reach out for any questions about the code or data)*

---

## 🗂️ Repository Structure

```
epu_shrinkage_combination/
├── BBG Tickers.xlsx                      # Bloomberg tickers for the 50 state-level stock market indices
├── README.md
├── README.html
├── comb_functions.R                      # core functions
└── main_reproducibility.R                # Replicates the entire workflow

```

---

## 💻 Computing Environment

- **Language:** R version 4.6.1
- **Operating system:** Windows 11 (code is platform-independent)

### Required R packages

| Package         | Version (used) | Purpose |
|-----------------|----------------|---------|
| `rumidas`       | 0.1.3          | GARCH-MIDAS estimation via `ugmfit` |
| `xts`           | 0.14.2         | Time series management |
| `zoo`           | 1.8.13         | Time series infrastructure |
| `fBasics`       | 4052.98        | Summary statistics |
| `rugarch`   | 1.5.6              | Model Confidence Set (MCS) procedure via `mcsTest` and GJR estimation |
| `np`            | 0.70.5         | Automatic block length selection via `b.star` for the MCS bootstrap |
| `DT`            | 0.34           | Interactive HTML tables |
| `roll`          | 1.2.1          | Rolling-window computations |
| `maxLik`        | 1.5.2.2        | Maximum likelihood estimation utilities |
| `highfrequency` | 1.0.3          | High-frequency financial data utilities |
| `xtable`        | 1.8.8          | LaTeX table generation |
| `datasets`      | 4.6.1          | Built-in R datasets, including U.S. state names |
| `lubridate`     | 1.9.5          | Date and time manipulation |
| `DescTools`     | 0.99.60        | Winsorization via `Winsorize` |
| `readxl`        | 1.5.0          | Import of Excel data files |
| `splm`          | 1.6.5          | U.S. state spatial adjacency data |
| `caret`         | 7.0.1          | Time-series validation for Elastic Net tuning |
| `glmnet`        | 5.0            | Elastic Net estimation |                        |



Install all packages at once with:

```r
install.packages(
  c(
    "rumidas",
    "xts",
    "zoo",
    "fBasics",
    "rugarch",
    "np",
    "DT",
    "roll",
    "maxLik",
    "highfrequency",
    "xtable",
    "lubridate",
    "DescTools",
    "readxl",
    "splm",
    "caret",
    "glmnet"
  )
)

```

No environment manager (e.g., `renv`) was used. To record a snapshot of your installed versions, run `sessionInfo()` in R after installing the packages.


---

## 💾 Data

The state-level stock market data used in this reproducibility check were obtained from Bloomberg (https://www.bloomberg.com/professional/) and are subject to licensing and copyright restrictions. Therefore, these data cannot be redistributed as part of the replication package and must be obtained independently by users with access to Bloomberg.

The dataset contains daily closing prices for the state-level capitalization-weighted stock market indices used in the paper, covering the 50 U.S. states. Due to differences in data availability, the starting date of the series varies across states, while the sample ends in December 2024.

The file `BBG Tickers.xlsx`, included in this repository, reports the Bloomberg tickers for the 50 state-level stock market indices used in the paper. Users with access to Bloomberg can use these tickers to retrieve the corresponding daily closing-price series.

For full replication, the downloaded series should be organized in a file named `close_new.csv` and placed in the working directory. The expected structure is:

- the first column contains dates in `dd/mm/yyyy` format;
- the remaining 50 columns contain the daily closing prices of the state-level indices, ordered alphabetically by state name.

The empirical analysis reported in the paper is based on the Bloomberg data downloaded in March 2025.

Users without access to the proprietary Bloomberg data can still reproduce the parts of the analysis based exclusively on publicly available data, as described below.

---

## 🔄 Workflow

The replication workflow is organized into two main stages:

1. **Replication of Table A.1, Table A.2 (columns 1–14), and Table A.3**
   - These results can be replicated using exclusively publicly available data. All required datasets, together with their sources and download instructions, are reported in `main_reproducibility.R`.

2. **Replication of the full empirical analysis (estimation, forecast combination, and evaluation)**
   - Replication of the full empirical analysis additionally requires the proprietary Bloomberg state-level stock market data described above. The file `BBG Tickers.xlsx`, included in this repository, provides the Bloomberg tickers required to retrieve these data. Users should organize the downloaded data as described in the **Data** section above and place `close_new.csv` in the working directory.
   - The script then performs the complete workflow, including the rolling estimation of the individual models, construction of the forecast combinations, and out-of-sample forecast evaluation.
   - This stage is computationally intensive. Since models are estimated independently across states, users may split the estimation across subsets of states and/or multiple R sessions or machines to reduce the overall computation time.
---

## ⏱️ Hardware and Expected Runtime

| Script         | Hardware used               | Approximate runtime |Note                 |
|:---------------|:----------------------------|:--------------------|:--------------------|
| `main_reproducibility.R` (estimation) | Dedicated workstation (Intel(R) Core(TM) i7-14700K @ 3.40 GHz, 32 GB RAM) | ~4 minutes per rolling step | A state with approximately 200 rolling steps requires about 13 hours of computation. The actual runtime varies across states depending on the number of rolling steps. |
| `main_reproducibility.R` (forecast combination) | Dedicated workstation (Intel(R) Core(TM) i7-14700K @ 3.40 GHz, 32 GB RAM) | ~15–30 seconds per rolling combination step | For a state with approximately 200 rolling steps, the Elastic Net combination stage may require about 1–2 hours of computation. |
| `main_reproducibility.R` (forecast evaluation) | Dedicated workstation (Intel(R) Core(TM) i7-14700K @ 3.40 GHz, 32 GB RAM) | ~30 seconds per MCS and state | A state requires four MCS procedures and therefore approximately 2 minutes of computation. |


---

## 📜 Licence

The code in this repository is shared for academic reproducibility purposes. If you use it, please cite the paper:

> Candila V., Cepni O., Gallo G.M. & R. Gupta (2026), "Forecasting U.S. State-Level Volatilities with Local and Global EPUs: A Mixed-Frequency Shrinkage Combination Approach", *International Journal of Forecasting*.
