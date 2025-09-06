# Single-Trait Random Regression (RRM) with BLUPF90

This repository shows how to fit a **single-trait Random Regression Model (RRM)** for a longitudinal secondary trait (e.g., canopy coverage across days) using **BLUPF90** from R.

The tutorial script:

* downloads BLUPF90 binaries,
* sources helper functions,
* **loads and prepares phenotypes (including a nearest-neighbor yield covariate, `covYLD`)**,
* defines the RRM (Legendre basis),
* runs the model,
* reads heritability and additive variance profiles over time.

---

## Quick start
Download the repository
```bash
git clone https://github.com/leonardosgloria/random_regression_blupf90.git
```

## Repository layout

```
.
├── R/                      # helper R functions (blup(), download_BLUPF90(), etc.)
├── data/
│   ├── pheno_RRM.txt       # phenotypes (see Section 4)
│   └── genotype_BLUPF90.txt
├── output/                 # model outputs (h2.txt, VC_Time_var.txt, solutions, logs)
└── scripts/
    └── run_single_rrm.R    # tutorial script
```

---

## Requirements

* **R ≥ 4.1**
  R packages: `pacman`, `orthopolynom`, `splines`, `dplyr`, `tidyr`, `stringr`, `data.table`, `readr`, `purrr`, `gtools`, **`NAM`**
* **BLUPF90** binaries (downloaded automatically by the script)
* A CPU with multiple threads (optional but recommended)

Install R deps (first run will auto-install via `pacman`):

```r
if (!require("pacman")) install.packages("pacman")
pacman::p_load(orthopolynom, splines, dplyr, tidyr, stringr,
               data.table, readr, purrr, gtools, NAM)
```

---

## Input data expected

* **Phenotypes** (`data/pheno_RRM.txt`): a long table with at least

  * `Geno`  — genotype/line ID
  * `Day`   — time index (e.g., days after planting)
  * `AdjCC` — response for RRM (e.g., adjusted canopy coverage)
  * `YLD`   — grain yield (used to build `covYLD`)
  * `Block`, `Row`, `Col` — field layout (for nearest-neighbor map)
  * Optional fixed covariates, e.g., `IntBlk` and an environment indicator `int`

* **Genotypes**: `data/genotype_BLUPF90.txt` in BLUPF90 format for building the genomic relationship matrix.

---

## What the model fits

### Model formula (R wrapper)

```r
model <- AdjCC ~ int/RRM + IntBlk + ped|RRM|Geno
```

**Interpretation in this wrapper:**

* `int/RRM` — random regression for the **population mean curve** over time (Legendre polynomials) with environment factor `int`.
* `IntBlk` — (optional) fixed effect / covariate (e.g., spatial or block).
* `ped|RRM|Geno` — **additive genetic random regression** for each genotype (coefficients linked via genomic K).

### Random regression basis & time window

```r
RRM_option1 <- list(
  poly   = 2,       # Legendre order (2=quadratic; try 3 or 4 if needed)
  Timevar = "Day",  # time column
  Pmin   = 17,      # min Day for summaries
  Pmax   = 73       # max Day for summaries
)
```

---

## Preparing the data (`datarenum1`) — **includes your missing block**

Add this block near the top of your script (before model fitting):

```r
# --- Load phenotypes and create nearest-neighbor yield covariate ---
pheno_RRM <- data.table::fread("data/pheno_RRM.txt", header = TRUE, data.table = FALSE)

# Use as working dataset for the wrapper
datarenum1 <- pheno_RRM

# Factors expected by the model
datarenum1$Geno  <- factor(datarenum1$Geno)
datarenum1$IntBlk <- factor(datarenum1$IntBlk)
datarenum1$Block <- factor(datarenum1$Block)
datarenum1$int   <- factor(datarenum1$int)
```

---

## Running the script (core lines)

```r
# Get BLUPF90 locally
download_BLUPF90(update = TRUE) # download the binaries to the folder: paste0(.libPaths()[1], "/blupf90"), you can use the binary in the repository, copy and paste to this folder

# Source helpers
R_script_folder <- "R"
files <- list.files(R_script_folder, pattern = "[.][Rr]$", full.names = TRUE, recursive = TRUE)
invisible(lapply(files, function(f) try(source(f, chdir = TRUE), silent = TRUE)))

# Model pieces
model <- AdjCC ~ int/RRM + IntBlk + ped|RRM|Geno

residual_start1 <- matrix(c(0.047))    # residual start
VC_start1 <- matrix(c(
  0.2404,  0.05354,  -0.2300,
  0.05354, 0.01589,  -0.05394,
 -0.2300, -0.05394,   0.2220
), nrow = 3, byrow = TRUE) # additive genetic effect start

RRM_option1 <- list(poly = 2, Timevar = "Day", Pmin = 17, Pmax = 73)

genotype_file1   <- "data/genotype_BLUPF90.txt"
ped_name1        <- NULL
PED_DEPTH1       <- 3
missing_values1  <- -99

# Fit single-trait RRM
model_single <- blup(
  datarenum        = datarenum1,
  formula          = model,
  fields_output    = NULL,
  weights_object   = NULL,
  residual_start   = residual_start1,
  VCA_RRM          = VC_start1,
  ped_name         = ped_name1,
  PED_DEPTH        = PED_DEPTH1,
  genotype_file    = genotype_file1,
  missing_values   = missing_values1,
  RRM_option       = RRM_option1,
  het_res_variance = NULL,
  fit_option = list(
    yams           = TRUE,
    solution_mean  = TRUE,
    VCE            = TRUE,
    sol_se         = TRUE,
    Inbreeding     = TRUE,
    alpha_size     = 30,
    EM_REML        = 1,
    maxrounds      = 3000,
    alpha_beta     = c(0.95, 0.05),
    tunedG         = 0,
    conv_crit      = 1e-10
  ),
  run_model = list(n_threads = 20),
  keep_files = TRUE
)

# Read key outputs
library(data.table)
h2 <- fread("output/h2.txt") %>%
  dplyr::filter(Time_var >= RRM_option1$Pmin & Time_var <= RRM_option1$Pmax)

vg <- fread("output/VC_Time_var.txt") %>%
  dplyr::filter(Time_var >= RRM_option1$Pmin & Time_var <= RRM_option1$Pmax) %>%
  dplyr::select(Time_var, Additive_variance = Geno)

```

---
## Bi-trait model — MT-RRM: longitudinal **AdjCC** + single-time **YLD**

This section shows how to fit a **multi-trait random regression model** where:

* **AdjCC** is modeled with **random regression** (Legendre polynomials) across `Day`;
* **YLD** is modeled as a **single-time** additive effect (no random regression).

> ⚠️ **Important:** Do **not** duplicate YLD across days. Keep **one YLD per Genotype **.

### Data prep (nearest-neighbor covariate + factor setup)

```r
# --- Load phenotypes ---
pheno_RRM <- data.table::fread("data/pheno_RRM.txt", header = TRUE, data.table = FALSE)
datarenum1 <- pheno_RRM

# Factors expected by the model
datarenum1$Geno   <- factor(datarenum1$Geno)
datarenum1$IntBlk <- factor(datarenum1$IntBlk)
datarenum1$Block  <- factor(datarenum1$Block)
datarenum1$int    <- factor(datarenum1$int)

# ---- Nearest-neighbor map (optional helper covariate; not used in the formulas below) ----
SP  <- data.frame(Block = datarenum1$Block, Row = datarenum1$Row, Col = datarenum1$Col)
MAP <- NAM::NNsrc(SP, 2.2, 1)      # tune radii as needed

# Clean YLD, compute smoothed covariate if you want to use it elsewhere
datarenum1$YLD <- replace(datarenum1$YLD, which(datarenum1$YLD < 0), NA)
datarenum1$covYLD <- NAM::NNcov(MAP, datarenum1$YLD)

# ---- Sanity check: YLD should not be replicated across days ----
chk <- dplyr::count(datarenum1, Geno, int, wt = !is.na(YLD), name = "nYLD_nonNA")
stopifnot(all(chk$nYLD_nonNA <= 1))   # 1 per Geno×int (or 0 if missing), never >1
```

### Model formulas (as a list)

* First formula = AdjCC with **RRM**;
* Second formula = YLD as **single-time**.

```r
model <- list(AdjCC ~   int/RRM + IntBlk + ped|RRM|Geno,
              YLD ~   IntBlk +covYLD+ ped|Geno
              )

)
```

### Residual and genetic start values

These are **starting values** 

```r
# Residual start: diag for [AdjCC, YLD]
residual_start1 <- matrix(c(0.0400, 0.0000,
                            0.0000, 90.0000),
                          nrow = 2, byrow = TRUE)

# Genetic start matrix (example provided by you; keep order consistent with your wrapper)
VC_start1 <- matrix(c(
  0.2528,    4.0250,   0.004715, 0.0000,  -0.2369, 0.0000,
  4.0250,   64.4200,   0.758600, 0.0000,  -3.7720, 0.0000,
  0.004715,  0.758600, 0.001352, 0.0000,  -0.04739,0.0000,
  0.0000,    0.0000,   0.0000,   0.0000,   0.0000, 0.0000,
 -0.2369,   -3.7720,  -0.04739,  0.0000,   0.2246, 0.0000,
  0.0000,    0.0000,   0.0000,   0.0000,   0.0000, 0.0000
), nrow = 6, byrow = TRUE)

genotype_file1 <- "genotype_BLUPF90.txt"
ped_name1 <- NULL
PED_DEPTH1 <- 3
fiels_output_pos1 <- NULL
weights_object1<- NULL
```

### Fit the bi-trait model

```r
RRM_option1 <- list(poly = 2, Timevar = "Day", Pmin = 17, Pmax = 73)

model_bi <- blup(datarenum=datarenum1,formula = model,fields_output=fiels_output_pos1,weights_object=weights_object1,
               residual_start=residual_start1,VCA_RRM=VC_start1,ped_name=ped_name1,
               PED_DEPTH=PED_DEPTH1,genotype_file=genotype_file1,
               missing_values=missing_values1,
               RRM_option=RRM_option1,het_res_variance=NULL,
               fit_option=list(yams=T,solution_mean=T,
                               VCE=T,sol_se=T,Inbreeding=T,
                               alpha_size=30,
                               EM_REML=1000,
                               maxrounds=3000,
                               alpha_beta=c(0.95,0.05),
                               tunedG=0,
                               conv_crit=1e-10),
              run_model=list(n_threads=20)
              ,keep_files=T)
```

### Reading multi-trait outputs


```r
# Heritability over time for SLT (AdjCC)
h2_slt <- data.table::fread("output/h2.txt") |>
  dplyr::filter(Time_var >= RRM_option1$Pmin & Time_var <= RRM_option1$Pmax)

# Additive genetic variance over time for SLT
vg_slt <- data.table::fread("output/VC_Time_var.txt") |>
  dplyr::filter(Time_var >= RRM_option1$Pmin & Time_var <= RRM_option1$Pmax) |>
  dplyr::select(Time_var, Additive_variance = Geno)
# Correlation

```

