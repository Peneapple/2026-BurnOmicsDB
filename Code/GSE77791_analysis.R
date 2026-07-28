# ==============================================================================
# BurnOmicsDB: GSE77791 microarray analysis
# BurnOmicsDB：GSE77791微阵列分析
#
# GEO accession / GEO编号:
#   GSE77791
#
# Study / 研究:
#   Prospective randomized double-blind assessment of transcriptome modulation
#   by hydrocortisone in severe burn shock.
#   严重烧伤休克中氢化可的松对转录组影响的前瞻性随机双盲研究。
#
# Data type / 数据类型:
#   Affymetrix Human Genome U133 Plus 2.0 Array (GPL570)
#   Raw CEL files contained in GSE77791_RAW.tar
#   GPL570 Affymetrix原始CEL文件。
#
# Study structure / 研究结构:
#   - 117 arrays in total.
#   - 13 healthy-volunteer arrays from 13 independent volunteers.
#   - 104 burn-patient arrays from 30 randomized patients.
#   - Hydrocortisone: 15 patients, 48 arrays.
#   - Placebo: 15 patients, 56 arrays.
#   - Burn time points:
#       S1: before treatment, n = 30
#       S2: approximately 24 h after treatment initiation, n = 27
#       S3: approximately 120 h after treatment initiation, n = 29
#       S4: approximately 168 h after treatment initiation, n = 18
#
#   共117张芯片：13名健康志愿者，以及30名随机分组烧伤患者的104张纵向芯片。
#
# Seven prespecified BurnOmicsDB contrasts / 七个预设比较:
#   Baseline burn-shock effect:
#   1. Burn shock S1 versus healthy volunteers
#
#   Natural longitudinal change in the placebo group:
#   2. Placebo S2 versus Placebo S1
#   3. Placebo S3 versus Placebo S1
#   4. Placebo S4 versus Placebo S1
#
#   Hydrocortisone effect at matched study time points:
#   5. Hydrocortisone versus Placebo at S2
#   6. Hydrocortisone versus Placebo at S3
#   7. Hydrocortisone versus Placebo at S4
#
# Contrast direction / 比较方向:
#   Positive log2FC always means higher expression in Case_Group than in
#   Control_Group.
#   正log2FC始终表示Case_Group相对于Control_Group表达更高。
#
# Two-model rationale / 两个模型的理由:
#   Model A uses only 30 pretreatment S1 arrays and 13 healthy-volunteer arrays.
#   Healthy-volunteer age and sex are unavailable in GEO, so this comparison is
#   intentionally not adjusted for age, sex, or TBSA.
#
#   Model B uses all 104 burn-patient arrays. It adjusts for age, sex, and TBSA,
#   and uses duplicateCorrelation with Patient_ID blocking to account for
#   repeated samples from the same patient.
#
#   模型A比较治疗前S1与健康志愿者；由于健康志愿者人口学信息缺失，不校正协变量。
#   模型B使用104张患者芯片，校正年龄、性别和TBSA，并处理患者内重复测量。
#
# Treatment-effect interpretation / 治疗效应解释:
#   Hydrocortisone-versus-placebo contrasts are adjusted comparisons at S2, S3,
#   and S4. They are not difference-in-differences contrasts.
#   S1 hydrocortisone-versus-placebo is calculated only as a randomization
#   baseline-balance QC check and is not exported as a core website contrast.
#   氢化可的松与安慰剂在S2、S3、S4进行同时间点校正比较；
#   S1组间差异仅作为随机化基线QC，不进入数据库核心结果。
#
# Sampling-time interpretation / 时间解释:
#   S1 is the onset of burn shock and occurs before trial treatment, generally
#   24-72 h after injury. S2-S4 are defined relative to treatment initiation,
#   not directly relative to the burn event.
#   S1为休克入组和治疗前时间点；S2-S4相对于治疗开始定义，而非直接相对于烧伤发生。
#
# Analysis workflow / 分析流程:
#   GSE77791_RAW.tar
#   -> validate and extract 117 CEL.gz files
#   -> decompress CEL files
#   -> batch-wise raw-array QC
#   -> RMA background correction, quantile normalization, and summarization
#   -> checkpoint the RMA probe-level matrix
#   -> GPL570 probe-to-NCBI-Gene-ID mapping
#   -> remove probes mapping to zero or multiple NCBI Gene IDs
#   -> select one representative probe per NCBI Gene ID using the highest mean
#      RMA expression across all 117 arrays
#   -> Model A: limma, Burn S1 versus healthy volunteer
#   -> Model B: limma + duplicateCorrelation, adjusted for age, sex, and TBSA
#   -> seven prespecified core contrasts
#
# Probe-selection rule / 代表探针规则:
#   For each unambiguous NCBI Gene ID, select the probe with the highest mean
#   RMA expression across all 117 arrays before differential testing.
#   The rule does not use group labels, p-values, or FDR.
#   对每个明确Gene ID选择全部117张芯片平均表达最高的探针，不依据统计显著性。
#
# Relation to the original publications / 与原论文的关系:
#   The primary paper used GCRMA, MAS5 detection filtering, COMBAT, and several
#   longitudinal tests. BurnOmicsDB uses a uniform RMA + limma workflow and
#   therefore does not claim exact numerical reproduction.
#
#   A secondary publication studied HERV-targeting GPL570 probes. BurnOmicsDB is
#   a standard human-gene resource; probes that cannot be mapped unambiguously
#   to one NCBI Gene ID are excluded from the core gene table.
#   原论文使用GCRMA、MAS5和COMBAT；本数据库使用统一RMA流程，不声称精确复现。
#   HERV专用解释不进入标准Gene ID核心表。
#
# Batch-effect policy / 批次效应策略:
#   The original paper used COMBAT, but the uploaded SOFT metadata does not
#   provide a reliable sample-level Batch_ID. This script does not infer a batch
#   variable or apply COMBAT without documented batch labels. PCA, correlation,
#   raw intensity, RMA distribution, and RLE metrics are exported for QC.
#   原论文使用COMBAT，但SOFT没有可靠Batch_ID，因此本代码不推测批次变量。
#
# Missing-data and evidence limitations / 缺失与证据限制:
#   - Healthy-volunteer age and sex are unavailable.
#   - S4 includes only 7 hydrocortisone and 11 placebo arrays.
#   - Longitudinal missingness may be related to death, illness severity, or
#     technical availability and may not be completely random.
#   - Whole-blood bulk expression can reflect changes in leukocyte composition.
#   - No sample is deleted automatically by this script.
#
# How to run in RStudio / 如何在RStudio中运行:
#   - Save this script in:
#     /Users/peter/Downloads/Project-2026-BurnOmicsDB/GSE77791/
#   - Keep these local files in the project folder:
#       GSE77791_RAW.tar
#       GSE77791_family.soft.gz
#       13054_2017_Article_1743.pdf
#       fimmu-09-03091.pdf
#   - Press Cmd + Shift + O on macOS to open the section outline.
#   - Run one section at a time using Cmd + Enter.
#   - During the first run, do not Source the entire script.
#
# Output-language rule / 输出语言规则:
#   Comments are bilingual. Figures, CSV files, TXT files, console messages,
#   and error messages are English only.
#   注释使用中英文；图片、CSV、TXT、Console信息和报错全部只使用英文。
# ==============================================================================


# ---- 00. Install required packages once / 首次安装所需软件包 ----

options(timeout = 1800)
options(download.file.method = "libcurl")
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Change this mirror if necessary.
# 如当前网络不可用，可以替换Bioconductor镜像。
options(
  BioC_mirror =
    "https://mirrors.tuna.tsinghua.edu.cn/bioconductor"
)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages(
    "BiocManager",
    repos = "https://cloud.r-project.org"
  )
}

cran_packages <- c(
  "ggplot2",
  "ggrepel",
  "pheatmap",
  "R.utils",
  "matrixStats",
  "statmod"
)

missing_cran <- cran_packages[
  !vapply(
    cran_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_cran) > 0) {
  install.packages(
    missing_cran,
    repos = "https://cloud.r-project.org",
    dependencies = TRUE
  )
}

bioc_packages <- c(
  "affy",
  "limma",
  "AnnotationDbi",
  "hgu133plus2cdf",
  "hgu133plus2.db",
  "org.Hs.eg.db"
)

missing_bioc <- bioc_packages[
  !vapply(
    bioc_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_bioc) > 0) {
  BiocManager::install(
    missing_bioc,
    ask = FALSE,
    update = FALSE
  )
}

all_packages <- c(
  cran_packages,
  bioc_packages
)

package_check <- vapply(
  all_packages,
  requireNamespace,
  logical(1),
  quietly = TRUE
)

print(package_check)

if (!all(package_check)) {
  stop(
    paste0(
      "Package installation is incomplete. Missing packages: ",
      paste(
        names(package_check)[!package_check],
        collapse = ", "
      )
    )
  )
}

cat("\nAll required packages are installed successfully.\n\n")


# ---- 01. Project settings and package loading / 项目设置与软件包加载 ----

GEO_ID <- "GSE77791"

PROJECT_DIR <-
  "/Users/peter/Downloads/Project-2026-BurnOmicsDB/GSE77791"

FDR_CUTOFF <- 0.05
LOG2FC_CUTOFF <- 1
TOP_VARIABLE_GENES_FOR_PCA <- 500
TOP_GENES_FOR_HEATMAP <- 30
VOLCANO_LABEL_GENE_N <- 12
RAW_QC_BATCH_SIZE <- 10
RLE_BATCH_SIZE <- 40
WRITE_PROBE_LEVEL_CSV <- TRUE
FORCE_RMA_RECOMPUTE <- FALSE

RESULTS_DIR <- file.path(PROJECT_DIR, "Results")
FIGURES_DIR <- file.path(PROJECT_DIR, "Figures")
QC_DIR <- file.path(PROJECT_DIR, "QC")
OBJECTS_DIR <- file.path(PROJECT_DIR, "R_objects")
INPUT_DIR <- file.path(PROJECT_DIR, "Input")
CEL_GZ_DIR <- file.path(INPUT_DIR, "CEL_gz")
CEL_DIR <- file.path(INPUT_DIR, "CEL")

dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(QC_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OBJECTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(INPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(CEL_GZ_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(CEL_DIR, recursive = TRUE, showWarnings = FALSE)

required_packages <- c(
  "ggplot2",
  "ggrepel",
  "pheatmap",
  "R.utils",
  "matrixStats",
  "statmod",
  "affy",
  "limma",
  "AnnotationDbi",
  "hgu133plus2cdf",
  "hgu133plus2.db",
  "org.Hs.eg.db"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  stop(
    paste0(
      "Missing packages: ",
      paste(missing_packages, collapse = ", "),
      ". Run Section 00 before continuing."
    )
  )
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(R.utils)
  library(matrixStats)
  library(statmod)
  library(affy)
  library(limma)
  library(AnnotationDbi)
  library(hgu133plus2cdf)
  library(hgu133plus2.db)
  library(org.Hs.eg.db)
})

set.seed(2026)

STUDY_GROUP_LEVELS <- c(
  "Healthy_volunteer",
  "Placebo",
  "Hydrocortisone"
)

STUDY_GROUP_COLORS <- c(
  "Healthy_volunteer" = "#0072B2",
  "Placebo" = "#009E73",
  "Hydrocortisone" = "#D55E00"
)

SAMPLE_TIME_LEVELS <- c(
  "HV",
  "S1",
  "S2",
  "S3",
  "S4"
)

SAMPLE_TIME_COLORS <- c(
  "HV" = "#0072B2",
  "S1" = "#D55E00",
  "S2" = "#009E73",
  "S3" = "#CC79A7",
  "S4" = "#E69F00"
)

SAMPLE_TIME_SHAPES <- c(
  "HV" = 4,
  "S1" = 16,
  "S2" = 17,
  "S3" = 15,
  "S4" = 18
)

VOLCANO_COLORS <- c(
  "Up_significant" = "#D55E00",
  "Down_significant" = "#0072B2",
  "Not_significant" = "#BDBDBD"
)

CONTRAST_GROUP_COLORS <- c(
  "Control" = "#0072B2",
  "Case" = "#D55E00"
)

EXPRESSION_HEATMAP_COLORS <-
  grDevices::colorRampPalette(
    c("#0072B2", "#F7F7F7", "#D55E00")
  )(101)

cat("Project directory:\n", PROJECT_DIR, "\n\n")


# ---- 02. Locate input files and extract CEL files / 定位输入文件并解压CEL ----

input_search_dirs <- c(
  INPUT_DIR,
  PROJECT_DIR
)

find_one_file <- function(
  pattern,
  search_dirs,
  required = TRUE
) {
  hits <- unlist(
    lapply(
      search_dirs,
      function(directory) {
        if (!dir.exists(directory)) {
          return(character(0))
        }

        list.files(
          path = directory,
          pattern = pattern,
          full.names = TRUE,
          recursive = FALSE,
          ignore.case = TRUE
        )
      }
    ),
    use.names = FALSE
  )

  hits <- unique(hits)

  if (length(hits) == 0) {
    if (required) {
      stop(
        paste0(
          "No input file matched this pattern: ",
          pattern
        )
      )
    }

    return(NA_character_)
  }

  if (length(hits) > 1) {
    message(
      "Multiple files matched. The first file will be used:\n",
      paste(hits, collapse = "\n")
    )
  }

  normalizePath(
    hits[1],
    mustWork = TRUE
  )
}

raw_tar_file <- find_one_file(
  "^GSE77791_RAW.*\\.tar$",
  input_search_dirs
)

soft_file <- find_one_file(
  "^GSE77791_family\\.soft.*\\.gz$",
  input_search_dirs
)

primary_paper_file <- find_one_file(
  "13054_2017_Article_1743.*\\.pdf$",
  input_search_dirs,
  required = FALSE
)

secondary_paper_file <- find_one_file(
  "fimmu-09-03091.*\\.pdf$",
  input_search_dirs,
  required = FALSE
)

cat("Raw CEL archive:\n", raw_tar_file, "\n\n")
cat("SOFT metadata:\n", soft_file, "\n\n")
cat("Primary paper PDF:\n", primary_paper_file, "\n\n")
cat("Secondary HERV paper PDF:\n", secondary_paper_file, "\n\n")

tar_members <- utils::untar(
  raw_tar_file,
  list = TRUE
)

cel_gz_members <- tar_members[
  grepl(
    "\\.CEL\\.gz$",
    tar_members,
    ignore.case = TRUE
  )
]

if (length(cel_gz_members) != 117) {
  stop(
    paste0(
      "The raw archive should contain 117 CEL.gz files, but ",
      length(cel_gz_members),
      " were found."
    )
  )
}

existing_cel_gz <- list.files(
  CEL_GZ_DIR,
  pattern = "\\.CEL\\.gz$",
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
)

if (length(existing_cel_gz) != 117) {
  cat("Extracting 117 CEL.gz files from the raw archive.\n")

  utils::untar(
    raw_tar_file,
    exdir = CEL_GZ_DIR
  )
}

cel_gz_files <- list.files(
  CEL_GZ_DIR,
  pattern = "\\.CEL\\.gz$",
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
)

if (length(cel_gz_files) != 117) {
  stop(
    paste0(
      "Expected 117 extracted CEL.gz files, but ",
      length(cel_gz_files),
      " were found."
    )
  )
}

cat("Decompressing CEL files when necessary.\n")

for (cel_gz_file in cel_gz_files) {
  cel_filename <- sub(
    "\\.gz$",
    "",
    basename(cel_gz_file),
    ignore.case = TRUE
  )

  destination_file <- file.path(
    CEL_DIR,
    cel_filename
  )

  if (!file.exists(destination_file)) {
    R.utils::gunzip(
      filename = cel_gz_file,
      destname = destination_file,
      remove = FALSE,
      overwrite = FALSE
    )
  }
}

cel_files <- list.files(
  CEL_DIR,
  pattern = "\\.CEL$",
  full.names = TRUE,
  recursive = FALSE,
  ignore.case = TRUE
)

if (length(cel_files) != 117) {
  stop(
    paste0(
      "Expected 117 decompressed CEL files, but ",
      length(cel_files),
      " were found."
    )
  )
}

# Raw filenames begin with the GSM accession followed by the study sample name.
# 原始文件名以GSM编号开头，后面带研究样本名。
cel_sample_ids <- sub(
  "^(GSM[0-9]+).*",
  "\\1",
  basename(cel_files),
  ignore.case = TRUE
)

if (anyDuplicated(cel_sample_ids) > 0) {
  stop(
    "Duplicated GSM identifiers were detected among the CEL filenames."
  )
}

cat(
  "CEL archive validation completed. ",
  length(cel_files),
  " CEL files are available.\n\n",
  sep = ""
)


# ---- 03. Stream-parse GEO SOFT metadata and validate study design / 流式解析SOFT元数据并检查设计 ----

# The SOFT file embeds a complete probe table for every sample. The parser
# reads sequentially and skips the probe tables.
# SOFT包含每个样本的完整探针表，本函数流式读取并跳过表达矩阵行。
parse_geo_soft_samples_streaming <- function(
  soft_gz_file,
  chunk_size = 50000
) {
  connection <- gzfile(
    soft_gz_file,
    open = "rt"
  )

  on.exit(
    close(connection),
    add = TRUE
  )

  records <- list()
  current <- NULL
  in_sample_table <- FALSE

  finalize_current <- function(current_record) {
    if (is.null(current_record)) {
      return(NULL)
    }

    current_record
  }

  repeat {
    chunk <- readLines(
      connection,
      n = chunk_size,
      warn = FALSE,
      encoding = "UTF-8"
    )

    if (length(chunk) == 0) {
      break
    }

    for (line in chunk) {
      if (grepl("^\\^SAMPLE\\s*=", line)) {
        finalized <- finalize_current(current)

        if (!is.null(finalized)) {
          records[[length(records) + 1]] <- finalized
        }

        current <- list(
          Sample_ID = sub(
            "^\\^SAMPLE\\s*=\\s*",
            "",
            line
          ),
          Characteristics = character(0),
          Relations = character(0)
        )

        in_sample_table <- FALSE
        next
      }

      if (is.null(current)) {
        next
      }

      if (grepl(
        "^!sample_table_begin",
        line,
        ignore.case = TRUE
      )) {
        in_sample_table <- TRUE
        next
      }

      if (grepl(
        "^!sample_table_end",
        line,
        ignore.case = TRUE
      )) {
        in_sample_table <- FALSE
        next
      }

      if (in_sample_table) {
        next
      }

      if (grepl("^!Sample_title\\s*=", line)) {
        current$Original_Title <- sub(
          "^!Sample_title\\s*=\\s*",
          "",
          line
        )
      } else if (grepl(
        "^!Sample_source_name_ch1\\s*=",
        line
      )) {
        current$Original_Source_Name <- sub(
          "^!Sample_source_name_ch1\\s*=\\s*",
          "",
          line
        )
      } else if (grepl(
        "^!Sample_characteristics_ch1\\s*=",
        line
      )) {
        current$Characteristics <- c(
          current$Characteristics,
          sub(
            "^!Sample_characteristics_ch1\\s*=\\s*",
            "",
            line
          )
        )
      } else if (grepl(
        "^!Sample_description\\s*=",
        line
      )) {
        current$Description <- sub(
          "^!Sample_description\\s*=\\s*",
          "",
          line
        )
      } else if (grepl(
        "^!Sample_platform_id\\s*=",
        line
      )) {
        current$Platform_ID <- sub(
          "^!Sample_platform_id\\s*=\\s*",
          "",
          line
        )
      } else if (grepl(
        "^!Sample_supplementary_file\\s*=",
        line
      )) {
        current$Supplementary_File <- sub(
          "^!Sample_supplementary_file\\s*=\\s*",
          "",
          line
        )
      } else if (grepl(
        "^!Sample_relation\\s*=",
        line
      )) {
        current$Relations <- c(
          current$Relations,
          sub(
            "^!Sample_relation\\s*=\\s*",
            "",
            line
          )
        )
      }
    }
  }

  finalized <- finalize_current(current)

  if (!is.null(finalized)) {
    records[[length(records) + 1]] <- finalized
  }

  if (length(records) == 0) {
    stop(
      "No SAMPLE records were parsed from the SOFT file."
    )
  }

  extract_characteristic <- function(
    characteristics,
    key
  ) {
    prefix <- paste0(
      tolower(key),
      ":"
    )

    match_index <- which(
      startsWith(
        tolower(characteristics),
        prefix
      )
    )

    if (length(match_index) == 0) {
      return(NA_character_)
    }

    trimws(
      sub(
        "^[^:]+:\\s*",
        "",
        characteristics[match_index[1]]
      )
    )
  }

  first_or_na <- function(value) {
    if (is.null(value) || length(value) == 0) {
      return(NA_character_)
    }

    value[1]
  }

  sample_list <- lapply(
    records,
    function(record) {
      characteristics <- record$Characteristics
      relation_text <- paste(
        record$Relations,
        collapse = " | "
      )

      data.frame(
        GEO_ID = GEO_ID,
        Sample_ID = first_or_na(record$Sample_ID),
        Sample_Name = first_or_na(record$Sample_ID),
        Original_Title = first_or_na(
          record$Original_Title
        ),
        Original_Source_Name = first_or_na(
          record$Original_Source_Name
        ),
        Tissue_Original = extract_characteristic(
          characteristics,
          "tissue"
        ),
        Patient_ID_Original = extract_characteristic(
          characteristics,
          "patient id"
        ),
        Sample_Time_Original = extract_characteristic(
          characteristics,
          "sample type"
        ),
        Treatment_Original = extract_characteristic(
          characteristics,
          "treatment"
        ),
        Survival_Original = extract_characteristic(
          characteristics,
          "survival (d28)"
        ),
        Sex_Original = extract_characteristic(
          characteristics,
          "Sex"
        ),
        Age_Original = extract_characteristic(
          characteristics,
          "age"
        ),
        TBSA_Original = extract_characteristic(
          characteristics,
          "tbsa"
        ),
        Platform_ID = first_or_na(
          record$Platform_ID
        ),
        Original_Characteristics = paste(
          characteristics,
          collapse = " | "
        ),
        BioSample_ID = NA_character_,
        SRA_Experiment = NA_character_,
        CEL_GZ_URL = first_or_na(
          record$Supplementary_File
        ),
        Relation_Text = relation_text,
        stringsAsFactors = FALSE
      )
    }
  )

  do.call(rbind, sample_list)
}

sample_metadata <-
  parse_geo_soft_samples_streaming(
    soft_file
  )

if (nrow(sample_metadata) != 117) {
  stop(
    paste0(
      "The SOFT file should contain 117 samples, but ",
      nrow(sample_metadata),
      " were parsed."
    )
  )
}

sample_metadata$Patient_ID <-
  trimws(sample_metadata$Patient_ID_Original)

if (any(
  is.na(sample_metadata$Patient_ID) |
    sample_metadata$Patient_ID == ""
)) {
  stop(
    "At least one sample lacks a valid Patient_ID."
  )
}

sample_metadata$Sample_Time <- toupper(
  trimws(sample_metadata$Sample_Time_Original)
)

if (!all(
  sample_metadata$Sample_Time %in%
    SAMPLE_TIME_LEVELS
)) {
  stop(
    paste0(
      "Unexpected sample-time labels were found: ",
      paste(
        unique(
          sample_metadata$Sample_Time[
            !sample_metadata$Sample_Time %in%
              SAMPLE_TIME_LEVELS
          ]
        ),
        collapse = ", "
      )
    )
  )
}

sample_metadata$Sample_Time <- factor(
  sample_metadata$Sample_Time,
  levels = SAMPLE_TIME_LEVELS
)

treatment_clean <- trimws(
  sample_metadata$Treatment_Original
)

sample_metadata$Treatment <- ifelse(
  sample_metadata$Sample_Time == "HV",
  "Healthy_volunteer",
  ifelse(
    treatment_clean == "Placebo",
    "Placebo",
    ifelse(
      treatment_clean == "Hydrocortisone",
      "Hydrocortisone",
      NA_character_
    )
  )
)

if (any(is.na(sample_metadata$Treatment))) {
  stop(
    "At least one sample has an unsupported treatment label."
  )
}

sample_metadata$Treatment <- factor(
  sample_metadata$Treatment,
  levels = STUDY_GROUP_LEVELS
)

strict_numeric <- function(
  character_values,
  field_name,
  allow_missing = TRUE
) {
  cleaned <- trimws(character_values)

  cleaned[
    tolower(cleaned) %in% c(
      "",
      "--",
      "unknown",
      "na",
      "n/a",
      "not reported"
    )
  ] <- NA_character_

  numeric_values <- suppressWarnings(
    as.numeric(cleaned)
  )

  invalid <- !is.na(cleaned) &
    is.na(numeric_values)

  if (any(invalid)) {
    stop(
      paste0(
        "Unsupported non-numeric values were found in ",
        field_name,
        ": ",
        paste(
          unique(character_values[invalid]),
          collapse = ", "
        )
      )
    )
  }

  if (!allow_missing &&
      any(is.na(numeric_values))) {
    stop(
      paste0(
        field_name,
        " contains missing values."
      )
    )
  }

  numeric_values
}

sample_metadata$Age <- strict_numeric(
  sample_metadata$Age_Original,
  "Age",
  allow_missing = TRUE
)

sample_metadata$TBSA_Percent <- strict_numeric(
  sample_metadata$TBSA_Original,
  "TBSA",
  allow_missing = TRUE
)

sex_lower <- tolower(
  trimws(sample_metadata$Sex_Original)
)

sample_metadata$Sex <- ifelse(
  sex_lower == "m",
  "Male",
  ifelse(
    sex_lower == "f",
    "Female",
    NA_character_
  )
)

burn_index <- sample_metadata$Sample_Time != "HV"
hv_index <- sample_metadata$Sample_Time == "HV"

if (any(is.na(sample_metadata$Age[burn_index]))) {
  stop(
    "At least one burn-patient sample has missing age."
  )
}

if (any(is.na(sample_metadata$Sex[burn_index]))) {
  stop(
    "At least one burn-patient sample has missing sex."
  )
}

if (any(is.na(sample_metadata$TBSA_Percent[burn_index]))) {
  stop(
    "At least one burn-patient sample has missing TBSA."
  )
}

if (!all(is.na(sample_metadata$Age[hv_index]))) {
  stop(
    "Healthy-volunteer age was expected to be unavailable in the uploaded SOFT metadata."
  )
}

if (!all(is.na(sample_metadata$Sex[hv_index]))) {
  stop(
    "Healthy-volunteer sex was expected to be unavailable in the uploaded SOFT metadata."
  )
}

if (!all(is.na(sample_metadata$TBSA_Percent[hv_index]))) {
  stop(
    "Healthy-volunteer TBSA was expected to be unavailable."
  )
}

sample_metadata$Outcome <- ifelse(
  sample_metadata$Sample_Time == "HV",
  NA_character_,
  ifelse(
    tolower(
      trimws(
        sample_metadata$Survival_Original
      )
    ) == "survivor",
    "D28_survivor",
    ifelse(
      tolower(
        trimws(
          sample_metadata$Survival_Original
        )
      ) == "non survivor",
      "D28_non_survivor",
      NA_character_
    )
  )
)

if (any(is.na(sample_metadata$Outcome[burn_index]))) {
  stop(
    "At least one burn-patient sample has an unsupported D28 survival label."
  )
}

sample_metadata$Group <-
  as.character(sample_metadata$Treatment)

sample_metadata$Group <- factor(
  sample_metadata$Group,
  levels = STUDY_GROUP_LEVELS
)

sample_metadata$Treatment_Time <- ifelse(
  sample_metadata$Sample_Time == "HV",
  "Healthy_volunteer",
  paste0(
    as.character(sample_metadata$Treatment),
    "_",
    as.character(sample_metadata$Sample_Time)
  )
)

TREATMENT_TIME_LEVELS <- c(
  "Healthy_volunteer",
  "Placebo_S1",
  "Hydrocortisone_S1",
  "Placebo_S2",
  "Hydrocortisone_S2",
  "Placebo_S3",
  "Hydrocortisone_S3",
  "Placebo_S4",
  "Hydrocortisone_S4"
)

sample_metadata$Treatment_Time <- factor(
  sample_metadata$Treatment_Time,
  levels = TREATMENT_TIME_LEVELS
)

sample_metadata$Plot_Group_Time <- factor(
  sample_metadata$Treatment_Time,
  levels = TREATMENT_TIME_LEVELS,
  labels = c(
    "Healthy volunteer",
    "Placebo S1",
    "Hydrocortisone S1",
    "Placebo S2",
    "Hydrocortisone S2",
    "Placebo S3",
    "Hydrocortisone S3",
    "Placebo S4",
    "Hydrocortisone S4"
  )
)

sample_metadata$Time_or_Stage <- ifelse(
  sample_metadata$Sample_Time == "HV",
  "Healthy volunteer",
  ifelse(
    sample_metadata$Sample_Time == "S1",
    "S1: burn-shock inclusion, before treatment",
    ifelse(
      sample_metadata$Sample_Time == "S2",
      "S2: approximately 24 h after treatment initiation",
      ifelse(
        sample_metadata$Sample_Time == "S3",
        "S3: approximately 120 h after treatment initiation",
        "S4: approximately 168 h after treatment initiation"
      )
    )
  )
)

sample_metadata$Tissue <- "Whole blood"
sample_metadata$Data_Type <- "Affymetrix raw CEL"
sample_metadata$Is_Pooled <- FALSE
sample_metadata$Metadata_Confidence <-
  "Direct_from_GEO_SOFT"

patient_array_count <- table(
  sample_metadata$Patient_ID
)

sample_metadata$Arrays_Per_Patient <- as.integer(
  patient_array_count[
    sample_metadata$Patient_ID
  ]
)

sample_metadata$Is_Repeated_Measure <-
  sample_metadata$Arrays_Per_Patient > 1

sample_metadata$Is_Paired <-
  sample_metadata$Is_Repeated_Measure

sample_metadata$Quality_Notes <- paste0(
  "Whole-blood bulk expression; 30 burn patients were randomized to ",
  "hydrocortisone or placebo and sampled longitudinally; repeated samples ",
  "from the same patient require Patient_ID blocking; healthy-volunteer age ",
  "and sex are unavailable, so Burn S1 versus healthy volunteer is unadjusted; ",
  "the burn-only longitudinal model adjusts for age, sex, and TBSA; S2-S4 are ",
  "defined relative to treatment initiation, not directly relative to burn ",
  "injury; missing later samples may not be completely random; the original ",
  "publication used GCRMA, MAS5 filtering, and COMBAT, whereas BurnOmicsDB ",
  "uses uniform RMA without inferred batch correction; whole-blood differences ",
  "may reflect both intracellular regulation and leukocyte-composition changes"
)

sample_metadata$Quality_Notes[
  sample_metadata$Sample_Time == "S4"
] <- paste0(
  sample_metadata$Quality_Notes[
    sample_metadata$Sample_Time == "S4"
  ],
  "; S4 has limited sample availability, especially in the hydrocortisone group"
)

sample_metadata$Expected_CEL_GZ_Filename <-
  basename(sample_metadata$CEL_GZ_URL)

# Validate sample counts and patient structure.
# 检查样本数与患者结构。
actual_time_counts <- table(
  sample_metadata$Sample_Time
)

expected_time_counts <- c(
  HV = 13,
  S1 = 30,
  S2 = 27,
  S3 = 29,
  S4 = 18
)

if (!all(
  actual_time_counts[
    names(expected_time_counts)
  ] == expected_time_counts
)) {
  stop(
    paste0(
      "Unexpected sample-time counts: ",
      paste(
        names(actual_time_counts),
        actual_time_counts,
        sep = "=",
        collapse = ", "
      )
    )
  )
}

treatment_time_counts <- table(
  sample_metadata$Treatment,
  sample_metadata$Sample_Time
)

expected_treatment_time_counts <- matrix(
  c(
    13, 0, 0, 0, 0,
    0, 15, 15, 15, 11,
    0, 15, 12, 14, 7
  ),
  nrow = 3,
  byrow = TRUE,
  dimnames = list(
    STUDY_GROUP_LEVELS,
    SAMPLE_TIME_LEVELS
  )
)

if (!all(
  treatment_time_counts[
    rownames(expected_treatment_time_counts),
    colnames(expected_treatment_time_counts)
  ] == expected_treatment_time_counts
)) {
  stop(
    "Treatment-by-time sample counts do not match the expected GSE77791 design."
  )
}

burn_patient_ids <- unique(
  sample_metadata$Patient_ID[burn_index]
)

healthy_patient_ids <- unique(
  sample_metadata$Patient_ID[hv_index]
)

if (length(burn_patient_ids) != 30) {
  stop(
    paste0(
      "Expected 30 burn patients, but ",
      length(burn_patient_ids),
      " were identified."
    )
  )
}

if (length(healthy_patient_ids) != 13) {
  stop(
    paste0(
      "Expected 13 healthy volunteers, but ",
      length(healthy_patient_ids),
      " were identified."
    )
  )
}

burn_treatment_by_patient <- tapply(
  as.character(
    sample_metadata$Treatment[burn_index]
  ),
  sample_metadata$Patient_ID[burn_index],
  function(values) {
    length(unique(values))
  }
)

if (any(burn_treatment_by_patient != 1)) {
  stop(
    "At least one burn patient has inconsistent treatment labels across time."
  )
}

s1_patient_n <- length(
  unique(
    sample_metadata$Patient_ID[
      sample_metadata$Sample_Time == "S1"
    ]
  )
)

if (s1_patient_n != 30) {
  stop(
    "All 30 burn patients should have an S1 sample."
  )
}

# Verify CEL and SOFT sample identifiers.
# 检查CEL与SOFT样本ID。
if (!setequal(
  cel_sample_ids,
  sample_metadata$Sample_ID
)) {
  missing_from_soft <- setdiff(
    cel_sample_ids,
    sample_metadata$Sample_ID
  )

  missing_from_cel <- setdiff(
    sample_metadata$Sample_ID,
    cel_sample_ids
  )

  stop(
    paste0(
      "CEL and SOFT sample identifiers do not match. Missing from SOFT: ",
      paste(missing_from_soft, collapse = ", "),
      ". Missing from CEL: ",
      paste(missing_from_cel, collapse = ", "),
      "."
    )
  )
}

cel_file_map <- setNames(
  cel_files,
  cel_sample_ids
)

cel_files <- unname(
  cel_file_map[
    sample_metadata$Sample_ID
  ]
)

cel_sample_ids <- sample_metadata$Sample_ID

if (any(is.na(cel_files))) {
  stop(
    "At least one CEL path is missing after sample-order alignment."
  )
}

metadata_audit <- data.frame(
  Check = c(
    "Raw CEL file count",
    "SOFT sample count",
    "Healthy-volunteer arrays",
    "Burn-patient arrays",
    "Healthy volunteers",
    "Burn patients",
    "Hydrocortisone patients",
    "Placebo patients",
    "S1 arrays",
    "S2 arrays",
    "S3 arrays",
    "S4 arrays",
    "Hydrocortisone S4 arrays",
    "Placebo S4 arrays",
    "Healthy-volunteer age available",
    "Healthy-volunteer sex available",
    "Burn-patient age complete",
    "Burn-patient sex complete",
    "Burn-patient TBSA complete",
    "Pooled samples",
    "Patient-level repeated measures",
    "Explicit sample-level batch ID available",
    "Primary core contrasts"
  ),
  Value = c(
    length(cel_files),
    nrow(sample_metadata),
    sum(hv_index),
    sum(burn_index),
    length(healthy_patient_ids),
    length(burn_patient_ids),
    length(
      unique(
        sample_metadata$Patient_ID[
          sample_metadata$Treatment ==
            "Hydrocortisone"
        ]
      )
    ),
    length(
      unique(
        sample_metadata$Patient_ID[
          sample_metadata$Treatment ==
            "Placebo"
        ]
      )
    ),
    actual_time_counts["S1"],
    actual_time_counts["S2"],
    actual_time_counts["S3"],
    actual_time_counts["S4"],
    treatment_time_counts[
      "Hydrocortisone",
      "S4"
    ],
    treatment_time_counts[
      "Placebo",
      "S4"
    ],
    "No",
    "No",
    sum(
      is.na(
        sample_metadata$Age[burn_index]
      )
    ) == 0,
    sum(
      is.na(
        sample_metadata$Sex[burn_index]
      )
    ) == 0,
    sum(
      is.na(
        sample_metadata$TBSA_Percent[
          burn_index
        ]
      )
    ) == 0,
    sum(sample_metadata$Is_Pooled),
    "Yes",
    "No",
    7
  ),
  Interpretation = c(
    "Directly observed in GSE77791_RAW.tar.",
    "Directly parsed from the streaming SOFT parser.",
    "One array per healthy volunteer.",
    "Longitudinal arrays from randomized burn patients.",
    "Independent healthy volunteers.",
    "Thirty burn patients remained after original study QC exclusions.",
    "Randomized treatment arm.",
    "Randomized placebo arm.",
    "Pretreatment inclusion samples.",
    "Approximately 24 h after treatment initiation.",
    "Approximately 120 h after treatment initiation.",
    "Approximately 168 h after treatment initiation.",
    "Lower-confidence treatment comparison due to limited sample count.",
    "Lower-confidence treatment comparison due to limited sample count.",
    "Unavailable in the uploaded SOFT metadata.",
    "Unavailable in the uploaded SOFT metadata.",
    "Complete for the burn-only longitudinal model.",
    "Complete for the burn-only longitudinal model.",
    "Complete for the burn-only longitudinal model.",
    "No pooling is reported.",
    "duplicateCorrelation and Patient_ID blocking are required.",
    "COMBAT is not applied without a documented sample-level batch variable.",
    "One baseline burn-shock contrast, three placebo longitudinal contrasts, and three treatment contrasts."
  ),
  stringsAsFactors = FALSE
)

write.csv(
  metadata_audit,
  file = file.path(
    QC_DIR,
    "GSE77791_metadata_consistency_check.csv"
  ),
  row.names = FALSE
)

cat("Sample-time counts:\n")
print(actual_time_counts)

cat("\nTreatment-by-time counts:\n")
print(treatment_time_counts)

cat("\nMetadata audit:\n")
print(metadata_audit)


# ---- 04. Batch-wise raw CEL QC / 分批进行原始CEL质量控制 ----

raw_qc_metrics <- data.frame(
  Sample_ID = sample_metadata$Sample_ID,
  Patient_ID = sample_metadata$Patient_ID,
  Group = sample_metadata$Group,
  Sample_Time = sample_metadata$Sample_Time,
  Treatment_Time = sample_metadata$Treatment_Time,
  Raw_Median_Log2_Intensity = NA_real_,
  Raw_IQR_Log2_Intensity = NA_real_,
  Raw_Min_Log2_Intensity = NA_real_,
  Raw_Max_Log2_Intensity = NA_real_,
  stringsAsFactors = FALSE
)

batch_starts <- seq(
  1,
  length(cel_files),
  by = RAW_QC_BATCH_SIZE
)

for (batch_start in batch_starts) {
  batch_end <- min(
    batch_start + RAW_QC_BATCH_SIZE - 1,
    length(cel_files)
  )

  batch_index <- batch_start:batch_end

  cat(
    "Reading raw CEL QC batch ",
    batch_start,
    "-",
    batch_end,
    " of ",
    length(cel_files),
    ".\n",
    sep = ""
  )

  raw_batch <- affy::ReadAffy(
    filenames = cel_files[batch_index]
  )

  sampleNames(raw_batch) <-
    sample_metadata$Sample_ID[batch_index]

  raw_log2_batch <- log2(
    exprs(raw_batch) + 1
  )

  raw_qc_metrics$Raw_Median_Log2_Intensity[
    batch_index
  ] <- matrixStats::colMedians(
    raw_log2_batch
  )

  raw_qc_metrics$Raw_IQR_Log2_Intensity[
    batch_index
  ] <- matrixStats::colIQRs(
    raw_log2_batch
  )

  raw_qc_metrics$Raw_Min_Log2_Intensity[
    batch_index
  ] <- matrixStats::colMins(
    raw_log2_batch
  )

  raw_qc_metrics$Raw_Max_Log2_Intensity[
    batch_index
  ] <- matrixStats::colMaxs(
    raw_log2_batch
  )

  rm(raw_batch, raw_log2_batch)
  invisible(gc())
}

if (any(
  is.na(
    raw_qc_metrics$Raw_Median_Log2_Intensity
  )
)) {
  stop(
    "Raw CEL QC metrics are incomplete."
  )
}

raw_median_center <- median(
  raw_qc_metrics$Raw_Median_Log2_Intensity
)

raw_median_mad <- mad(
  raw_qc_metrics$Raw_Median_Log2_Intensity
)

raw_qc_metrics$Raw_Median_Robust_Z <-
  if (raw_median_mad > 0) {
    (
      raw_qc_metrics$Raw_Median_Log2_Intensity -
        raw_median_center
    ) / raw_median_mad
  } else {
    0
  }

raw_qc_metrics$Raw_QC_Flag <- abs(
  raw_qc_metrics$Raw_Median_Robust_Z
) > 4

write.csv(
  raw_qc_metrics,
  file = file.path(
    QC_DIR,
    "GSE77791_raw_array_QC_metrics.csv"
  ),
  row.names = FALSE
)

p_raw_qc <- ggplot(
  raw_qc_metrics,
  aes(
    x = Treatment_Time,
    y = Raw_Median_Log2_Intensity,
    color = Group,
    fill = Group
  )
) +
  geom_boxplot(
    alpha = 0.18,
    outlier.shape = NA,
    width = 0.65
  ) +
  geom_jitter(
    width = 0.14,
    height = 0,
    size = 1.6,
    alpha = 0.72
  ) +
  scale_color_manual(
    values = STUDY_GROUP_COLORS,
    drop = FALSE
  ) +
  scale_fill_manual(
    values = STUDY_GROUP_COLORS,
    drop = FALSE
  ) +
  scale_x_discrete(
    labels = levels(
      sample_metadata$Plot_Group_Time
    )
  ) +
  labs(
    title = "GSE77791 raw CEL intensity QC",
    subtitle =
      "Each point is the median raw log2 probe intensity of one array",
    x = "Study group and sample time",
    y = "Median raw log2 probe intensity",
    color = "Study group",
    fill = "Study group"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(
      angle = 35,
      hjust = 1
    )
  )

print(p_raw_qc)

ggsave(
  filename = file.path(
    FIGURES_DIR,
    "01_GSE77791_raw_data_QC.png"
  ),
  plot = p_raw_qc,
  width = 11,
  height = 6.2,
  dpi = 300
)


# ---- 05. RMA normalization with checkpoint and normalized QC / 带断点的RMA标准化与QC ----

RMA_CHECKPOINT_FILE <- file.path(
  OBJECTS_DIR,
  "GSE77791_RMA_probe_expression_checkpoint.rds"
)

if (
  file.exists(RMA_CHECKPOINT_FILE) &&
    !FORCE_RMA_RECOMPUTE
) {
  cat(
    "Loading the existing RMA checkpoint instead of repeating RMA.\n"
  )

  rma_probe_expression <- readRDS(
    RMA_CHECKPOINT_FILE
  )
} else {
  cat(
    "Starting RMA normalization for 117 CEL files.\n"
  )

  rma_cel_filenames <- basename(cel_files)
  rma_cel_paths <- file.path(
    CEL_DIR,
    rma_cel_filenames
  )

  if (any(!file.exists(rma_cel_paths))) {
    stop(
      "RMA input validation failed because one or more CEL files are missing."
    )
  }

  # justRMA prepends celfile.path to filenames; pass basenames only.
  # justRMA会拼接celfile.path，因此filenames仅传入文件名。
  rma_eset <- affy::justRMA(
    filenames = rma_cel_filenames,
    celfile.path = CEL_DIR,
    sampleNames = sample_metadata$Sample_ID,
    compress = FALSE,
    destructive = TRUE,
    verbose = TRUE
  )

  rma_probe_expression <- exprs(
    rma_eset
  )

  saveRDS(
    rma_probe_expression,
    file = RMA_CHECKPOINT_FILE,
    compress = FALSE
  )

  rm(rma_eset)
  invisible(gc())

  cat(
    "RMA normalization completed and checkpointed.\n"
  )
}

if (ncol(rma_probe_expression) != 117) {
  stop(
    "The RMA expression matrix does not contain 117 arrays."
  )
}

if (!identical(
  colnames(rma_probe_expression),
  sample_metadata$Sample_ID
)) {
  stop(
    "RMA expression columns do not match the sample metadata."
  )
}

if (any(is.na(rma_probe_expression))) {
  stop(
    "The RMA expression matrix contains missing values."
  )
}

if (any(!is.finite(rma_probe_expression))) {
  stop(
    "The RMA expression matrix contains non-finite values."
  )
}

if (WRITE_PROBE_LEVEL_CSV) {
  probe_csv_file <- file.path(
    RESULTS_DIR,
    "GSE77791_RMA_probe_level_expression.csv.gz"
  )

  if (!file.exists(probe_csv_file)) {
    probe_connection <- gzfile(
      probe_csv_file,
      open = "wt"
    )

    write.csv(
      data.frame(
        Probe_ID = rownames(rma_probe_expression),
        rma_probe_expression,
        check.names = FALSE
      ),
      file = probe_connection,
      row.names = FALSE
    )

    close(probe_connection)
    invisible(gc())
  }
}

rma_qc_metrics <- data.frame(
  Sample_ID = sample_metadata$Sample_ID,
  Patient_ID = sample_metadata$Patient_ID,
  Group = sample_metadata$Group,
  Sample_Time = sample_metadata$Sample_Time,
  Treatment_Time = sample_metadata$Treatment_Time,
  RMA_Median = matrixStats::colMedians(
    rma_probe_expression
  ),
  RMA_IQR = matrixStats::colIQRs(
    rma_probe_expression
  ),
  RLE_Median = NA_real_,
  RLE_IQR = NA_real_,
  stringsAsFactors = FALSE
)

probe_row_medians <- matrixStats::rowMedians(
  rma_probe_expression
)

rle_batch_starts <- seq(
  1,
  ncol(rma_probe_expression),
  by = RLE_BATCH_SIZE
)

for (batch_start in rle_batch_starts) {
  batch_end <- min(
    batch_start + RLE_BATCH_SIZE - 1,
    ncol(rma_probe_expression)
  )

  batch_index <- batch_start:batch_end

  rle_batch <- sweep(
    rma_probe_expression[
      ,
      batch_index,
      drop = FALSE
    ],
    1,
    probe_row_medians,
    "-"
  )

  rma_qc_metrics$RLE_Median[
    batch_index
  ] <- matrixStats::colMedians(
    rle_batch
  )

  rma_qc_metrics$RLE_IQR[
    batch_index
  ] <- matrixStats::colIQRs(
    rle_batch
  )

  rm(rle_batch)
  invisible(gc())
}

write.csv(
  rma_qc_metrics,
  file = file.path(
    QC_DIR,
    "GSE77791_RMA_array_QC_metrics.csv"
  ),
  row.names = FALSE
)

p_normalized <- ggplot(
  rma_qc_metrics,
  aes(
    x = Treatment_Time,
    y = RMA_Median,
    color = Group,
    fill = Group
  )
) +
  geom_boxplot(
    alpha = 0.18,
    outlier.shape = NA,
    width = 0.65
  ) +
  geom_jitter(
    width = 0.14,
    height = 0,
    size = 1.6,
    alpha = 0.72
  ) +
  scale_color_manual(
    values = STUDY_GROUP_COLORS,
    drop = FALSE
  ) +
  scale_fill_manual(
    values = STUDY_GROUP_COLORS,
    drop = FALSE
  ) +
  scale_x_discrete(
    labels = levels(
      sample_metadata$Plot_Group_Time
    )
  ) +
  labs(
    title =
      "GSE77791 RMA-normalized expression distribution",
    subtitle =
      "Each point is the median RMA log2 expression of one array",
    x = "Study group and sample time",
    y = "Median RMA log2 expression",
    color = "Study group",
    fill = "Study group"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(
      angle = 35,
      hjust = 1
    )
  )

print(p_normalized)

ggsave(
  filename = file.path(
    FIGURES_DIR,
    "02_GSE77791_normalized_expression_distribution.png"
  ),
  plot = p_normalized,
  width = 11,
  height = 6.2,
  dpi = 300
)


# ---- 06. Probe annotation and gene-level representative selection / 探针注释与基因级代表探针选择 ----

probe_ids <- rownames(
  rma_probe_expression
)

probe_annotation_raw <- AnnotationDbi::select(
  hgu133plus2.db,
  keys = probe_ids,
  columns = "ENTREZID",
  keytype = "PROBEID"
)

probe_annotation_raw$PROBEID <- as.character(
  probe_annotation_raw$PROBEID
)

probe_annotation_raw$ENTREZID <- as.character(
  probe_annotation_raw$ENTREZID
)

probe_entrez_count <- vapply(
  split(
    probe_annotation_raw$ENTREZID,
    probe_annotation_raw$PROBEID
  ),
  function(entrez_values) {
    length(
      unique(
        entrez_values[
          !is.na(entrez_values) &
            entrez_values != ""
        ]
      )
    )
  },
  integer(1)
)

probe_mapping_status <- data.frame(
  Probe_ID = probe_ids,
  Entrez_Mapping_Count = unname(
    probe_entrez_count[probe_ids]
  ),
  stringsAsFactors = FALSE
)

probe_mapping_status$Entrez_Mapping_Count[
  is.na(
    probe_mapping_status$Entrez_Mapping_Count
  )
] <- 0

probe_mapping_status$Mapping_Status <- ifelse(
  probe_mapping_status$Entrez_Mapping_Count == 0,
  "Unmapped",
  ifelse(
    probe_mapping_status$Entrez_Mapping_Count == 1,
    "Unambiguous",
    "Ambiguous"
  )
)

unambiguous_probe_ids <-
  probe_mapping_status$Probe_ID[
    probe_mapping_status$Mapping_Status ==
      "Unambiguous"
  ]

unambiguous_annotation_rows <-
  probe_annotation_raw[
    probe_annotation_raw$PROBEID %in%
      unambiguous_probe_ids &
      !is.na(probe_annotation_raw$ENTREZID) &
      probe_annotation_raw$ENTREZID != "",
  ]

unambiguous_probe_annotation <- unique(
  unambiguous_annotation_rows[
    ,
    c("PROBEID", "ENTREZID")
  ]
)

if (anyDuplicated(
  unambiguous_probe_annotation$PROBEID
) > 0) {
  stop(
    "A probe classified as unambiguous still maps to multiple Entrez Gene IDs."
  )
}

probe_mean_expression <- matrixStats::rowMeans2(
  rma_probe_expression
)

unambiguous_probe_annotation$Mean_RMA_Expression <-
  probe_mean_expression[
    unambiguous_probe_annotation$PROBEID
  ]

unambiguous_probe_annotation <-
  unambiguous_probe_annotation[
    order(
      unambiguous_probe_annotation$ENTREZID,
      -unambiguous_probe_annotation$Mean_RMA_Expression,
      unambiguous_probe_annotation$PROBEID
    ),
  ]

representative_probe_table <-
  unambiguous_probe_annotation[
    !duplicated(
      unambiguous_probe_annotation$ENTREZID
    ),
  ]

colnames(representative_probe_table)[
  colnames(representative_probe_table) ==
    "PROBEID"
] <- "Representative_Probe_ID"

colnames(representative_probe_table)[
  colnames(representative_probe_table) ==
    "ENTREZID"
] <- "NCBI_Gene_ID"

gene_expression <- rma_probe_expression[
  representative_probe_table$Representative_Probe_ID,
  ,
  drop = FALSE
]

rownames(gene_expression) <-
  representative_probe_table$NCBI_Gene_ID

if (anyDuplicated(
  rownames(gene_expression)
) > 0) {
  stop(
    "Duplicated NCBI Gene IDs remain after representative-probe selection."
  )
}

selected_entrez_ids <- rownames(
  gene_expression
)

gene_symbol <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = selected_entrez_ids,
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
)

ensembl_id <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = selected_entrez_ids,
  column = "ENSEMBL",
  keytype = "ENTREZID",
  multiVals = "first"
)

gene_name <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = selected_entrez_ids,
  column = "GENENAME",
  keytype = "ENTREZID",
  multiVals = "first"
)

gene_annotation <- data.frame(
  NCBI_Gene_ID = selected_entrez_ids,
  Gene_Symbol = unname(
    gene_symbol[selected_entrez_ids]
  ),
  Ensembl_ID = unname(
    ensembl_id[selected_entrez_ids]
  ),
  Gene_Name = unname(
    gene_name[selected_entrez_ids]
  ),
  Representative_Probe_ID =
    representative_probe_table$Representative_Probe_ID[
      match(
        selected_entrez_ids,
        representative_probe_table$NCBI_Gene_ID
      )
    ],
  Mean_RMA_Expression =
    representative_probe_table$Mean_RMA_Expression[
      match(
        selected_entrez_ids,
        representative_probe_table$NCBI_Gene_ID
      )
    ],
  stringsAsFactors = FALSE
)

gene_annotation$Mapping_Status <- ifelse(
  is.na(gene_annotation$Gene_Symbol),
  "Entrez_mapped_symbol_unavailable",
  "Mapped"
)

probe_mapping_audit <- merge(
  probe_mapping_status,
  unambiguous_probe_annotation,
  by.x = "Probe_ID",
  by.y = "PROBEID",
  all.x = TRUE,
  sort = FALSE
)

probe_mapping_audit$Selected_as_Representative <-
  FALSE

selected_match <- match(
  representative_probe_table$Representative_Probe_ID,
  probe_mapping_audit$Probe_ID
)

probe_mapping_audit$Selected_as_Representative[
  selected_match[!is.na(selected_match)]
] <- TRUE

probe_mapping_connection <- gzfile(
  file.path(
    QC_DIR,
    "GSE77791_probe_mapping_and_selection.csv.gz"
  ),
  open = "wt"
)

write.csv(
  probe_mapping_audit,
  file = probe_mapping_connection,
  row.names = FALSE
)

close(probe_mapping_connection)

mapping_summary <- data.frame(
  Metric = c(
    "Input probe sets",
    "Unmapped probe sets",
    "Ambiguous probe sets",
    "Unambiguous probe sets",
    "Unique NCBI Gene IDs represented",
    "Representative probes selected"
  ),
  Count = c(
    length(probe_ids),
    sum(
      probe_mapping_status$Mapping_Status ==
        "Unmapped"
    ),
    sum(
      probe_mapping_status$Mapping_Status ==
        "Ambiguous"
    ),
    sum(
      probe_mapping_status$Mapping_Status ==
        "Unambiguous"
    ),
    length(
      unique(
        representative_probe_table$NCBI_Gene_ID
      )
    ),
    nrow(representative_probe_table)
  ),
  stringsAsFactors = FALSE
)

write.csv(
  mapping_summary,
  file = file.path(
    QC_DIR,
    "GSE77791_probe_mapping_summary.csv"
  ),
  row.names = FALSE
)

cat("Probe-mapping summary:\n")
print(mapping_summary)

normalized_gene_connection <- gzfile(
  file.path(
    RESULTS_DIR,
    "GSE77791_normalized_expression.csv.gz"
  ),
  open = "wt"
)

write.csv(
  data.frame(
    NCBI_Gene_ID = rownames(gene_expression),
    gene_expression,
    check.names = FALSE
  ),
  file = normalized_gene_connection,
  row.names = FALSE
)

close(normalized_gene_connection)


# ---- 07. PCA and sample correlation / PCA与样本相关性 ----

gene_variance <- matrixStats::rowVars(
  gene_expression
)

n_pca_genes <- min(
  TOP_VARIABLE_GENES_FOR_PCA,
  length(gene_variance)
)

top_variable_gene_ids <- rownames(
  gene_expression
)[
  order(
    gene_variance,
    decreasing = TRUE
  )[seq_len(n_pca_genes)]
]

pca_result <- prcomp(
  t(
    gene_expression[
      top_variable_gene_ids,
      ,
      drop = FALSE
    ]
  ),
  center = TRUE,
  scale. = FALSE
)

pca_variance <- 100 * (
  pca_result$sdev^2 /
    sum(pca_result$sdev^2)
)

pca_table <- data.frame(
  Sample_ID = rownames(pca_result$x),
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  stringsAsFactors = FALSE
)

pca_table$Patient_ID <- sample_metadata$Patient_ID[
  match(
    pca_table$Sample_ID,
    sample_metadata$Sample_ID
  )
]

pca_table$Group <- sample_metadata$Group[
  match(
    pca_table$Sample_ID,
    sample_metadata$Sample_ID
  )
]

pca_table$Sample_Time <- sample_metadata$Sample_Time[
  match(
    pca_table$Sample_ID,
    sample_metadata$Sample_ID
  )
]

write.csv(
  pca_table,
  file = file.path(
    RESULTS_DIR,
    "GSE77791_PCA_scores.csv"
  ),
  row.names = FALSE
)

p_pca <- ggplot(
  pca_table,
  aes(
    x = PC1,
    y = PC2,
    color = Group,
    shape = Sample_Time
  )
) +
  geom_point(
    size = 2.5,
    alpha = 0.78
  ) +
  scale_color_manual(
    values = STUDY_GROUP_COLORS,
    drop = FALSE
  ) +
  scale_shape_manual(
    values = SAMPLE_TIME_SHAPES,
    drop = FALSE
  ) +
  labs(
    title = "GSE77791 PCA",
    subtitle = paste0(
      "Top ",
      n_pca_genes,
      " variable gene-level RMA features; 117 arrays"
    ),
    x = paste0(
      "PC1 (",
      round(pca_variance[1], 1),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(pca_variance[2], 1),
      "%)"
    ),
    color = "Study group",
    shape = "Sample time"
  ) +
  theme_classic(base_size = 12)

print(p_pca)

ggsave(
  filename = file.path(
    FIGURES_DIR,
    "03_GSE77791_PCA.png"
  ),
  plot = p_pca,
  width = 8.5,
  height = 6.3,
  dpi = 300
)

sample_correlation <- cor(
  gene_expression,
  method = "pearson"
)

correlation_connection <- gzfile(
  file.path(
    QC_DIR,
    "GSE77791_sample_correlation.csv.gz"
  ),
  open = "wt"
)

write.csv(
  data.frame(
    Sample_ID = rownames(sample_correlation),
    sample_correlation,
    check.names = FALSE
  ),
  file = correlation_connection,
  row.names = FALSE
)

close(correlation_connection)

correlation_annotation <- data.frame(
  Group = sample_metadata$Group,
  Sample_Time = sample_metadata$Sample_Time,
  row.names = sample_metadata$Sample_ID
)

correlation_annotation_colors <- list(
  Group = STUDY_GROUP_COLORS,
  Sample_Time = SAMPLE_TIME_COLORS
)

pheatmap::pheatmap(
  sample_correlation,
  annotation_col = correlation_annotation,
  annotation_row = correlation_annotation,
  annotation_colors = correlation_annotation_colors,
  show_colnames = FALSE,
  show_rownames = FALSE,
  border_color = NA,
  color = EXPRESSION_HEATMAP_COLORS,
  main = "GSE77791 sample correlation",
  silent = FALSE
)

pheatmap::pheatmap(
  sample_correlation,
  annotation_col = correlation_annotation,
  annotation_row = correlation_annotation,
  annotation_colors = correlation_annotation_colors,
  show_colnames = FALSE,
  show_rownames = FALSE,
  border_color = NA,
  color = EXPRESSION_HEATMAP_COLORS,
  main = "GSE77791 sample correlation",
  filename = file.path(
    FIGURES_DIR,
    "04_GSE77791_sample_correlation_heatmap.png"
  ),
  width = 11,
  height = 10,
  silent = TRUE
)


# ---- 08. Differential-expression models / 差异表达模型 ----

# Model A: pretreatment burn shock versus healthy volunteers.
# 模型A：治疗前烧伤休克与健康志愿者。
baseline_index <- sample_metadata$Sample_Time %in%
  c("HV", "S1")

baseline_metadata <- sample_metadata[
  baseline_index,
]

baseline_expression <- gene_expression[
  ,
  baseline_metadata$Sample_ID,
  drop = FALSE
]

baseline_metadata$Baseline_Group <- factor(
  ifelse(
    baseline_metadata$Sample_Time == "HV",
    "Healthy_volunteer",
    "Burn_shock_S1"
  ),
  levels = c(
    "Healthy_volunteer",
    "Burn_shock_S1"
  )
)

if (nrow(baseline_metadata) != 43) {
  stop(
    "The baseline model should contain 43 arrays."
  )
}

if (length(
  unique(
    baseline_metadata$Patient_ID
  )
) != 43) {
  stop(
    "The baseline model should contain one array per individual."
  )
}

baseline_design <- model.matrix(
  ~ 0 + Baseline_Group,
  data = baseline_metadata
)

colnames(baseline_design) <- sub(
  "^Baseline_Group",
  "",
  colnames(baseline_design)
)

rownames(baseline_design) <-
  baseline_metadata$Sample_ID

if (qr(baseline_design)$rank <
    ncol(baseline_design)) {
  stop(
    "The baseline design matrix is not full rank."
  )
}

baseline_fit <- limma::lmFit(
  baseline_expression,
  design = baseline_design
)

baseline_contrast_matrix <- limma::makeContrasts(
  GSE77791_BurnShockS1_vs_HealthyVolunteer =
    Burn_shock_S1 - Healthy_volunteer,
  levels = baseline_design
)

baseline_fit_contrasts <- limma::contrasts.fit(
  baseline_fit,
  contrasts = baseline_contrast_matrix
)

baseline_fit_contrasts <- limma::eBayes(
  baseline_fit_contrasts,
  robust = TRUE,
  trend = TRUE
)

# Model B: burn-only longitudinal treatment model.
# 模型B：患者内纵向治疗模型。
burn_metadata <- sample_metadata[
  burn_index,
]

burn_expression <- gene_expression[
  ,
  burn_metadata$Sample_ID,
  drop = FALSE
]

if (nrow(burn_metadata) != 104) {
  stop(
    "The burn-only longitudinal model should contain 104 arrays."
  )
}

burn_metadata$Treatment_Time <- factor(
  as.character(
    burn_metadata$Treatment_Time
  ),
  levels = c(
    "Placebo_S1",
    "Placebo_S2",
    "Placebo_S3",
    "Placebo_S4",
    "Hydrocortisone_S1",
    "Hydrocortisone_S2",
    "Hydrocortisone_S3",
    "Hydrocortisone_S4"
  )
)

burn_metadata$Sex <- factor(
  burn_metadata$Sex,
  levels = c("Female", "Male")
)

burn_metadata$Age_Centered <-
  burn_metadata$Age -
    mean(burn_metadata$Age)

burn_metadata$TBSA_Centered <-
  burn_metadata$TBSA_Percent -
    mean(burn_metadata$TBSA_Percent)

longitudinal_design <- model.matrix(
  ~ 0 + Treatment_Time +
    Age_Centered +
    Sex +
    TBSA_Centered,
  data = burn_metadata
)

colnames(longitudinal_design) <- sub(
  "^Treatment_Time",
  "",
  colnames(longitudinal_design)
)

rownames(longitudinal_design) <-
  burn_metadata$Sample_ID

if (qr(longitudinal_design)$rank <
    ncol(longitudinal_design)) {
  stop(
    "The burn-only longitudinal design matrix is not full rank."
  )
}

patient_block <- factor(
  burn_metadata$Patient_ID
)

correlation_fit <- limma::duplicateCorrelation(
  burn_expression,
  design = longitudinal_design,
  block = patient_block
)

within_patient_correlation <-
  correlation_fit$consensus.correlation

if (!is.finite(within_patient_correlation)) {
  stop(
    "Within-patient correlation could not be estimated."
  )
}

cat(
  "Estimated within-patient correlation: ",
  round(within_patient_correlation, 6),
  "\n",
  sep = ""
)

longitudinal_fit <- limma::lmFit(
  burn_expression,
  design = longitudinal_design,
  block = patient_block,
  correlation = within_patient_correlation
)

longitudinal_contrast_matrix <- limma::makeContrasts(
  GSE77791_PlaceboS2_vs_PlaceboS1 =
    Placebo_S2 - Placebo_S1,
  GSE77791_PlaceboS3_vs_PlaceboS1 =
    Placebo_S3 - Placebo_S1,
  GSE77791_PlaceboS4_vs_PlaceboS1 =
    Placebo_S4 - Placebo_S1,
  GSE77791_Hydrocortisone_vs_Placebo_at_S2 =
    Hydrocortisone_S2 - Placebo_S2,
  GSE77791_Hydrocortisone_vs_Placebo_at_S3 =
    Hydrocortisone_S3 - Placebo_S3,
  GSE77791_Hydrocortisone_vs_Placebo_at_S4 =
    Hydrocortisone_S4 - Placebo_S4,
  GSE77791_RandomizationBaseline_Hydrocortisone_vs_Placebo_S1 =
    Hydrocortisone_S1 - Placebo_S1,
  levels = longitudinal_design
)

longitudinal_fit_contrasts <- limma::contrasts.fit(
  longitudinal_fit,
  contrasts = longitudinal_contrast_matrix
)

longitudinal_fit_contrasts <- limma::eBayes(
  longitudinal_fit_contrasts,
  robust = TRUE,
  trend = TRUE
)

cat("Baseline design-matrix columns:\n")
print(colnames(baseline_design))

cat("\nLongitudinal design-matrix columns:\n")
print(colnames(longitudinal_design))


# ---- 09. Gene annotation and complete result tables / Gene注释与完整结果表 ----

extract_limma_results <- function(
  fitted_object,
  contrast_id,
  annotation_data,
  contrast_info
) {
  contrast_result <- limma::topTable(
    fitted_object,
    coef = contrast_id,
    number = Inf,
    sort.by = "P"
  )

  contrast_result$NCBI_Gene_ID <-
    rownames(contrast_result)

  annotation_match <- match(
    contrast_result$NCBI_Gene_ID,
    annotation_data$NCBI_Gene_ID
  )

  output <- data.frame(
    NCBI_Gene_ID =
      contrast_result$NCBI_Gene_ID,
    Gene_Symbol =
      annotation_data$Gene_Symbol[
        annotation_match
      ],
    Ensembl_ID =
      annotation_data$Ensembl_ID[
        annotation_match
      ],
    Gene_Name =
      annotation_data$Gene_Name[
        annotation_match
      ],
    Mapping_Status =
      annotation_data$Mapping_Status[
        annotation_match
      ],
    Representative_Probe_ID =
      annotation_data$Representative_Probe_ID[
        annotation_match
      ],
    Model_ID =
      contrast_info$Model_ID,
    Contrast_ID =
      contrast_id,
    Contrast_Label =
      contrast_info$Contrast_Label,
    Case_Group =
      contrast_info$Case_Group,
    Control_Group =
      contrast_info$Control_Group,
    Case_N =
      contrast_info$Case_N,
    Control_N =
      contrast_info$Control_N,
    Case_Patient_N =
      contrast_info$Case_Patient_N,
    Control_Patient_N =
      contrast_info$Control_Patient_N,
    Complete_Pair_N =
      contrast_info$Complete_Pair_N,
    Is_Paired_Contrast =
      contrast_info$Is_Paired_Contrast,
    log2FC =
      contrast_result$logFC,
    Fold_Change =
      2^contrast_result$logFC,
    Mean_Normalized_Expression =
      contrast_result$AveExpr,
    Statistic =
      contrast_result$t,
    P_value =
      contrast_result$P.Value,
    FDR =
      contrast_result$adj.P.Val,
    B_statistic =
      contrast_result$B,
    Age_Adjusted =
      contrast_info$Age_Adjusted,
    Sex_Adjusted =
      contrast_info$Sex_Adjusted,
    TBSA_Adjusted =
      contrast_info$TBSA_Adjusted,
    Within_Patient_Correlation =
      contrast_info$Within_Patient_Correlation,
    Evidence_Confidence =
      contrast_info$Evidence_Confidence,
    stringsAsFactors = FALSE
  )

  output$Direction <- ifelse(
    output$log2FC > 0,
    "Up",
    ifelse(
      output$log2FC < 0,
      "Down",
      "No_change"
    )
  )

  output$DE_Status <- ifelse(
    output$FDR < FDR_CUTOFF &
      output$log2FC >= LOG2FC_CUTOFF,
    "Up_significant",
    ifelse(
      output$FDR < FDR_CUTOFF &
        output$log2FC <= -LOG2FC_CUTOFF,
      "Down_significant",
      "Not_significant"
    )
  )

  output$NegLog10_FDR <- -log10(
    pmax(
      output$FDR,
      .Machine$double.xmin
    )
  )

  output
}

complete_pair_count <- function(
  treatment_name,
  case_time,
  control_time
) {
  patient_time_list <- split(
    as.character(
      burn_metadata$Sample_Time[
        burn_metadata$Treatment ==
          treatment_name
      ]
    ),
    burn_metadata$Patient_ID[
      burn_metadata$Treatment ==
        treatment_name
    ]
  )

  sum(
    vapply(
      patient_time_list,
      function(values) {
        case_time %in% values &&
          control_time %in% values
      },
      logical(1)
    )
  )
}

core_contrast_definitions <- data.frame(
  Contrast_ID = c(
    "GSE77791_BurnShockS1_vs_HealthyVolunteer",
    "GSE77791_PlaceboS2_vs_PlaceboS1",
    "GSE77791_PlaceboS3_vs_PlaceboS1",
    "GSE77791_PlaceboS4_vs_PlaceboS1",
    "GSE77791_Hydrocortisone_vs_Placebo_at_S2",
    "GSE77791_Hydrocortisone_vs_Placebo_at_S3",
    "GSE77791_Hydrocortisone_vs_Placebo_at_S4"
  ),
  Model_ID = c(
    "Baseline_unadjusted_model",
    rep(
      "Burn_longitudinal_adjusted_model",
      6
    )
  ),
  Contrast_Label = c(
    "Burn shock S1 vs healthy volunteer",
    "Placebo S2 vs Placebo S1",
    "Placebo S3 vs Placebo S1",
    "Placebo S4 vs Placebo S1",
    "Hydrocortisone vs Placebo at S2",
    "Hydrocortisone vs Placebo at S3",
    "Hydrocortisone vs Placebo at S4"
  ),
  Case_Group = c(
    "Burn shock S1",
    "Placebo S2",
    "Placebo S3",
    "Placebo S4",
    "Hydrocortisone S2",
    "Hydrocortisone S3",
    "Hydrocortisone S4"
  ),
  Control_Group = c(
    "Healthy volunteer",
    "Placebo S1",
    "Placebo S1",
    "Placebo S1",
    "Placebo S2",
    "Placebo S3",
    "Placebo S4"
  ),
  Case_N = c(
    30,
    15,
    15,
    11,
    12,
    14,
    7
  ),
  Control_N = c(
    13,
    15,
    15,
    15,
    15,
    15,
    11
  ),
  Case_Patient_N = c(
    30,
    15,
    15,
    11,
    12,
    14,
    7
  ),
  Control_Patient_N = c(
    13,
    15,
    15,
    15,
    15,
    15,
    11
  ),
  Complete_Pair_N = c(
    NA,
    complete_pair_count(
      "Placebo",
      "S2",
      "S1"
    ),
    complete_pair_count(
      "Placebo",
      "S3",
      "S1"
    ),
    complete_pair_count(
      "Placebo",
      "S4",
      "S1"
    ),
    NA,
    NA,
    NA
  ),
  Is_Paired_Contrast = c(
    FALSE,
    TRUE,
    TRUE,
    TRUE,
    FALSE,
    FALSE,
    FALSE
  ),
  Sample_Context = c(
    "Whole blood at burn-shock inclusion before treatment",
    "Placebo-treated whole blood longitudinal response",
    "Placebo-treated whole blood longitudinal response",
    "Placebo-treated whole blood longitudinal response",
    "Whole blood during randomized hydrocortisone treatment",
    "Whole blood during randomized hydrocortisone treatment",
    "Whole blood during randomized hydrocortisone treatment"
  ),
  Time_or_Stage = c(
    "S1: burn-shock inclusion before treatment",
    "S2: approximately 24 h after treatment initiation",
    "S3: approximately 120 h after treatment initiation",
    "S4: approximately 168 h after treatment initiation",
    "S2: approximately 24 h after treatment initiation",
    "S3: approximately 120 h after treatment initiation",
    "S4: approximately 168 h after treatment initiation"
  ),
  Age_Adjusted = c(
    FALSE,
    rep(TRUE, 6)
  ),
  Sex_Adjusted = c(
    FALSE,
    rep(TRUE, 6)
  ),
  TBSA_Adjusted = c(
    FALSE,
    rep(TRUE, 6)
  ),
  Within_Patient_Correlation = c(
    NA_real_,
    rep(
      within_patient_correlation,
      6
    )
  ),
  Evidence_Confidence = c(
    "Core",
    "Core",
    "Core",
    "Moderate_due_to_missing_S4_samples",
    "Core",
    "Core",
    "Lower_due_to_small_S4_groups"
  ),
  stringsAsFactors = FALSE
)

baseline_info <- core_contrast_definitions[
  core_contrast_definitions$Model_ID ==
    "Baseline_unadjusted_model",
]

baseline_result <- extract_limma_results(
  fitted_object =
    baseline_fit_contrasts,
  contrast_id =
    baseline_info$Contrast_ID,
  annotation_data =
    gene_annotation,
  contrast_info =
    baseline_info
)

longitudinal_result_list <- lapply(
  core_contrast_definitions$Contrast_ID[
    core_contrast_definitions$Model_ID ==
      "Burn_longitudinal_adjusted_model"
  ],
  function(contrast_id) {
    contrast_info <-
      core_contrast_definitions[
        core_contrast_definitions$Contrast_ID ==
          contrast_id,
      ]

    extract_limma_results(
      fitted_object =
        longitudinal_fit_contrasts,
      contrast_id =
        contrast_id,
      annotation_data =
        gene_annotation,
      contrast_info =
        contrast_info
    )
  }
)

all_gene_results <- do.call(
  rbind,
  c(
    list(baseline_result),
    longitudinal_result_list
  )
)

row.names(all_gene_results) <- NULL

all_results_connection <- gzfile(
  file.path(
    RESULTS_DIR,
    "GSE77791_all_gene_results.csv.gz"
  ),
  open = "wt"
)

write.csv(
  all_gene_results,
  file = all_results_connection,
  row.names = FALSE
)

close(all_results_connection)

cat("Differential-expression result counts:\n")
print(
  with(
    all_gene_results,
    table(
      Contrast_ID,
      DE_Status
    )
  )
)

# Randomization baseline QC is deliberately separate from the seven core results.
# 随机化基线QC单独保存，不进入七个核心比较。
baseline_randomization_qc_id <-
  "GSE77791_RandomizationBaseline_Hydrocortisone_vs_Placebo_S1"

baseline_randomization_qc_info <- data.frame(
  Model_ID =
    "Burn_longitudinal_adjusted_model",
  Contrast_Label =
    "Randomization baseline: Hydrocortisone vs Placebo at S1",
  Case_Group =
    "Hydrocortisone S1",
  Control_Group =
    "Placebo S1",
  Case_N = 15,
  Control_N = 15,
  Case_Patient_N = 15,
  Control_Patient_N = 15,
  Complete_Pair_N = NA,
  Is_Paired_Contrast = FALSE,
  Age_Adjusted = TRUE,
  Sex_Adjusted = TRUE,
  TBSA_Adjusted = TRUE,
  Within_Patient_Correlation =
    within_patient_correlation,
  Evidence_Confidence =
    "Randomization_baseline_QC_only",
  stringsAsFactors = FALSE
)

baseline_randomization_qc <- extract_limma_results(
  fitted_object =
    longitudinal_fit_contrasts,
  contrast_id =
    baseline_randomization_qc_id,
  annotation_data =
    gene_annotation,
  contrast_info =
    baseline_randomization_qc_info
)

write.csv(
  baseline_randomization_qc,
  file = file.path(
    QC_DIR,
    "GSE77791_randomization_baseline_gene_QC.csv"
  ),
  row.names = FALSE
)


# ---- 10. Volcano plots / 火山图 ----

for (
  contrast_id in
    core_contrast_definitions$Contrast_ID
) {
  volcano_data <- all_gene_results[
    all_gene_results$Contrast_ID ==
      contrast_id,
  ]

  contrast_info <-
    core_contrast_definitions[
      core_contrast_definitions$Contrast_ID ==
        contrast_id,
    ]

  label_candidates <- volcano_data[
    volcano_data$DE_Status !=
      "Not_significant" &
      !is.na(volcano_data$Gene_Symbol),
  ]

  label_candidates <- label_candidates[
    order(
      label_candidates$FDR,
      -abs(label_candidates$log2FC)
    ),
  ]

  label_candidates <- head(
    label_candidates,
    VOLCANO_LABEL_GENE_N
  )

  model_label <- ifelse(
    contrast_info$Model_ID ==
      "Baseline_unadjusted_model",
    "RMA + limma baseline model",
    "RMA + limma repeated-measures model"
  )

  p_volcano <- ggplot(
    volcano_data,
    aes(
      x = log2FC,
      y = NegLog10_FDR,
      color = DE_Status
    )
  ) +
    geom_point(
      alpha = 0.64,
      size = 1.1
    ) +
    geom_vline(
      xintercept = c(
        -LOG2FC_CUTOFF,
        LOG2FC_CUTOFF
      ),
      linetype = "dashed"
    ) +
    geom_hline(
      yintercept = -log10(FDR_CUTOFF),
      linetype = "dashed"
    ) +
    ggrepel::geom_text_repel(
      data = label_candidates,
      aes(label = Gene_Symbol),
      size = 3,
      max.overlaps = Inf
    ) +
    scale_color_manual(
      values = VOLCANO_COLORS,
      breaks = c(
        "Down_significant",
        "Not_significant",
        "Up_significant"
      ),
      drop = FALSE
    ) +
    labs(
      title = paste0(
        "GSE77791: ",
        contrast_info$Contrast_Label
      ),
      subtitle = paste0(
        model_label,
        "; |log2FC| >= ",
        LOG2FC_CUTOFF,
        ", FDR < ",
        FDR_CUTOFF
      ),
      x = "log2 fold change",
      y = "-log10(FDR)",
      color = "DE status"
    ) +
    theme_classic(base_size = 12)

  print(p_volcano)

  ggsave(
    filename = file.path(
      FIGURES_DIR,
      paste0(
        "05_GSE77791_volcano_",
        contrast_id,
        ".png"
      )
    ),
    plot = p_volcano,
    width = 8,
    height = 6,
    dpi = 300
  )
}


# ---- 11. Top differential-gene heatmaps / 主要差异基因热图 ----

get_contrast_sample_ids <- function(
  contrast_id
) {
  if (
    contrast_id ==
      "GSE77791_BurnShockS1_vs_HealthyVolunteer"
  ) {
    sample_metadata$Sample_ID[
      sample_metadata$Sample_Time %in%
        c("HV", "S1")
    ]
  } else if (
    contrast_id ==
      "GSE77791_PlaceboS2_vs_PlaceboS1"
  ) {
    sample_metadata$Sample_ID[
      sample_metadata$Treatment ==
        "Placebo" &
        sample_metadata$Sample_Time %in%
        c("S1", "S2")
    ]
  } else if (
    contrast_id ==
      "GSE77791_PlaceboS3_vs_PlaceboS1"
  ) {
    sample_metadata$Sample_ID[
      sample_metadata$Treatment ==
        "Placebo" &
        sample_metadata$Sample_Time %in%
        c("S1", "S3")
    ]
  } else if (
    contrast_id ==
      "GSE77791_PlaceboS4_vs_PlaceboS1"
  ) {
    sample_metadata$Sample_ID[
      sample_metadata$Treatment ==
        "Placebo" &
        sample_metadata$Sample_Time %in%
        c("S1", "S4")
    ]
  } else if (
    contrast_id ==
      "GSE77791_Hydrocortisone_vs_Placebo_at_S2"
  ) {
    sample_metadata$Sample_ID[
      sample_metadata$Sample_Time == "S2"
    ]
  } else if (
    contrast_id ==
      "GSE77791_Hydrocortisone_vs_Placebo_at_S3"
  ) {
    sample_metadata$Sample_ID[
      sample_metadata$Sample_Time == "S3"
    ]
  } else if (
    contrast_id ==
      "GSE77791_Hydrocortisone_vs_Placebo_at_S4"
  ) {
    sample_metadata$Sample_ID[
      sample_metadata$Sample_Time == "S4"
    ]
  } else {
    stop(
      paste0(
        "No heatmap sample rule is defined for ",
        contrast_id,
        "."
      )
    )
  }
}

get_contrast_annotation <- function(
  contrast_id,
  sample_ids
) {
  metadata_subset <- sample_metadata[
    match(
      sample_ids,
      sample_metadata$Sample_ID
    ),
  ]

  contrast_group <- if (
    contrast_id ==
      "GSE77791_BurnShockS1_vs_HealthyVolunteer"
  ) {
    ifelse(
      metadata_subset$Sample_Time == "HV",
      "Control",
      "Case"
    )
  } else if (
    contrast_id %in%
      c(
        "GSE77791_PlaceboS2_vs_PlaceboS1",
        "GSE77791_PlaceboS3_vs_PlaceboS1",
        "GSE77791_PlaceboS4_vs_PlaceboS1"
      )
  ) {
    case_time <- ifelse(
      grepl("S2", contrast_id),
      "S2",
      ifelse(
        grepl("S3", contrast_id),
        "S3",
        "S4"
      )
    )

    ifelse(
      metadata_subset$Sample_Time == "S1",
      "Control",
      ifelse(
        metadata_subset$Sample_Time ==
          case_time,
        "Case",
        NA_character_
      )
    )
  } else {
    ifelse(
      metadata_subset$Treatment ==
        "Placebo",
      "Control",
      "Case"
    )
  }

  data.frame(
    Contrast_Group = factor(
      contrast_group,
      levels = c(
        "Control",
        "Case"
      )
    ),
    row.names = sample_ids
  )
}

for (
  contrast_id in
    core_contrast_definitions$Contrast_ID
) {
  contrast_info <-
    core_contrast_definitions[
      core_contrast_definitions$Contrast_ID ==
        contrast_id,
    ]

  contrast_results <- all_gene_results[
    all_gene_results$Contrast_ID ==
      contrast_id,
  ]

  ranked_gene_ids <-
    contrast_results$NCBI_Gene_ID[
      order(
        contrast_results$FDR,
        -abs(contrast_results$log2FC)
      )
    ]

  n_heatmap_genes <- min(
    TOP_GENES_FOR_HEATMAP,
    length(ranked_gene_ids)
  )

  top_heatmap_ids <- ranked_gene_ids[
    seq_len(n_heatmap_genes)
  ]

  contrast_sample_ids <-
    get_contrast_sample_ids(
      contrast_id
    )

  heatmap_matrix <- gene_expression[
    top_heatmap_ids,
    contrast_sample_ids,
    drop = FALSE
  ]

  heatmap_z <- t(
    scale(
      t(heatmap_matrix)
    )
  )

  gene_label_match <- match(
    top_heatmap_ids,
    contrast_results$NCBI_Gene_ID
  )

  heatmap_labels <-
    contrast_results$Gene_Symbol[
      gene_label_match
    ]

  missing_heatmap_labels <-
    is.na(heatmap_labels) |
      heatmap_labels == ""

  heatmap_labels[
    missing_heatmap_labels
  ] <- top_heatmap_ids[
    missing_heatmap_labels
  ]

  rownames(heatmap_z) <-
    make.unique(heatmap_labels)

  heatmap_annotation <-
    get_contrast_annotation(
      contrast_id,
      contrast_sample_ids
    )

  heatmap_annotation_colors <- list(
    Contrast_Group =
      CONTRAST_GROUP_COLORS
  )

  show_sample_names <-
    length(contrast_sample_ids) <= 30

  heatmap_title <- paste0(
    "Top ",
    n_heatmap_genes,
    " genes: ",
    contrast_info$Contrast_Label
  )

  pheatmap::pheatmap(
    heatmap_z,
    annotation_col =
      heatmap_annotation,
    annotation_colors =
      heatmap_annotation_colors,
    show_colnames =
      show_sample_names,
    show_rownames = TRUE,
    border_color = NA,
    color =
      EXPRESSION_HEATMAP_COLORS,
    main = heatmap_title,
    silent = FALSE
  )

  pheatmap::pheatmap(
    heatmap_z,
    annotation_col =
      heatmap_annotation,
    annotation_colors =
      heatmap_annotation_colors,
    show_colnames =
      show_sample_names,
    show_rownames = TRUE,
    border_color = NA,
    color =
      EXPRESSION_HEATMAP_COLORS,
    main = heatmap_title,
    filename = file.path(
      FIGURES_DIR,
      paste0(
        "06_GSE77791_top_differential_genes_heatmap_",
        contrast_id,
        ".png"
      )
    ),
    width = 10,
    height = 10,
    silent = TRUE
  )
}


# ---- 12. Publication-context and HERV-probe checks / 论文背景与HERV探针核对 ----

# Primary-paper feature genes. The publication used GCRMA, MAS5 filtering,
# COMBAT, and longitudinal analyses, so this is not an exact replication.
# 原论文使用不同预处理和统计流程，因此这里只作定性核对。
primary_publication_gene_context <- data.frame(
  Gene_Symbol = c(
    "NR3C1",
    "ADRB2",
    "IL6",
    "HLA-DRA",
    "CD3D",
    "CD3E",
    "CD55",
    "CD300LF",
    "SLC2A3",
    "MMP8"
  ),
  Published_Context = c(
    "Glucocorticoid-receptor expression was similar at inclusion.",
    "Adrenergic receptor discussed in hydrocortisone-response analysis.",
    "Innate immune pathway context.",
    "Antigen-presentation and immunosuppression context.",
    "Adaptive T-cell response context.",
    "Adaptive T-cell response context.",
    "Immune and HERV-adjacent gene in the secondary publication.",
    "Myeloid inhibitory receptor and HERV-adjacent gene.",
    "Burn-responsive metabolic gene.",
    "Burn-associated neutrophil protease."
  ),
  Direct_Replication = FALSE,
  stringsAsFactors = FALSE
)

calculated_primary_genes <- all_gene_results[
  all_gene_results$Gene_Symbol %in%
    primary_publication_gene_context$Gene_Symbol,
  c(
    "NCBI_Gene_ID",
    "Gene_Symbol",
    "Representative_Probe_ID",
    "Contrast_ID",
    "Contrast_Label",
    "log2FC",
    "Fold_Change",
    "P_value",
    "FDR",
    "Direction",
    "DE_Status"
  )
]

primary_publication_feature_check <- merge(
  primary_publication_gene_context,
  calculated_primary_genes,
  by = "Gene_Symbol",
  all.x = TRUE,
  sort = FALSE
)

write.csv(
  primary_publication_feature_check,
  file = file.path(
    RESULTS_DIR,
    "GSE77791_publication_feature_gene_check.csv"
  ),
  row.names = FALSE
)

# Context for six HERV-targeting probes reported in the secondary paper.
# These probes are audited at probe level but do not automatically enter the
# standard gene-level database table.
herv_probe_context <- data.frame(
  Probe_ID = c(
    "1556107_at",
    "230354_at",
    "1553043_a_at",
    "1560527_at",
    "1559777_at",
    "236982_at"
  ),
  HERV_Family = c(
    "MLT1H",
    "LTR33",
    "MLT1D",
    "LTR101_Mam",
    "LTR16B2",
    "ERV24B_Prim-int"
  ),
  Nearby_Gene = c(
    "CD55",
    "SLC8A1",
    "CD300LF",
    "NFE4",
    "MIR3945HG",
    "PTTG1IP"
  ),
  Published_Burn_log2FC = c(
    1.13,
    1.73,
    1.48,
    -0.55,
    1.05,
    0.77
  ),
  stringsAsFactors = FALSE
)

herv_probe_context$Probe_Mapping_Status <-
  probe_mapping_audit$Mapping_Status[
    match(
      herv_probe_context$Probe_ID,
      probe_mapping_audit$Probe_ID
    )
  ]

herv_probe_context$Mapped_NCBI_Gene_ID <-
  probe_mapping_audit$ENTREZID[
    match(
      herv_probe_context$Probe_ID,
      probe_mapping_audit$Probe_ID
    )
  ]

herv_probe_context$Selected_as_Representative <-
  probe_mapping_audit$Selected_as_Representative[
    match(
      herv_probe_context$Probe_ID,
      probe_mapping_audit$Probe_ID
    )
  ]

herv_probe_context$Core_Database_Interpretation <- ifelse(
  herv_probe_context$Probe_Mapping_Status ==
    "Unambiguous" &
    herv_probe_context$Selected_as_Representative,
  "Included only through its unambiguous standard-gene mapping.",
  ifelse(
    herv_probe_context$Probe_Mapping_Status ==
      "Unambiguous",
    "Mapped to a standard gene but not selected as the representative probe.",
    "Excluded from the core standard-gene table."
  )
)

write.csv(
  herv_probe_context,
  file = file.path(
    QC_DIR,
    "GSE77791_HERV_probe_mapping_context.csv"
  ),
  row.names = FALSE
)

# Compare method-dependent counts without claiming exact replication.
baseline_fdr_only_n <- sum(
  baseline_result$FDR < 0.05
)

baseline_unified_threshold_n <- sum(
  baseline_result$DE_Status !=
    "Not_significant"
)

treatment_core_gene_n <- length(
  unique(
    all_gene_results$NCBI_Gene_ID[
      all_gene_results$Contrast_ID %in%
        c(
          "GSE77791_Hydrocortisone_vs_Placebo_at_S2",
          "GSE77791_Hydrocortisone_vs_Placebo_at_S3",
          "GSE77791_Hydrocortisone_vs_Placebo_at_S4"
        ) &
        all_gene_results$DE_Status !=
          "Not_significant"
    ]
  )
)

publication_method_context <- data.frame(
  Metric = c(
    "BurnOmicsDB baseline genes with FDR < 0.05",
    "BurnOmicsDB baseline genes meeting FDR < 0.05 and absolute log2FC >= 1",
    "Primary paper inclusion comparison reported genes",
    "Unique BurnOmicsDB genes significant in at least one pointwise hydrocortisone contrast",
    "Primary paper reported hydrocortisone-modulated genes across its longitudinal analysis"
  ),
  Count = c(
    baseline_fdr_only_n,
    baseline_unified_threshold_n,
    967,
    treatment_core_gene_n,
    175
  ),
  Directly_Comparable = FALSE,
  Explanation = c(
    "Uniform RMA gene-level analysis without healthy-volunteer covariate adjustment.",
    "BurnOmicsDB unified database threshold.",
    "Publication used GCRMA, MAS5 filtering, COMBAT, and probe-level analysis.",
    "Three pointwise treatment contrasts with the BurnOmicsDB unified threshold.",
    "Publication used its own longitudinal treatment-analysis strategy."
  ),
  stringsAsFactors = FALSE
)

write.csv(
  publication_method_context,
  file = file.path(
    RESULTS_DIR,
    "GSE77791_publication_method_context.csv"
  ),
  row.names = FALSE
)


# ---- 13. Create the BurnOmicsDB-ready result table / 创建BurnOmicsDB标准结果表 ----

contrast_match <- match(
  all_gene_results$Contrast_ID,
  core_contrast_definitions$Contrast_ID
)

database_quality_notes <- ifelse(
  all_gene_results$Contrast_ID ==
    "GSE77791_BurnShockS1_vs_HealthyVolunteer",
  paste0(
    "Core pretreatment burn-shock comparison; healthy-volunteer age and sex ",
    "are unavailable in GEO, so this contrast is not adjusted for age, sex, ",
    "or TBSA; S1 occurs at shock inclusion before treatment and generally ",
    "24-72 h after burn injury; the original publication used GCRMA, MAS5 ",
    "filtering, and COMBAT, whereas BurnOmicsDB uses uniform RMA without an ",
    "inferred batch variable; whole-blood differences may reflect both ",
    "intracellular regulation and leukocyte-composition changes"
  ),
  paste0(
    "Randomized longitudinal burn-shock cohort; the burn-only model adjusts ",
    "for age, sex, and TBSA and accounts for repeated samples with Patient_ID ",
    "blocking; S2-S4 are defined relative to treatment initiation; later ",
    "sample missingness may not be completely random; the original publication ",
    "used GCRMA, MAS5 filtering, and COMBAT, whereas BurnOmicsDB uses uniform ",
    "RMA without an inferred batch variable; whole-blood differences may ",
    "reflect both intracellular regulation and leukocyte-composition changes"
  )
)

database_quality_notes[
  all_gene_results$Contrast_ID %in%
    c(
      "GSE77791_PlaceboS4_vs_PlaceboS1",
      "GSE77791_Hydrocortisone_vs_Placebo_at_S4"
    )
] <- paste0(
  database_quality_notes[
    all_gene_results$Contrast_ID %in%
      c(
        "GSE77791_PlaceboS4_vs_PlaceboS1",
        "GSE77791_Hydrocortisone_vs_Placebo_at_S4"
      )
  ],
  "; S4 evidence has reduced precision because only 7 hydrocortisone and ",
  "11 placebo arrays are available"
)

database_ready <- data.frame(
  NCBI_Gene_ID =
    all_gene_results$NCBI_Gene_ID,
  Gene_Symbol =
    all_gene_results$Gene_Symbol,
  Ensembl_ID =
    all_gene_results$Ensembl_ID,
  Gene_Name =
    all_gene_results$Gene_Name,
  GEO_ID = GEO_ID,
  Organism = "Homo sapiens",
  Study_Population =
    "Adults with severe burn shock randomized to hydrocortisone or placebo, with healthy volunteers",
  Tissue = "Whole blood",
  Sample_Context =
    core_contrast_definitions$Sample_Context[
      contrast_match
    ],
  Time_or_Stage =
    core_contrast_definitions$Time_or_Stage[
      contrast_match
    ],
  Contrast_ID =
    all_gene_results$Contrast_ID,
  Contrast_Label =
    all_gene_results$Contrast_Label,
  Case_Group =
    all_gene_results$Case_Group,
  Control_Group =
    all_gene_results$Control_Group,
  Case_N =
    all_gene_results$Case_N,
  Control_N =
    all_gene_results$Control_N,
  log2FC =
    all_gene_results$log2FC,
  Fold_Change =
    all_gene_results$Fold_Change,
  Direction =
    all_gene_results$Direction,
  Mean_log2CPM = NA_real_,
  P_value =
    all_gene_results$P_value,
  FDR =
    all_gene_results$FDR,
  DE_Status =
    all_gene_results$DE_Status,
  Platform =
    "Affymetrix Human Genome U133 Plus 2.0 Array (GPL570)",
  Input_Data =
    "Raw CEL files from GSE77791_RAW.tar",
  Normalization = paste0(
    "RMA: background correction, quantile normalization, ",
    "and median-polish probe-set summarization"
  ),
  Analysis_Method = ifelse(
    all_gene_results$Model_ID ==
      "Baseline_unadjusted_model",
    paste0(
      "limma baseline linear model with robust trend-aware empirical Bayes; ",
      "one pretreatment S1 array per burn patient versus one array per healthy ",
      "volunteer; no age, sex, or TBSA adjustment because healthy-volunteer ",
      "covariates are unavailable; one representative probe per NCBI Gene ID ",
      "selected by highest mean RMA expression across all 117 arrays; positive ",
      "log2FC represents Case_Group minus Control_Group"
    ),
    paste0(
      "limma longitudinal linear model with robust trend-aware empirical Bayes; ",
      "duplicateCorrelation and Patient_ID blocking; adjusted for continuous ",
      "age, sex, and TBSA; one representative probe per NCBI Gene ID selected ",
      "by highest mean RMA expression across all 117 arrays; positive log2FC ",
      "represents Case_Group minus Control_Group"
    )
  ),
  Annotation_Method = paste0(
    "hgu133plus2.db ",
    as.character(
      packageVersion("hgu133plus2.db")
    ),
    " used for probe-to-Entrez mapping; probes mapping to zero or multiple ",
    "Entrez Gene IDs excluded; org.Hs.eg.db ",
    as.character(
      packageVersion("org.Hs.eg.db")
    ),
    " used for current official symbol, Ensembl ID, and gene name"
  ),
  Quality_Notes =
    database_quality_notes,
  Mean_Normalized_Expression =
    all_gene_results$Mean_Normalized_Expression,
  Expression_Scale =
    "RMA log2 expression",
  Representative_Probe_ID =
    all_gene_results$Representative_Probe_ID,
  Case_Patient_N =
    all_gene_results$Case_Patient_N,
  Control_Patient_N =
    all_gene_results$Control_Patient_N,
  Complete_Pair_N =
    all_gene_results$Complete_Pair_N,
  Is_Paired_Contrast =
    all_gene_results$Is_Paired_Contrast,
  Within_Patient_Correlation =
    all_gene_results$Within_Patient_Correlation,
  Age_Adjusted =
    all_gene_results$Age_Adjusted,
  Sex_Adjusted =
    all_gene_results$Sex_Adjusted,
  TBSA_Adjusted =
    all_gene_results$TBSA_Adjusted,
  Is_Pooled = FALSE,
  Evidence_Confidence =
    all_gene_results$Evidence_Confidence,
  Quality_Grade =
    "Core_randomized_longitudinal",
  stringsAsFactors = FALSE
)

database_connection <- gzfile(
  file.path(
    RESULTS_DIR,
    "GSE77791_database_ready_all_genes.csv.gz"
  ),
  open = "wt"
)

write.csv(
  database_ready,
  file = database_connection,
  row.names = FALSE
)

close(database_connection)

cat(
  "Database-ready rows: ",
  nrow(database_ready),
  "\n",
  sep = ""
)


# ---- 14. Save metadata, analysis objects, summary, and session information / 保存元数据、对象、摘要与环境信息 ----

sample_metadata$Raw_Median_Log2_Intensity <-
  raw_qc_metrics$Raw_Median_Log2_Intensity[
    match(
      sample_metadata$Sample_ID,
      raw_qc_metrics$Sample_ID
    )
  ]

sample_metadata$Raw_IQR_Log2_Intensity <-
  raw_qc_metrics$Raw_IQR_Log2_Intensity[
    match(
      sample_metadata$Sample_ID,
      raw_qc_metrics$Sample_ID
    )
  ]

sample_metadata$Raw_QC_Flag <-
  raw_qc_metrics$Raw_QC_Flag[
    match(
      sample_metadata$Sample_ID,
      raw_qc_metrics$Sample_ID
    )
  ]

sample_metadata$RMA_Median <-
  rma_qc_metrics$RMA_Median[
    match(
      sample_metadata$Sample_ID,
      rma_qc_metrics$Sample_ID
    )
  ]

sample_metadata$RMA_IQR <-
  rma_qc_metrics$RMA_IQR[
    match(
      sample_metadata$Sample_ID,
      rma_qc_metrics$Sample_ID
    )
  ]

sample_metadata$RLE_Median <-
  rma_qc_metrics$RLE_Median[
    match(
      sample_metadata$Sample_ID,
      rma_qc_metrics$Sample_ID
    )
  ]

sample_metadata$RLE_IQR <-
  rma_qc_metrics$RLE_IQR[
    match(
      sample_metadata$Sample_ID,
      rma_qc_metrics$Sample_ID
    )
  ]

sample_metadata_output <- sample_metadata[
  ,
  c(
    "GEO_ID",
    "Sample_ID",
    "Sample_Name",
    "Patient_ID",
    "Group",
    "Tissue",
    "Time_or_Stage",
    "Treatment",
    "Outcome",
    "Platform_ID",
    "Data_Type",
    "Is_Paired",
    "Is_Pooled",
    "Is_Repeated_Measure",
    "Metadata_Confidence",
    "Original_Title",
    "Original_Source_Name",
    "Original_Characteristics",
    "BioSample_ID",
    "SRA_Experiment",
    "Quality_Notes",
    "Sample_Time",
    "Treatment_Time",
    "Age",
    "Sex",
    "TBSA_Percent",
    "Survival_Original",
    "Arrays_Per_Patient",
    "Tissue_Original",
    "CEL_GZ_URL",
    "Expected_CEL_GZ_Filename",
    "Raw_Median_Log2_Intensity",
    "Raw_IQR_Log2_Intensity",
    "Raw_QC_Flag",
    "RMA_Median",
    "RMA_IQR",
    "RLE_Median",
    "RLE_IQR"
  )
]

write.csv(
  sample_metadata_output,
  file = file.path(
    RESULTS_DIR,
    "GSE77791_sample_metadata.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

significant_count_table <- as.data.frame(
  with(
    all_gene_results,
    table(
      Contrast_ID,
      DE_Status
    )
  ),
  stringsAsFactors = FALSE
)

randomization_baseline_summary <- data.frame(
  Metric = c(
    "Genes with FDR < 0.05 at randomized pretreatment baseline",
    "Genes meeting FDR < 0.05 and absolute log2FC >= 1 at randomized pretreatment baseline",
    "Maximum absolute baseline log2FC",
    "Minimum baseline FDR"
  ),
  Value = c(
    sum(
      baseline_randomization_qc$FDR <
        0.05
    ),
    sum(
      baseline_randomization_qc$DE_Status !=
        "Not_significant"
    ),
    max(
      abs(
        baseline_randomization_qc$log2FC
      )
    ),
    min(
      baseline_randomization_qc$FDR
    )
  ),
  stringsAsFactors = FALSE
)

write.csv(
  randomization_baseline_summary,
  file = file.path(
    QC_DIR,
    "GSE77791_randomization_baseline_summary.csv"
  ),
  row.names = FALSE
)

# The full probe-level matrix has its own checkpoint and optional CSV and is not
# duplicated in the final analysis RDS.
# 完整探针矩阵已有独立断点文件，不重复写入最终RDS。
analysis_objects <- list(
  project = "BurnOmicsDB",
  GEO_ID = GEO_ID,
  raw_tar_file = raw_tar_file,
  soft_file = soft_file,
  rma_checkpoint_file =
    RMA_CHECKPOINT_FILE,
  sample_metadata =
    sample_metadata_output,
  metadata_audit =
    metadata_audit,
  raw_qc_metrics =
    raw_qc_metrics,
  rma_qc_metrics =
    rma_qc_metrics,
  probe_mapping_summary =
    mapping_summary,
  representative_probe_table =
    representative_probe_table,
  gene_annotation =
    gene_annotation,
  gene_expression =
    gene_expression,
  baseline_design =
    baseline_design,
  baseline_fitted_model =
    baseline_fit_contrasts,
  longitudinal_design =
    longitudinal_design,
  within_patient_correlation =
    within_patient_correlation,
  longitudinal_fitted_model =
    longitudinal_fit_contrasts,
  contrast_definitions =
    core_contrast_definitions,
  all_gene_results =
    all_gene_results,
  randomization_baseline_qc =
    baseline_randomization_qc,
  publication_feature_check =
    primary_publication_feature_check,
  publication_method_context =
    publication_method_context,
  herv_probe_context =
    herv_probe_context,
  thresholds = list(
    FDR_CUTOFF = FDR_CUTOFF,
    LOG2FC_CUTOFF =
      LOG2FC_CUTOFF
  )
)

saveRDS(
  analysis_objects,
  file = file.path(
    OBJECTS_DIR,
    "GSE77791_analysis_objects.rds"
  ),
  compress = "gzip"
)

summary_lines <- c(
  "BurnOmicsDB - GSE77791 analysis summary",
  "",
  paste0(
    "Analysis date: ",
    Sys.Date()
  ),
  "Data type: Affymetrix Human Genome U133 Plus 2.0 raw CEL",
  paste0(
    "Input CEL arrays: ",
    length(cel_files)
  ),
  "Healthy volunteers: 13",
  "Burn patients: 30",
  "Hydrocortisone patients: 15",
  "Placebo patients: 15",
  "Pooling: No",
  "",
  "Sample counts:",
  paste0(
    "  Healthy volunteers: ",
    actual_time_counts["HV"]
  ),
  paste0(
    "  S1 before treatment: ",
    actual_time_counts["S1"]
  ),
  paste0(
    "  S2 approximately 24 h after treatment initiation: ",
    actual_time_counts["S2"]
  ),
  paste0(
    "  S3 approximately 120 h after treatment initiation: ",
    actual_time_counts["S3"]
  ),
  paste0(
    "  S4 approximately 168 h after treatment initiation: ",
    actual_time_counts["S4"]
  ),
  "",
  "Treatment-by-time counts:",
  paste0(
    "  Placebo S1/S2/S3/S4: ",
    paste(
      treatment_time_counts[
        "Placebo",
        c("S1", "S2", "S3", "S4")
      ],
      collapse = "/"
    )
  ),
  paste0(
    "  Hydrocortisone S1/S2/S3/S4: ",
    paste(
      treatment_time_counts[
        "Hydrocortisone",
        c("S1", "S2", "S3", "S4")
      ],
      collapse = "/"
    )
  ),
  "",
  paste0(
    "Gene-level representatives tested: ",
    nrow(gene_expression)
  ),
  paste0(
    "Unmapped probe sets: ",
    mapping_summary$Count[
      mapping_summary$Metric ==
        "Unmapped probe sets"
    ]
  ),
  paste0(
    "Ambiguous probe sets excluded: ",
    mapping_summary$Count[
      mapping_summary$Metric ==
        "Ambiguous probe sets"
    ]
  ),
  paste0(
    "Representative-probe rule: highest mean RMA expression across all ",
    "117 arrays for each NCBI Gene ID"
  ),
  "",
  paste0(
    "Normalization: RMA background correction, quantile normalization, ",
    "and median-polish summarization"
  ),
  paste0(
    "Baseline model: unadjusted limma comparison of 30 pretreatment S1 arrays ",
    "and 13 healthy-volunteer arrays because healthy-volunteer age and sex ",
    "are unavailable"
  ),
  paste0(
    "Longitudinal model: limma with robust trend-aware empirical Bayes, ",
    "duplicateCorrelation blocking Patient_ID, and adjustment for age, sex, ",
    "and TBSA"
  ),
  paste0(
    "Estimated within-patient correlation: ",
    round(
      within_patient_correlation,
      6
    )
  ),
  paste0(
    "Threshold: |log2FC| >= ",
    LOG2FC_CUTOFF,
    " and FDR < ",
    FDR_CUTOFF
  ),
  "",
  "Core contrasts:",
  paste0(
    "  ",
    core_contrast_definitions$Contrast_ID,
    ": ",
    core_contrast_definitions$Contrast_Label,
    " (case arrays = ",
    core_contrast_definitions$Case_N,
    ", control arrays = ",
    core_contrast_definitions$Control_N,
    ", complete pairs = ",
    ifelse(
      is.na(
        core_contrast_definitions$Complete_Pair_N
      ),
      "NA",
      core_contrast_definitions$Complete_Pair_N
    ),
    ", evidence confidence = ",
    core_contrast_definitions$Evidence_Confidence,
    ")"
  ),
  "",
  "Differential-expression counts:",
  apply(
    significant_count_table,
    1,
    function(row_values) {
      paste0(
        "  ",
        row_values["Contrast_ID"],
        " | ",
        row_values["DE_Status"],
        ": ",
        row_values["Freq"]
      )
    }
  ),
  "",
  paste0(
    "Positive log2FC means higher expression in Case_Group than in ",
    "Control_Group."
  ),
  paste0(
    "No separate low-expression filter was applied after RMA and ",
    "representative-probe selection; all selected unambiguous gene ",
    "representatives were statistically tested and exported."
  ),
  "",
  "Major limitations:",
  paste0(
    "  Healthy-volunteer age and sex are unavailable, so Burn S1 versus ",
    "healthy volunteer is not covariate-adjusted."
  ),
  paste0(
    "  S2-S4 are defined relative to treatment initiation, not directly ",
    "relative to burn injury."
  ),
  paste0(
    "  S4 includes only 7 hydrocortisone and 11 placebo arrays and therefore ",
    "has reduced precision."
  ),
  paste0(
    "  Longitudinal missingness may be related to clinical severity, death, ",
    "RNA quality, or technical availability."
  ),
  paste0(
    "  The original paper used GCRMA, MAS5 detection filtering, and COMBAT; ",
    "the current uniform reanalysis uses RMA and does not infer an undocumented ",
    "batch variable."
  ),
  paste0(
    "  Standard gene-level mapping excludes HERV probes that do not map ",
    "unambiguously to one NCBI Gene ID."
  ),
  paste0(
    "  Whole-blood expression differences may reflect both intracellular ",
    "regulation and changes in leukocyte composition."
  ),
  "Quality grade: Core_randomized_longitudinal"
)

writeLines(
  summary_lines,
  con = file.path(
    RESULTS_DIR,
    "GSE77791_analysis_summary.txt"
  ),
  useBytes = TRUE
)

capture.output(
  sessionInfo(),
  file = file.path(
    RESULTS_DIR,
    "GSE77791_R_sessionInfo.txt"
  )
)

cat("Analysis completed.\n")
cat("Results directory:\n", RESULTS_DIR, "\n")
cat("Figures directory:\n", FIGURES_DIR, "\n")
cat("QC directory:\n", QC_DIR, "\n")
cat("R objects directory:\n", OBJECTS_DIR, "\n")


# ---- 15. Optional individual-gene checking / 可选的单基因检查 ----

# Run after Section 09.
# 在第09部分完成后运行。
genes_to_check <- c(
  "NR3C1",
  "ADRB2",
  "IL6",
  "HLA-DRA",
  "CD3D",
  "CD3E",
  "CD55",
  "CD300LF",
  "SLC2A3",
  "MMP8",
  "S100A8",
  "S100A9"
)

gene_check_table <- all_gene_results[
  all_gene_results$Gene_Symbol %in%
    genes_to_check,
  c(
    "NCBI_Gene_ID",
    "Gene_Symbol",
    "Representative_Probe_ID",
    "Contrast_ID",
    "Contrast_Label",
    "log2FC",
    "Fold_Change",
    "P_value",
    "FDR",
    "Direction",
    "DE_Status",
    "Evidence_Confidence"
  )
]

gene_check_table <- gene_check_table[
  order(
    match(
      gene_check_table$Gene_Symbol,
      genes_to_check
    ),
    match(
      gene_check_table$Contrast_ID,
      core_contrast_definitions$Contrast_ID
    )
  ),
]

print(gene_check_table)

# End of script / 代码结束

# 【需要记录的问题：SOFT与论文人口学统计不完全一致，20270728】
# 根据SOFT逐患者汇总，Hydrocortisone组为：
# Female = 3
# Median TBSA = 64%
# 但原论文Table 1报告：
# Female = 2
# Median TBSA = 75%
# Placebo组的性别、年龄和TBSA汇总则基本匹配论文。原论文确实报告Hydrocortisone组2名女性、TBSA中位数75%。
# 这说明论文汇总值与GEO样本级元数据之间存在一处无法从现有文件解决的不一致。当前分析使用SOFT提供的逐患者数据，做法是合理的；没有可靠的逐患者修正版，因此不建议凭论文汇总数字修改患者记录，也不需要重跑。
# 在后续study_manifest或Dataset页面加入一句即可：
# Sample-level sex and TBSA metadata were obtained from GEO SOFT. The SOFT-derived hydrocortisone-group summary differs slightly from the aggregate values reported in the primary publication.
