# ==============================================================================
# BurnOmicsDB: GSE37069 microarray analysis
# BurnOmicsDB：GSE37069微阵列分析
#
# GEO accession / GEO编号:
#   GSE37069
#
# Study / 研究:
#   Gene response to major burn injuries
#   严重烧伤后的长期外周血基因表达反应
#
# Data type / 数据类型:
#   Affymetrix Human Genome U133 Plus 2.0 Array (GPL570)
#   Raw CEL files contained in GSE37069_RAW.tar
#   GPL570 Affymetrix原始CEL文件。
#
# Study structure / 研究结构:
#   - 590 expression arrays in total.
#   - 37 healthy-control arrays representing 35 control subjects.
#   - 553 burn arrays collected longitudinally.
#   - GEO reports 244 burn patients, whereas sample titles in the uploaded
#     SOFT file contain 248 distinct burn Subject IDs. These IDs are retained
#     exactly as reported and are not merged by inference.
#   - Burn subjects contributed between 1 and 6 arrays.
#
#   共590张芯片：37张健康对照芯片代表35名对照，553张烧伤芯片为纵向采样。
#   GEO总体描述报告244名烧伤患者，但SOFT样本标题可解析出248个不同Subject ID。
#   本代码保留SOFT原始ID，不进行推测性合并。
#
# Tissue definition / 组织定义:
#   GEO and SOFT identify the samples as white blood cells from blood.
#   BurnOmicsDB therefore classifies this dataset as peripheral-blood
#   leukocyte expression, not adipose-tissue expression.
#   GEO和SOFT均将样本定义为血液白细胞，因此本数据库按外周血白细胞处理。
#
# Five prespecified time windows / 五个预设时间窗口:
#   1. 0-3 days:    hours since injury <= 72
#   2. 4-7 days:    72 < hours since injury <= 168
#   3. 8-14 days:   168 < hours since injury <= 336
#   4. 15-28 days:  336 < hours since injury <= 672
#   5. >28 days:    hours since injury > 672
#
# Prespecified contrasts / 预设比较:
#   1. Burn 0-3 days vs healthy control
#   2. Burn 4-7 days vs healthy control
#   3. Burn 8-14 days vs healthy control
#   4. Burn 15-28 days vs healthy control
#   5. Burn >28 days vs healthy control
#
# Positive log2FC always means higher expression in the burn time-window group
# than in healthy controls.
# 正log2FC始终表示烧伤时间窗组相对于健康对照表达更高。
#
# Analysis workflow / 分析流程:
#   GSE37069_RAW.tar
#   -> validate and extract 590 CEL.gz files
#   -> decompress CEL files
#   -> batch-wise raw CEL QC
#   -> memory-conscious RMA using justRMA()
#   -> RMA checkpoint saved to avoid repeating the expensive step
#   -> GPL570 probe-to-NCBI-Gene-ID mapping
#   -> remove unmapped and ambiguous probes
#   -> select one representative probe per NCBI Gene ID using the highest
#      mean RMA expression across all 590 arrays
#   -> limma linear model
#   -> duplicateCorrelation blocking Patient_ID
#   -> adjustment for continuous age and sex
#   -> five prespecified contrasts
#
# Repeated-measures rationale / 重复测量理由:
#   Burn subjects have 1-6 arrays, and two healthy controls have two arrays.
#   All arrays from the same derived subject ID are treated as correlated.
#   duplicateCorrelation estimates a common within-subject correlation.
#   烧伤患者有1-6张芯片，另有两名对照各检测两次。
#   同一Subject ID的芯片按重复测量处理。
#
# Probe-selection rule / 代表探针规则:
#   Probes mapping to zero or multiple NCBI Gene IDs are excluded.
#   For each remaining NCBI Gene ID, the probe with the highest mean RMA
#   expression across all 590 arrays is selected before differential testing.
#   The rule does not use group labels, p-values, or FDR.
#   无Gene ID或对应多个Gene ID的探针被排除；
#   每个Gene ID选择全部590张芯片平均RMA表达最高的探针。
#
# Cohort overlap / 队列重叠:
#   GSE19743 is a substantially overlapping subset of this burn cohort.
#   The two GEO Series may be displayed separately because their predefined
#   contrasts differ, but they must not be counted as independent replication.
#   GSE19743与本队列高度重叠，不能作为独立验证证据重复计数。
#
# Memory and disk strategy / 内存与磁盘策略:
#   - Raw CEL QC reads small batches rather than all raw intensities at once.
#   - justRMA() is used with destructive = TRUE.
#   - RMA probe-level expression is checkpointed as an RDS file.
#   - If a later section fails, rerunning Section 05 loads the checkpoint
#     instead of repeating RMA.
#   - RLE metrics are computed in column batches.
#   - The final analysis-object RDS excludes the full probe-level matrix.
#
# Hardware note / 硬件说明:
#   A 32-GB Intel Mac can attempt this analysis, but RMA on 590 GPL570 arrays
#   can approach the available memory. Close memory-intensive applications,
#   keep at least 20 GB of free disk space, and do not run the entire script
#   with Source on the first attempt.
#   590张GPL570芯片的RMA可能接近32GB机器的内存上限，请关闭其他大型程序。
#
# How to run in RStudio / 如何在RStudio中运行:
#   - Save this script in:
#     /Users/peter/Downloads/Project-2026-BurnOmicsDB/GSE37069/
#   - Keep these local files in the project folder:
#       GSE37069_RAW.tar
#       GSE37069_family.soft.gz
#       pnas.201222878.pdf
#       Dataset Note.rtf   (optional documentation)
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

options(timeout = 3600)
options(download.file.method = "libcurl")
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Change the mirror if it is unavailable in the current network.
# 如果当前网络无法使用，可以替换镜像。
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

GEO_ID <- "GSE37069"

PROJECT_DIR <-
  "/Users/peter/Downloads/Project-2026-BurnOmicsDB/GSE37069"

FDR_CUTOFF <- 0.05
LOG2FC_CUTOFF <- 1
TOP_VARIABLE_GENES_FOR_PCA <- 500
TOP_GENES_FOR_HEATMAP <- 30
VOLCANO_LABEL_GENE_N <- 12
RAW_QC_BATCH_SIZE <- 10
RLE_BATCH_SIZE <- 50
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

GROUP_LEVELS <- c(
  "Healthy_control",
  "Burn_0_3d",
  "Burn_4_7d",
  "Burn_8_14d",
  "Burn_15_28d",
  "Burn_gt_28d"
)

GROUP_COLORS <- c(
  "Healthy_control" = "#0072B2",
  "Burn_0_3d" = "#D55E00",
  "Burn_4_7d" = "#009E73",
  "Burn_8_14d" = "#CC79A7",
  "Burn_15_28d" = "#E69F00",
  "Burn_gt_28d" = "#56B4E9"
)

VOLCANO_COLORS <- c(
  "Up_significant" = "#D55E00",
  "Down_significant" = "#0072B2",
  "Not_significant" = "#BDBDBD"
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
  "^GSE37069_RAW.*\\.tar$",
  input_search_dirs
)

soft_file <- find_one_file(
  "^GSE37069_family\\.soft.*\\.gz$",
  input_search_dirs
)

paper_file <- find_one_file(
  "pnas.*\\.pdf$",
  input_search_dirs,
  required = FALSE
)

dataset_note_file <- find_one_file(
  "Dataset Note.*\\.rtf$",
  input_search_dirs,
  required = FALSE
)

cat("Raw CEL archive:\n", raw_tar_file, "\n\n")
cat("SOFT metadata:\n", soft_file, "\n\n")
cat("Primary paper PDF:\n", paper_file, "\n\n")
cat("Dataset note:\n", dataset_note_file, "\n\n")

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

if (length(cel_gz_members) != 590) {
  stop(
    paste0(
      "The raw archive should contain 590 CEL.gz files, but ",
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

if (length(existing_cel_gz) != 590) {
  cat("Extracting 590 CEL.gz files from the raw archive.\n")

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

if (length(cel_gz_files) != 590) {
  stop(
    paste0(
      "Expected 590 extracted CEL.gz files, but ",
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

if (length(cel_files) != 590) {
  stop(
    paste0(
      "Expected 590 decompressed CEL files, but ",
      length(cel_files),
      " were found."
    )
  )
}

# Filenames contain strings such as GSM909644_c12003273.CEL.
# Only the leading GSM accession is used as the sample identifier.
# 文件名包含额外内部编号，仅提取开头GSM号。
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

# The SOFT file embeds a processed probe table for each sample. The parser
# reads the file sequentially and skips sample-table rows.
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
        Sex_Original = extract_characteristic(
          characteristics,
          "Sex"
        ),
        Age_Original = extract_characteristic(
          characteristics,
          "age"
        ),
        Hours_Since_Injury_Original =
          extract_characteristic(
            characteristics,
            "hours_since_injury"
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
        Reanalysis_Relation = relation_text,
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

if (nrow(sample_metadata) != 590) {
  stop(
    paste0(
      "The SOFT file should contain 590 samples, but ",
      nrow(sample_metadata),
      " were parsed."
    )
  )
}

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
  allow_missing = FALSE
)

sample_metadata$Hours_Since_Injury <-
  strict_numeric(
    sample_metadata$Hours_Since_Injury_Original,
    "Hours since injury",
    allow_missing = TRUE
  )

sample_metadata$Days_Since_Injury <-
  sample_metadata$Hours_Since_Injury / 24

sex_lower <- tolower(
  trimws(sample_metadata$Sex_Original)
)

sample_metadata$Sex <- ifelse(
  sex_lower == "f",
  "Female",
  ifelse(
    sex_lower == "m",
    "Male",
    NA_character_
  )
)

if (any(is.na(sample_metadata$Sex))) {
  stop(
    "At least one sample has an unsupported sex label."
  )
}

source_lower <- tolower(
  trimws(sample_metadata$Original_Source_Name)
)

sample_metadata$Is_Control <-
  source_lower == "control"

if (sum(sample_metadata$Is_Control) != 37) {
  stop(
    paste0(
      "Expected 37 control arrays, but ",
      sum(sample_metadata$Is_Control),
      " were identified."
    )
  )
}

if (any(
  is.na(
    sample_metadata$Hours_Since_Injury[
      !sample_metadata$Is_Control
    ]
  )
)) {
  stop(
    "At least one burn sample lacks a numeric hours-since-injury value."
  )
}

if (any(
  !is.na(
    sample_metadata$Hours_Since_Injury[
      sample_metadata$Is_Control
    ]
  )
)) {
  stop(
    "At least one control sample unexpectedly has a numeric injury time."
  )
}

# Derive subject IDs from explicit sample titles.
# 根据样本标题提取Subject ID。
control_raw_id <- sub(
  "^.*control\\s+",
  "",
  sample_metadata$Original_Title,
  ignore.case = TRUE
)

control_subject_id <- sub(
  "-[12]$",
  "",
  control_raw_id
)

burn_subject_id <- sub(
  "^.*Subject\\s+([^,]+).*$",
  "\\1",
  sample_metadata$Original_Title,
  ignore.case = TRUE
)

sample_metadata$Patient_ID <- ifelse(
  sample_metadata$Is_Control,
  paste0("Control_", control_subject_id),
  paste0("Burn_", burn_subject_id)
)

invalid_burn_patient <- !sample_metadata$Is_Control &
  !grepl(
    "^Blood, Subject\\s+[^,]+,",
    sample_metadata$Original_Title,
    ignore.case = TRUE
  )

if (any(invalid_burn_patient)) {
  stop(
    paste0(
      "At least one burn Patient_ID could not be safely parsed: ",
      paste(
        sample_metadata$Original_Title[
          invalid_burn_patient
        ],
        collapse = " | "
      )
    )
  )
}

# Define continuous, non-overlapping time windows.
# 定义连续且互不重叠的时间窗口。
sample_metadata$Group <- ifelse(
  sample_metadata$Is_Control,
  "Healthy_control",
  ifelse(
    sample_metadata$Hours_Since_Injury <= 72,
    "Burn_0_3d",
    ifelse(
      sample_metadata$Hours_Since_Injury <= 168,
      "Burn_4_7d",
      ifelse(
        sample_metadata$Hours_Since_Injury <= 336,
        "Burn_8_14d",
        ifelse(
          sample_metadata$Hours_Since_Injury <= 672,
          "Burn_15_28d",
          "Burn_gt_28d"
        )
      )
    )
  )
)

sample_metadata$Group <- factor(
  sample_metadata$Group,
  levels = GROUP_LEVELS
)

sample_metadata$Time_or_Stage <- ifelse(
  sample_metadata$Group == "Healthy_control",
  "Healthy control",
  ifelse(
    sample_metadata$Group == "Burn_0_3d",
    "0-3 days after injury",
    ifelse(
      sample_metadata$Group == "Burn_4_7d",
      "4-7 days after injury",
      ifelse(
        sample_metadata$Group == "Burn_8_14d",
        "8-14 days after injury",
        ifelse(
          sample_metadata$Group == "Burn_15_28d",
          "15-28 days after injury",
          ">28 days after injury"
        )
      )
    )
  )
)

sample_metadata$Tissue <-
  "Peripheral blood leukocytes"

sample_metadata$Treatment <- NA_character_
sample_metadata$Outcome <- NA_character_
sample_metadata$Data_Type <- "Affymetrix raw CEL"
sample_metadata$Is_Pooled <- FALSE
sample_metadata$Metadata_Confidence <-
  "Direct_from_GEO_SOFT"

subject_array_count <- table(
  sample_metadata$Patient_ID
)

sample_metadata$Arrays_Per_Subject <- as.integer(
  subject_array_count[
    sample_metadata$Patient_ID
  ]
)

sample_metadata$Is_Repeated_Measure <-
  sample_metadata$Arrays_Per_Subject > 1

sample_metadata$Is_Paired <-
  sample_metadata$Is_Repeated_Measure

sample_metadata$Quality_Notes <- paste0(
  "Peripheral-blood leukocyte bulk expression; longitudinal arrays from the ",
  "same SOFT-derived subject ID are modeled with Patient_ID blocking; the ",
  "main model adjusts for continuous age and sex; age ranges differ between ",
  "burn patients and controls; clinical variables such as TBSA, inhalation ",
  "injury, treatment, infection, and survival are not available in the GEO ",
  "sample metadata used here; bulk-blood differences may reflect both ",
  "intracellular regulation and leukocyte-composition changes; GSE19743 is ",
  "a substantially overlapping cohort and is not independent replication."
)

sample_metadata$Quality_Notes[
  sample_metadata$Group == "Burn_gt_28d"
] <- paste0(
  sample_metadata$Quality_Notes[
    sample_metadata$Group == "Burn_gt_28d"
  ],
  " The >28-day group spans a broad range of late follow-up times."
)

sample_metadata$Expected_CEL_GZ_Filename <-
  basename(sample_metadata$CEL_GZ_URL)

# Verify CEL and SOFT identifiers.
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

# Align all files and metadata to the SOFT order.
# 按SOFT顺序统一排列文件和元数据。
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

actual_group_counts <- table(
  sample_metadata$Group
)

expected_group_counts <- c(
  Healthy_control = 37,
  Burn_0_3d = 209,
  Burn_4_7d = 62,
  Burn_8_14d = 71,
  Burn_15_28d = 98,
  Burn_gt_28d = 113
)

if (!all(
  actual_group_counts[
    names(expected_group_counts)
  ] == expected_group_counts
)) {
  stop(
    paste0(
      "Unexpected time-window counts: ",
      paste(
        names(actual_group_counts),
        actual_group_counts,
        sep = "=",
        collapse = ", "
      )
    )
  )
}

burn_patient_ids <- unique(
  sample_metadata$Patient_ID[
    !sample_metadata$Is_Control
  ]
)

control_patient_ids <- unique(
  sample_metadata$Patient_ID[
    sample_metadata$Is_Control
  ]
)

burn_patient_n <- length(burn_patient_ids)
control_patient_n <- length(control_patient_ids)

# The SOFT sample titles produce 248 burn IDs although GEO reports 244.
# Preserve the explicit IDs and report the inconsistency.
if (burn_patient_n != 248) {
  stop(
    paste0(
      "Expected 248 distinct burn Subject IDs from the uploaded SOFT titles, ",
      "but ",
      burn_patient_n,
      " were derived."
    )
  )
}

if (control_patient_n != 35) {
  stop(
    paste0(
      "Expected 35 healthy-control subject IDs, but ",
      control_patient_n,
      " were derived."
    )
  )
}

time_window_patient_n <- vapply(
  GROUP_LEVELS,
  function(group_code) {
    length(
      unique(
        sample_metadata$Patient_ID[
          sample_metadata$Group == group_code
        ]
      )
    )
  },
  numeric(1)
)

metadata_audit <- data.frame(
  Check = c(
    "Raw CEL file count",
    "SOFT sample count",
    "Healthy-control array count",
    "Burn array count",
    "Healthy-control subject IDs derived from SOFT titles",
    "Burn subject IDs derived from SOFT titles",
    "Burn patients reported by GEO",
    "Total SOFT-derived subject IDs",
    "Minimum burn arrays per subject",
    "Maximum burn arrays per subject",
    "Burn 0-3-day arrays",
    "Burn 4-7-day arrays",
    "Burn 8-14-day arrays",
    "Burn 15-28-day arrays",
    "Burn >28-day arrays",
    "Missing age values",
    "Missing sex values",
    "Pooled samples",
    "Related GEO dataset",
    "Independent replication relative to GSE19743"
  ),
  Value = c(
    length(cel_files),
    nrow(sample_metadata),
    actual_group_counts["Healthy_control"],
    sum(!sample_metadata$Is_Control),
    control_patient_n,
    burn_patient_n,
    244,
    length(unique(sample_metadata$Patient_ID)),
    min(
      sample_metadata$Arrays_Per_Subject[
        !sample_metadata$Is_Control
      ]
    ),
    max(
      sample_metadata$Arrays_Per_Subject[
        !sample_metadata$Is_Control
      ]
    ),
    actual_group_counts["Burn_0_3d"],
    actual_group_counts["Burn_4_7d"],
    actual_group_counts["Burn_8_14d"],
    actual_group_counts["Burn_15_28d"],
    actual_group_counts["Burn_gt_28d"],
    sum(is.na(sample_metadata$Age)),
    sum(is.na(sample_metadata$Sex)),
    sum(sample_metadata$Is_Pooled),
    "GSE19743",
    "No"
  ),
  Interpretation = c(
    "Directly observed in GSE37069_RAW.tar.",
    "Directly parsed from the streaming SOFT parser.",
    "Thirty-seven arrays represent thirty-five control subjects.",
    "Longitudinal burn arrays.",
    "Two control subjects have two arrays.",
    "Subject IDs are retained exactly as parsed and are not inferred or merged.",
    "The GEO Series overall design reports 244 burn patients.",
    "Thirty-five control IDs plus 248 burn IDs.",
    "Observed among SOFT-derived burn subject IDs.",
    "Observed among SOFT-derived burn subject IDs.",
    "Hours since injury <=72.",
    "Hours since injury >72 and <=168.",
    "Hours since injury >168 and <=336.",
    "Hours since injury >336 and <=672.",
    "Hours since injury >672.",
    "Age is complete for the primary model.",
    "Sex is complete for the primary model.",
    "No pooling is reported.",
    "Substantial cohort overlap is documented.",
    "GSE19743 and GSE37069 must not be counted as independent evidence."
  ),
  stringsAsFactors = FALSE
)

write.csv(
  metadata_audit,
  file = file.path(
    QC_DIR,
    "GSE37069_metadata_consistency_check.csv"
  ),
  row.names = FALSE
)

cat("Sample-group counts:\n")
print(actual_group_counts)

cat(
  "\nSOFT-derived subjects: ",
  burn_patient_n,
  " burn IDs and ",
  control_patient_n,
  " control IDs.\n\n",
  sep = ""
)

cat("Metadata audit:\n")
print(metadata_audit)


# ---- 04. Batch-wise raw CEL QC / 分批进行原始CEL质量控制 ----

raw_qc_metrics <- data.frame(
  Sample_ID = sample_metadata$Sample_ID,
  Patient_ID = sample_metadata$Patient_ID,
  Group = sample_metadata$Group,
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
    "GSE37069_raw_array_QC_metrics.csv"
  ),
  row.names = FALSE
)

p_raw_qc <- ggplot(
  raw_qc_metrics,
  aes(
    x = Group,
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
    size = 1.15,
    alpha = 0.58
  ) +
  scale_color_manual(
    values = GROUP_COLORS,
    drop = FALSE
  ) +
  scale_fill_manual(
    values = GROUP_COLORS,
    drop = FALSE
  ) +
  labs(
    title = "GSE37069 raw CEL intensity QC",
    subtitle =
      "Each point is the median raw log2 probe intensity of one array",
    x = "Sample group",
    y = "Median raw log2 probe intensity"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      angle = 25,
      hjust = 1
    )
  )

print(p_raw_qc)

ggsave(
  filename = file.path(
    FIGURES_DIR,
    "01_GSE37069_raw_data_QC.png"
  ),
  plot = p_raw_qc,
  width = 10,
  height = 5.8,
  dpi = 300
)


# ---- 05. Memory-conscious RMA normalization with checkpoint / 带断点的内存友好RMA ----

RMA_CHECKPOINT_FILE <- file.path(
  OBJECTS_DIR,
  "GSE37069_RMA_probe_expression_checkpoint.rds"
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
    "Starting RMA normalization for 590 CEL files. ",
    "This can take several hours and may use substantial memory.\n"
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

if (ncol(rma_probe_expression) != 590) {
  stop(
    "The RMA expression matrix does not contain 590 arrays."
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
    "GSE37069_RMA_probe_level_expression.csv.gz"
  )

  if (!file.exists(probe_csv_file)) {
    cat(
      "Writing the large probe-level RMA CSV file.\n"
    )

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
    "GSE37069_RMA_array_QC_metrics.csv"
  ),
  row.names = FALSE
)

p_normalized <- ggplot(
  rma_qc_metrics,
  aes(
    x = Group,
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
    size = 1.15,
    alpha = 0.58
  ) +
  scale_color_manual(
    values = GROUP_COLORS,
    drop = FALSE
  ) +
  scale_fill_manual(
    values = GROUP_COLORS,
    drop = FALSE
  ) +
  labs(
    title =
      "GSE37069 RMA-normalized expression distribution",
    subtitle =
      "Each point is the median RMA log2 expression of one array",
    x = "Sample group",
    y = "Median RMA log2 expression"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      angle = 25,
      hjust = 1
    )
  )

print(p_normalized)

ggsave(
  filename = file.path(
    FIGURES_DIR,
    "02_GSE37069_normalized_expression_distribution.png"
  ),
  plot = p_normalized,
  width = 10,
  height = 5.8,
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
    "GSE37069_probe_mapping_and_selection.csv.gz"
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
    "GSE37069_probe_mapping_summary.csv"
  ),
  row.names = FALSE
)

cat("Probe-mapping summary:\n")
print(mapping_summary)

normalized_gene_connection <- gzfile(
  file.path(
    RESULTS_DIR,
    "GSE37069_normalized_expression.csv.gz"
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

# Release the large probe-level matrix after gene-level selection.
# 基因级矩阵建立后释放大探针矩阵。
rm(
  rma_probe_expression,
  probe_row_medians
)

invisible(gc())


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

write.csv(
  pca_table,
  file = file.path(
    RESULTS_DIR,
    "GSE37069_PCA_scores.csv"
  ),
  row.names = FALSE
)

p_pca <- ggplot(
  pca_table,
  aes(
    x = PC1,
    y = PC2,
    color = Group
  )
) +
  geom_point(
    size = 1.7,
    alpha = 0.62
  ) +
  scale_color_manual(
    values = GROUP_COLORS,
    drop = FALSE
  ) +
  labs(
    title = "GSE37069 PCA",
    subtitle = paste0(
      "Top ",
      n_pca_genes,
      " variable gene-level RMA features; 590 arrays"
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
    color = "Group"
  ) +
  theme_classic(base_size = 12)

print(p_pca)

ggsave(
  filename = file.path(
    FIGURES_DIR,
    "03_GSE37069_PCA.png"
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
    "GSE37069_sample_correlation.csv.gz"
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
  row.names = sample_metadata$Sample_ID
)

correlation_annotation_colors <- list(
  Group = GROUP_COLORS
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
  main = "GSE37069 sample correlation",
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
  main = "GSE37069 sample correlation",
  filename = file.path(
    FIGURES_DIR,
    "04_GSE37069_sample_correlation_heatmap.png"
  ),
  width = 12,
  height = 11,
  silent = TRUE
)


# ---- 08. Repeated-measures differential-expression model / 重复测量差异表达模型 ----

sample_metadata$Group <- factor(
  sample_metadata$Group,
  levels = GROUP_LEVELS
)

sample_metadata$Sex <- factor(
  sample_metadata$Sex,
  levels = c("Female", "Male")
)

sample_metadata$Age_Centered <-
  sample_metadata$Age -
    mean(sample_metadata$Age)

design <- model.matrix(
  ~ 0 + Group + Age_Centered + Sex,
  data = sample_metadata
)

colnames(design) <- sub(
  "^Group",
  "",
  colnames(design)
)

rownames(design) <-
  sample_metadata$Sample_ID

if (qr(design)$rank < ncol(design)) {
  stop(
    "The differential-expression design matrix is not full rank."
  )
}

patient_block <- factor(
  sample_metadata$Patient_ID
)

cat("Design-matrix columns:\n")
print(colnames(design))

correlation_fit <- limma::duplicateCorrelation(
  gene_expression,
  design = design,
  block = patient_block
)

within_patient_correlation <-
  correlation_fit$consensus.correlation

if (!is.finite(within_patient_correlation)) {
  stop(
    "Within-subject correlation could not be estimated."
  )
}

cat(
  "Estimated within-subject correlation: ",
  round(within_patient_correlation, 6),
  "\n",
  sep = ""
)

fit <- limma::lmFit(
  gene_expression,
  design = design,
  block = patient_block,
  correlation = within_patient_correlation
)

contrast_matrix <- limma::makeContrasts(
  GSE37069_Burn0to3d_vs_HealthyControl =
    Burn_0_3d - Healthy_control,
  GSE37069_Burn4to7d_vs_HealthyControl =
    Burn_4_7d - Healthy_control,
  GSE37069_Burn8to14d_vs_HealthyControl =
    Burn_8_14d - Healthy_control,
  GSE37069_Burn15to28d_vs_HealthyControl =
    Burn_15_28d - Healthy_control,
  GSE37069_BurnGt28d_vs_HealthyControl =
    Burn_gt_28d - Healthy_control,
  levels = design
)

fit_contrasts <- limma::contrasts.fit(
  fit,
  contrasts = contrast_matrix
)

fit_contrasts <- limma::eBayes(
  fit_contrasts,
  robust = TRUE,
  trend = TRUE
)

contrast_definitions <- data.frame(
  Contrast_ID = colnames(contrast_matrix),
  Contrast_Label = c(
    "Burn 0-3 days vs healthy control",
    "Burn 4-7 days vs healthy control",
    "Burn 8-14 days vs healthy control",
    "Burn 15-28 days vs healthy control",
    "Burn >28 days vs healthy control"
  ),
  Case_Group = c(
    "Burn 0-3 days",
    "Burn 4-7 days",
    "Burn 8-14 days",
    "Burn 15-28 days",
    "Burn >28 days"
  ),
  Case_Group_Code = c(
    "Burn_0_3d",
    "Burn_4_7d",
    "Burn_8_14d",
    "Burn_15_28d",
    "Burn_gt_28d"
  ),
  Control_Group = "Healthy control",
  Control_Group_Code = "Healthy_control",
  Sample_Context =
    "Peripheral blood leukocytes after severe burn injury",
  Time_or_Stage = c(
    "0-3 days after injury",
    "4-7 days after injury",
    "8-14 days after injury",
    "15-28 days after injury",
    ">28 days after injury"
  ),
  stringsAsFactors = FALSE
)

contrast_definitions$Case_N <- vapply(
  contrast_definitions$Case_Group_Code,
  function(group_code) {
    sum(sample_metadata$Group == group_code)
  },
  numeric(1)
)

contrast_definitions$Control_N <- vapply(
  contrast_definitions$Control_Group_Code,
  function(group_code) {
    sum(sample_metadata$Group == group_code)
  },
  numeric(1)
)

contrast_definitions$Case_Patient_N <- vapply(
  contrast_definitions$Case_Group_Code,
  function(group_code) {
    length(
      unique(
        sample_metadata$Patient_ID[
          sample_metadata$Group == group_code
        ]
      )
    )
  },
  numeric(1)
)

contrast_definitions$Control_Patient_N <- vapply(
  contrast_definitions$Control_Group_Code,
  function(group_code) {
    length(
      unique(
        sample_metadata$Patient_ID[
          sample_metadata$Group == group_code
        ]
      )
    )
  },
  numeric(1)
)

cat("Contrast definitions:\n")
print(contrast_definitions)


# ---- 09. Create complete gene-level result tables / 创建完整基因级结果表 ----

extract_contrast_results <- function(
  contrast_id,
  fitted_object,
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
    Contrast_ID = contrast_id,
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

contrast_result_list <- lapply(
  contrast_definitions$Contrast_ID,
  function(contrast_id) {
    contrast_info <- contrast_definitions[
      contrast_definitions$Contrast_ID ==
        contrast_id,
    ]

    extract_contrast_results(
      contrast_id = contrast_id,
      fitted_object = fit_contrasts,
      annotation_data = gene_annotation,
      contrast_info = contrast_info
    )
  }
)

all_gene_results <- do.call(
  rbind,
  contrast_result_list
)

row.names(all_gene_results) <- NULL

all_results_connection <- gzfile(
  file.path(
    RESULTS_DIR,
    "GSE37069_all_gene_results.csv.gz"
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


# ---- 10. Volcano plots / 火山图 ----

for (
  contrast_id in
    contrast_definitions$Contrast_ID
) {
  volcano_data <- all_gene_results[
    all_gene_results$Contrast_ID ==
      contrast_id,
  ]

  contrast_info <- contrast_definitions[
    contrast_definitions$Contrast_ID ==
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

  p_volcano <- ggplot(
    volcano_data,
    aes(
      x = log2FC,
      y = NegLog10_FDR,
      color = DE_Status
    )
  ) +
    geom_point(
      alpha = 0.61,
      size = 1.05
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
        "GSE37069: ",
        contrast_info$Contrast_Label
      ),
      subtitle = paste0(
        "RMA + limma repeated-measures model; ",
        "|log2FC| >= ",
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
        "05_GSE37069_volcano_",
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

for (
  contrast_id in
    contrast_definitions$Contrast_ID
) {
  contrast_info <- contrast_definitions[
    contrast_definitions$Contrast_ID ==
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

  contrast_sample_ids <- sample_metadata$Sample_ID[
    sample_metadata$Group %in%
      c(
        contrast_info$Control_Group_Code,
        contrast_info$Case_Group_Code
      )
  ]

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

  heatmap_annotation <- data.frame(
    Group = sample_metadata$Group[
      match(
        contrast_sample_ids,
        sample_metadata$Sample_ID
      )
    ],
    row.names = contrast_sample_ids
  )

  heatmap_annotation$Group <- factor(
    heatmap_annotation$Group,
    levels = c(
      contrast_info$Control_Group_Code,
      contrast_info$Case_Group_Code
    )
  )

  heatmap_annotation_colors <- list(
    Group = GROUP_COLORS[
      c(
        contrast_info$Control_Group_Code,
        contrast_info$Case_Group_Code
      )
    ]
  )

  heatmap_title <- paste0(
    "Top ",
    n_heatmap_genes,
    " genes: ",
    contrast_info$Contrast_Label
  )

  pheatmap::pheatmap(
    heatmap_z,
    annotation_col = heatmap_annotation,
    annotation_colors =
      heatmap_annotation_colors,
    show_colnames = FALSE,
    show_rownames = TRUE,
    border_color = NA,
    color = EXPRESSION_HEATMAP_COLORS,
    main = heatmap_title,
    silent = FALSE
  )

  pheatmap::pheatmap(
    heatmap_z,
    annotation_col = heatmap_annotation,
    annotation_colors =
      heatmap_annotation_colors,
    show_colnames = FALSE,
    show_rownames = TRUE,
    border_color = NA,
    color = EXPRESSION_HEATMAP_COLORS,
    main = heatmap_title,
    filename = file.path(
      FIGURES_DIR,
      paste0(
        "06_GSE37069_top_differential_genes_heatmap_",
        contrast_id,
        ".png"
      )
    ),
    width = 13,
    height = 10,
    silent = TRUE
  )
}


# ---- 12. Qualitative comparison with the primary publication / 与主要论文进行定性核对 ----

# The primary publication used longitudinal spline trajectories and EDGE,
# whereas BurnOmicsDB uses five fixed windows. The checks below are therefore
# contextual and not direct reproduction.
# 原论文使用纵向样条轨迹和EDGE，本数据库使用固定时间窗，因此仅作定性核对。

publication_gene_context <- data.frame(
  Gene_Symbol = c(
    "HLA-DRA",
    "S100A8",
    "S100A9",
    "CD177",
    "MMP8",
    "SLC2A3",
    "IL10",
    "IGKC",
    "CD3D",
    "CD3E"
  ),
  Published_Context = c(
    "Example of a prolonged human injury response in the primary publication.",
    "Innate-myeloid inflammatory marker.",
    "Innate-myeloid inflammatory marker.",
    "Neutrophil-associated marker.",
    "Neutrophil-associated protease.",
    "Burn-responsive metabolic gene.",
    "Inflammatory signaling context.",
    "Adaptive-humoral immune context.",
    "Adaptive T-cell context.",
    "Adaptive T-cell context."
  ),
  Direct_Replication = FALSE,
  stringsAsFactors = FALSE
)

calculated_publication_genes <- all_gene_results[
  all_gene_results$Gene_Symbol %in%
    publication_gene_context$Gene_Symbol,
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

publication_feature_check <- merge(
  publication_gene_context,
  calculated_publication_genes,
  by = "Gene_Symbol",
  all.x = TRUE,
  sort = FALSE
)

write.csv(
  publication_feature_check,
  file = file.path(
    RESULTS_DIR,
    "GSE37069_publication_feature_gene_check.csv"
  ),
  row.names = FALSE
)

strict_publication_like <- all_gene_results[
  all_gene_results$FDR < 0.001 &
    abs(all_gene_results$log2FC) >= 1,
]

strict_gene_summary <- data.frame(
  Metric = c(
    "Unique genes meeting FDR < 0.001 and absolute log2FC >= 1 in at least one BurnOmicsDB window",
    "Genes reported as burn-responsive in the primary publication"
  ),
  Count = c(
    length(
      unique(
        strict_publication_like$NCBI_Gene_ID
      )
    ),
    3250
  ),
  Directly_Comparable = c(
    FALSE,
    FALSE
  ),
  Explanation = c(
    paste0(
      "BurnOmicsDB uses five fixed time-window contrasts with age and sex ",
      "adjustment and repeated-measures modeling."
    ),
    paste0(
      "The publication used longitudinal spline trajectories, maximum ",
      "deviation from controls, and EDGE permutation testing."
    )
  ),
  stringsAsFactors = FALSE
)

write.csv(
  strict_gene_summary,
  file = file.path(
    RESULTS_DIR,
    "GSE37069_publication_method_context.csv"
  ),
  row.names = FALSE
)

cat(
  "Publication-context files created. These are not direct replication tests.\n"
)


# ---- 13. Create the BurnOmicsDB-ready result table / 创建BurnOmicsDB标准结果表 ----

contrast_match <- match(
  all_gene_results$Contrast_ID,
  contrast_definitions$Contrast_ID
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
    "Pediatric and adult patients with severe burns and healthy controls",
  Tissue = "Peripheral blood",
  Sample_Context =
    contrast_definitions$Sample_Context[
      contrast_match
    ],
  Time_or_Stage =
    contrast_definitions$Time_or_Stage[
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
    "Raw CEL files from GSE37069_RAW.tar",
  Normalization = paste0(
    "RMA: background correction, quantile normalization, ",
    "and median-polish probe-set summarization"
  ),
  Analysis_Method = paste0(
    "limma linear model with robust trend-aware empirical Bayes moderation; ",
    "duplicateCorrelation and Patient_ID blocking for longitudinal and ",
    "repeated arrays; adjusted for continuous age and sex; one representative ",
    "probe per NCBI Gene ID selected by highest mean RMA expression across ",
    "all 590 arrays; positive log2FC represents Case_Group minus Control_Group"
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
  Quality_Notes = paste0(
    "Core longitudinal human blood dataset; GEO reports 244 burn patients, ",
    "whereas SOFT titles yield 248 distinct burn Subject IDs that are retained ",
    "without inferred merging; 37 control arrays represent 35 control subjects; ",
    "the model adjusts for age and sex and accounts for within-subject ",
    "correlation; burn ages extend beyond the control age range, limiting ",
    "covariate overlap; the >28-day window is biologically broad; clinical ",
    "severity, treatment, infection, and outcome variables are not present in ",
    "the GEO sample metadata used here; GSE19743 is a substantially overlapping ",
    "subset and is not independent replication; GEO/SOFT define the tissue as ",
    "white blood cells; bulk-blood differences may reflect both intracellular ",
    "regulation and leukocyte-composition changes"
  ),
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
  Within_Patient_Correlation =
    within_patient_correlation,
  Age_Adjusted = TRUE,
  Sex_Adjusted = TRUE,
  Is_Pooled = FALSE,
  Quality_Grade =
    "Core_high_quality",
  Related_GEO_ID =
    "GSE19743",
  Cohort_Overlap =
    "Substantial",
  Independent_Replication =
    FALSE,
  stringsAsFactors = FALSE
)

database_connection <- gzfile(
  file.path(
    RESULTS_DIR,
    "GSE37069_database_ready_all_genes.csv.gz"
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
    "Age",
    "Sex",
    "Hours_Since_Injury",
    "Days_Since_Injury",
    "Arrays_Per_Subject",
    "Tissue_Original",
    "Reanalysis_Relation",
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
    "GSE37069_sample_metadata.csv"
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

# The full probe-level matrix is intentionally excluded because it has a
# separate checkpoint and optional CSV file.
# 完整探针矩阵已有独立断点文件，不重复写入最终分析RDS。
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
  design_matrix =
    design,
  within_patient_correlation =
    within_patient_correlation,
  contrast_matrix =
    contrast_matrix,
  contrast_definitions =
    contrast_definitions,
  fitted_model =
    fit_contrasts,
  all_gene_results =
    all_gene_results,
  publication_feature_check =
    publication_feature_check,
  publication_method_context =
    strict_gene_summary,
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
    "GSE37069_analysis_objects.rds"
  ),
  compress = "gzip"
)

summary_lines <- c(
  "BurnOmicsDB - GSE37069 analysis summary",
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
  paste0(
    "Healthy-control arrays: ",
    actual_group_counts["Healthy_control"]
  ),
  paste0(
    "Burn arrays: ",
    sum(!sample_metadata$Is_Control)
  ),
  paste0(
    "SOFT-derived healthy-control subject IDs: ",
    control_patient_n
  ),
  paste0(
    "SOFT-derived burn subject IDs: ",
    burn_patient_n
  ),
  "GEO-reported burn patients: 244",
  paste0(
    "Total SOFT-derived subject IDs used for blocking: ",
    length(
      unique(
        sample_metadata$Patient_ID
      )
    )
  ),
  "Pooling: No",
  "",
  "Five time windows:",
  paste0(
    "  Healthy control: ",
    actual_group_counts["Healthy_control"],
    " arrays, ",
    time_window_patient_n["Healthy_control"],
    " subject IDs"
  ),
  paste0(
    "  Burn 0-3 days: ",
    actual_group_counts["Burn_0_3d"],
    " arrays, ",
    time_window_patient_n["Burn_0_3d"],
    " subject IDs"
  ),
  paste0(
    "  Burn 4-7 days: ",
    actual_group_counts["Burn_4_7d"],
    " arrays, ",
    time_window_patient_n["Burn_4_7d"],
    " subject IDs"
  ),
  paste0(
    "  Burn 8-14 days: ",
    actual_group_counts["Burn_8_14d"],
    " arrays, ",
    time_window_patient_n["Burn_8_14d"],
    " subject IDs"
  ),
  paste0(
    "  Burn 15-28 days: ",
    actual_group_counts["Burn_15_28d"],
    " arrays, ",
    time_window_patient_n["Burn_15_28d"],
    " subject IDs"
  ),
  paste0(
    "  Burn >28 days: ",
    actual_group_counts["Burn_gt_28d"],
    " arrays, ",
    time_window_patient_n["Burn_gt_28d"],
    " subject IDs"
  ),
  "",
  paste0(
    "Observed burn hours since injury: ",
    round(
      min(
        sample_metadata$Hours_Since_Injury[
          !sample_metadata$Is_Control
        ]
      ),
      1
    ),
    "-",
    round(
      max(
        sample_metadata$Hours_Since_Injury[
          !sample_metadata$Is_Control
        ]
      ),
      1
    )
  ),
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
    "590 arrays for each NCBI Gene ID"
  ),
  "",
  paste0(
    "Normalization: RMA background correction, quantile normalization, ",
    "and median-polish summarization"
  ),
  paste0(
    "Statistical model: limma with robust trend-aware empirical Bayes; ",
    "duplicateCorrelation blocking Patient_ID; adjusted for continuous age ",
    "and sex"
  ),
  paste0(
    "Estimated within-subject correlation: ",
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
  "Contrasts:",
  paste0(
    "  ",
    contrast_definitions$Contrast_ID,
    ": ",
    contrast_definitions$Contrast_Label,
    " (case arrays = ",
    contrast_definitions$Case_N,
    ", case subject IDs = ",
    contrast_definitions$Case_Patient_N,
    ", control arrays = ",
    contrast_definitions$Control_N,
    ", control subject IDs = ",
    contrast_definitions$Control_Patient_N,
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
    "Positive log2FC means higher expression in the burn time-window group ",
    "than in healthy controls."
  ),
  paste0(
    "No separate low-expression filter was applied after RMA and ",
    "representative-probe selection; all selected unambiguous gene ",
    "representatives were statistically tested and exported."
  ),
  "",
  "Major limitations:",
  paste0(
    "  GEO reports 244 burn patients, but SOFT titles yield 248 distinct ",
    "burn Subject IDs; explicit IDs were retained without inferred merging."
  ),
  paste0(
    "  Burn ages extend beyond the healthy-control age range, so linear age ",
    "adjustment cannot fully correct the limited covariate overlap."
  ),
  paste0(
    "  The >28-day group combines a broad range of late follow-up times."
  ),
  paste0(
    "  TBSA, inhalation injury, treatment, infection, and outcome variables ",
    "are not present in the GEO sample metadata used in this analysis."
  ),
  paste0(
    "  GSE19743 is a substantially overlapping cohort and must not be counted ",
    "as independent replication."
  ),
  paste0(
    "  The primary publication used longitudinal spline trajectories and ",
    "EDGE rather than the five fixed-window contrasts used here."
  ),
  paste0(
    "  Bulk-blood expression differences may reflect both intracellular ",
    "regulation and changes in leukocyte composition."
  ),
  "Quality grade: Core_high_quality"
)

writeLines(
  summary_lines,
  con = file.path(
    RESULTS_DIR,
    "GSE37069_analysis_summary.txt"
  ),
  useBytes = TRUE
)

capture.output(
  sessionInfo(),
  file = file.path(
    RESULTS_DIR,
    "GSE37069_R_sessionInfo.txt"
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
  "HLA-DRA",
  "S100A8",
  "S100A9",
  "CD177",
  "MMP8",
  "SLC2A3",
  "IL6",
  "CXCL8",
  "IL10",
  "IGKC",
  "CD3D",
  "CD3E"
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
    "DE_Status"
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
      contrast_definitions$Contrast_ID
    )
  ),
]

print(gene_check_table)

# End of script / 代码结束
