# ==============================================================================
# BurnOmicsDB: GSE19743 microarray analysis
# BurnOmicsDB：GSE19743微阵列分析
#
# GEO accession / GEO编号:
#   GSE19743
#
# Study / 研究:
#   A large-scale clinical study of gene expression response to severe burn
#   injury.
#   严重烧伤后基因表达反应的大规模临床研究。
#
# Data type / 数据类型:
#   Affymetrix Human Genome U133 Plus 2.0 Array (GPL570)
#   Raw CEL files contained in GSE19743_RAW.tar
#   GPL570 Affymetrix原始CEL文件。
#
# Study structure / 研究结构:
#   Healthy controls:
#     63 arrays from 63 independent subjects.
#   Severe-burn patients:
#     57 patients, each measured at two stages:
#       Early: <11 days after injury
#       Middle: 11-49 days after injury
#     114 burn arrays in total.
#   Overall:
#     177 arrays representing 120 individuals.
#
#   健康对照63人，每人1张芯片；严重烧伤患者57人，
#   每人具有Early和Mid两个重复时间点；共177张芯片、120名受试者。
#
# Prespecified contrasts / 预设比较:
#   1. Early burn versus healthy control
#   2. Middle burn versus healthy control
#   3. Middle burn versus early burn
#
# Positive log2FC always means higher expression in the case group than in the
# control/reference group.
# 正log2FC始终表示Case_Group相对于Control_Group表达更高。
#
# Analysis workflow / 分析流程:
#   Raw CEL
#   -> batch-wise raw-array QC to limit memory use
#   -> RMA background correction
#   -> quantile normalization
#   -> median-polish probe-set summarization
#   -> probe-to-NCBI-Gene-ID mapping
#   -> one representative probe per NCBI Gene ID selected by the highest mean
#      RMA expression across all 177 arrays
#   -> limma linear model
#   -> duplicateCorrelation blocking Patient_ID
#   -> adjustment for continuous age and sex
#   -> three prespecified contrasts
#
# Repeated-measures rationale / 重复测量处理:
#   Each burn patient contributes an Early and a Mid array. These two arrays
#   must not be treated as independent. duplicateCorrelation estimates a common
#   within-patient correlation and lmFit applies it during model fitting.
#   每名烧伤患者提供Early和Mid两张芯片，二者并非独立。
#   duplicateCorrelation估计患者内共同相关性，并在lmFit中应用。
#
# Healthy-control handling / 健康对照处理:
#   Healthy controls have one array each and enter the model as singleton
#   blocks. Unlike the original TANOVA analysis, this script does not duplicate
#   control arrays to create pseudo time courses.
#   健康对照每人只有一张芯片，以单例block进入模型。
#   本代码不复制对照芯片构造伪时间序列。
#
# Covariates / 协变量:
#   The main model adjusts for continuous age and sex.
#   TBSA, inhalation injury, survival, and hospital length of stay are retained
#   in metadata but are not included in the three main contrasts because they
#   are unavailable for healthy controls and are not the primary study question.
#   主模型校正连续年龄和性别。TBSA、吸入伤、生存和住院时长保留在元数据中，
#   但不进入三个核心比较。
#
# Relation to the associated paper / 与原论文的关系:
#   The associated paper developed TANOVA to study age-dependent temporal
#   responses and used specially balanced subsets. BurnOmicsDB instead estimates
#   three simpler overall contrasts using all 57 burn patients and 63 controls.
#   Therefore, this analysis is not a direct reproduction of TANOVA or its
#   age-by-burn interaction findings.
#   原论文重点研究年龄依赖的时间反应并使用平衡子集；
#   本数据库分析使用全部样本估计三个简化总体比较，不直接复现TANOVA。
#
# Probe-selection rule / 代表探针规则:
#   Probes mapping to zero or multiple NCBI Gene IDs are excluded.
#   For each remaining NCBI Gene ID, the probe with the highest mean RMA
#   expression across all 177 arrays is selected before differential testing.
#   This rule is independent of contrast labels and p-values.
#   无Gene ID或对应多个Gene ID的探针被排除；
#   每个Gene ID选择全部177张芯片平均RMA表达最高的探针。
#
# Memory strategy / 内存策略:
#   - Raw QC reads CEL files in small batches instead of loading all raw
#     probe intensities at once.
#   - justRMA() is used to reduce memory compared with retaining a full
#     raw AffyBatch object.
#   - The raw AffyBatch is not saved in the final RDS.
#   原始QC采用分批读取，RMA使用justRMA，并且不保存完整原始AffyBatch。
#
# How to run in RStudio / 如何在RStudio中运行:
#   - Save this script in:
#     /Users/peter/Downloads/Project-2026-BurnOmicsDB/GSE19743/
#   - Keep these local files in the project folder:
#       GSE19743_RAW.tar
#       GSE19743_family.soft.gz
#       pnas.1002757107.pdf
#   - Press Cmd + Shift + O on macOS to open the section outline.
#     在macOS中按Cmd + Shift + O打开代码分区目录。
#   - Run one section at a time using Cmd + Enter.
#     使用Cmd + Enter逐段运行。
#   - During the first run, do not Source the entire script at once.
#     第一次运行时不要直接Source全文。
#
# Output-language rule / 输出语言规则:
#   Comments are bilingual. Figures, CSV files, TXT files, console messages,
#   and error messages are English only.
#   注释使用中英文；图片、CSV、TXT、Console信息和报错全部只使用英文。
#
#【代码运行后笔记】
# GSE19743已经成功完成，可以作为BurnOmicsDB的核心血液数据集纳入。
# 目前不需要修改代码或重跑。唯一建议是在最终Dataset页面的QC说明中标记一个RLE离散度较高但仍被保留的样本，同时保留下面这条解释：
# This database reports overall age- and sex-adjusted burn-stage contrasts and does not reproduce the age-dependent TANOVA interaction analysis of the original publication.
# ==============================================================================


# ---- 00. Install required packages once / 首次安装所需软件包 ----

# Run this section only if packages are missing.
# 仅在软件包缺失时运行本节。

options(
  timeout = 1800
)

options(
  download.file.method = "libcurl"
)

options(
  repos = c(
    CRAN = "https://cloud.r-project.org"
  )
)

# This mirror can be changed if it is slow or unavailable.
# 如果该镜像速度慢或不可用，可以自行替换。
options(
  BioC_mirror =
    "https://mirrors.tuna.tsinghua.edu.cn/bioconductor"
)

if (!requireNamespace(
  "BiocManager",
  quietly = TRUE
)) {
  install.packages(
    "BiocManager",
    repos = "https://cloud.r-project.org"
  )
}

cran_packages <- c(
  "ggplot2",
  "ggrepel",
  "pheatmap",
  "R.utils"
)

missing_cran_packages <- cran_packages[
  !vapply(
    cran_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(
  missing_cran_packages
) > 0) {
  install.packages(
    missing_cran_packages,
    repos = "https://cloud.r-project.org",
    dependencies = TRUE
  )
}

bioconductor_packages <- c(
  "affy",
  "limma",
  "AnnotationDbi",
  "hgu133plus2cdf",
  "hgu133plus2.db",
  "org.Hs.eg.db"
)

missing_bioconductor_packages <- bioconductor_packages[
  !vapply(
    bioconductor_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(
  missing_bioconductor_packages
) > 0) {
  BiocManager::install(
    missing_bioconductor_packages,
    ask = FALSE,
    update = FALSE
  )
}

all_required_packages <- c(
  cran_packages,
  bioconductor_packages
)

package_check <- vapply(
  all_required_packages,
  requireNamespace,
  logical(1),
  quietly = TRUE
)

print(
  package_check
)

if (!all(
  package_check
)) {
  stop(
    paste0(
      "Package installation is incomplete. Missing packages: ",
      paste(
        names(
          package_check
        )[
          !package_check
        ],
        collapse = ", "
      )
    )
  )
}

cat(
  "\nAll required packages are installed successfully.\n\n"
)


# ---- 01. Project settings and package loading / 项目设置与软件包加载 ----

GEO_ID <- "GSE19743"

PROJECT_DIR <-
  "/Users/peter/Downloads/Project-2026-BurnOmicsDB/GSE19743"

FDR_CUTOFF <- 0.05
LOG2FC_CUTOFF <- 1
TOP_VARIABLE_GENES_FOR_PCA <- 500
TOP_GENES_FOR_HEATMAP <- 30
VOLCANO_LABEL_GENE_N <- 12
RAW_QC_BATCH_SIZE <- 10

RESULTS_DIR <- file.path(
  PROJECT_DIR,
  "Results"
)

FIGURES_DIR <- file.path(
  PROJECT_DIR,
  "Figures"
)

QC_DIR <- file.path(
  PROJECT_DIR,
  "QC"
)

OBJECTS_DIR <- file.path(
  PROJECT_DIR,
  "R_objects"
)

INPUT_DIR <- file.path(
  PROJECT_DIR,
  "Input"
)

CEL_GZ_DIR <- file.path(
  INPUT_DIR,
  "CEL_gz"
)

CEL_DIR <- file.path(
  INPUT_DIR,
  "CEL"
)

dir.create(
  RESULTS_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  FIGURES_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  QC_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  OBJECTS_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  INPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  CEL_GZ_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  CEL_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

required_packages <- c(
  "ggplot2",
  "ggrepel",
  "pheatmap",
  "R.utils",
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

if (length(
  missing_packages
) > 0) {
  stop(
    paste0(
      "Missing packages: ",
      paste(
        missing_packages,
        collapse = ", "
      ),
      ". Run Section 00 before continuing."
    )
  )
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(R.utils)
  library(affy)
  library(limma)
  library(AnnotationDbi)
  library(hgu133plus2cdf)
  library(hgu133plus2.db)
  library(org.Hs.eg.db)
})

set.seed(
  2026
)

GROUP_COLORS <- c(
  "Healthy_control" = "#0072B2",
  "Early_burn" = "#D55E00",
  "Mid_burn" = "#009E73"
)

VOLCANO_COLORS <- c(
  "Up_significant" = "#D55E00",
  "Down_significant" = "#0072B2",
  "Not_significant" = "#BDBDBD"
)

EXPRESSION_HEATMAP_COLORS <-
  grDevices::colorRampPalette(
    c(
      "#0072B2",
      "#F7F7F7",
      "#D55E00"
    )
  )(101)

GROUP_LEVELS <- c(
  "Healthy_control",
  "Early_burn",
  "Mid_burn"
)

cat(
  "Project directory:\n",
  PROJECT_DIR,
  "\n\n"
)


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
        if (!dir.exists(
          directory
        )) {
          return(
            character(0)
          )
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

  hits <- unique(
    hits
  )

  if (length(
    hits
  ) == 0) {
    if (required) {
      stop(
        paste0(
          "No input file matched this pattern: ",
          pattern
        )
      )
    }

    return(
      NA_character_
    )
  }

  if (length(
    hits
  ) > 1) {
    message(
      "Multiple files matched. The first file will be used:\n",
      paste(
        hits,
        collapse = "\n"
      )
    )
  }

  normalizePath(
    hits[1],
    mustWork = TRUE
  )
}

raw_tar_file <- find_one_file(
  "^GSE19743_RAW.*\\.tar$",
  input_search_dirs
)

soft_file <- find_one_file(
  "^GSE19743_family\\.soft.*\\.gz$",
  input_search_dirs
)

paper_file <- find_one_file(
  ".*\\.pdf$",
  input_search_dirs,
  required = FALSE
)

cat(
  "Raw CEL archive:\n",
  raw_tar_file,
  "\n\n"
)

cat(
  "SOFT metadata:\n",
  soft_file,
  "\n\n"
)

cat(
  "Paper PDF:\n",
  paper_file,
  "\n\n"
)

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

if (length(
  cel_gz_members
) != 177) {
  stop(
    paste0(
      "The raw archive should contain 177 CEL.gz files, but ",
      length(
        cel_gz_members
      ),
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

if (length(
  existing_cel_gz
) != 177) {
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

if (length(
  cel_gz_files
) != 177) {
  stop(
    paste0(
      "Expected 177 extracted CEL.gz files, but ",
      length(
        cel_gz_files
      ),
      " were found."
    )
  )
}

cat(
  "Decompressing CEL files when necessary.\n"
)

for (
  cel_gz_file in cel_gz_files
) {
  cel_filename <- sub(
    "\\.gz$",
    "",
    basename(
      cel_gz_file
    ),
    ignore.case = TRUE
  )

  destination_file <- file.path(
    CEL_DIR,
    cel_filename
  )

  if (!file.exists(
    destination_file
  )) {
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

if (length(
  cel_files
) != 177) {
  stop(
    paste0(
      "Expected 177 decompressed CEL files, but ",
      length(
        cel_files
      ),
      " were found."
    )
  )
}

cel_sample_ids <- sub(
  "\\.CEL$",
  "",
  basename(
    cel_files
  ),
  ignore.case = TRUE
)

if (anyDuplicated(
  cel_sample_ids
) > 0) {
  stop(
    "Duplicated GSM identifiers were detected among the CEL filenames."
  )
}

cat(
  "CEL archive validation completed. ",
  length(
    cel_files
  ),
  " CEL files are available.\n\n",
  sep = ""
)


# ---- 03. Stream-parse GEO SOFT metadata and validate study design / 流式解析SOFT元数据并检查设计 ----

# The SOFT file contains approximately 9.7 million lines because processed
# probe tables are embedded for all 177 samples. This parser reads the file in
# chunks and ignores probe-table rows to avoid loading the complete SOFT file
# into memory.
# SOFT文件因包含177个样本的处理后探针表而接近970万行。
# 本函数采用分块流式读取并跳过探针表，避免整体载入内存。

parse_geo_soft_samples_streaming <- function(
  soft_gz_file,
  chunk_size = 50000
) {
  connection <- gzfile(
    soft_gz_file,
    open = "rt"
  )

  on.exit(
    close(
      connection
    ),
    add = TRUE
  )

  records <- list()
  current <- NULL
  in_sample_table <- FALSE

  finalize_current <- function(
    current_record
  ) {
    if (is.null(
      current_record
    )) {
      return(
        NULL
      )
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

    if (length(
      chunk
    ) == 0) {
      break
    }

    for (
      line in chunk
    ) {
      if (grepl(
        "^\\^SAMPLE\\s*=",
        line
      )) {
        finalized <- finalize_current(
          current
        )

        if (!is.null(
          finalized
        )) {
          records[[
            length(
              records
            ) + 1
          ]] <- finalized
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

      if (is.null(
        current
      )) {
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

      if (grepl(
        "^!Sample_title\\s*=",
        line
      )) {
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

  finalized <- finalize_current(
    current
  )

  if (!is.null(
    finalized
  )) {
    records[[
      length(
        records
      ) + 1
    ]] <- finalized
  }

  if (length(
    records
  ) == 0) {
    stop(
      paste0(
        "No SAMPLE records were parsed from the SOFT file. ",
        "Verify that soft_file points to the unmodified ",
        "GSE19743_family.soft.gz file."
      )
    )
  }

  extract_characteristic <- function(
    characteristics,
    key
  ) {
    key_prefix <- paste0(
      tolower(
        key
      ),
      ":"
    )

    matching_index <- which(
      startsWith(
        tolower(
          characteristics
        ),
        key_prefix
      )
    )

    if (length(
      matching_index
    ) == 0) {
      return(
        NA_character_
      )
    }

    trimws(
      sub(
        "^[^:]+:\\s*",
        "",
        characteristics[
          matching_index[1]
        ]
      )
    )
  }

  first_or_na <- function(
    value
  ) {
    if (is.null(
      value
    ) ||
        length(
          value
        ) == 0) {
      return(
        NA_character_
      )
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

      biosample_id <- ifelse(
        grepl(
          "biosample/",
          relation_text,
          ignore.case = TRUE
        ),
        sub(
          ".*biosample/([^| ]+).*",
          "\\1",
          relation_text,
          ignore.case = TRUE
        ),
        NA_character_
      )

      sra_experiment <- ifelse(
        grepl(
          "sra\\?term=",
          relation_text,
          ignore.case = TRUE
        ),
        sub(
          ".*sra\\?term=([^| ]+).*",
          "\\1",
          relation_text,
          ignore.case = TRUE
        ),
        NA_character_
      )

      data.frame(
        GEO_ID = GEO_ID,
        Sample_ID = first_or_na(
          record$Sample_ID
        ),
        Sample_Name = first_or_na(
          record$Sample_ID
        ),
        Original_Title = first_or_na(
          record$Original_Title
        ),
        Original_Source_Name = first_or_na(
          record$Original_Source_Name
        ),
        Sample_Type_Original =
          extract_characteristic(
            characteristics,
            "sample type"
          ),
        Age_Original =
          extract_characteristic(
            characteristics,
            "age (yrs)"
          ),
        Age_Group_Original =
          extract_characteristic(
            characteristics,
            "age group"
          ),
        Sex_Original =
          extract_characteristic(
            characteristics,
            "sex"
          ),
        Hours_Post_Injury_Original =
          extract_characteristic(
            characteristics,
            "hours post injury"
          ),
        Sampling_Time_Group_Original =
          extract_characteristic(
            characteristics,
            "sampling time group"
          ),
        TBSA_Original =
          extract_characteristic(
            characteristics,
            "total burn surface area (tbsa, %)"
          ),
        Inhalation_Injury_Original =
          extract_characteristic(
            characteristics,
            "injury inhalation"
          ),
        Survival_Original =
          extract_characteristic(
            characteristics,
            "survival"
          ),
        Hospital_Length_of_Stay_Original =
          extract_characteristic(
            characteristics,
            "hospital length of stay (days)"
          ),
        Platform_ID = first_or_na(
          record$Platform_ID
        ),
        Original_Characteristics = paste(
          characteristics,
          collapse = " | "
        ),
        BioSample_ID = biosample_id,
        SRA_Experiment = sra_experiment,
        CEL_GZ_URL = first_or_na(
          record$Supplementary_File
        ),
        stringsAsFactors = FALSE
      )
    }
  )

  do.call(
    rbind,
    sample_list
  )
}

sample_metadata <-
  parse_geo_soft_samples_streaming(
    soft_file
  )

if (nrow(
  sample_metadata
) != 177) {
  stop(
    paste0(
      "The SOFT file should contain 177 samples, but ",
      nrow(
        sample_metadata
      ),
      " were parsed."
    )
  )
}

sample_type_lower <- tolower(
  trimws(
    sample_metadata$Sample_Type_Original
  )
)

time_group_lower <- tolower(
  trimws(
    sample_metadata$Sampling_Time_Group_Original
  )
)

sample_metadata$Group <- ifelse(
  sample_type_lower == "control",
  "Healthy_control",
  ifelse(
    sample_type_lower == "burn" &
      time_group_lower == "early",
    "Early_burn",
    ifelse(
      sample_type_lower == "burn" &
        time_group_lower == "mid",
      "Mid_burn",
      NA_character_
    )
  )
)

if (any(
  is.na(
    sample_metadata$Group
  )
)) {
  stop(
    paste0(
      "Unexpected sample type or time-group labels were found: ",
      paste(
        unique(
          paste(
            sample_metadata$Sample_Type_Original[
              is.na(
                sample_metadata$Group
              )
            ],
            sample_metadata$Sampling_Time_Group_Original[
              is.na(
                sample_metadata$Group
              )
            ],
            sep = " / "
          )
        ),
        collapse = ", "
      )
    )
  )
}

sample_metadata$Group <- factor(
  sample_metadata$Group,
  levels = GROUP_LEVELS
)

sample_metadata$Patient_ID <- ifelse(
  sample_metadata$Group ==
    "Healthy_control",
  sample_metadata$Original_Title,
  sub(
    "_(Early|Mid)$",
    "",
    sample_metadata$Original_Title
  )
)

if (any(
  is.na(
    sample_metadata$Patient_ID
  )
) ||
    any(
      sample_metadata$Patient_ID == ""
    )) {
  stop(
    "At least one Patient_ID could not be derived from the sample title."
  )
}

strict_numeric <- function(
  character_values,
  field_name,
  allow_missing = TRUE
) {
  cleaned <- trimws(
    character_values
  )

  cleaned[
    tolower(
      cleaned
    ) %in% c(
      "",
      "--",
      "unknown",
      "na",
      "n/a",
      "not reported"
    )
  ] <- NA_character_

  numeric_values <- suppressWarnings(
    as.numeric(
      cleaned
    )
  )

  invalid <- !is.na(
    cleaned
  ) &
    is.na(
      numeric_values
    )

  if (any(
    invalid
  )) {
    stop(
      paste0(
        "Unsupported non-numeric values were found in ",
        field_name,
        ": ",
        paste(
          unique(
            character_values[
              invalid
            ]
          ),
          collapse = ", "
        )
      )
    )
  }

  if (!allow_missing &&
      any(
        is.na(
          numeric_values
        )
      )) {
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

sample_metadata$Hours_Post_Injury <-
  strict_numeric(
    sample_metadata$Hours_Post_Injury_Original,
    "Hours post injury",
    allow_missing = TRUE
  )

sample_metadata$Days_Post_Injury <-
  sample_metadata$Hours_Post_Injury /
    24

sample_metadata$TBSA_Percent <-
  strict_numeric(
    sample_metadata$TBSA_Original,
    "TBSA",
    allow_missing = TRUE
  )

sample_metadata$Hospital_Length_of_Stay_Days <-
  strict_numeric(
    sample_metadata$Hospital_Length_of_Stay_Original,
    "Hospital length of stay",
    allow_missing = TRUE
  )

sex_lower <- tolower(
  trimws(
    sample_metadata$Sex_Original
  )
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

if (any(
  is.na(
    sample_metadata$Sex
  )
)) {
  stop(
    "At least one sample has an unsupported sex label."
  )
}

age_group_lower <- tolower(
  trimws(
    sample_metadata$Age_Group_Original
  )
)

sample_metadata$Age_Group <- ifelse(
  age_group_lower == "ped",
  "Pediatric",
  ifelse(
    age_group_lower == "adult",
    "Adult",
    NA_character_
  )
)

if (any(
  is.na(
    sample_metadata$Age_Group
  )
)) {
  stop(
    "At least one sample has an unsupported age-group label."
  )
}

sample_metadata$Time_or_Stage <- ifelse(
  sample_metadata$Group ==
    "Healthy_control",
  "Healthy control",
  ifelse(
    sample_metadata$Group ==
      "Early_burn",
    "Early stage (<11 days after injury)",
    "Middle stage (11-49 days after injury)"
  )
)

sample_metadata$Tissue <-
  "Peripheral blood leukocytes"

sample_metadata$Treatment <-
  NA_character_

sample_metadata$Outcome <- ifelse(
  tolower(
    sample_metadata$Survival_Original
  ) == "yes",
  "Survived",
  ifelse(
    tolower(
      sample_metadata$Survival_Original
    ) == "no",
    "Died",
    NA_character_
  )
)

sample_metadata$Inhalation_Injury <- ifelse(
  tolower(
    sample_metadata$Inhalation_Injury_Original
  ) == "yes",
  "Yes",
  ifelse(
    tolower(
      sample_metadata$Inhalation_Injury_Original
    ) == "no",
    "No",
    NA_character_
  )
)

sample_metadata$Data_Type <-
  "Affymetrix raw CEL"

sample_metadata$Is_Pooled <-
  FALSE

sample_metadata$Is_Paired <-
  sample_metadata$Group !=
    "Healthy_control"

sample_metadata$Is_Repeated_Measure <-
  sample_metadata$Group !=
    "Healthy_control"

sample_metadata$Metadata_Confidence <-
  "Direct_from_GEO_SOFT"

sample_metadata$Expected_CEL_GZ_Filename <-
  basename(
    sample_metadata$CEL_GZ_URL
  )

sample_metadata$Quality_Notes <- paste0(
  "Peripheral-blood leukocyte bulk expression; Early and Mid burn samples ",
  "are repeated measurements from the same patient; healthy controls have ",
  "one array each; the main model adjusts for age and sex and blocks ",
  "Patient_ID; TBSA, inhalation injury, survival, and hospital stay are ",
  "retained as metadata but are not covariates in the three main contrasts; ",
  "bulk-blood differences may reflect both intracellular regulation and ",
  "changes in leukocyte composition."
)

sample_metadata$Quality_Notes[
  sample_metadata$Group ==
    "Healthy_control"
] <- paste0(
  sample_metadata$Quality_Notes[
    sample_metadata$Group ==
      "Healthy_control"
  ],
  " Healthy-control arrays are not duplicated to create pseudo time courses ",
  "in the current BurnOmicsDB analysis."
)

# Verify repeated-measure structure.
# 检查重复测量结构。
burn_patient_counts <- table(
  sample_metadata$Patient_ID[
    sample_metadata$Group !=
      "Healthy_control"
  ]
)

control_patient_counts <- table(
  sample_metadata$Patient_ID[
    sample_metadata$Group ==
      "Healthy_control"
  ]
)

if (length(
  burn_patient_counts
) != 57 ||
    !all(
      burn_patient_counts == 2
    )) {
  stop(
    "The burn cohort should contain 57 patients with exactly two arrays each."
  )
}

if (length(
  control_patient_counts
) != 63 ||
    !all(
      control_patient_counts == 1
    )) {
  stop(
    "The healthy-control cohort should contain 63 subjects with one array each."
  )
}

patient_group_list <- split(
  as.character(
    sample_metadata$Group
  ),
  sample_metadata$Patient_ID
)

burn_pair_check <- vapply(
  patient_group_list[
    grepl(
      "^BurnPatient",
      names(
        patient_group_list
      )
    )
  ],
  function(group_values) {
    setequal(
      group_values,
      c(
        "Early_burn",
        "Mid_burn"
      )
    )
  },
  logical(1)
)

if (!all(
  burn_pair_check
)) {
  stop(
    "At least one burn patient does not have one Early and one Mid sample."
  )
}

actual_group_counts <- table(
  sample_metadata$Group
)

expected_group_counts <- c(
  Healthy_control = 63,
  Early_burn = 57,
  Mid_burn = 57
)

if (!all(
  actual_group_counts[
    names(
      expected_group_counts
    )
  ] == expected_group_counts
)) {
  stop(
    paste0(
      "Unexpected group counts: ",
      paste(
        names(
          actual_group_counts
        ),
        actual_group_counts,
        sep = "=",
        collapse = ", "
      )
    )
  )
}

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
      paste(
        missing_from_soft,
        collapse = ", "
      ),
      ". Missing from CEL: ",
      paste(
        missing_from_cel,
        collapse = ", "
      ),
      "."
    )
  )
}

sample_metadata <- sample_metadata[
  match(
    cel_sample_ids,
    sample_metadata$Sample_ID
  ),
]

if (!identical(
  sample_metadata$Sample_ID,
  cel_sample_ids
)) {
  stop(
    "Sample metadata could not be reordered to match the CEL files."
  )
}

# Reorder CEL paths to the validated metadata order.
# 按已核对的元数据顺序排列CEL文件。
cel_file_map <- setNames(
  cel_files,
  cel_sample_ids
)

cel_files <- unname(
  cel_file_map[
    sample_metadata$Sample_ID
  ]
)

cel_sample_ids <-
  sample_metadata$Sample_ID

if (any(
  is.na(
    cel_files
  )
)) {
  stop(
    "At least one CEL path is missing after sample-order alignment."
  )
}

early_hours <- sample_metadata$Hours_Post_Injury[
  sample_metadata$Group ==
    "Early_burn"
]

mid_hours <- sample_metadata$Hours_Post_Injury[
  sample_metadata$Group ==
    "Mid_burn"
]

metadata_audit <- data.frame(
  Check = c(
    "Raw CEL file count",
    "SOFT sample count",
    "Healthy-control array count",
    "Early-burn array count",
    "Mid-burn array count",
    "Unique burn-patient count",
    "Unique healthy-control count",
    "Total unique individuals",
    "Burn arrays per patient",
    "Healthy-control arrays per subject",
    "Early-stage observed hour range",
    "Mid-stage observed hour range",
    "Missing age values",
    "Missing sex values",
    "Pooled samples",
    "Patient-level repeated measures",
    "Primary model covariates"
  ),
  Value = c(
    length(
      cel_files
    ),
    nrow(
      sample_metadata
    ),
    actual_group_counts[
      "Healthy_control"
    ],
    actual_group_counts[
      "Early_burn"
    ],
    actual_group_counts[
      "Mid_burn"
    ],
    length(
      burn_patient_counts
    ),
    length(
      control_patient_counts
    ),
    length(
      unique(
        sample_metadata$Patient_ID
      )
    ),
    "2",
    "1",
    paste0(
      round(
        min(
          early_hours
        ),
        1
      ),
      "-",
      round(
        max(
          early_hours
        ),
        1
      )
    ),
    paste0(
      round(
        min(
          mid_hours
        ),
        1
      ),
      "-",
      round(
        max(
          mid_hours
        ),
        1
      )
    ),
    sum(
      is.na(
        sample_metadata$Age
      )
    ),
    sum(
      is.na(
        sample_metadata$Sex
      )
    ),
    sum(
      sample_metadata$Is_Pooled
    ),
    "Yes",
    "Continuous age and sex"
  ),
  Interpretation = c(
    "Directly observed in GSE19743_RAW.tar.",
    "Directly parsed from the streaming SOFT parser.",
    "One independent array per healthy subject.",
    "One Early array per burn patient.",
    "One Mid array per burn patient.",
    "Each burn patient contributes Early and Mid arrays.",
    "Each healthy control contributes one array.",
    "Fifty-seven burn patients plus sixty-three healthy controls.",
    "Repeated longitudinal measurements.",
    "Singleton block in the repeated-measures model.",
    "Observed values in the uploaded SOFT metadata.",
    "Observed values in the uploaded SOFT metadata.",
    "Age is complete for the primary model.",
    "Sex is complete for the primary model.",
    "No pooling is reported.",
    "duplicateCorrelation and Patient_ID blocking are required.",
    "TBSA and outcome variables are retained in metadata but not modeled."
  ),
  stringsAsFactors = FALSE
)

write.csv(
  metadata_audit,
  file = file.path(
    QC_DIR,
    "GSE19743_metadata_consistency_check.csv"
  ),
  row.names = FALSE
)

cat(
  "Sample-group counts:\n"
)

print(
  actual_group_counts
)

cat(
  "\nUnique individuals: ",
  length(
    unique(
      sample_metadata$Patient_ID
    )
  ),
  "\n\n",
  sep = ""
)

cat(
  "Metadata audit:\n"
)

print(
  metadata_audit
)


# ---- 04. Batch-wise raw CEL QC / 分批进行原始CEL质量控制 ----

# Reading all 177 raw probe-intensity matrices at once is unnecessary for the
# requested QC. CEL files are read in batches and only per-array summary
# metrics are retained.
# 无需一次载入177张芯片的全部原始强度。本节分批读取，
# 仅保留每张芯片的摘要指标。

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
  length(
    cel_files
  ),
  by = RAW_QC_BATCH_SIZE
)

for (
  batch_start in batch_starts
) {
  batch_end <- min(
    batch_start +
      RAW_QC_BATCH_SIZE -
      1,
    length(
      cel_files
    )
  )

  batch_index <- batch_start:
    batch_end

  cat(
    "Reading raw CEL QC batch ",
    batch_start,
    "-",
    batch_end,
    " of ",
    length(
      cel_files
    ),
    ".\n",
    sep = ""
  )

  raw_batch <- affy::ReadAffy(
    filenames = cel_files[
      batch_index
    ]
  )

  sampleNames(
    raw_batch
  ) <- sample_metadata$Sample_ID[
    batch_index
  ]

  raw_log2_batch <- log2(
    exprs(
      raw_batch
    ) + 1
  )

  raw_qc_metrics$Raw_Median_Log2_Intensity[
    batch_index
  ] <- apply(
    raw_log2_batch,
    2,
    median
  )

  raw_qc_metrics$Raw_IQR_Log2_Intensity[
    batch_index
  ] <- apply(
    raw_log2_batch,
    2,
    IQR
  )

  raw_qc_metrics$Raw_Min_Log2_Intensity[
    batch_index
  ] <- apply(
    raw_log2_batch,
    2,
    min
  )

  raw_qc_metrics$Raw_Max_Log2_Intensity[
    batch_index
  ] <- apply(
    raw_log2_batch,
    2,
    max
  )

  rm(
    raw_batch,
    raw_log2_batch
  )

  invisible(
    gc()
  )
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

raw_qc_metrics$Raw_Median_Robust_Z <- if (
  raw_median_mad > 0
) {
  (
    raw_qc_metrics$Raw_Median_Log2_Intensity -
      raw_median_center
  ) /
    raw_median_mad
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
    "GSE19743_raw_array_QC_metrics.csv"
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
    width = 0.15,
    height = 0,
    size = 1.5,
    alpha = 0.70
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
    title = "GSE19743 raw CEL intensity QC",
    subtitle = "Each point is the median raw log2 probe intensity of one array",
    x = "Sample group",
    y = "Median raw log2 probe intensity"
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      angle = 25,
      hjust = 1
    )
  )

print(
  p_raw_qc
)

ggsave(
  filename = file.path(
    FIGURES_DIR,
    "01_GSE19743_raw_data_QC.png"
  ),
  plot = p_raw_qc,
  width = 9,
  height = 5.5,
  dpi = 300
)


# ---- 05. Memory-conscious RMA normalization / 内存友好的RMA标准化 ----

cat(
  "Starting RMA normalization for 177 CEL files. ",
  "This step can take substantial time on an Intel Mac.\n"
)

# justRMA avoids retaining the complete raw AffyBatch after normalization.
# justRMA避免在RMA完成后继续保留完整原始AffyBatch。
rma_cel_filenames <- basename(
  cel_files
)

rma_cel_paths <- file.path(
  CEL_DIR,
  rma_cel_filenames
)

if (any(
  !file.exists(
    rma_cel_paths
  )
)) {
  missing_rma_files <- rma_cel_paths[
    !file.exists(
      rma_cel_paths
    )
  ]

  stop(
    paste0(
      "RMA input validation failed. Missing CEL files: ",
      paste(
        missing_rma_files,
        collapse = ", "
      )
    )
  )
}

# justRMA prepends celfile.path to filenames. Pass basenames here to prevent
# an absolute path from being prefixed a second time.
# justRMA会将celfile.path拼接到filenames之前，因此此处仅传入CEL文件名，
# 避免完整绝对路径被重复拼接。
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

if (ncol(
  rma_probe_expression
) != 177) {
  stop(
    "The RMA expression matrix does not contain 177 arrays."
  )
}

if (!identical(
  colnames(
    rma_probe_expression
  ),
  sample_metadata$Sample_ID
)) {
  stop(
    "RMA expression columns do not match the sample metadata."
  )
}

if (any(
  is.na(
    rma_probe_expression
  )
)) {
  stop(
    "The RMA expression matrix contains missing values."
  )
}

if (any(
  !is.finite(
    rma_probe_expression
  )
)) {
  stop(
    "The RMA expression matrix contains non-finite values."
  )
}

cat(
  "RMA normalization completed.\n"
)

# Save probe-level RMA expression for reproducibility.
# 保存探针集层面的RMA表达用于复现。
rma_probe_output <- data.frame(
  Probe_ID = rownames(
    rma_probe_expression
  ),
  rma_probe_expression,
  check.names = FALSE
)

rma_probe_connection <- gzfile(
  file.path(
    RESULTS_DIR,
    "GSE19743_RMA_probe_level_expression.csv.gz"
  ),
  open = "wt"
)

write.csv(
  rma_probe_output,
  file = rma_probe_connection,
  row.names = FALSE
)

close(
  rma_probe_connection
)

rm(
  rma_probe_output
)

invisible(
  gc()
)

rma_qc_metrics <- data.frame(
  Sample_ID = sample_metadata$Sample_ID,
  Patient_ID = sample_metadata$Patient_ID,
  Group = sample_metadata$Group,
  RMA_Median = apply(
    rma_probe_expression,
    2,
    median
  ),
  RMA_IQR = apply(
    rma_probe_expression,
    2,
    IQR
  ),
  stringsAsFactors = FALSE
)

rma_qc_metrics$RLE_Median <- apply(
  sweep(
    rma_probe_expression,
    1,
    apply(
      rma_probe_expression,
      1,
      median
    ),
    "-"
  ),
  2,
  median
)

rma_qc_metrics$RLE_IQR <- apply(
  sweep(
    rma_probe_expression,
    1,
    apply(
      rma_probe_expression,
      1,
      median
    ),
    "-"
  ),
  2,
  IQR
)

write.csv(
  rma_qc_metrics,
  file = file.path(
    QC_DIR,
    "GSE19743_RMA_array_QC_metrics.csv"
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
    width = 0.15,
    height = 0,
    size = 1.5,
    alpha = 0.70
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
    title = "GSE19743 RMA-normalized expression distribution",
    subtitle = "Each point is the median RMA log2 expression of one array",
    x = "Sample group",
    y = "Median RMA log2 expression"
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      angle = 25,
      hjust = 1
    )
  )

print(
  p_normalized
)

ggsave(
  filename = file.path(
    FIGURES_DIR,
    "02_GSE19743_normalized_expression_distribution.png"
  ),
  plot = p_normalized,
  width = 9,
  height = 5.5,
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
          !is.na(
            entrez_values
          ) &
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
    probe_entrez_count[
      probe_ids
    ]
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

unambiguous_probe_ids <- probe_mapping_status$Probe_ID[
  probe_mapping_status$Mapping_Status ==
    "Unambiguous"
]

unambiguous_annotation_rows <- probe_annotation_raw[
  probe_annotation_raw$PROBEID %in%
    unambiguous_probe_ids &
    !is.na(
      probe_annotation_raw$ENTREZID
    ) &
    probe_annotation_raw$ENTREZID != "",
]

unambiguous_probe_annotation <- unique(
  unambiguous_annotation_rows[
    ,
    c(
      "PROBEID",
      "ENTREZID"
    )
  ]
)

if (anyDuplicated(
  unambiguous_probe_annotation$PROBEID
) > 0) {
  stop(
    "A probe classified as unambiguous still maps to multiple Entrez Gene IDs."
  )
}

probe_mean_expression <- rowMeans(
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

colnames(
  representative_probe_table
)[
  colnames(
    representative_probe_table
  ) == "PROBEID"
] <- "Representative_Probe_ID"

colnames(
  representative_probe_table
)[
  colnames(
    representative_probe_table
  ) == "ENTREZID"
] <- "NCBI_Gene_ID"

gene_expression <- rma_probe_expression[
  representative_probe_table$Representative_Probe_ID,
  ,
  drop = FALSE
]

rownames(
  gene_expression
) <- representative_probe_table$NCBI_Gene_ID

if (anyDuplicated(
  rownames(
    gene_expression
  )
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
    gene_symbol[
      selected_entrez_ids
    ]
  ),
  Ensembl_ID = unname(
    ensembl_id[
      selected_entrez_ids
    ]
  ),
  Gene_Name = unname(
    gene_name[
      selected_entrez_ids
    ]
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
  is.na(
    gene_annotation$Gene_Symbol
  ),
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
  selected_match[
    !is.na(
      selected_match
    )
  ]
] <- TRUE

probe_mapping_connection <- gzfile(
  file.path(
    QC_DIR,
    "GSE19743_probe_mapping_and_selection.csv.gz"
  ),
  open = "wt"
)

write.csv(
  probe_mapping_audit,
  file = probe_mapping_connection,
  row.names = FALSE
)

close(
  probe_mapping_connection
)

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
    length(
      probe_ids
    ),
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
    nrow(
      representative_probe_table
    )
  ),
  stringsAsFactors = FALSE
)

write.csv(
  mapping_summary,
  file = file.path(
    QC_DIR,
    "GSE19743_probe_mapping_summary.csv"
  ),
  row.names = FALSE
)

cat(
  "Probe-mapping summary:\n"
)

print(
  mapping_summary
)

normalized_gene_output <- data.frame(
  NCBI_Gene_ID = rownames(
    gene_expression
  ),
  gene_expression,
  check.names = FALSE
)

normalized_gene_connection <- gzfile(
  file.path(
    RESULTS_DIR,
    "GSE19743_normalized_expression.csv.gz"
  ),
  open = "wt"
)

write.csv(
  normalized_gene_output,
  file = normalized_gene_connection,
  row.names = FALSE
)

close(
  normalized_gene_connection
)

rm(
  normalized_gene_output
)

invisible(
  gc()
)


# ---- 07. PCA and sample correlation / PCA与样本相关性 ----

gene_variance <- apply(
  gene_expression,
  1,
  var
)

n_pca_genes <- min(
  TOP_VARIABLE_GENES_FOR_PCA,
  length(
    gene_variance
  )
)

top_variable_gene_ids <- names(
  sort(
    gene_variance,
    decreasing = TRUE
  )
)[
  seq_len(
    n_pca_genes
  )
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
    sum(
      pca_result$sdev^2
    )
)

pca_table <- data.frame(
  Sample_ID = rownames(
    pca_result$x
  ),
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

pca_table$Age_Group <- sample_metadata$Age_Group[
  match(
    pca_table$Sample_ID,
    sample_metadata$Sample_ID
  )
]

write.csv(
  pca_table,
  file = file.path(
    RESULTS_DIR,
    "GSE19743_PCA_scores.csv"
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
    size = 2.2,
    alpha = 0.72
  ) +
  scale_color_manual(
    values = GROUP_COLORS,
    drop = FALSE
  ) +
  labs(
    title = "GSE19743 PCA",
    subtitle = paste0(
      "Top ",
      n_pca_genes,
      " variable gene-level RMA features; 177 arrays"
    ),
    x = paste0(
      "PC1 (",
      round(
        pca_variance[1],
        1
      ),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(
        pca_variance[2],
        1
      ),
      "%)"
    ),
    color = "Group"
  ) +
  theme_classic(
    base_size = 12
  )

print(
  p_pca
)

ggsave(
  filename = file.path(
    FIGURES_DIR,
    "03_GSE19743_PCA.png"
  ),
  plot = p_pca,
  width = 8,
  height = 6,
  dpi = 300
)

sample_correlation <- cor(
  gene_expression,
  method = "pearson"
)

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
  main = "GSE19743 sample correlation"
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
  main = "GSE19743 sample correlation",
  filename = file.path(
    FIGURES_DIR,
    "04_GSE19743_sample_correlation_heatmap.png"
  ),
  width = 11,
  height = 10
)


# ---- 08. Repeated-measures differential-expression model / 重复测量差异表达模型 ----

sample_metadata$Group <- factor(
  sample_metadata$Group,
  levels = GROUP_LEVELS
)

sample_metadata$Sex <- factor(
  sample_metadata$Sex,
  levels = c(
    "Female",
    "Male"
  )
)

sample_metadata$Age_Centered <-
  sample_metadata$Age -
    mean(
      sample_metadata$Age
    )

design <- model.matrix(
  ~ 0 + Group + Age_Centered + Sex,
  data = sample_metadata
)

colnames(
  design
) <- sub(
  "^Group",
  "",
  colnames(
    design
  )
)

rownames(
  design
) <- sample_metadata$Sample_ID

if (qr(
  design
)$rank < ncol(
  design
)) {
  stop(
    "The differential-expression design matrix is not full rank."
  )
}

patient_block <- factor(
  sample_metadata$Patient_ID
)

cat(
  "Design-matrix columns:\n"
)

print(
  colnames(
    design
  )
)

correlation_fit <- limma::duplicateCorrelation(
  gene_expression,
  design = design,
  block = patient_block
)

within_patient_correlation <-
  correlation_fit$consensus.correlation

if (!is.finite(
  within_patient_correlation
)) {
  stop(
    "Within-patient correlation could not be estimated."
  )
}

cat(
  "Estimated within-patient correlation: ",
  round(
    within_patient_correlation,
    6
  ),
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
  GSE19743_EarlyBurn_vs_HealthyControl =
    Early_burn - Healthy_control,
  GSE19743_MidBurn_vs_HealthyControl =
    Mid_burn - Healthy_control,
  GSE19743_MidBurn_vs_EarlyBurn =
    Mid_burn - Early_burn,
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
  Contrast_ID = colnames(
    contrast_matrix
  ),
  Contrast_Label = c(
    "Early burn vs healthy control",
    "Middle burn vs healthy control",
    "Middle burn vs early burn"
  ),
  Case_Group = c(
    "Early burn",
    "Middle burn",
    "Middle burn"
  ),
  Case_Group_Code = c(
    "Early_burn",
    "Mid_burn",
    "Mid_burn"
  ),
  Control_Group = c(
    "Healthy control",
    "Healthy control",
    "Early burn"
  ),
  Control_Group_Code = c(
    "Healthy_control",
    "Healthy_control",
    "Early_burn"
  ),
  Sample_Context = c(
    "Peripheral blood leukocytes after severe burn injury",
    "Peripheral blood leukocytes after severe burn injury",
    "Longitudinal peripheral blood leukocyte response after severe burn injury"
  ),
  Time_or_Stage = c(
    "Early stage (<11 days after injury)",
    "Middle stage (11-49 days after injury; observed 11.1-31.7 days)",
    "Longitudinal change from early to middle stage"
  ),
  Is_Paired_Contrast = c(
    FALSE,
    FALSE,
    TRUE
  ),
  stringsAsFactors = FALSE
)

contrast_definitions$Case_N <- vapply(
  contrast_definitions$Case_Group_Code,
  function(group_code) {
    sum(
      sample_metadata$Group ==
        group_code
    )
  },
  numeric(1)
)

contrast_definitions$Control_N <- vapply(
  contrast_definitions$Control_Group_Code,
  function(group_code) {
    sum(
      sample_metadata$Group ==
        group_code
    )
  },
  numeric(1)
)

contrast_definitions$Case_Patient_N <- vapply(
  contrast_definitions$Case_Group_Code,
  function(group_code) {
    length(
      unique(
        sample_metadata$Patient_ID[
          sample_metadata$Group ==
            group_code
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
          sample_metadata$Group ==
            group_code
        ]
      )
    )
  },
  numeric(1)
)

cat(
  "Contrast definitions:\n"
)

print(
  contrast_definitions
)


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
    rownames(
      contrast_result
    )

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
        output$log2FC <=
          -LOG2FC_CUTOFF,
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

row.names(
  all_gene_results
) <- NULL

all_results_connection <- gzfile(
  file.path(
    RESULTS_DIR,
    "GSE19743_all_gene_results.csv.gz"
  ),
  open = "wt"
)

write.csv(
  all_gene_results,
  file = all_results_connection,
  row.names = FALSE
)

close(
  all_results_connection
)

cat(
  "Differential-expression result counts:\n"
)

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
      !is.na(
        volcano_data$Gene_Symbol
      ),
  ]

  label_candidates <- label_candidates[
    order(
      label_candidates$FDR,
      -abs(
        label_candidates$log2FC
      )
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
      alpha = 0.63,
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
      yintercept = -log10(
        FDR_CUTOFF
      ),
      linetype = "dashed"
    ) +
    ggrepel::geom_text_repel(
      data = label_candidates,
      aes(
        label = Gene_Symbol
      ),
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
        "GSE19743: ",
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
    theme_classic(
      base_size = 12
    )

  print(
    p_volcano
  )

  ggsave(
    filename = file.path(
      FIGURES_DIR,
      paste0(
        "05_GSE19743_volcano_",
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

  ranked_gene_ids <- contrast_results$NCBI_Gene_ID[
    order(
      contrast_results$FDR,
      -abs(
        contrast_results$log2FC
      )
    )
  ]

  n_heatmap_genes <- min(
    TOP_GENES_FOR_HEATMAP,
    length(
      ranked_gene_ids
    )
  )

  top_heatmap_ids <- ranked_gene_ids[
    seq_len(
      n_heatmap_genes
    )
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
      t(
        heatmap_matrix
      )
    )
  )

  gene_label_match <- match(
    top_heatmap_ids,
    contrast_results$NCBI_Gene_ID
  )

  heatmap_labels <- contrast_results$Gene_Symbol[
    gene_label_match
  ]

  missing_heatmap_labels <- is.na(
    heatmap_labels
  ) |
    heatmap_labels == ""

  heatmap_labels[
    missing_heatmap_labels
  ] <- top_heatmap_ids[
    missing_heatmap_labels
  ]

  rownames(
    heatmap_z
  ) <- make.unique(
    heatmap_labels
  )

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

  show_sample_names <- length(
    contrast_sample_ids
  ) <= 30

  heatmap_title <- paste0(
    "Top ",
    n_heatmap_genes,
    " genes: ",
    contrast_info$Contrast_Label
  )

  pheatmap::pheatmap(
    heatmap_z,
    annotation_col = heatmap_annotation,
    annotation_colors = heatmap_annotation_colors,
    show_colnames = show_sample_names,
    show_rownames = TRUE,
    border_color = NA,
    color = EXPRESSION_HEATMAP_COLORS,
    main = heatmap_title
  )

  pheatmap::pheatmap(
    heatmap_z,
    annotation_col = heatmap_annotation,
    annotation_colors = heatmap_annotation_colors,
    show_colnames = show_sample_names,
    show_rownames = TRUE,
    border_color = NA,
    color = EXPRESSION_HEATMAP_COLORS,
    main = heatmap_title,
    filename = file.path(
      FIGURES_DIR,
      paste0(
        "06_GSE19743_top_differential_genes_heatmap_",
        contrast_id,
        ".png"
      )
    ),
    width = 12,
    height = 10
  )
}


# ---- 12. Qualitative check of genes emphasized in the associated paper / 定性核对论文重点基因 ----

# The paper reports age-dependent TANOVA findings rather than simple overall
# fold changes. This table only surfaces the corresponding BurnOmicsDB results;
# it is not a direct reproduction of the published TANOVA analysis.
# 论文报告年龄依赖的TANOVA结果，而非简单总体fold change。
# 本表仅列出相关基因在数据库三个contrast中的结果，不构成直接复现。

paper_gene_context <- data.frame(
  Gene_Symbol = c(
    "SLC2A3",
    "IGLJ3",
    "NDUFA3",
    "NDUFA7",
    "COX11",
    "NDUFA11",
    "NDUFA13",
    "UCP2",
    "OGDH",
    "NDUFB7",
    "NDUFS7",
    "HTRA2",
    "MAPK8",
    "BACE1",
    "IGKC",
    "IGHD"
  ),
  Published_Context = c(
    "Example of a burn-responsive gene across both stages.",
    "Example of an age-dependent interaction concentrated at the middle stage.",
    "Mitochondrial gene highlighted in age-dependent burn response.",
    "Mitochondrial gene highlighted in age-dependent burn response.",
    "Mitochondrial gene highlighted in age-dependent burn response.",
    "Mitochondrial gene highlighted in age-dependent burn response.",
    "Mitochondrial gene highlighted in age-dependent burn response.",
    "Mitochondrial gene highlighted in age-dependent burn response.",
    "Mitochondrial gene highlighted in age-dependent burn response.",
    "Mitochondrial gene highlighted in age-dependent burn response.",
    "Mitochondrial gene highlighted in age-dependent burn response.",
    "Mitochondrial gene highlighted in age-dependent burn response.",
    "Mitochondrial gene highlighted in age-dependent burn response.",
    "Mitochondrial gene highlighted in age-dependent burn response.",
    "Immunoglobulin gene highlighted in age-specific middle-stage response.",
    "Immunoglobulin gene highlighted in age-specific middle-stage response."
  ),
  Direct_Replication = FALSE,
  stringsAsFactors = FALSE
)

calculated_paper_genes <- all_gene_results[
  all_gene_results$Gene_Symbol %in%
    paper_gene_context$Gene_Symbol,
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

paper_feature_check <- merge(
  paper_gene_context,
  calculated_paper_genes,
  by = "Gene_Symbol",
  all.x = TRUE,
  sort = FALSE
)

write.csv(
  paper_feature_check,
  file = file.path(
    RESULTS_DIR,
    "GSE19743_paper_feature_gene_check.csv"
  ),
  row.names = FALSE
)

cat(
  "Paper-feature gene check created. ",
  "This is a qualitative context table, not a TANOVA replication.\n"
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
    "Pediatric and adult severe-burn patients with healthy controls",
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
    "Raw CEL files from GSE19743_RAW.tar",
  Normalization = paste0(
    "RMA: background correction, quantile normalization, ",
    "and median-polish probe-set summarization"
  ),
  Analysis_Method = paste0(
    "limma linear model with robust trend-aware empirical Bayes ",
    "moderation; duplicateCorrelation and Patient_ID blocking for Early/Mid ",
    "repeated measurements; adjusted for continuous age and sex; one ",
    "representative probe per NCBI Gene ID selected by highest mean RMA ",
    "expression across all 177 arrays; positive log2FC represents Case_Group ",
    "minus Control_Group"
  ),
  Annotation_Method = paste0(
    "hgu133plus2.db ",
    as.character(
      packageVersion(
        "hgu133plus2.db"
      )
    ),
    " used for probe-to-Entrez mapping; probes mapping to zero or multiple ",
    "Entrez Gene IDs excluded; org.Hs.eg.db ",
    as.character(
      packageVersion(
        "org.Hs.eg.db"
      )
    ),
    " used for current official symbol, Ensembl ID, and gene name"
  ),
  Quality_Notes = paste0(
    "Core human longitudinal blood dataset; 57 burn patients have paired ",
    "Early and Mid arrays and 63 healthy controls have one array each; ",
    "the model adjusts for age and sex and accounts for patient-level ",
    "correlation; the associated paper focused on age-dependent TANOVA using ",
    "balanced subsets, whereas the current analysis estimates three overall ",
    "BurnOmicsDB contrasts using all samples; healthy controls are not ",
    "duplicated; the Mid category is defined as 11-49 days, but uploaded SOFT ",
    "samples are observed through approximately 31.7 days; TBSA, inhalation ",
    "injury, survival, and hospital stay are not covariates in the main model; ",
    "bulk-blood differences may reflect both intracellular regulation and ",
    "changes in leukocyte composition"
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
  Is_Paired_Contrast =
    all_gene_results$Is_Paired_Contrast,
  Within_Patient_Correlation =
    within_patient_correlation,
  Age_Adjusted = TRUE,
  Sex_Adjusted = TRUE,
  Is_Pooled = FALSE,
  Quality_Grade =
    "Core_high_quality",
  stringsAsFactors = FALSE
)

database_connection <- gzfile(
  file.path(
    RESULTS_DIR,
    "GSE19743_database_ready_all_genes.csv.gz"
  ),
  open = "wt"
)

write.csv(
  database_ready,
  file = database_connection,
  row.names = FALSE
)

close(
  database_connection
)

cat(
  "Database-ready rows: ",
  nrow(
    database_ready
  ),
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
    "Age_Group",
    "Sex",
    "Hours_Post_Injury",
    "Days_Post_Injury",
    "TBSA_Percent",
    "Inhalation_Injury",
    "Hospital_Length_of_Stay_Days",
    "Survival_Original",
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
    "GSE19743_sample_metadata.csv"
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

analysis_objects <- list(
  project = "BurnOmicsDB",
  GEO_ID = GEO_ID,
  raw_tar_file = raw_tar_file,
  soft_file = soft_file,
  sample_metadata = sample_metadata_output,
  metadata_audit = metadata_audit,
  raw_qc_metrics = raw_qc_metrics,
  rma_eset = rma_eset,
  rma_qc_metrics = rma_qc_metrics,
  probe_mapping_summary = mapping_summary,
  representative_probe_table = representative_probe_table,
  gene_annotation = gene_annotation,
  gene_expression = gene_expression,
  design_matrix = design,
  within_patient_correlation =
    within_patient_correlation,
  contrast_matrix = contrast_matrix,
  contrast_definitions = contrast_definitions,
  fitted_model = fit_contrasts,
  all_gene_results = all_gene_results,
  paper_feature_check = paper_feature_check,
  database_ready = database_ready,
  thresholds = list(
    FDR_CUTOFF = FDR_CUTOFF,
    LOG2FC_CUTOFF = LOG2FC_CUTOFF
  )
)

saveRDS(
  analysis_objects,
  file = file.path(
    OBJECTS_DIR,
    "GSE19743_analysis_objects.rds"
  ),
  compress = "xz"
)

summary_lines <- c(
  "BurnOmicsDB - GSE19743 analysis summary",
  "",
  paste0(
    "Analysis date: ",
    Sys.Date()
  ),
  "Data type: Affymetrix Human Genome U133 Plus 2.0 raw CEL",
  paste0(
    "Input CEL arrays: ",
    length(
      cel_files
    )
  ),
  paste0(
    "Unique individuals: ",
    length(
      unique(
        sample_metadata$Patient_ID
      )
    )
  ),
  "Burn patients: 57",
  "Healthy controls: 63",
  "Pooling: No",
  "",
  "Groups:",
  paste0(
    "  Healthy control: ",
    actual_group_counts[
      "Healthy_control"
    ],
    " arrays from 63 independent controls"
  ),
  paste0(
    "  Early burn: ",
    actual_group_counts[
      "Early_burn"
    ],
    " arrays from 57 burn patients"
  ),
  paste0(
    "  Middle burn: ",
    actual_group_counts[
      "Mid_burn"
    ],
    " arrays from the same 57 burn patients"
  ),
  "",
  paste0(
    "Observed Early hours post injury: ",
    round(
      min(
        early_hours
      ),
      1
    ),
    "-",
    round(
      max(
        early_hours
      ),
      1
    )
  ),
  paste0(
    "Observed Middle hours post injury: ",
    round(
      min(
        mid_hours
      ),
      1
    ),
    "-",
    round(
      max(
        mid_hours
      ),
      1
    )
  ),
  paste0(
    "Input Affymetrix probe sets: ",
    nrow(
      rma_probe_expression
    )
  ),
  paste0(
    "Gene-level representatives tested: ",
    nrow(
      gene_expression
    )
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
    "177 arrays for each NCBI Gene ID"
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
  "Contrasts:",
  paste0(
    "  ",
    contrast_definitions$Contrast_ID,
    ": ",
    contrast_definitions$Contrast_Label,
    " (case arrays = ",
    contrast_definitions$Case_N,
    ", case patients = ",
    contrast_definitions$Case_Patient_N,
    ", control arrays = ",
    contrast_definitions$Control_N,
    ", control patients = ",
    contrast_definitions$Control_Patient_N,
    ", paired contrast = ",
    contrast_definitions$Is_Paired_Contrast,
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
    "Positive log2FC means higher expression in the case group than in the ",
    "control/reference group."
  ),
  paste0(
    "No separate low-expression filter was applied after RMA and ",
    "representative-probe selection; all selected unambiguous gene ",
    "representatives were statistically tested and exported."
  ),
  "",
  "Major limitations:",
  paste0(
    "  The associated paper focused on age-dependent TANOVA using balanced ",
    "subsets; the current analysis estimates simpler overall contrasts using ",
    "all 177 arrays."
  ),
  paste0(
    "  Continuous age and sex are adjusted, but age-by-burn interactions are ",
    "not estimated."
  ),
  paste0(
    "  Healthy controls have one array each and are not duplicated into ",
    "pseudo time courses."
  ),
  paste0(
    "  The Middle category is defined as 11-49 days, but the uploaded SOFT ",
    "samples are observed through approximately ",
    round(
      max(
        mid_hours
      ) /
        24,
      1
    ),
    " days."
  ),
  paste0(
    "  TBSA, inhalation injury, survival, and hospital length of stay are ",
    "not covariates in the primary model."
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
    "GSE19743_analysis_summary.txt"
  ),
  useBytes = TRUE
)

capture.output(
  sessionInfo(),
  file = file.path(
    RESULTS_DIR,
    "GSE19743_R_sessionInfo.txt"
  )
)

cat(
  "Analysis completed.\n"
)

cat(
  "Results directory:\n",
  RESULTS_DIR,
  "\n"
)

cat(
  "Figures directory:\n",
  FIGURES_DIR,
  "\n"
)

cat(
  "QC directory:\n",
  QC_DIR,
  "\n"
)

cat(
  "R objects directory:\n",
  OBJECTS_DIR,
  "\n"
)


# ---- 15. Optional individual-gene checking / 可选的单基因检查 ----

# Run this section after Section 09.
# 在完成第09部分后运行本节。

genes_to_check <- c(
  "SLC2A3",
  "IL6",
  "CXCL8",
  "S100A8",
  "S100A9",
  "NDUFA3",
  "NDUFA7",
  "NDUFA11",
  "NDUFA13",
  "UCP2",
  "IGKC",
  "IGHD"
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

print(
  gene_check_table
)

# End of script / 代码结束
