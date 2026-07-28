# ==============================================================================
# BurnOmicsDB: GSE178411 RNA-seq analysis
# BurnOmicsDB：GSE178411 RNA-seq分析
#
# Script revision / 代码修订:
#   Version 1.1 fixes GEO age values recorded as "unknown".
#   Version 1.1处理GEO中年龄标记为unknown的样本。
#
# GEO accession / GEO编号:
#   GSE178411
#
# Project / 项目:
#   Whole-transcriptome analysis of human uninjured skin, acute burn wounds,
#   and hypertrophic scars.
#   人体未损伤皮肤、急性烧伤创面和增生性瘢痕的全转录组分析。
#
# Input files / 输入文件:
#   1. GSE178411_counts.txt.gz
#      Author-provided NCBI Gene ID raw-count matrix.
#      作者提供的NCBI Gene ID原始计数矩阵。
#   2. GSE178411_family.soft.gz
#      GEO Series and sample-level metadata.
#      GEO研究及样本级元数据。
#
# No linked paper / 未关联正式论文:
#   GEO currently reports "Citation missing". This script therefore does not
#   attempt to reproduce a peer-reviewed paper's results.
#   GEO目前显示Citation missing，因此本代码不尝试复现正式论文结果。
#
# Differential-expression contrasts / 差异表达比较:
#   1. Early wound versus uninjured skin
#      早期创面 vs 未损伤皮肤
#   2. Late wound versus uninjured skin
#      晚期创面 vs 未损伤皮肤
#   3. Hypertrophic scar versus uninjured skin
#      增生性瘢痕 vs 未损伤皮肤
#
# Positive log2FC always means higher expression in the case group than in
# uninjured skin.
# 正log2FC始终表示case组相对于未损伤皮肤表达更高。
#
# Important design features / 重要实验设计:
#   - The uploaded count matrix contains 28,395 NCBI Gene IDs and 108 samples.
#     计数矩阵包含28,395个NCBI Gene ID和108个样本。
#   - SOFT metadata identifies 75 patients. Some patients contributed multiple
#     tissues or wound samples.
#     SOFT元数据包含75名患者，部分患者提供多个组织或创面样本。
#   - The four eligible groups contain 103 samples, but two samples have age
#     recorded as "unknown" in GEO:
#       GSM5390619 / E728-W / Patient 97 / Late wound
#       GSM5390626 / E742-S / Patient 104 / Hypertrophic scar
#     Because the model adjusts for age, these two samples are retained in
#     metadata and raw-data QC but excluded from differential-expression
#     modeling. The final model therefore uses 101 samples:
#       24 uninjured skin
#       22 early wound
#       28 late wound
#       27 hypertrophic scar
#     四个目标组原有103个样本，但两个样本的年龄在GEO中记录为unknown。
#     由于模型需要校正年龄，这两个样本保留在元数据和原始QC中，
#     但不进入差异表达模型。最终模型使用101个样本。
#   - Three chronic-wound samples and two normal-scar samples are retained in
#     metadata but excluded because their groups are not part of the three
#     prespecified contrasts.
#     3个慢性创面和2个正常瘢痕样本保留在元数据中，
#     但其分组不属于三个预设比较，因此不进入差异分析。
#   - TMM normalization and low-expression filtering are performed with edgeR.
#     TMM标准化和低表达过滤使用edgeR完成。
#   - Because some patients contributed repeated samples, differential
#     expression is modeled with limma-voom and duplicateCorrelation using
#     Patient_ID as the blocking factor. Age and sex are included as covariates.
#     由于部分患者存在重复样本，差异分析使用limma-voom和
#     duplicateCorrelation，并以Patient_ID作为阻断因素，同时校正年龄和性别。
#
# Why not use an ordinary unpaired edgeR GLM? / 为什么不使用普通非配对edgeR模型?
#   An ordinary model would treat repeated tissues from the same patient as
#   independent observations. The voom + duplicateCorrelation workflow retains
#   raw-count-aware mean-variance modeling while accounting for within-patient
#   correlation.
#   普通模型会错误地把同一患者的多个组织当作独立样本。voom结合
#   duplicateCorrelation可以在处理RNA-seq均值-方差关系的同时校正患者内相关性。
#
# How to run in RStudio / 如何在RStudio中运行:
#   - Save this script in:
#     /Users/peter/Downloads/Project-2026-BurnOmicsDB/GSE178411/
#   - Press Cmd + Shift + O on macOS to open the section outline.
#     在macOS中按Cmd + Shift + O打开代码分区目录。
#   - Run one section at a time with Cmd + Enter.
#     使用Cmd + Enter逐段运行。
#   - During the first run, do not Source the entire script at once.
#     第一次运行时不要直接Source全文。
#
# Output-language rule / 输出语言规则:
#   Comments are bilingual, but figures, CSV files, TXT files, console messages,
#   and error messages are English only.
#   注释使用中英文，但图片、CSV、TXT、Console信息和报错全部只使用英文。
# ==============================================================================


# ---- 00. Install required packages once / 首次安装所需软件包 ----

# Run this section only when packages are not installed.
# 仅在相应软件包尚未安装时运行本节。

# CRAN packages / CRAN软件包:
# install.packages(c(
#   "ggplot2",
#   "ggrepel",
#   "pheatmap"
# ))

# Bioconductor packages / Bioconductor软件包:
# if (!requireNamespace("BiocManager", quietly = TRUE)) {
#   install.packages("BiocManager")
# }
#
# BiocManager::install(c(
#   "edgeR",
#   "limma",
#   "AnnotationDbi",
#   "org.Hs.eg.db"
# ))


# ---- 01. Project settings and package loading / 项目设置与软件包加载 ----

GEO_ID <- "GSE178411"

PROJECT_DIR <- "/Users/peter/Downloads/Project-2026-BurnOmicsDB/GSE178411"

FDR_CUTOFF <- 0.05
LOG2FC_CUTOFF <- 1
TOP_VARIABLE_GENES_FOR_PCA <- 500
TOP_GENES_FOR_HEATMAP <- 30
VOLCANO_LABEL_GENE_N <- 12

RESULTS_DIR <- file.path(PROJECT_DIR, "Results")
FIGURES_DIR <- file.path(PROJECT_DIR, "Figures")
QC_DIR <- file.path(PROJECT_DIR, "QC")
OBJECTS_DIR <- file.path(PROJECT_DIR, "R_objects")

dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(QC_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OBJECTS_DIR, recursive = TRUE, showWarnings = FALSE)

required_packages <- c(
  "ggplot2",
  "ggrepel",
  "pheatmap",
  "edgeR",
  "limma",
  "AnnotationDbi",
  "org.Hs.eg.db"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
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
  library(edgeR)
  library(limma)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

set.seed(2026)

# Fixed semantic colors used across BurnOmicsDB.
# BurnOmicsDB所有项目使用的固定语义配色。
GROUP_COLORS <- c(
  "Uninjured_skin" = "#0072B2",
  "Early_wound" = "#D55E00",
  "Late_wound" = "#009E73",
  "Hypertrophic_scar" = "#CC79A7",
  "Chronic_wound" = "#E69F00",
  "Normal_scar" = "#56B4E9"
)

VOLCANO_COLORS <- c(
  "Up_significant" = "#D55E00",
  "Down_significant" = "#0072B2",
  "Not_significant" = "#BDBDBD"
)

EXPRESSION_HEATMAP_COLORS <- grDevices::colorRampPalette(
  c("#0072B2", "#F7F7F7", "#D55E00")
)(101)

ANALYSIS_GROUP_LEVELS <- c(
  "Uninjured_skin",
  "Early_wound",
  "Late_wound",
  "Hypertrophic_scar"
)

ALL_GROUP_LEVELS <- c(
  ANALYSIS_GROUP_LEVELS,
  "Chronic_wound",
  "Normal_scar"
)

cat("Project directory:\n", PROJECT_DIR, "\n\n")


# ---- 02. Locate input files / 定位输入文件 ----

# Search both the project root and an optional Input folder.
# 同时搜索项目根目录和可选的Input子文件夹。
input_search_dirs <- c(
  file.path(PROJECT_DIR, "Input"),
  PROJECT_DIR
)

find_one_file <- function(pattern, search_dirs, required = TRUE) {
  hits <- unlist(
    lapply(search_dirs, function(directory) {
      if (!dir.exists(directory)) {
        return(character(0))
      }

      list.files(
        path = directory,
        pattern = pattern,
        full.names = TRUE,
        ignore.case = TRUE
      )
    }),
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

  normalizePath(hits[1], mustWork = TRUE)
}

counts_file <- find_one_file(
  "^GSE178411_counts\\.txt.*\\.gz$",
  input_search_dirs
)

soft_file <- find_one_file(
  "^GSE178411_family\\.soft.*\\.gz$",
  input_search_dirs
)

cat("Count matrix:\n", counts_file, "\n\n")
cat("SOFT metadata:\n", soft_file, "\n\n")


# ---- 03. Import and validate the raw count matrix / 导入并检查原始计数矩阵 ----

# The first data column contains NCBI Gene IDs, but its header is blank.
# 数据第一列为NCBI Gene ID，但原文件没有给该列设置列名。
counts_table <- read.delim(
  gzfile(counts_file),
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

count_matrix <- as.matrix(counts_table)

if (nrow(count_matrix) == 0 || ncol(count_matrix) == 0) {
  stop("The count matrix is empty.")
}

if (any(is.na(count_matrix))) {
  stop("The count matrix contains missing values.")
}

if (any(count_matrix < 0)) {
  stop("The count matrix contains negative values.")
}

if (any(abs(count_matrix - round(count_matrix)) > 1e-8)) {
  stop("The count matrix contains non-integer values.")
}

if (any(is.na(rownames(count_matrix))) || any(rownames(count_matrix) == "")) {
  stop("The count matrix contains missing Gene IDs.")
}

if (anyDuplicated(rownames(count_matrix)) > 0) {
  stop("Duplicated Gene IDs were detected in the count matrix.")
}

if (anyDuplicated(colnames(count_matrix)) > 0) {
  stop("Duplicated sample names were detected in the count matrix.")
}

storage.mode(count_matrix) <- "integer"

input_qc <- list(
  input_gene_n = nrow(count_matrix),
  input_sample_n = ncol(count_matrix),
  unique_gene_id_n = length(unique(rownames(count_matrix))),
  duplicated_gene_id_n = anyDuplicated(rownames(count_matrix)),
  all_zero_gene_n = sum(rowSums(count_matrix) == 0),
  minimum_library_size = min(colSums(count_matrix)),
  maximum_library_size = max(colSums(count_matrix)),
  median_library_size = median(colSums(count_matrix)),
  minimum_zero_fraction = min(colMeans(count_matrix == 0)),
  maximum_zero_fraction = max(colMeans(count_matrix == 0))
)

cat("Input count-matrix dimensions:\n")
print(dim(count_matrix))

cat("\nInput QC summary:\n")
print(input_qc)

capture.output(
  input_qc,
  file = file.path(
    QC_DIR,
    "GSE178411_input_QC_summary.txt"
  )
)


# ---- 04. Parse GEO SOFT metadata and audit sample structure / 解析SOFT元数据并审计样本结构 ----

# Parse all ^SAMPLE blocks from the local SOFT file.
# 从本地SOFT文件解析所有^SAMPLE区块。
parse_geo_soft_samples <- function(soft_gz_file) {
  connection <- gzfile(soft_gz_file, open = "rt")
  on.exit(close(connection), add = TRUE)

  lines <- readLines(
    connection,
    warn = FALSE,
    encoding = "UTF-8"
  )

  sample_start <- grep("^\\^SAMPLE =", lines)

  if (length(sample_start) == 0) {
    stop("No SAMPLE blocks were found in the SOFT file.")
  }

  sample_end <- c(
    sample_start[-1] - 1,
    length(lines)
  )

  get_values <- function(block, field_name) {
    pattern <- paste0(
      "^",
      field_name,
      "\\s*=\\s*"
    )

    values <- grep(
      pattern,
      block,
      value = TRUE
    )

    values <- sub(
      pattern,
      "",
      values
    )

    values
  }

  first_or_na <- function(values) {
    if (length(values) == 0) {
      return(NA_character_)
    }

    values[1]
  }

  extract_characteristic <- function(
    characteristics,
    key,
    occurrence = 1
  ) {
    pattern <- paste0(
      "^",
      key,
      "\\s*:\\s*"
    )

    matches <- grep(
      pattern,
      characteristics,
      value = TRUE,
      ignore.case = TRUE
    )

    if (length(matches) < occurrence) {
      return(NA_character_)
    }

    sub(
      pattern,
      "",
      matches[occurrence],
      ignore.case = TRUE
    )
  }

  sample_list <- lapply(
    seq_along(sample_start),
    function(index) {
      block <- lines[
        sample_start[index]:
          sample_end[index]
      ]

      gsm_id <- sub(
        "^\\^SAMPLE =\\s*",
        "",
        block[1]
      )

      original_title <- first_or_na(
        get_values(
          block,
          "!Sample_title"
        )
      )

      characteristics <- get_values(
        block,
        "!Sample_characteristics_ch1"
      )

      relation_values <- get_values(
        block,
        "!Sample_relation"
      )

      relation_text <- paste(
        relation_values,
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
        Sample_ID = gsm_id,
        Sample_Name = sub(
          ":.*$",
          "",
          original_title
        ),
        Patient_ID = extract_characteristic(
          characteristics,
          "subject"
        ),
        Specific_Group_Original = extract_characteristic(
          characteristics,
          "wound type",
          occurrence = 1
        ),
        Broad_Group_Original = extract_characteristic(
          characteristics,
          "wound type",
          occurrence = 2
        ),
        Tissue = extract_characteristic(
          characteristics,
          "tissue"
        ),
        Days_Since_Injury_Original = extract_characteristic(
          characteristics,
          "days since injury"
        ),
        Age_Original = extract_characteristic(
          characteristics,
          "age"
        ),
        Sex_Original = extract_characteristic(
          characteristics,
          "Sex"
        ),
        Hispanic_Original = extract_characteristic(
          characteristics,
          "hispanic"
        ),
        Race_Original = extract_characteristic(
          characteristics,
          "race"
        ),
        Burn_Type_Original = extract_characteristic(
          characteristics,
          "burn type"
        ),
        Location_Original = extract_characteristic(
          characteristics,
          "location"
        ),
        Platform_ID = first_or_na(
          get_values(
            block,
            "!Sample_platform_id"
          )
        ),
        Instrument = first_or_na(
          get_values(
            block,
            "!Sample_instrument_model"
          )
        ),
        Library_Strategy = first_or_na(
          get_values(
            block,
            "!Sample_library_strategy"
          )
        ),
        Original_Title = original_title,
        Original_Source_Name = first_or_na(
          get_values(
            block,
            "!Sample_source_name_ch1"
          )
        ),
        Original_Characteristics = paste(
          characteristics,
          collapse = " | "
        ),
        BioSample_ID = biosample_id,
        SRA_Experiment = sra_experiment,
        stringsAsFactors = FALSE
      )
    }
  )

  do.call(
    rbind,
    sample_list
  )
}

sample_metadata <- parse_geo_soft_samples(
  soft_file
)

group_mapping <- c(
  "normal skin" = "Uninjured_skin",
  "early wound" = "Early_wound",
  "late wound" = "Late_wound",
  "hts" = "Hypertrophic_scar",
  "chronic wound" = "Chronic_wound",
  "normal scar" = "Normal_scar"
)

sample_metadata$Group <- unname(
  group_mapping[
    tolower(
      sample_metadata$Specific_Group_Original
    )
  ]
)

if (any(is.na(sample_metadata$Group))) {
  unknown_groups <- unique(
    sample_metadata$Specific_Group_Original[
      is.na(sample_metadata$Group)
    ]
  )

  stop(
    paste0(
      "Unexpected sample-group labels were found in SOFT metadata: ",
      paste(
        unknown_groups,
        collapse = ", "
      )
    )
  )
}

sample_metadata$Group <- factor(
  sample_metadata$Group,
  levels = ALL_GROUP_LEVELS
)

# Convert age to numeric while preserving explicitly unknown values as NA.
# 将年龄转换为数值；GEO中明确写为unknown的年龄保留为NA。
age_clean <- trimws(
  sample_metadata$Age_Original
)

age_clean[
  tolower(age_clean) %in% c(
    "",
    "--",
    "unknown",
    "na",
    "n/a",
    "not reported"
  )
] <- NA_character_

sample_metadata$Age <- suppressWarnings(
  as.numeric(
    age_clean
  )
)

# Stop only when a non-missing age value cannot be parsed.
# 仅当某个非缺失年龄无法转换为数值时才停止。
invalid_age_format <- !is.na(age_clean) &
  is.na(sample_metadata$Age)

if (any(invalid_age_format)) {
  stop(
    paste0(
      "Unsupported non-numeric age values were found: ",
      paste(
        unique(
          sample_metadata$Age_Original[
            invalid_age_format
          ]
        ),
        collapse = ", "
      )
    )
  )
}

days_clean <- trimws(
  sample_metadata$Days_Since_Injury_Original
)

days_clean[
  tolower(days_clean) %in% c(
    "",
    "--",
    "unknown",
    "na",
    "n/a",
    "not reported"
  )
] <- NA_character_

sample_metadata$Days_Since_Injury <- suppressWarnings(
  as.numeric(
    days_clean
  )
)

invalid_day_format <- !is.na(days_clean) &
  is.na(sample_metadata$Days_Since_Injury)

if (any(invalid_day_format)) {
  stop(
    paste0(
      "Unsupported non-numeric days-since-injury values were found: ",
      paste(
        unique(
          sample_metadata$Days_Since_Injury_Original[
            invalid_day_format
          ]
        ),
        collapse = ", "
      )
    )
  )
}

sex_clean <- tolower(
  trimws(
    sample_metadata$Sex_Original
  )
)

sample_metadata$Sex <- ifelse(
  sex_clean == "male",
  "Male",
  ifelse(
    sex_clean == "female",
    "Female",
    NA_character_
  )
)

unsupported_sex <- !is.na(
  sample_metadata$Sex_Original
) &
  trimws(
    sample_metadata$Sex_Original
  ) != "" &
  is.na(
    sample_metadata$Sex
  )

if (any(unsupported_sex)) {
  stop(
    paste0(
      "Unsupported sex labels were found: ",
      paste(
        unique(
          sample_metadata$Sex_Original[
            unsupported_sex
          ]
        ),
        collapse = ", "
      )
    )
  )
}

missing_age_rows <- which(
  is.na(
    sample_metadata$Age
  )
)

if (length(missing_age_rows) > 0) {
  cat(
    "Samples with missing age will be retained in metadata but excluded ",
    "from the age-adjusted differential-expression model:\n",
    sep = ""
  )

  print(
    sample_metadata[
      missing_age_rows,
      c(
        "Sample_ID",
        "Sample_Name",
        "Patient_ID",
        "Specific_Group_Original",
        "Age_Original"
      )
    ]
  )
}

sample_metadata$Time_or_Stage <- ifelse(
  sample_metadata$Group == "Early_wound",
  "3-7 days after injury",
  ifelse(
    sample_metadata$Group == "Late_wound",
    "8-27 days after injury",
    ifelse(
      sample_metadata$Group == "Hypertrophic_scar",
      "147-5287 days after injury",
      ifelse(
        sample_metadata$Group == "Uninjured_skin",
        "Uninjured skin",
        ifelse(
          sample_metadata$Group == "Chronic_wound",
          "Chronic wound",
          "Normal scar"
        )
      )
    )
  )
)

sample_metadata$Treatment <- NA_character_
sample_metadata$Outcome <- NA_character_
sample_metadata$Data_Type <- "RNA-seq raw counts"
sample_metadata$Is_Pooled <- FALSE
sample_metadata$Metadata_Confidence <- "Direct_from_GEO_SOFT"

# Samples must belong to one of the four target groups and have complete
# age/sex covariates to enter the adjusted differential-expression model.
# 样本必须属于四个目标分组，且年龄和性别信息完整，才能进入校正模型。
sample_metadata$Eligible_Group <- sample_metadata$Group %in%
  ANALYSIS_GROUP_LEVELS

sample_metadata$Complete_Model_Covariates <- !is.na(
  sample_metadata$Age
) &
  !is.na(
    sample_metadata$Sex
  )

sample_metadata$Include_in_Analysis <- sample_metadata$Eligible_Group &
  sample_metadata$Complete_Model_Covariates

sample_metadata$Exclusion_Reason <- ifelse(
  !sample_metadata$Eligible_Group,
  "Group not included in the three prespecified contrasts.",
  ifelse(
    is.na(
      sample_metadata$Age
    ),
    "Missing age; excluded from the age-adjusted differential-expression model.",
    ifelse(
      is.na(
        sample_metadata$Sex
      ),
      "Missing sex; excluded from the sex-adjusted differential-expression model.",
      NA_character_
    )
  )
)

patient_sample_n <- table(
  sample_metadata$Patient_ID
)

sample_metadata$Is_Repeated_Measure <- as.integer(
  patient_sample_n[
    sample_metadata$Patient_ID
  ]
) > 1

patient_group_list <- split(
  as.character(
    sample_metadata$Group
  ),
  sample_metadata$Patient_ID
)

paired_patient_ids <- names(
  patient_group_list
)[
  vapply(
    patient_group_list,
    function(group_values) {
      "Uninjured_skin" %in% group_values &&
        any(
          c(
            "Early_wound",
            "Late_wound",
            "Hypertrophic_scar"
          ) %in% group_values
        )
    },
    logical(1)
  )
]

sample_metadata$Is_Paired <- sample_metadata$Patient_ID %in%
  paired_patient_ids

# Combine multiple transparent sample-level quality notes.
# 合并多个透明、可追溯的样本级质量说明。
sample_metadata$Quality_Notes <- NA_character_

append_quality_note <- function(
  existing_note,
  new_note
) {
  ifelse(
    is.na(existing_note) |
      existing_note == "",
    new_note,
    paste0(
      existing_note,
      " ",
      new_note
    )
  )
}

missing_exact_day <- sample_metadata$Group == "Late_wound" &
  is.na(
    sample_metadata$Days_Since_Injury
  )

sample_metadata$Quality_Notes[
  missing_exact_day
] <- append_quality_note(
  sample_metadata$Quality_Notes[
    missing_exact_day
  ],
  "Exact days since injury are missing in GEO metadata."
)

missing_age <- is.na(
  sample_metadata$Age
)

sample_metadata$Quality_Notes[
  missing_age
] <- append_quality_note(
  sample_metadata$Quality_Notes[
    missing_age
  ],
  paste0(
    "Age is recorded as unknown in GEO; sample retained in metadata ",
    "but excluded from the age-adjusted differential-expression model."
  )
)

# Verify sample names and reorder metadata to match the count matrix.
# 检查样本名并按计数矩阵列顺序重新排列元数据。
if (!setequal(
  sample_metadata$Sample_Name,
  colnames(count_matrix)
)) {
  missing_from_metadata <- setdiff(
    colnames(count_matrix),
    sample_metadata$Sample_Name
  )

  missing_from_matrix <- setdiff(
    sample_metadata$Sample_Name,
    colnames(count_matrix)
  )

  stop(
    paste0(
      "Count-matrix and SOFT sample names do not match. ",
      "Missing from metadata: ",
      paste(
        missing_from_metadata,
        collapse = ", "
      ),
      ". Missing from count matrix: ",
      paste(
        missing_from_matrix,
        collapse = ", "
      ),
      "."
    )
  )
}

sample_metadata <- sample_metadata[
  match(
    colnames(count_matrix),
    sample_metadata$Sample_Name
  ),
]

if (!identical(
  sample_metadata$Sample_Name,
  colnames(count_matrix)
)) {
  stop(
    "Sample metadata could not be reordered to match the count matrix."
  )
}

sample_metadata$Raw_Library_Size <- colSums(
  count_matrix
)

sample_metadata$TMM_Normalization_Factor <- NA_real_
sample_metadata$Effective_Library_Size <- NA_real_

# Keep the requested core metadata columns first.
# 将用户要求的核心元数据字段置于前面。
sample_metadata <- sample_metadata[, c(
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
  "Days_Since_Injury",
  "Age",
  "Sex",
  "Hispanic_Original",
  "Race_Original",
  "Burn_Type_Original",
  "Location_Original",
  "Instrument",
  "Library_Strategy",
  "Specific_Group_Original",
  "Broad_Group_Original",
  "Eligible_Group",
  "Complete_Model_Covariates",
  "Include_in_Analysis",
  "Exclusion_Reason",
  "Raw_Library_Size",
  "TMM_Normalization_Factor",
  "Effective_Library_Size"
)]

actual_group_counts <- table(
  sample_metadata$Group
)

actual_patient_n <- length(
  unique(
    sample_metadata$Patient_ID
  )
)

metadata_audit <- data.frame(
  Check = c(
    "Count-matrix sample count",
    "SOFT sample count",
    "Unique patient count",
    "Actual uninjured-skin samples",
    "Actual early-wound samples",
    "Actual late-wound samples",
    "Actual chronic-wound samples",
    "Actual normal-scar samples",
    "Actual hypertrophic-scar samples",
    "Actual total wound samples",
    "Actual total scar samples",
    "GEO summary reported uninjured skin",
    "GEO summary reported acute wounds",
    "GEO summary reported hypertrophic scars",
    "GEO summary subgroup sum",
    "Samples in the four eligible groups",
    "Samples excluded because age is unknown",
    "Samples included in the adjusted differential-expression model",
    "Samples excluded from the adjusted differential-expression model"
  ),
  Value = c(
    ncol(count_matrix),
    nrow(sample_metadata),
    actual_patient_n,
    actual_group_counts["Uninjured_skin"],
    actual_group_counts["Early_wound"],
    actual_group_counts["Late_wound"],
    actual_group_counts["Chronic_wound"],
    actual_group_counts["Normal_scar"],
    actual_group_counts["Hypertrophic_scar"],
    actual_group_counts["Early_wound"] +
      actual_group_counts["Late_wound"] +
      actual_group_counts["Chronic_wound"],
    actual_group_counts["Normal_scar"] +
      actual_group_counts["Hypertrophic_scar"],
    26,
    54,
    30,
    110,
    sum(sample_metadata$Eligible_Group),
    sum(
      sample_metadata$Eligible_Group &
        is.na(sample_metadata$Age)
    ),
    sum(sample_metadata$Include_in_Analysis),
    sum(!sample_metadata$Include_in_Analysis)
  ),
  Interpretation = c(
    "Directly observed in the uploaded count matrix.",
    "Directly parsed from the uploaded SOFT file.",
    "Calculated from the SOFT subject field.",
    "Directly parsed from sample-level SOFT metadata.",
    "Directly parsed from sample-level SOFT metadata.",
    "Directly parsed from sample-level SOFT metadata.",
    "Retained in metadata but excluded from differential analysis.",
    "Retained in metadata but excluded from differential analysis.",
    "Directly parsed from sample-level SOFT metadata.",
    "Early plus late plus chronic wound.",
    "Normal scar plus hypertrophic scar.",
    "Reported in the GEO Series summary.",
    "Reported in the GEO Series summary.",
    "Reported in the GEO Series summary.",
    "The reported subgroup counts sum to 110, not 108.",
    "Uninjured skin, early wound, late wound, and hypertrophic scar before covariate filtering.",
    "GSM5390619 and GSM5390626 have age recorded as unknown.",
    "Eligible groups with complete age and sex covariates.",
    "Five samples belong to excluded groups and two eligible samples have unknown age."
  ),
  stringsAsFactors = FALSE
)

write.csv(
  metadata_audit,
  file = file.path(
    QC_DIR,
    "GSE178411_metadata_consistency_check.csv"
  ),
  row.names = FALSE
)

write.csv(
  sample_metadata,
  file = file.path(
    RESULTS_DIR,
    "GSE178411_sample_metadata.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

cat("Sample-group counts:\n")
print(actual_group_counts)

cat("\nUnique patients:\n")
print(actual_patient_n)

cat("\nMetadata consistency audit:\n")
print(metadata_audit)


# ---- 05. Raw-data QC / 原始数据质量控制 ----

library_size_table <- data.frame(
  Sample_Name = sample_metadata$Sample_Name,
  Group = factor(
    sample_metadata$Group,
    levels = ALL_GROUP_LEVELS
  ),
  Library_Size = sample_metadata$Raw_Library_Size,
  Library_Size_Million = sample_metadata$Raw_Library_Size / 1e6,
  stringsAsFactors = FALSE
)

write.csv(
  library_size_table,
  file = file.path(
    QC_DIR,
    "GSE178411_library_sizes.csv"
  ),
  row.names = FALSE
)

p_raw_qc <- ggplot(
  library_size_table,
  aes(
    x = Group,
    y = Library_Size_Million,
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
    size = 1.7,
    alpha = 0.75
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
    title = "GSE178411 raw library-size QC",
    subtitle = "Each point represents one RNA-seq sample",
    x = "Sample group",
    y = "Library size (million counts)"
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
    "01_GSE178411_raw_data_QC.png"
  ),
  plot = p_raw_qc,
  width = 9,
  height = 5.5,
  dpi = 300
)


# ---- 06. Filtering, TMM normalization, and voom transformation / 过滤、TMM标准化与voom转换 ----

analysis_metadata <- sample_metadata[
  sample_metadata$Include_in_Analysis,
]

if (any(is.na(analysis_metadata$Age))) {
  stop(
    "The analysis metadata still contains missing age values after complete-case filtering."
  )
}

if (any(is.na(analysis_metadata$Sex))) {
  stop(
    "The analysis metadata still contains missing sex values after complete-case filtering."
  )
}

if (any(is.na(analysis_metadata$Patient_ID)) |
    any(analysis_metadata$Patient_ID == "")) {
  stop(
    "The analysis metadata contains missing Patient_ID values."
  )
}

cat(
  "Samples included in the adjusted model: ",
  nrow(analysis_metadata),
  "\n",
  sep = ""
)

cat("Adjusted-model group counts:\n")
print(
  table(
    analysis_metadata$Group
  )
)

analysis_count_matrix <- count_matrix[
  ,
  analysis_metadata$Sample_Name,
  drop = FALSE
]

analysis_group <- factor(
  as.character(
    analysis_metadata$Group
  ),
  levels = ANALYSIS_GROUP_LEVELS
)

analysis_metadata$Group <- analysis_group

dge <- edgeR::DGEList(
  counts = analysis_count_matrix,
  group = analysis_group,
  genes = data.frame(
    Entrez_ID = rownames(
      analysis_count_matrix
    ),
    stringsAsFactors = FALSE
  )
)

keep_gene <- edgeR::filterByExpr(
  dge,
  group = analysis_group
)

cat(
  "Genes before filtering: ",
  nrow(dge),
  "\nGenes retained after filtering: ",
  sum(keep_gene),
  "\nGenes removed by filterByExpr: ",
  sum(!keep_gene),
  "\n",
  sep = ""
)

dge_filtered <- dge[
  keep_gene,
  ,
  keep.lib.sizes = FALSE
]

dge_filtered <- edgeR::calcNormFactors(
  dge_filtered,
  method = "TMM"
)

analysis_metadata$Age_Centered <- analysis_metadata$Age -
  mean(
    analysis_metadata$Age
  )

analysis_metadata$Sex <- factor(
  analysis_metadata$Sex,
  levels = c(
    "Female",
    "Male"
  )
)

design <- model.matrix(
  ~ 0 + Group + Age_Centered + Sex,
  data = analysis_metadata
)

colnames(design) <- sub(
  "^Group",
  "",
  colnames(design)
)

rownames(design) <- analysis_metadata$Sample_Name

if (qr(design)$rank < ncol(design)) {
  stop(
    "The differential-expression design matrix is not full rank."
  )
}

cat("\nDifferential-expression design matrix columns:\n")
print(colnames(design))

# First voom pass for initial correlation estimation.
# 第一次voom用于初步估计患者内相关性。
voom_initial <- limma::voom(
  dge_filtered,
  design = design,
  plot = FALSE
)

correlation_initial <- limma::duplicateCorrelation(
  voom_initial,
  design = design,
  block = analysis_metadata$Patient_ID
)

if (!is.finite(
  correlation_initial$consensus.correlation
)) {
  stop(
    "Within-patient correlation could not be estimated."
  )
}

# Second voom pass using the initial within-patient correlation.
# 第二次voom使用初步患者内相关性，并保存均值-方差趋势图。
png(
  filename = file.path(
    QC_DIR,
    "GSE178411_voom_mean_variance.png"
  ),
  width = 1800,
  height = 1400,
  res = 220
)

voom_object <- limma::voom(
  dge_filtered,
  design = design,
  block = analysis_metadata$Patient_ID,
  correlation = correlation_initial$consensus.correlation,
  plot = TRUE
)

dev.off()

correlation_final <- limma::duplicateCorrelation(
  voom_object,
  design = design,
  block = analysis_metadata$Patient_ID
)

within_patient_correlation <- correlation_final$consensus.correlation

if (!is.finite(
  within_patient_correlation
)) {
  stop(
    "The final within-patient correlation estimate is not finite."
  )
}

cat(
  "\nEstimated within-patient correlation: ",
  round(
    within_patient_correlation,
    4
  ),
  "\n",
  sep = ""
)

normalized_expression <- voom_object$E

normalized_output <- data.frame(
  NCBI_Gene_ID = rownames(
    normalized_expression
  ),
  normalized_expression,
  check.names = FALSE
)

normalized_connection <- gzfile(
  file.path(
    RESULTS_DIR,
    "GSE178411_normalized_expression.csv.gz"
  ),
  open = "wt"
)

write.csv(
  normalized_output,
  file = normalized_connection,
  row.names = FALSE
)

close(normalized_connection)

# Add normalization metrics to the full sample metadata.
# 将标准化参数补充到完整样本元数据表。
analysis_index <- match(
  analysis_metadata$Sample_Name,
  sample_metadata$Sample_Name
)

sample_metadata$TMM_Normalization_Factor[
  analysis_index
] <- dge_filtered$samples$norm.factors

sample_metadata$Effective_Library_Size[
  analysis_index
] <- with(
  dge_filtered$samples,
  lib.size * norm.factors
)

write.csv(
  sample_metadata,
  file = file.path(
    RESULTS_DIR,
    "GSE178411_sample_metadata.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# Plot the distribution of sample-level median normalized expression.
# 绘制每个样本标准化表达中位数的分布。
normalized_distribution <- data.frame(
  Sample_Name = colnames(
    normalized_expression
  ),
  Median_log2CPM = apply(
    normalized_expression,
    2,
    median
  ),
  Group = analysis_metadata$Group,
  stringsAsFactors = FALSE
)

p_normalized <- ggplot(
  normalized_distribution,
  aes(
    x = Group,
    y = Median_log2CPM,
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
    size = 1.7,
    alpha = 0.75
  ) +
  scale_color_manual(
    values = GROUP_COLORS[
      ANALYSIS_GROUP_LEVELS
    ],
    drop = FALSE
  ) +
  scale_fill_manual(
    values = GROUP_COLORS[
      ANALYSIS_GROUP_LEVELS
    ],
    drop = FALSE
  ) +
  labs(
    title = "GSE178411 normalized expression distribution",
    subtitle = paste0(
      "Each point is the median voom log2 CPM across filtered genes ",
      "for one sample"
    ),
    x = "Sample group",
    y = "Median log2 CPM"
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
    "02_GSE178411_normalized_expression_distribution.png"
  ),
  plot = p_normalized,
  width = 9,
  height = 5.5,
  dpi = 300
)


# ---- 07. PCA and sample correlation / PCA与样本相关性 ----

gene_variance <- apply(
  normalized_expression,
  1,
  var
)

n_pca_genes <- min(
  TOP_VARIABLE_GENES_FOR_PCA,
  length(gene_variance)
)

top_variable_ids <- names(
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
    normalized_expression[
      top_variable_ids,
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
  Sample_Name = rownames(
    pca_result$x
  ),
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  stringsAsFactors = FALSE
)

pca_table$Sample_ID <- analysis_metadata$Sample_ID[
  match(
    pca_table$Sample_Name,
    analysis_metadata$Sample_Name
  )
]

pca_table$Patient_ID <- analysis_metadata$Patient_ID[
  match(
    pca_table$Sample_Name,
    analysis_metadata$Sample_Name
  )
]

pca_table$Group <- analysis_metadata$Group[
  match(
    pca_table$Sample_Name,
    analysis_metadata$Sample_Name
  )
]

write.csv(
  pca_table,
  file = file.path(
    RESULTS_DIR,
    "GSE178411_PCA_scores.csv"
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
    size = 2.5,
    alpha = 0.78
  ) +
  scale_color_manual(
    values = GROUP_COLORS[
      ANALYSIS_GROUP_LEVELS
    ],
    drop = FALSE
  ) +
  labs(
    title = "GSE178411 PCA",
    subtitle = paste0(
      "Top ",
      n_pca_genes,
      " variable genes; ",
      nrow(analysis_metadata),
      " samples included in the adjusted model"
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
  theme_classic(base_size = 12)

print(p_pca)

ggsave(
  filename = file.path(
    FIGURES_DIR,
    "03_GSE178411_PCA.png"
  ),
  plot = p_pca,
  width = 8,
  height = 6,
  dpi = 300
)

sample_correlation <- cor(
  normalized_expression,
  method = "pearson"
)

correlation_annotation <- data.frame(
  Group = analysis_metadata$Group,
  row.names = analysis_metadata$Sample_Name
)

correlation_annotation_colors <- list(
  Group = GROUP_COLORS[
    ANALYSIS_GROUP_LEVELS
  ]
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
  main = "GSE178411 sample correlation"
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
  main = "GSE178411 sample correlation",
  filename = file.path(
    FIGURES_DIR,
    "04_GSE178411_sample_correlation_heatmap.png"
  ),
  width = 11,
  height = 10
)


# ---- 08. Differential-expression model / 差异表达模型 ----

# Fit the repeated-measures model using the final common correlation estimate.
# 使用最终患者内共同相关系数拟合重复测量模型。
fit <- limma::lmFit(
  voom_object,
  design = design,
  block = analysis_metadata$Patient_ID,
  correlation = within_patient_correlation
)

contrast_matrix <- limma::makeContrasts(
  GSE178411_EarlyWound_vs_UninjuredSkin =
    Early_wound - Uninjured_skin,
  GSE178411_LateWound_vs_UninjuredSkin =
    Late_wound - Uninjured_skin,
  GSE178411_HypertrophicScar_vs_UninjuredSkin =
    Hypertrophic_scar - Uninjured_skin,
  levels = design
)

fit_contrasts <- limma::contrasts.fit(
  fit,
  contrasts = contrast_matrix
)

fit_contrasts <- limma::eBayes(
  fit_contrasts,
  robust = TRUE
)

contrast_definitions <- data.frame(
  Contrast_ID = colnames(
    contrast_matrix
  ),
  Contrast_Label = c(
    "Early wound vs uninjured skin",
    "Late wound vs uninjured skin",
    "Hypertrophic scar vs uninjured skin"
  ),
  Case_Group = c(
    "Early wound",
    "Late wound",
    "Hypertrophic scar"
  ),
  Case_Group_Code = c(
    "Early_wound",
    "Late_wound",
    "Hypertrophic_scar"
  ),
  Control_Group = "Uninjured skin",
  Control_Group_Code = "Uninjured_skin",
  Sample_Context = c(
    "Acute burn wound tissue",
    "Acute burn wound tissue",
    "Hypertrophic scar tissue"
  ),
  Time_or_Stage = c(
    "3-7 days after injury",
    "8-27 days after injury",
    "147-5287 days after injury"
  ),
  stringsAsFactors = FALSE
)

contrast_definitions$Case_N <- vapply(
  contrast_definitions$Case_Group_Code,
  function(group_code) {
    sum(
      analysis_metadata$Group == group_code
    )
  },
  numeric(1)
)

contrast_definitions$Control_N <- vapply(
  contrast_definitions$Control_Group_Code,
  function(group_code) {
    sum(
      analysis_metadata$Group == group_code
    )
  },
  numeric(1)
)

contrast_definitions$Case_Patient_N <- vapply(
  contrast_definitions$Case_Group_Code,
  function(group_code) {
    length(
      unique(
        analysis_metadata$Patient_ID[
          analysis_metadata$Group == group_code
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
        analysis_metadata$Patient_ID[
          analysis_metadata$Group == group_code
        ]
      )
    )
  },
  numeric(1)
)

cat("Contrast definitions:\n")
print(contrast_definitions)


# ---- 09. Gene-ID mapping and complete result table / Gene ID映射与完整结果表 ----

filtered_entrez_ids <- rownames(
  voom_object$E
)

gene_symbol <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = filtered_entrez_ids,
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
)

ensembl_id <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = filtered_entrez_ids,
  column = "ENSEMBL",
  keytype = "ENTREZID",
  multiVals = "first"
)

gene_name <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = filtered_entrez_ids,
  column = "GENENAME",
  keytype = "ENTREZID",
  multiVals = "first"
)

annotation_table <- data.frame(
  NCBI_Gene_ID = filtered_entrez_ids,
  Gene_Symbol = unname(
    gene_symbol[
      filtered_entrez_ids
    ]
  ),
  Ensembl_ID = unname(
    ensembl_id[
      filtered_entrez_ids
    ]
  ),
  Gene_Name = unname(
    gene_name[
      filtered_entrez_ids
    ]
  ),
  stringsAsFactors = FALSE
)

annotation_table$Mapping_Status <- ifelse(
  is.na(
    annotation_table$Gene_Symbol
  ),
  "Unmapped",
  "Mapped"
)

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

  contrast_result$NCBI_Gene_ID <- rownames(
    contrast_result
  )

  contrast_result$Gene_Symbol <- annotation_data$Gene_Symbol[
    match(
      contrast_result$NCBI_Gene_ID,
      annotation_data$NCBI_Gene_ID
    )
  ]

  contrast_result$Ensembl_ID <- annotation_data$Ensembl_ID[
    match(
      contrast_result$NCBI_Gene_ID,
      annotation_data$NCBI_Gene_ID
    )
  ]

  contrast_result$Gene_Name <- annotation_data$Gene_Name[
    match(
      contrast_result$NCBI_Gene_ID,
      annotation_data$NCBI_Gene_ID
    )
  ]

  contrast_result$Mapping_Status <- annotation_data$Mapping_Status[
    match(
      contrast_result$NCBI_Gene_ID,
      annotation_data$NCBI_Gene_ID
    )
  ]

  contrast_result$Contrast_ID <- contrast_id
  contrast_result$Contrast_Label <- contrast_info$Contrast_Label
  contrast_result$Case_Group <- contrast_info$Case_Group
  contrast_result$Control_Group <- contrast_info$Control_Group
  contrast_result$Case_N <- contrast_info$Case_N
  contrast_result$Control_N <- contrast_info$Control_N
  contrast_result$Case_Patient_N <- contrast_info$Case_Patient_N
  contrast_result$Control_Patient_N <- contrast_info$Control_Patient_N

  contrast_result$Fold_Change <- 2^contrast_result$logFC

  contrast_result$Direction <- ifelse(
    contrast_result$logFC > 0,
    "Up",
    ifelse(
      contrast_result$logFC < 0,
      "Down",
      "No_change"
    )
  )

  contrast_result$DE_Status <- ifelse(
    contrast_result$adj.P.Val < FDR_CUTOFF &
      contrast_result$logFC >= LOG2FC_CUTOFF,
    "Up_significant",
    ifelse(
      contrast_result$adj.P.Val < FDR_CUTOFF &
        contrast_result$logFC <= -LOG2FC_CUTOFF,
      "Down_significant",
      "Not_significant"
    )
  )

  contrast_result$NegLog10_FDR <- -log10(
    pmax(
      contrast_result$adj.P.Val,
      .Machine$double.xmin
    )
  )

  output <- data.frame(
    NCBI_Gene_ID = contrast_result$NCBI_Gene_ID,
    Gene_Symbol = contrast_result$Gene_Symbol,
    Ensembl_ID = contrast_result$Ensembl_ID,
    Gene_Name = contrast_result$Gene_Name,
    Mapping_Status = contrast_result$Mapping_Status,
    Contrast_ID = contrast_result$Contrast_ID,
    Contrast_Label = contrast_result$Contrast_Label,
    Case_Group = contrast_result$Case_Group,
    Control_Group = contrast_result$Control_Group,
    Case_N = contrast_result$Case_N,
    Control_N = contrast_result$Control_N,
    Case_Patient_N = contrast_result$Case_Patient_N,
    Control_Patient_N = contrast_result$Control_Patient_N,
    log2FC = contrast_result$logFC,
    Fold_Change = contrast_result$Fold_Change,
    Mean_log2CPM = contrast_result$AveExpr,
    Statistic = contrast_result$t,
    P_value = contrast_result$P.Value,
    FDR = contrast_result$adj.P.Val,
    B_statistic = contrast_result$B,
    Direction = contrast_result$Direction,
    DE_Status = contrast_result$DE_Status,
    NegLog10_FDR = contrast_result$NegLog10_FDR,
    stringsAsFactors = FALSE
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
      annotation_data = annotation_table,
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
    "GSE178411_all_gene_results.csv.gz"
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
  contrast_id in contrast_definitions$Contrast_ID
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
      alpha = 0.62,
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
        "GSE178411: ",
        contrast_info$Contrast_Label
      ),
      subtitle = paste0(
        "TMM + voom repeated-measures model; ",
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
        "05_GSE178411_volcano_",
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
  contrast_id in contrast_definitions$Contrast_ID
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

  contrast_sample_names <- analysis_metadata$Sample_Name[
    analysis_metadata$Group %in%
      c(
        contrast_info$Control_Group_Code,
        contrast_info$Case_Group_Code
      )
  ]

  heatmap_matrix <- normalized_expression[
    top_heatmap_ids,
    contrast_sample_names,
    drop = FALSE
  ]

  heatmap_z <- t(
    scale(
      t(
        heatmap_matrix
      )
    )
  )

  heatmap_labels <- contrast_results$Gene_Symbol[
    match(
      top_heatmap_ids,
      contrast_results$NCBI_Gene_ID
    )
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

  rownames(heatmap_z) <- make.unique(
    heatmap_labels
  )

  heatmap_annotation <- data.frame(
    Group = analysis_metadata$Group[
      match(
        contrast_sample_names,
        analysis_metadata$Sample_Name
      )
    ],
    row.names = contrast_sample_names
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

  show_sample_names <- ncol(
    heatmap_z
  ) <= 30

  pheatmap::pheatmap(
    heatmap_z,
    annotation_col = heatmap_annotation,
    annotation_colors = heatmap_annotation_colors,
    show_colnames = show_sample_names,
    show_rownames = TRUE,
    border_color = NA,
    color = EXPRESSION_HEATMAP_COLORS,
    main = heatmap_title,
    cluster_rows = TRUE,
    cluster_cols = TRUE
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
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    filename = file.path(
      FIGURES_DIR,
      paste0(
        "06_GSE178411_top_differential_genes_heatmap_",
        contrast_id,
        ".png"
      )
    ),
    width = 10,
    height = 10
  )
}


# ---- 12. Metadata consistency review / 元数据一致性复核 ----

# No peer-reviewed paper is linked to GEO, so this section summarizes
# internal consistency instead of reproducing published results.
# GEO未关联同行评议论文，因此本节复核内部一致性，而不复现论文结果。

cat("Metadata consistency review:\n")
cat(
  "- Count-matrix samples: ",
  ncol(count_matrix),
  "\n",
  sep = ""
)
cat(
  "- SOFT samples: ",
  nrow(sample_metadata),
  "\n",
  sep = ""
)
cat(
  "- Unique patients: ",
  actual_patient_n,
  "\n",
  sep = ""
)
cat(
  "- Actual sample groups: 24 uninjured skin, 22 early wound, ",
  "29 late wound, 3 chronic wound, 2 normal scar, ",
  "28 hypertrophic scar.\n",
  sep = ""
)
cat(
  "- Two eligible samples have age recorded as unknown and are excluded ",
  "from the age-adjusted model: GSM5390619 and GSM5390626.\n",
  sep = ""
)
cat(
  "- Adjusted-model group counts: 24 uninjured skin, 22 early wound, ",
  "28 late wound, and 27 hypertrophic scar.\n",
  sep = ""
)
cat(
  "- GEO summary subgroup counts sum to 110 although the Series contains ",
  "108 samples.\n",
  sep = ""
)
cat(
  "- GEO summary mentions Illumina Hi-Seq, whereas the registered platform ",
  "and sample metadata specify Illumina NovaSeq 6000.\n",
  sep = ""
)
cat(
  "- No linked peer-reviewed citation is available in the GEO record.\n",
  sep = ""
)


# ---- 13. Create the BurnOmicsDB-ready result table / 创建BurnOmicsDB标准结果表 ----

contrast_match <- match(
  all_gene_results$Contrast_ID,
  contrast_definitions$Contrast_ID
)

database_ready <- data.frame(
  NCBI_Gene_ID = all_gene_results$NCBI_Gene_ID,
  Gene_Symbol = all_gene_results$Gene_Symbol,
  Ensembl_ID = all_gene_results$Ensembl_ID,
  Gene_Name = all_gene_results$Gene_Name,
  GEO_ID = GEO_ID,
  Organism = "Homo sapiens",
  Study_Population = "Mixed-age human burn-center cohort",
  Tissue = "Skin",
  Sample_Context = contrast_definitions$Sample_Context[
    contrast_match
  ],
  Time_or_Stage = contrast_definitions$Time_or_Stage[
    contrast_match
  ],
  Contrast_ID = all_gene_results$Contrast_ID,
  Contrast_Label = all_gene_results$Contrast_Label,
  Case_Group = all_gene_results$Case_Group,
  Control_Group = all_gene_results$Control_Group,
  Case_N = all_gene_results$Case_N,
  Control_N = all_gene_results$Control_N,
  log2FC = all_gene_results$log2FC,
  Fold_Change = all_gene_results$Fold_Change,
  Direction = all_gene_results$Direction,
  Mean_log2CPM = all_gene_results$Mean_log2CPM,
  P_value = all_gene_results$P_value,
  FDR = all_gene_results$FDR,
  DE_Status = all_gene_results$DE_Status,
  Platform = "Illumina NovaSeq 6000 (GPL24676)",
  Input_Data = "Author-provided NCBI Gene ID raw counts",
  Normalization = "TMM followed by voom log2 CPM",
  Analysis_Method = paste0(
    "edgeR filterByExpr and TMM; limma-voom linear model with ",
    "duplicateCorrelation blocking Patient_ID; adjusted for age and sex; ",
    "positive log2FC represents Case_Group minus Control_Group"
  ),
  Annotation_Method = paste0(
    "org.Hs.eg.db ",
    as.character(
      packageVersion(
        "org.Hs.eg.db"
      )
    ),
    "; Entrez ID mapped with mapIds(multiVals='first')"
  ),
  Quality_Notes = paste0(
    "No linked peer-reviewed citation in GEO; actual sample-level SOFT ",
    "counts differ from the GEO summary subgroup counts; two eligible ",
    "samples with age recorded as unknown were excluded from the age-adjusted ",
    "model; repeated samples from the same patient were modeled with ",
    "duplicateCorrelation; ",
    "uninjured skin represents donor or uninjured skin from a burn-center ",
    "repository and is not necessarily healthy-volunteer skin; bulk-tissue ",
    "differences may reflect both expression regulation and cell-composition ",
    "changes"
  ),
  Case_Patient_N = all_gene_results$Case_Patient_N,
  Control_Patient_N = all_gene_results$Control_Patient_N,
  Within_Patient_Correlation = within_patient_correlation,
  Age_Adjusted = TRUE,
  Sex_Adjusted = TRUE,
  stringsAsFactors = FALSE
)

database_connection <- gzfile(
  file.path(
    RESULTS_DIR,
    "GSE178411_database_ready_all_genes.csv.gz"
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


# ---- 14. Save analysis objects, summary, and session information / 保存分析对象、摘要与环境信息 ----

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
  count_matrix_dimensions = dim(
    count_matrix
  ),
  sample_metadata = sample_metadata,
  analysis_metadata = analysis_metadata,
  metadata_audit = metadata_audit,
  filtered_edgeR_object = dge_filtered,
  design_matrix = design,
  initial_within_patient_correlation =
    correlation_initial$consensus.correlation,
  final_within_patient_correlation =
    within_patient_correlation,
  voom_object = voom_object,
  contrast_matrix = contrast_matrix,
  contrast_definitions = contrast_definitions,
  fitted_model = fit_contrasts,
  annotation_table = annotation_table,
  all_gene_results = all_gene_results,
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
    "GSE178411_analysis_objects.rds"
  ),
  compress = "xz"
)

summary_lines <- c(
  "BurnOmicsDB - GSE178411 analysis summary",
  "",
  paste0(
    "Analysis date: ",
    Sys.Date()
  ),
  paste0(
    "Input genes: ",
    nrow(count_matrix)
  ),
  paste0(
    "Input samples: ",
    ncol(count_matrix)
  ),
  paste0(
    "Unique patients in the full metadata: ",
    actual_patient_n
  ),
  paste0(
    "Samples included in the differential-expression model: ",
    nrow(analysis_metadata)
  ),
  paste0(
    "Patients represented in the differential-expression model: ",
    length(
      unique(
        analysis_metadata$Patient_ID
      )
    )
  ),
  paste0(
    "Samples excluded from differential-expression modeling: ",
    sum(
      !sample_metadata$Include_in_Analysis
    ),
    " (3 chronic-wound samples, 2 normal-scar samples, and 2 eligible ",
    "samples with age recorded as unknown)"
  ),
  "",
  "Actual sample groups parsed from SOFT:",
  paste0(
    "  Uninjured skin: ",
    actual_group_counts["Uninjured_skin"]
  ),
  paste0(
    "  Early wound: ",
    actual_group_counts["Early_wound"]
  ),
  paste0(
    "  Late wound: ",
    actual_group_counts["Late_wound"]
  ),
  paste0(
    "  Chronic wound: ",
    actual_group_counts["Chronic_wound"]
  ),
  paste0(
    "  Normal scar: ",
    actual_group_counts["Normal_scar"]
  ),
  paste0(
    "  Hypertrophic scar: ",
    actual_group_counts["Hypertrophic_scar"]
  ),
  "",
  paste0(
    "Genes retained after filterByExpr: ",
    nrow(dge_filtered)
  ),
  "Normalization: TMM followed by voom log2 CPM",
  paste0(
    "Final estimated within-patient correlation: ",
    round(
      within_patient_correlation,
      6
    )
  ),
  paste0(
    "Statistical model: limma-voom with duplicateCorrelation blocking ",
    "Patient_ID and adjustment for age and sex"
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
    " (case samples = ",
    contrast_definitions$Case_N,
    ", case patients = ",
    contrast_definitions$Case_Patient_N,
    ", control samples = ",
    contrast_definitions$Control_N,
    ", control patients = ",
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
    "Positive log2FC means higher expression in the case group than in ",
    "uninjured skin."
  ),
  paste0(
    "Only genes retained by filterByExpr were statistically tested and ",
    "exported."
  ),
  "",
  "Major limitations:",
  paste0(
    "  No peer-reviewed citation is linked to the GEO Series record."
  ),
  paste0(
    "  GEO summary subgroup counts sum to 110, whereas the count matrix ",
    "and SOFT file contain 108 samples."
  ),
  paste0(
    "  The GEO summary mentions Illumina Hi-Seq, whereas the registered ",
    "platform and sample metadata specify Illumina NovaSeq 6000."
  ),
  paste0(
    "  Two eligible samples have age recorded as unknown in GEO and were ",
    "excluded from the age-adjusted differential-expression model."
  ),
  paste0(
    "  Uninjured skin is donor or uninjured skin from a burn-center ",
    "repository and is not necessarily healthy-volunteer skin."
  ),
  paste0(
    "  Bulk-tissue expression differences can reflect both intracellular ",
    "regulation and changes in cell composition or tissue structure."
  )
)

writeLines(
  summary_lines,
  con = file.path(
    RESULTS_DIR,
    "GSE178411_analysis_summary.txt"
  ),
  useBytes = TRUE
)

capture.output(
  sessionInfo(),
  file = file.path(
    RESULTS_DIR,
    "GSE178411_R_sessionInfo.txt"
  )
)

cat("Analysis completed.\n")
cat("Results directory:\n", RESULTS_DIR, "\n")
cat("Figures directory:\n", FIGURES_DIR, "\n")
cat("QC directory:\n", QC_DIR, "\n")
cat("R objects directory:\n", OBJECTS_DIR, "\n")


# ---- 15. Optional individual-gene checking / 可选的单基因检查 ----

# Run this section after Section 09.
# 在完成第09部分后运行本节。

genes_to_check <- c(
  "IL6",
  "MMP1",
  "MMP8",
  "MMP9",
  "MMP13",
  "COL1A1",
  "COL3A1",
  "TGFB1",
  "KRT1",
  "KRT10"
)

gene_check_table <- all_gene_results[
  all_gene_results$Gene_Symbol %in%
    genes_to_check,
  c(
    "NCBI_Gene_ID",
    "Gene_Symbol",
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
