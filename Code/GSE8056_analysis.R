# ==============================================================================
# BurnOmicsDB: GSE8056 microarray analysis
# BurnOmicsDB：GSE8056微阵列分析
#
# GEO accession / GEO编号:
#   GSE8056
#
# Study / 研究:
#   Gene Expression Profiles in Thermally Injured Human Skin:
#   A Temporal Microarray Analysis
#   热损伤人体皮肤的时间序列基因表达微阵列分析
#
# Data type / 数据类型:
#   Affymetrix Human Genome U133 Plus 2.0 Array (GPL570)
#   Raw CEL files contained in GSE8056_RAW.tar
#   GPL570 Affymetrix原始CEL文件
#
# Experimental structure / 实验结构:
#   Normal skin:  3 pooled arrays, 5 patients per pool
#   Early wound:  3 pooled arrays, 5 patients per pool
#   Middle wound: 3 pooled arrays, 5 patients per pool
#   Late wound:   3 pooled arrays, 5 patients per pool
#
#   A total of 60 patients are represented, but the statistical unit is the
#   pooled array. Therefore, each group has n = 3 statistical replicates,
#   not n = 15 independent patient-level expression profiles.
#   本研究代表60名患者，但统计单位是pool芯片。因此每组统计学n=3，
#   不能将每组15名患者视为15个独立表达样本。
#
# Prespecified contrasts / 预设比较:
#   1. Early wound versus normal skin
#   2. Middle wound versus normal skin
#   3. Late wound versus normal skin
#
# Positive log2FC always means higher expression in the burn-wound case group
# than in normal skin.
# 正log2FC始终表示烧伤创缘组相对于正常皮肤组表达更高。
#
# Analysis workflow / 分析流程:
#   Raw CEL
#   -> RMA background correction
#   -> quantile normalization
#   -> median-polish probe-set summarization
#   -> probe-to-NCBI-Gene-ID mapping
#   -> one representative probe per NCBI Gene ID selected using the highest
#      mean RMA expression across all 12 arrays
#   -> limma linear model and empirical Bayes statistics
#
# Probe-selection rationale / 代表探针选择理由:
#   The representative probe is selected before differential testing and uses
#   all 12 arrays without considering case/control labels. This avoids choosing
#   a probe because it has the smallest p-value in a particular contrast.
#   代表探针在差异分析之前选择，并使用全部12张芯片的平均表达，
#   不依据任何contrast的p值，避免人为放大显著性。
#
# Pooling and covariates / Pooling与协变量:
#   Patient-level age, sex, TBSA, anatomical location, and post-burn day are
#   aggregated within pools. They cannot be included as independent patient-level
#   covariates. No patient-level covariate adjustment is performed.
#   年龄、性别、TBSA、部位和伤后时间均被pool聚合，
#   无法作为患者级独立协变量，因此不进行患者级协变量校正。
#
# Quality classification / 质量分类:
#   Auxiliary / lower-confidence evidence because:
#   - only 3 pooled arrays per group are available;
#   - individual-patient variability is not observable;
#   - normal and burn groups differ substantially in sex composition;
#   - control skin was collected during elective cosmetic surgery;
#   - one SOFT control-pool metadata record appears inconsistent with Table 1
#     of the associated paper.
#   由于每组仅3个pool、无法观察患者个体差异、性别构成不平衡、
#   正常组织来源不同以及一条SOFT元数据疑点，本数据作为辅助/较低置信证据。
#
# How to run in RStudio / 如何在RStudio中运行:
#   - Save this script in:
#     /Users/peter/Downloads/Project-2026-BurnOmicsDB/GSE8056/
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

# 【结果运行后注释】
# 纳入网站时需要明确显示：
# Pooled historical microarray study
# n = 3 pools per group
# 5 patients represented per pool
# Auxiliary / lower-confidence evidence
# 同时在后续生物学分析中，建议排除或单独标记明显受性别影响的Y染色体基因，不能把它们作为烧伤特异性发现。
# ==============================================================================


# ---- 00. Install required packages once / 首次安装所需软件包 ----

# This section checks all required packages and installs only missing ones.
# 本节检查全部依赖，只安装尚未安装的软件包。

# Increase the download timeout for large Bioconductor packages.
# 延长下载超时时间，避免大型Bioconductor软件包安装中断。
options(
  timeout = 1200
)

# Use libcurl for more reliable downloads.
# 使用libcurl提高下载稳定性。
options(
  download.file.method = "libcurl"
)

# Set the CRAN mirror.
# 设置CRAN镜像。
options(
  repos = c(
    CRAN = "https://cloud.r-project.org"
  )
)

# Use the Tsinghua University Bioconductor mirror.
# 使用清华大学Bioconductor镜像。
options(
  BioC_mirror =
    "https://mirrors.tuna.tsinghua.edu.cn/bioconductor"
)

# Install or update BiocManager when necessary.
# 如有必要，安装或更新BiocManager。
biocmanager_needs_installation <-
  !requireNamespace(
    "BiocManager",
    quietly = TRUE
  )

if (
  !biocmanager_needs_installation
) {
  biocmanager_needs_installation <-
    packageVersion(
      "BiocManager"
    ) <
    package_version(
      "1.30.12"
    )
}

if (
  biocmanager_needs_installation
) {
  install.packages(
    "BiocManager",
    repos = "https://cloud.r-project.org"
  )
}

# CRAN packages required by this analysis.
# 本分析所需的CRAN软件包。
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

if (
  length(
    missing_cran_packages
  ) > 0
) {
  install.packages(
    missing_cran_packages,
    repos = "https://cloud.r-project.org",
    dependencies = TRUE
  )
}

# Bioconductor packages required by this analysis.
# 本分析所需的Bioconductor软件包。
#
# hgu133plus2cdf:
#   Supplies the HG-U133_Plus_2 chip-layout environment required by affy::rma().
#
# hgu133plus2.db:
#   Supplies probe-set annotation used for gene-level mapping.
#
# hgu133plus2cdf提供RMA所需的芯片探针布局；
# hgu133plus2.db提供后续探针到基因的注释。
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

if (
  length(
    missing_bioconductor_packages
  ) > 0
) {
  BiocManager::install(
    missing_bioconductor_packages,
    ask = FALSE,
    update = FALSE
  )
}

# Verify that every required package is now available.
# 检查所有软件包是否均已成功安装。
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

if (
  !all(
    package_check
  )
) {
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

GEO_ID <- "GSE8056"

PROJECT_DIR <- "/Users/peter/Downloads/Project-2026-BurnOmicsDB/GSE8056"

FDR_CUTOFF <- 0.05
LOG2FC_CUTOFF <- 1
TOP_VARIABLE_GENES_FOR_PCA <- 500
TOP_GENES_FOR_HEATMAP <- 30
VOLCANO_LABEL_GENE_N <- 12

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
  library(affy)
  library(limma)
  library(AnnotationDbi)
  library(hgu133plus2cdf)
  library(hgu133plus2.db)
  library(org.Hs.eg.db)
})

set.seed(2026)

# Fixed semantic colors used across BurnOmicsDB.
# BurnOmicsDB各项目采用固定语义配色。
GROUP_COLORS <- c(
  "Normal_skin" = "#0072B2",
  "Early_wound" = "#D55E00",
  "Middle_wound" = "#009E73",
  "Late_wound" = "#CC79A7"
)

VOLCANO_COLORS <- c(
  "Up_significant" = "#D55E00",
  "Down_significant" = "#0072B2",
  "Not_significant" = "#BDBDBD"
)

EXPRESSION_HEATMAP_COLORS <- grDevices::colorRampPalette(
  c("#0072B2", "#F7F7F7", "#D55E00")
)(101)

GROUP_LEVELS <- c(
  "Normal_skin",
  "Early_wound",
  "Middle_wound",
  "Late_wound"
)

cat("Project directory:\n", PROJECT_DIR, "\n\n")


# ---- 02. Locate input files and extract CEL files / 定位输入文件并解压CEL ----

# Search both the project root and the optional Input folder.
# 同时搜索项目根目录和Input子文件夹。
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
  "^GSE8056_RAW.*\\.tar$",
  input_search_dirs
)

soft_file <- find_one_file(
  "^GSE8056_family\\.soft.*\\.gz$",
  input_search_dirs
)

paper_file <- find_one_file(
  ".*\\.pdf$",
  input_search_dirs,
  required = FALSE
)

cat("Raw CEL archive:\n", raw_tar_file, "\n\n")
cat("SOFT metadata:\n", soft_file, "\n\n")
cat("Paper PDF:\n", paper_file, "\n\n")

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

if (length(cel_gz_members) != 12) {
  stop(
    paste0(
      "The raw archive should contain 12 CEL.gz files, but ",
      length(cel_gz_members),
      " were found."
    )
  )
}

# Extract the compressed CEL files when necessary.
# 如尚未解包，则将CEL.gz文件提取到本地Input目录。
existing_cel_gz <- list.files(
  CEL_GZ_DIR,
  pattern = "\\.CEL\\.gz$",
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
)

if (length(existing_cel_gz) != 12) {
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

if (length(cel_gz_files) != 12) {
  stop(
    paste0(
      "Expected 12 extracted CEL.gz files, but ",
      length(cel_gz_files),
      " were found."
    )
  )
}

# Decompress CEL.gz files while retaining the original compressed files.
# 解压CEL.gz，同时保留原始压缩文件。
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

if (length(cel_files) != 12) {
  stop(
    paste0(
      "Expected 12 decompressed CEL files, but ",
      length(cel_files),
      " were found."
    )
  )
}

cel_sample_ids <- sub(
  "\\.CEL$",
  "",
  basename(cel_files),
  ignore.case = TRUE
)

if (anyDuplicated(cel_sample_ids) > 0) {
  stop(
    "Duplicated GSM identifiers were detected among the CEL filenames."
  )
}

cat("CEL files detected:\n")
print(
  data.frame(
    Sample_ID = cel_sample_ids,
    CEL_File = basename(cel_files),
    stringsAsFactors = FALSE
  )
)


# ---- 03. Parse GEO SOFT metadata and validate sample structure / 解析SOFT元数据并检查样本结构 ----

parse_geo_soft_samples <- function(
  soft_gz_file
) {
  connection <- gzfile(
    soft_gz_file,
    open = "rt"
  )

  on.exit(
    close(connection),
    add = TRUE
  )

  lines <- readLines(
    connection,
    warn = FALSE,
    encoding = "UTF-8"
  )

  sample_start <- grep(
    "^\\^SAMPLE =",
    lines
  )

  if (length(sample_start) == 0) {
    stop(
      "No SAMPLE blocks were found in the SOFT file."
    )
  }

  sample_end <- c(
    sample_start[-1] - 1,
    length(lines)
  )

  get_values <- function(
    block,
    field_name
  ) {
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
    key
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

    if (length(matches) == 0) {
      return(NA_character_)
    }

    sub(
      pattern,
      "",
      matches[1],
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

      sample_id <- sub(
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

      original_source <- first_or_na(
        get_values(
          block,
          "!Sample_source_name_ch1"
        )
      )

      characteristics <- get_values(
        block,
        "!Sample_characteristics_ch1"
      )

      data.frame(
        GEO_ID = GEO_ID,
        Sample_ID = sample_id,
        Sample_Name = sample_id,
        Original_Title = original_title,
        Original_Source_Name = original_source,
        Median_Age_Original = extract_characteristic(
          characteristics,
          "median age"
        ),
        Mean_Age_Original = extract_characteristic(
          characteristics,
          "mean age"
        ),
        Procedures_Original = extract_characteristic(
          characteristics,
          "Procedure"
        ),
        TBSA_Values_Original = extract_characteristic(
          characteristics,
          "Total Body Surface Are Burned"
        ),
        Post_Burn_Day_Values_Original = extract_characteristic(
          characteristics,
          "Post Burn Day"
        ),
        Sex_Values_Original = extract_characteristic(
          characteristics,
          "Sex"
        ),
        Platform_ID = first_or_na(
          get_values(
            block,
            "!Sample_platform_id"
          )
        ),
        Original_Characteristics = paste(
          characteristics,
          collapse = " | "
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

sample_metadata <- parse_geo_soft_samples(
  soft_file
)

if (nrow(sample_metadata) != 12) {
  stop(
    paste0(
      "The SOFT file should contain 12 samples, but ",
      nrow(sample_metadata),
      " were parsed."
    )
  )
}

# Assign groups using directly reported sample titles.
# 根据SOFT中直接报告的样本标题分组。
sample_metadata$Group <- ifelse(
  grepl(
    "normal human",
    sample_metadata$Original_Title,
    ignore.case = TRUE
  ),
  "Normal_skin",
  ifelse(
    grepl(
      "0-3 days",
      sample_metadata$Original_Title,
      fixed = TRUE
    ),
    "Early_wound",
    ifelse(
      grepl(
        "4-7 days",
        sample_metadata$Original_Title,
        fixed = TRUE
      ),
      "Middle_wound",
      ifelse(
        grepl(
          ">7 days",
          sample_metadata$Original_Title,
          fixed = TRUE
        ),
        "Late_wound",
        NA_character_
      )
    )
  )
)

if (any(is.na(sample_metadata$Group))) {
  stop(
    paste0(
      "Unexpected sample titles prevented group assignment: ",
      paste(
        sample_metadata$Original_Title[
          is.na(sample_metadata$Group)
        ],
        collapse = " | "
      )
    )
  )
}

sample_metadata$Group <- factor(
  sample_metadata$Group,
  levels = GROUP_LEVELS
)

sample_metadata$Pool_Replicate <- suppressWarnings(
  as.integer(
    sub(
      ".*-",
      "",
      sample_metadata$Original_Title
    )
  )
)

if (any(is.na(sample_metadata$Pool_Replicate))) {
  stop(
    "At least one pool replicate number could not be parsed."
  )
}

pool_prefix <- c(
  "Normal_skin" = "NP",
  "Early_wound" = "EP",
  "Middle_wound" = "MP",
  "Late_wound" = "LP"
)

sample_metadata$Pool_ID <- paste0(
  unname(
    pool_prefix[
      as.character(
        sample_metadata$Group
      )
    ]
  ),
  sample_metadata$Pool_Replicate
)

sample_metadata$Patient_ID <- NA_character_
sample_metadata$Tissue <- ifelse(
  sample_metadata$Group == "Normal_skin",
  "Uninjured human skin",
  "Partial-thickness burn wound margin"
)

sample_metadata$Time_or_Stage <- ifelse(
  sample_metadata$Group == "Normal_skin",
  "Uninjured skin",
  ifelse(
    sample_metadata$Group == "Early_wound",
    "1-3 days after injury",
    ifelse(
      sample_metadata$Group == "Middle_wound",
      "4-7 days after injury",
      "8-17 days after injury"
    )
  )
)

sample_metadata$Treatment <- NA_character_
sample_metadata$Outcome <- NA_character_
sample_metadata$Data_Type <- "Affymetrix raw CEL"
sample_metadata$Is_Paired <- FALSE
sample_metadata$Is_Pooled <- TRUE
sample_metadata$Is_Repeated_Measure <- FALSE
sample_metadata$Biological_Specimens_Per_Assay <- 5
sample_metadata$Patient_N_Represented <- 5
sample_metadata$Metadata_Confidence <- "Direct_from_GEO_SOFT_and_paper"
sample_metadata$BioSample_ID <- NA_character_
sample_metadata$SRA_Experiment <- NA_character_

base_quality_note <- paste0(
  "Each array contains an equal-mass RNA pool from five patients; ",
  "individual-patient expression variability and patient-level covariate ",
  "adjustment are unavailable; normal and burn groups differ in sex ",
  "composition; normal skin was collected during elective cosmetic surgery; ",
  "this dataset is classified as auxiliary lower-confidence evidence."
)

sample_metadata$Quality_Notes <- base_quality_note

# The SOFT metadata for GSM198877 appears inconsistent with the NP1 patient
# composition shown in Table 1 of the paper. This does not change the control
# group assignment, but it prevents demographic adjustment.
# GSM198877的SOFT人口学描述与论文Table 1中的NP1组成疑似不一致。
sample_metadata$Quality_Notes[
  sample_metadata$Sample_ID == "GSM198877"
] <- paste0(
  base_quality_note,
  " The GSM198877 demographic summary in SOFT appears inconsistent with ",
  "the NP1 composition reported in Table 1 of the associated paper; group ",
  "assignment remains Normal_skin, but demographic fields are not used ",
  "for modeling."
)

# Verify that CEL files and SOFT records contain exactly the same GSM IDs.
# 检查CEL文件与SOFT样本是否完全对应。
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
      "CEL and SOFT sample identifiers do not match. ",
      "Missing from SOFT: ",
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
    "Sample metadata could not be reordered to match CEL files."
  )
}

actual_group_counts <- table(
  sample_metadata$Group
)

if (!all(
  actual_group_counts[
    GROUP_LEVELS
  ] == 3
)) {
  stop(
    paste0(
      "Each group should contain three pooled arrays. Observed counts: ",
      paste(
        names(actual_group_counts),
        actual_group_counts,
        sep = "=",
        collapse = ", "
      )
    )
  )
}

metadata_audit <- data.frame(
  Check = c(
    "Raw CEL file count",
    "SOFT sample count",
    "Normal-skin pool count",
    "Early-wound pool count",
    "Middle-wound pool count",
    "Late-wound pool count",
    "Patients represented per pool",
    "Total patients represented",
    "Statistical replicates per group",
    "Patient-level covariate adjustment",
    "Paired or repeated samples",
    "SOFT demographic inconsistency flag"
  ),
  Value = c(
    length(cel_files),
    nrow(sample_metadata),
    actual_group_counts["Normal_skin"],
    actual_group_counts["Early_wound"],
    actual_group_counts["Middle_wound"],
    actual_group_counts["Late_wound"],
    5,
    60,
    3,
    "No",
    "No",
    "GSM198877"
  ),
  Interpretation = c(
    "Directly observed in GSE8056_RAW.tar.",
    "Directly parsed from the GSE family SOFT file.",
    "Three independent RNA pools; five patients per pool.",
    "Three independent RNA pools; five patients per pool.",
    "Three independent RNA pools; five patients per pool.",
    "Three independent RNA pools; five patients per pool.",
    "Equal-mass RNA from five tissue specimens was combined per array.",
    "Fifteen normal-skin and forty-five burn patients are represented.",
    "The pool, not the individual patient, is the statistical unit.",
    "Unavailable because patient-level expression profiles were pooled.",
    "Each patient contributes to only one pool in the study design.",
    "The SOFT control-pool summary appears inconsistent with Table 1 of the paper."
  ),
  stringsAsFactors = FALSE
)

write.csv(
  metadata_audit,
  file = file.path(
    QC_DIR,
    "GSE8056_metadata_consistency_check.csv"
  ),
  row.names = FALSE
)

cat("Sample-group counts:\n")
print(actual_group_counts)

cat("\nMetadata audit:\n")
print(metadata_audit)


# ---- 04. Read raw Affymetrix CEL files / 读取Affymetrix原始CEL文件 ----

# ReadAffy imports the probe-level intensity values from the 12 CEL files.
# ReadAffy读取12个CEL文件中的探针级原始强度。
raw_affy <- affy::ReadAffy(
  filenames = cel_files
)

if (length(
  sampleNames(raw_affy)
) != 12) {
  stop(
    "ReadAffy did not return 12 arrays."
  )
}

sampleNames(raw_affy) <- cel_sample_ids

# Reorder the AffyBatch and metadata to the same GSM order.
# 将AffyBatch和元数据统一为同一GSM顺序。
raw_affy <- raw_affy[
  ,
  sample_metadata$Sample_ID
]

if (!identical(
  sampleNames(raw_affy),
  sample_metadata$Sample_ID
)) {
  stop(
    "AffyBatch sample order does not match the metadata."
  )
}

raw_intensity_matrix <- exprs(
  raw_affy
)

if (any(is.na(raw_intensity_matrix))) {
  stop(
    "Raw CEL intensities contain missing values."
  )
}

if (any(raw_intensity_matrix < 0)) {
  stop(
    "Raw CEL intensities contain negative values."
  )
}

raw_log2_intensity <- log2(
  raw_intensity_matrix + 1
)

raw_qc_metrics <- data.frame(
  Sample_ID = sample_metadata$Sample_ID,
  Pool_ID = sample_metadata$Pool_ID,
  Group = sample_metadata$Group,
  Raw_Median_Log2_Intensity = apply(
    raw_log2_intensity,
    2,
    median
  ),
  Raw_IQR_Log2_Intensity = apply(
    raw_log2_intensity,
    2,
    IQR
  ),
  Raw_Min_Log2_Intensity = apply(
    raw_log2_intensity,
    2,
    min
  ),
  Raw_Max_Log2_Intensity = apply(
    raw_log2_intensity,
    2,
    max
  ),
  stringsAsFactors = FALSE
)

write.csv(
  raw_qc_metrics,
  file = file.path(
    QC_DIR,
    "GSE8056_raw_array_QC_metrics.csv"
  ),
  row.names = FALSE
)

# Main raw-data QC figure.
# 主要原始数据QC图。
png(
  filename = file.path(
    FIGURES_DIR,
    "01_GSE8056_raw_data_QC.png"
  ),
  width = 2200,
  height = 1500,
  res = 240
)

par(
  mar = c(
    8,
    5,
    4,
    2
  )
)

boxplot(
  raw_log2_intensity,
  names = sample_metadata$Pool_ID,
  col = GROUP_COLORS[
    as.character(
      sample_metadata$Group
    )
  ],
  las = 2,
  outline = FALSE,
  main = "GSE8056 raw CEL intensity distribution",
  sub = "Each array is a pool of five independent human skin specimens",
  xlab = "",
  ylab = "Raw log2 probe intensity",
  cex.axis = 0.9
)

legend(
  "topright",
  legend = GROUP_LEVELS,
  fill = GROUP_COLORS[
    GROUP_LEVELS
  ],
  border = NA,
  bty = "n",
  title = "Group"
)

dev.off()


# ---- 05. RMA normalization and platform-specific QC / RMA标准化与平台质控 ----

# RMA performs background correction, quantile normalization, and median-polish
# summarization at the Affymetrix probe-set level.
# RMA执行背景校正、分位数标准化和探针集median-polish汇总。
rma_eset <- affy::rma(
  raw_affy
)

rma_probe_expression <- exprs(
  rma_eset
)

if (ncol(rma_probe_expression) != 12) {
  stop(
    "The RMA expression matrix does not contain 12 arrays."
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

# Save the probe-level matrix for reproducibility.
# 保存探针集层面的RMA矩阵用于复现。
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
    "GSE8056_RMA_probe_level_expression.csv.gz"
  ),
  open = "wt"
)

write.csv(
  rma_probe_output,
  file = rma_probe_connection,
  row.names = FALSE
)

close(rma_probe_connection)

normalized_qc_metrics <- data.frame(
  Sample_ID = sample_metadata$Sample_ID,
  Pool_ID = sample_metadata$Pool_ID,
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

write.csv(
  normalized_qc_metrics,
  file = file.path(
    QC_DIR,
    "GSE8056_RMA_array_QC_metrics.csv"
  ),
  row.names = FALSE
)

png(
  filename = file.path(
    FIGURES_DIR,
    "02_GSE8056_normalized_expression_distribution.png"
  ),
  width = 2200,
  height = 1500,
  res = 240
)

par(
  mar = c(
    8,
    5,
    4,
    2
  )
)

boxplot(
  rma_probe_expression,
  names = sample_metadata$Pool_ID,
  col = GROUP_COLORS[
    as.character(
      sample_metadata$Group
    )
  ],
  las = 2,
  outline = FALSE,
  main = "GSE8056 RMA-normalized expression distribution",
  sub = "RMA log2 expression at the Affymetrix probe-set level",
  xlab = "",
  ylab = "RMA log2 expression",
  cex.axis = 0.9
)

legend(
  "topright",
  legend = GROUP_LEVELS,
  fill = GROUP_COLORS[
    GROUP_LEVELS
  ],
  border = NA,
  bty = "n",
  title = "Group"
)

dev.off()


# ---- 06. Probe annotation and gene-level representative selection / 探针注释与基因级代表探针选择 ----

probe_ids <- rownames(
  rma_probe_expression
)

probe_annotation_raw <- AnnotationDbi::select(
  hgu133plus2.db,
  keys = probe_ids,
  columns = c(
    "ENTREZID",
    "SYMBOL",
    "GENENAME"
  ),
  keytype = "PROBEID"
)

probe_annotation_raw$PROBEID <- as.character(
  probe_annotation_raw$PROBEID
)

probe_annotation_raw$ENTREZID <- as.character(
  probe_annotation_raw$ENTREZID
)

# Count the number of unique Entrez Gene IDs assigned to each probe.
# 统计每个探针对应的唯一Entrez Gene ID数量。
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

# Reduce duplicated annotation rows so each unambiguous probe has one Entrez ID.
# 去除重复注释行，使每个明确探针对应一个Entrez ID。
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

unambiguous_probe_annotation$Mean_RMA_Expression <- probe_mean_expression[
  unambiguous_probe_annotation$PROBEID
]

# Select the highest-average-expression probe for each Entrez Gene ID.
# 对每个Entrez Gene ID选择全体12张芯片平均表达最高的探针。
unambiguous_probe_annotation <- unambiguous_probe_annotation[
  order(
    unambiguous_probe_annotation$ENTREZID,
    -unambiguous_probe_annotation$Mean_RMA_Expression,
    unambiguous_probe_annotation$PROBEID
  ),
]

representative_probe_table <- unambiguous_probe_annotation[
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

# Map current official gene symbols, Ensembl IDs, and names using org.Hs.eg.db.
# 使用org.Hs.eg.db映射当前官方symbol、Ensembl ID和基因名称。
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

# Build a full probe-mapping audit table.
# 建立完整探针映射审计表。
probe_mapping_audit <- merge(
  probe_mapping_status,
  unambiguous_probe_annotation,
  by.x = "Probe_ID",
  by.y = "PROBEID",
  all.x = TRUE,
  sort = FALSE
)

probe_mapping_audit$Selected_as_Representative <- FALSE

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
    "GSE8056_probe_mapping_and_selection.csv.gz"
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
    "GSE8056_probe_mapping_summary.csv"
  ),
  row.names = FALSE
)

cat("Probe-mapping summary:\n")
print(mapping_summary)

# Save the gene-level normalized expression matrix.
# 保存基因级RMA标准化表达矩阵。
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
    "GSE8056_normalized_expression.csv.gz"
  ),
  open = "wt"
)

write.csv(
  normalized_gene_output,
  file = normalized_gene_connection,
  row.names = FALSE
)

close(normalized_gene_connection)


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

pca_table$Pool_ID <- sample_metadata$Pool_ID[
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
    "GSE8056_PCA_scores.csv"
  ),
  row.names = FALSE
)

p_pca <- ggplot(
  pca_table,
  aes(
    x = PC1,
    y = PC2,
    color = Group,
    label = Pool_ID
  )
) +
  geom_point(
    size = 3.2
  ) +
  ggrepel::geom_text_repel(
    max.overlaps = Inf,
    size = 3.5
  ) +
  scale_color_manual(
    values = GROUP_COLORS,
    drop = FALSE
  ) +
  labs(
    title = "GSE8056 PCA",
    subtitle = paste0(
      "Top ",
      n_pca_genes,
      " variable gene-level RMA features"
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

print(p_pca)

ggsave(
  filename = file.path(
    FIGURES_DIR,
    "03_GSE8056_PCA.png"
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
  labels_col = sample_metadata$Pool_ID,
  labels_row = sample_metadata$Pool_ID,
  border_color = NA,
  color = EXPRESSION_HEATMAP_COLORS,
  main = "GSE8056 sample correlation"
)

pheatmap::pheatmap(
  sample_correlation,
  annotation_col = correlation_annotation,
  annotation_row = correlation_annotation,
  annotation_colors = correlation_annotation_colors,
  labels_col = sample_metadata$Pool_ID,
  labels_row = sample_metadata$Pool_ID,
  border_color = NA,
  color = EXPRESSION_HEATMAP_COLORS,
  main = "GSE8056 sample correlation",
  filename = file.path(
    FIGURES_DIR,
    "04_GSE8056_sample_correlation_heatmap.png"
  ),
  width = 8,
  height = 7
)


# ---- 08. Differential-expression model / 差异表达模型 ----

sample_metadata$Group <- factor(
  sample_metadata$Group,
  levels = GROUP_LEVELS
)

design <- model.matrix(
  ~ 0 + Group,
  data = sample_metadata
)

colnames(design) <- sub(
  "^Group",
  "",
  colnames(design)
)

rownames(design) <- sample_metadata$Sample_ID

if (qr(design)$rank < ncol(design)) {
  stop(
    "The limma design matrix is not full rank."
  )
}

cat("Design-matrix columns:\n")
print(
  colnames(design)
)

fit <- limma::lmFit(
  gene_expression,
  design = design
)

contrast_matrix <- limma::makeContrasts(
  GSE8056_EarlyWound_vs_NormalSkin =
    Early_wound - Normal_skin,
  GSE8056_MiddleWound_vs_NormalSkin =
    Middle_wound - Normal_skin,
  GSE8056_LateWound_vs_NormalSkin =
    Late_wound - Normal_skin,
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
    "Early wound vs normal skin",
    "Middle wound vs normal skin",
    "Late wound vs normal skin"
  ),
  Case_Group = c(
    "Early wound",
    "Middle wound",
    "Late wound"
  ),
  Case_Group_Code = c(
    "Early_wound",
    "Middle_wound",
    "Late_wound"
  ),
  Control_Group = "Normal skin",
  Control_Group_Code = "Normal_skin",
  Sample_Context = "Pooled partial-thickness burn wound margin",
  Time_or_Stage = c(
    "1-3 days after injury",
    "4-7 days after injury",
    "8-17 days after injury"
  ),
  Case_N = 3,
  Control_N = 3,
  Case_Patient_N = 15,
  Control_Patient_N = 15,
  stringsAsFactors = FALSE
)

cat("Contrast definitions:\n")
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

  contrast_result$NCBI_Gene_ID <- rownames(
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

row.names(
  all_gene_results
) <- NULL

all_results_connection <- gzfile(
  file.path(
    RESULTS_DIR,
    "GSE8056_all_gene_results.csv.gz"
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
      alpha = 0.65,
      size = 1.15
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
        "GSE8056: ",
        contrast_info$Contrast_Label
      ),
      subtitle = paste0(
        "RMA + limma pooled-array model; ",
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
        "05_GSE8056_volcano_",
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

  heatmap_column_labels <- sample_metadata$Pool_ID[
    match(
      contrast_sample_ids,
      sample_metadata$Sample_ID
    )
  ]

  pheatmap::pheatmap(
    heatmap_z,
    annotation_col = heatmap_annotation,
    annotation_colors = heatmap_annotation_colors,
    labels_col = heatmap_column_labels,
    show_colnames = TRUE,
    show_rownames = TRUE,
    border_color = NA,
    color = EXPRESSION_HEATMAP_COLORS,
    main = paste0(
      "Top ",
      n_heatmap_genes,
      " genes: ",
      contrast_info$Contrast_Label
    )
  )

  pheatmap::pheatmap(
    heatmap_z,
    annotation_col = heatmap_annotation,
    annotation_colors = heatmap_annotation_colors,
    labels_col = heatmap_column_labels,
    show_colnames = TRUE,
    show_rownames = TRUE,
    border_color = NA,
    color = EXPRESSION_HEATMAP_COLORS,
    main = paste0(
      "Top ",
      n_heatmap_genes,
      " genes: ",
      contrast_info$Contrast_Label
    ),
    filename = file.path(
      FIGURES_DIR,
      paste0(
        "06_GSE8056_top_differential_genes_heatmap_",
        contrast_id,
        ".png"
      )
    ),
    width = 8,
    height = 10
  )
}


# ---- 12. Cross-check selected genes reported in the paper / 核对论文报告的代表基因 ----

# Published microarray fold changes from Table 3 of the associated paper.
# 以下数值来自论文Table 3的microarray fold change。
published_gene_fc <- data.frame(
  Gene_Symbol = rep(
    c(
      "TNFRSF10B",
      "THBS1",
      "IL6",
      "CXCL8",
      "SPP1"
    ),
    each = 3
  ),
  Contrast_ID = rep(
    contrast_definitions$Contrast_ID,
    times = 5
  ),
  Published_Fold_Change = c(
    2.27,
    2.20,
    1.99,
    4.80,
    8.90,
    6.43,
    19.39,
    24.56,
    12.80,
    121.72,
    146.85,
    66.31,
    22.06,
    78.38,
    72.05
  ),
  stringsAsFactors = FALSE
)

calculated_crosscheck <- all_gene_results[
  all_gene_results$Gene_Symbol %in%
    unique(
      published_gene_fc$Gene_Symbol
    ),
  c(
    "Gene_Symbol",
    "NCBI_Gene_ID",
    "Representative_Probe_ID",
    "Contrast_ID",
    "log2FC",
    "Fold_Change",
    "P_value",
    "FDR",
    "Direction",
    "DE_Status"
  )
]

published_crosscheck <- merge(
  published_gene_fc,
  calculated_crosscheck,
  by = c(
    "Gene_Symbol",
    "Contrast_ID"
  ),
  all.x = TRUE,
  sort = FALSE
)

published_crosscheck$Direction_Consistent <- sign(
  log2(
    published_crosscheck$Published_Fold_Change
  )
) == sign(
  published_crosscheck$log2FC
)

published_crosscheck$Calculated_to_Published_FC_Ratio <-
  published_crosscheck$Fold_Change /
    published_crosscheck$Published_Fold_Change

write.csv(
  published_crosscheck,
  file = file.path(
    RESULTS_DIR,
    "GSE8056_published_gene_crosscheck.csv"
  ),
  row.names = FALSE
)

cat("Published-gene cross-check:\n")
print(
  published_crosscheck
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
    "Mixed-age pooled human skin specimens",
  Tissue = "Skin",
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
    "Raw CEL files from GSE8056_RAW.tar",
  Normalization = paste0(
    "RMA: background correction, quantile normalization, ",
    "and median-polish probe-set summarization"
  ),
  Analysis_Method = paste0(
    "limma linear model with empirical Bayes moderation; one ",
    "representative probe per NCBI Gene ID selected using the highest ",
    "mean RMA expression across all 12 arrays; no patient-level ",
    "covariate adjustment; positive log2FC represents Case_Group minus ",
    "Control_Group"
  ),
  Annotation_Method = paste0(
    "hgu133plus2.db ",
    as.character(
      packageVersion(
        "hgu133plus2.db"
      )
    ),
    " used for probe-to-Entrez mapping; probes mapping to zero or ",
    "multiple Entrez Gene IDs excluded; org.Hs.eg.db ",
    as.character(
      packageVersion(
        "org.Hs.eg.db"
      )
    ),
    " used for current official symbol, Ensembl ID, and gene name"
  ),
  Quality_Notes = paste0(
    "Auxiliary lower-confidence evidence; each statistical sample is an ",
    "RNA pool from five patients; n=3 pools per group; individual-patient ",
    "variability and covariate adjustment unavailable; normal and burn ",
    "groups differ in sex composition; normal skin was obtained during ",
    "elective cosmetic surgery; GSM198877 SOFT demographic summary appears ",
    "inconsistent with Table 1 of the associated paper; bulk-tissue ",
    "differences may reflect intracellular regulation, cell-composition ",
    "changes, and tissue-structure differences"
  ),
  Mean_Normalized_Expression =
    all_gene_results$Mean_Normalized_Expression,
  Expression_Scale = "RMA log2 expression",
  Representative_Probe_ID =
    all_gene_results$Representative_Probe_ID,
  Case_Patient_N =
    all_gene_results$Case_Patient_N,
  Control_Patient_N =
    all_gene_results$Control_Patient_N,
  Is_Pooled = TRUE,
  Biological_Specimens_Per_Assay = 5,
  Quality_Grade =
    "Auxiliary_lower_confidence",
  stringsAsFactors = FALSE
)

database_connection <- gzfile(
  file.path(
    RESULTS_DIR,
    "GSE8056_database_ready_all_genes.csv.gz"
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
  nrow(
    database_ready
  ),
  "\n",
  sep = ""
)


# ---- 14. Save metadata, analysis objects, summary, and session information / 保存元数据、对象、摘要与环境信息 ----

# Place the required core metadata columns first.
# 将要求的核心元数据字段置于前面。
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
    "Pool_ID",
    "Pool_Replicate",
    "Biological_Specimens_Per_Assay",
    "Patient_N_Represented",
    "Median_Age_Original",
    "Mean_Age_Original",
    "Procedures_Original",
    "TBSA_Values_Original",
    "Post_Burn_Day_Values_Original",
    "Sex_Values_Original"
  )
]

write.csv(
  sample_metadata_output,
  file = file.path(
    RESULTS_DIR,
    "GSE8056_sample_metadata.csv"
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
  raw_affy = raw_affy,
  raw_qc_metrics = raw_qc_metrics,
  rma_eset = rma_eset,
  rma_probe_expression = rma_probe_expression,
  normalized_qc_metrics = normalized_qc_metrics,
  probe_mapping_summary = mapping_summary,
  probe_mapping_audit = probe_mapping_audit,
  representative_probe_table = representative_probe_table,
  gene_annotation = gene_annotation,
  gene_expression = gene_expression,
  design_matrix = design,
  contrast_matrix = contrast_matrix,
  contrast_definitions = contrast_definitions,
  fitted_model = fit_contrasts,
  all_gene_results = all_gene_results,
  published_crosscheck = published_crosscheck,
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
    "GSE8056_analysis_objects.rds"
  ),
  compress = "xz"
)

summary_lines <- c(
  "BurnOmicsDB - GSE8056 analysis summary",
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
  "Patients represented: 60",
  "Patients represented per pool: 5",
  "Statistical unit: pooled array",
  "Statistical replicates per group: 3",
  "",
  "Groups:",
  "  Normal skin: 3 pools representing 15 patients",
  "  Early wound: 3 pools representing 15 patients",
  "  Middle wound: 3 pools representing 15 patients",
  "  Late wound: 3 pools representing 15 patients",
  "",
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
    "12 arrays for each NCBI Gene ID"
  ),
  "",
  paste0(
    "Normalization: RMA background correction, quantile normalization, ",
    "and median-polish summarization"
  ),
  paste0(
    "Statistical model: limma with robust, trend-aware empirical Bayes ",
    "moderation; no patient-level covariate adjustment"
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
    " (case pools = ",
    contrast_definitions$Case_N,
    ", case patients represented = ",
    contrast_definitions$Case_Patient_N,
    ", control pools = ",
    contrast_definitions$Control_N,
    ", control patients represented = ",
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
    "Positive log2FC means higher expression in the burn-wound case group ",
    "than in normal skin."
  ),
  paste0(
    "No separate low-expression filter was applied after RMA and ",
    "representative-probe selection; all selected unambiguous gene ",
    "representatives were statistically tested and exported."
  ),
  "",
  "Major limitations:",
  paste0(
    "  Each array is an RNA pool from five patients; patient-level ",
    "expression variability is unavailable."
  ),
  paste0(
    "  Each group contains only three statistical replicates even though ",
    "fifteen patients are represented."
  ),
  paste0(
    "  Patient-level age, sex, TBSA, anatomical location, and post-burn ",
    "day cannot be adjusted in the pooled-array model."
  ),
  paste0(
    "  Normal and burn groups differ substantially in sex composition."
  ),
  paste0(
    "  Normal skin was collected during elective cosmetic surgery, whereas ",
    "burn samples were partial-thickness wound-margin tissue."
  ),
  paste0(
    "  The GSM198877 SOFT demographic summary appears inconsistent with ",
    "the NP1 composition reported in Table 1 of the paper."
  ),
  paste0(
    "  Bulk-tissue expression differences may reflect intracellular ",
    "regulation, cell-composition changes, and tissue-structure differences."
  ),
  "Quality grade: Auxiliary_lower_confidence"
)

writeLines(
  summary_lines,
  con = file.path(
    RESULTS_DIR,
    "GSE8056_analysis_summary.txt"
  ),
  useBytes = TRUE
)

capture.output(
  sessionInfo(),
  file = file.path(
    RESULTS_DIR,
    "GSE8056_R_sessionInfo.txt"
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
  "TNFRSF10B",
  "THBS1",
  "IL6",
  "CXCL8",
  "SPP1",
  "MMP1",
  "MMP3",
  "MMP9",
  "COL4A1",
  "COL4A2",
  "KRT6A",
  "KRT16",
  "KRT17"
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
