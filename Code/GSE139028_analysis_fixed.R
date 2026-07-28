# ==============================================================================
# BurnOmicsDB: GSE139028 RNA-seq analysis
# BurnOmicsDB：GSE139028 RNA-seq分析
#
# GEO accession / GEO编号:
#   GSE139028
#
# Study / 研究:
#   RNA-seq analysis of pediatric burn-eschar tissue and normal skin.
#   儿童烧伤焦痂组织与正常皮肤的RNA-seq分析。
#
# Data type / 数据类型:
#   Author-provided raw RNA-seq counts summarized at NCBI Entrez Gene level.
#   作者提供的NCBI Entrez Gene层面RNA-seq原始计数矩阵。
#
# Study structure / 研究结构:
#   - 9 RNA-seq expression samples.
#   - 3 normal-skin samples.
#   - 6 burn-eschar samples.
#   - The associated publication identifies 3 normal-skin tissue donors and
#     6 burn-eschar tissue donors for RNA-seq.
#   - The GEO SOFT metadata does not link AE1-AE9 to the individual
#     demographic rows in the publication; sample-level age, sex, TBSA,
#     burn type, wound depth, anatomical site, and exact post-burn day are
#     therefore not assigned.
#
#   共9个RNA-seq样本：3个正常皮肤和6个烧伤焦痂。
#   论文可确认RNA-seq涉及3名正常皮肤供体和6名烧伤组织供体，
#   但AE1-AE9与论文人口学表之间没有明确逐一对应，因此不推测样本级临床信息。
#
# Prespecified contrast / 预设比较:
#   GSE139028_BurnEschar_vs_NormalSkin
#   Burn eschar versus normal skin
#
# Contrast direction / 比较方向:
#   Positive log2FC always means higher expression in burn eschar than in
#   normal skin.
#   正log2FC始终表示烧伤焦痂中的表达高于正常皮肤。
#
# Analysis workflow / 分析流程:
#   Author-provided raw counts
#   -> strict count-matrix and sample-order validation
#   -> edgeR filterByExpr low-expression filtering
#   -> TMM normalization
#   -> normalized log2 CPM for PCA, correlation, and heatmaps
#   -> edgeR robust quasi-likelihood negative-binomial GLM
#   -> Entrez Gene ID annotation with org.Hs.eg.db
#   -> one standardized BurnOmicsDB contrast
#
# Statistical unit and design / 统计单位与设计:
#   Each expression column is treated as one independent biological tissue
#   sample. The comparison is unpaired because GEO does not provide a
#   patient-linkage or matched-control structure.
#   每列表达数据作为一个独立组织样本；由于没有患者配对信息，采用非配对设计。
#
# Filtering policy / 过滤策略:
#   Only genes retained by edgeR::filterByExpr() are statistically tested and
#   exported. Low-expression genes removed before testing are not added back
#   to the result tables. The website can display:
#   "This gene was not found in the differential-expression result list."
#   仅输出通过filterByExpr并完成统计检验的基因，不补充低表达过滤基因。
#
# Relation to the publication / 与论文的关系:
#   The paper reports an RNA-seq significance criterion of FDR < 0.01 and
#   emphasizes protease genes. BurnOmicsDB uses its unified database threshold:
#   FDR < 0.05 and absolute log2FC >= 1.
#   A separate cross-check file compares published protease fold changes with
#   the current uniform reanalysis; it is QC only and is not a website input.
#   论文重点报告蛋白酶并使用更严格的FDR标准；本数据库使用统一阈值，
#   论文结果核对文件仅用于QC。
#
# Bulk-tissue interpretation / Bulk组织解释:
#   Differences can reflect intracellular regulation, leukocyte infiltration,
#   cell-composition changes, and major tissue-structure differences between
#   burn eschar and elective-surgery normal skin.
#   表达差异可能同时反映细胞内调控、炎症细胞浸润、细胞组成及组织结构差异。
#
# How to run in RStudio / 如何在RStudio中运行:
#   - Save this script in:
#     /Users/peter/Downloads/Project-2026-BurnOmicsDB/GSE139028/
#   - Keep the count matrix, SOFT file, and paper PDF in the project root or
#     in an optional Input folder.
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
  "readxl",
  "ggplot2",
  "ggrepel",
  "pheatmap",
  "matrixStats"
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
  "edgeR",
  "limma",
  "AnnotationDbi",
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

GEO_ID <- "GSE139028"

PROJECT_DIR <-
  "/Users/peter/Downloads/Project-2026-BurnOmicsDB/GSE139028"

FDR_CUTOFF <- 0.05
LOG2FC_CUTOFF <- 1
PAPER_FDR_CUTOFF <- 0.01
TOP_VARIABLE_GENES_FOR_PCA <- 500
TOP_GENES_FOR_HEATMAP <- 30
VOLCANO_LABEL_GENE_N <- 12

CONTRAST_ID <-
  "GSE139028_BurnEschar_vs_NormalSkin"

CONTRAST_LABEL <-
  "Burn eschar vs normal skin"

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

required_packages <- c(
  "readxl",
  "ggplot2",
  "ggrepel",
  "pheatmap",
  "matrixStats",
  "edgeR",
  "limma",
  "AnnotationDbi",
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
  library(readxl)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(matrixStats)
  library(edgeR)
  library(limma)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

set.seed(2026)

GROUP_LEVELS <- c(
  "Normal_skin",
  "Burn_eschar"
)

GROUP_COLORS <- c(
  "Normal_skin" = "#0072B2",
  "Burn_eschar" = "#D55E00"
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

cat(
  "Project directory:\n",
  PROJECT_DIR,
  "\n\n"
)


# ---- 02. Locate input files / 定位输入文件 ----

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

counts_file <- find_one_file(
  "^GSE139028_Counts_Table_for_RNA-Seq_Data.*\\.xlsx$",
  input_search_dirs
)

soft_file <- find_one_file(
  "^GSE139028_family\\.soft.*\\.gz$",
  input_search_dirs
)

paper_file <- find_one_file(
  ".*139028.*\\.pdf$|.*reference.*\\.pdf$",
  input_search_dirs,
  required = FALSE
)

cat("Count matrix:\n", counts_file, "\n\n")
cat("SOFT metadata:\n", soft_file, "\n\n")
cat("Paper PDF:\n", paper_file, "\n\n")


# ---- 03. Import and validate the raw count matrix / 导入并检查原始计数矩阵 ----

sheet_names <- readxl::excel_sheets(
  counts_file
)

if (length(sheet_names) != 1) {
  stop(
    paste0(
      "The count workbook should contain one worksheet, but ",
      length(sheet_names),
      " were found."
    )
  )
}

counts_table <- readxl::read_excel(
  path = counts_file,
  sheet = sheet_names[1]
)

if (ncol(counts_table) != 10) {
  stop(
    paste0(
      "The count table should contain one Gene ID column and nine sample ",
      "columns, but ",
      ncol(counts_table),
      " columns were found."
    )
  )
}

colnames(counts_table)[1] <-
  "NCBI_Gene_ID"

expected_sample_names <- paste0(
  "AE",
  1:9
)

observed_sample_names <- colnames(
  counts_table
)[-1]

if (!identical(
  observed_sample_names,
  expected_sample_names
)) {
  stop(
    paste0(
      "Unexpected count-matrix sample columns. Observed: ",
      paste(
        observed_sample_names,
        collapse = ", "
      ),
      ". Expected: ",
      paste(
        expected_sample_names,
        collapse = ", "
      ),
      "."
    )
  )
}

if (any(
  is.na(
    counts_table$NCBI_Gene_ID
  )
)) {
  stop(
    "The count matrix contains missing NCBI Gene IDs."
  )
}

if (anyDuplicated(
  counts_table$NCBI_Gene_ID
) > 0) {
  stop(
    "Duplicated NCBI Gene IDs were detected in the count matrix."
  )
}

count_matrix <- as.matrix(
  counts_table[
    ,
    -1
  ]
)

rownames(count_matrix) <- as.character(
  counts_table$NCBI_Gene_ID
)

if (any(is.na(count_matrix))) {
  stop(
    "The count matrix contains missing expression values."
  )
}

if (any(count_matrix < 0)) {
  stop(
    "The count matrix contains negative values."
  )
}

if (any(
  abs(
    count_matrix -
      round(count_matrix)
  ) > 1e-8
)) {
  stop(
    "The count matrix contains non-integer values."
  )
}

if (anyDuplicated(
  colnames(count_matrix)
) > 0) {
  stop(
    "Duplicated sample names were detected in the count matrix."
  )
}

storage.mode(count_matrix) <-
  "integer"

raw_library_sizes <- colSums(
  count_matrix
)

zero_count_fraction <- colMeans(
  count_matrix == 0
)

input_qc <- list(
  input_gene_n =
    nrow(count_matrix),
  input_sample_n =
    ncol(count_matrix),
  unique_gene_id_n =
    length(
      unique(
        rownames(count_matrix)
      )
    ),
  duplicated_gene_id_n =
    anyDuplicated(
      rownames(count_matrix)
    ),
  all_zero_gene_n =
    sum(
      rowSums(count_matrix) == 0
    ),
  minimum_library_size =
    min(raw_library_sizes),
  maximum_library_size =
    max(raw_library_sizes),
  median_library_size =
    median(raw_library_sizes),
  minimum_zero_fraction =
    min(zero_count_fraction),
  maximum_zero_fraction =
    max(zero_count_fraction)
)

cat(
  "Input count-matrix dimensions:\n"
)

print(
  dim(count_matrix)
)

cat(
  "\nInput QC summary:\n"
)

print(
  input_qc
)

capture.output(
  input_qc,
  file = file.path(
    QC_DIR,
    "GSE139028_input_QC_summary.txt"
  )
)


# ---- 04. Parse GEO SOFT metadata and audit sample structure / 解析SOFT元数据并审计样本结构 ----

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
    "^\\^SAMPLE\\s*=",
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

  collapse_or_na <- function(values) {
    if (length(values) == 0) {
      return(NA_character_)
    }

    paste(
      values,
      collapse = " | "
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
        "^\\^SAMPLE\\s*=\\s*",
        "",
        block[1]
      )

      relation_text <- collapse_or_na(
        get_values(
          block,
          "!Sample_relation"
        )
      )

      biosample_id <- ifelse(
        !is.na(relation_text) &&
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
        !is.na(relation_text) &&
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
        Sample_ID = sample_id,
        Sample_Name = first_or_na(
          get_values(
            block,
            "!Sample_title"
          )
        ),
        Original_Title = first_or_na(
          get_values(
            block,
            "!Sample_title"
          )
        ),
        Original_Source_Name = first_or_na(
          get_values(
            block,
            "!Sample_source_name_ch1"
          )
        ),
        Original_Characteristics = collapse_or_na(
          get_values(
            block,
            "!Sample_characteristics_ch1"
          )
        ),
        Organism = first_or_na(
          get_values(
            block,
            "!Sample_organism_ch1"
          )
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

if (nrow(sample_metadata) != 9) {
  stop(
    paste0(
      "The SOFT file should contain nine samples, but ",
      nrow(sample_metadata),
      " were parsed."
    )
  )
}

source_lower <- tolower(
  trimws(
    sample_metadata$Original_Source_Name
  )
)

sample_metadata$Group <- ifelse(
  grepl(
    "normal skin",
    source_lower
  ),
  "Normal_skin",
  ifelse(
    grepl(
      "burn skin|eschar",
      source_lower
    ),
    "Burn_eschar",
    NA_character_
  )
)

if (any(
  is.na(
    sample_metadata$Group
  )
)) {
  stop(
    paste0(
      "Unexpected SOFT source-name labels were found: ",
      paste(
        unique(
          sample_metadata$Original_Source_Name[
            is.na(
              sample_metadata$Group
            )
          ]
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

sample_metadata$Patient_ID <-
  NA_character_

sample_metadata$Tissue <-
  "Skin"

sample_metadata$Time_or_Stage <- ifelse(
  sample_metadata$Group ==
    "Burn_eschar",
  "Within one week after injury; exact AE-level day unavailable",
  "Uninjured skin"
)

sample_metadata$Treatment <-
  NA_character_

sample_metadata$Outcome <-
  NA_character_

sample_metadata$Data_Type <-
  "RNA-seq raw counts"

sample_metadata$Is_Paired <-
  FALSE

sample_metadata$Is_Pooled <-
  FALSE

sample_metadata$Is_Repeated_Measure <-
  FALSE

sample_metadata$Metadata_Confidence <-
  "Direct_from_GEO_SOFT_for_group_and_platform; cohort_context_from_publication"

sample_metadata$Quality_Notes <- ifelse(
  sample_metadata$Group ==
    "Burn_eschar",
  paste0(
    "Burn-eschar tissue excised within one week after injury; exact AE-level ",
    "patient identity, age, sex, race, TBSA, burn type, wound depth, ",
    "anatomical site, and post-burn day cannot be linked from GEO to the ",
    "publication demographic table."
  ),
  paste0(
    "Normal skin was collected during elective plastic-surgery procedures; ",
    "exact AE-level patient identity, age, sex, race, and anatomical site ",
    "cannot be linked from GEO to the publication demographic table."
  )
)

if (!setequal(
  sample_metadata$Sample_Name,
  colnames(count_matrix)
)) {
  stop(
    paste0(
      "Count-matrix and SOFT sample names do not match. Missing from SOFT: ",
      paste(
        setdiff(
          colnames(count_matrix),
          sample_metadata$Sample_Name
        ),
        collapse = ", "
      ),
      ". Missing from count matrix: ",
      paste(
        setdiff(
          sample_metadata$Sample_Name,
          colnames(count_matrix)
        ),
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

actual_group_counts <- table(
  sample_metadata$Group
)

expected_group_counts <- c(
  Normal_skin = 3,
  Burn_eschar = 6
)

if (!all(
  actual_group_counts[
    names(expected_group_counts)
  ] == expected_group_counts
)) {
  stop(
    paste0(
      "Unexpected group counts: ",
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
    "Count-matrix sample count",
    "SOFT sample count",
    "Normal-skin sample count",
    "Burn-eschar sample count",
    "Count-matrix columns match SOFT titles",
    "Patient-level identifiers available",
    "Cohort-level RNA-seq donor count reported in publication",
    "Paired design",
    "Repeated-measure design",
    "Pooled samples",
    "Exact AE-level post-burn day available",
    "Exact AE-to-publication-demographic mapping available"
  ),
  Value = c(
    ncol(count_matrix),
    nrow(sample_metadata),
    actual_group_counts["Normal_skin"],
    actual_group_counts["Burn_eschar"],
    identical(
      sample_metadata$Sample_Name,
      colnames(count_matrix)
    ),
    "No",
    9,
    "No",
    "No",
    "No",
    "No",
    "No"
  ),
  Interpretation = c(
    "Directly observed in the uploaded raw-count workbook.",
    "Directly parsed from the uploaded SOFT file.",
    "Normal-skin biological tissue samples.",
    "Burn-eschar biological tissue samples.",
    "Sample order was validated before analysis.",
    "Patient_ID is unavailable at the AE-sample level.",
    "The publication reports three normal-skin and six burn-eschar RNA-seq tissue donors.",
    "The comparison is unpaired.",
    "No longitudinal or repeated-measure structure is reported.",
    "No RNA pooling is reported.",
    "Burn tissue was collected within one week, but the exact AE-level day is unavailable.",
    "Clinical rows in the paper cannot be assigned to AE1-AE9 without inference."
  ),
  stringsAsFactors = FALSE
)

write.csv(
  metadata_audit,
  file = file.path(
    QC_DIR,
    "GSE139028_metadata_consistency_check.csv"
  ),
  row.names = FALSE
)

cat("Sample metadata:\n")

print(
  sample_metadata[
    ,
    c(
      "Sample_ID",
      "Sample_Name",
      "Original_Source_Name",
      "Group"
    )
  ]
)

cat("\nSample-group counts:\n")

print(
  actual_group_counts
)

cat("\nMetadata audit:\n")

print(
  metadata_audit
)


# ---- 05. Raw-data QC / 原始数据质量控制 ----

raw_qc_metrics <- data.frame(
  Sample_ID =
    sample_metadata$Sample_ID,
  Sample_Name =
    sample_metadata$Sample_Name,
  Group =
    sample_metadata$Group,
  Raw_Library_Size =
    as.numeric(
      raw_library_sizes[
        sample_metadata$Sample_Name
      ]
    ),
  Library_Size_Million =
    as.numeric(
      raw_library_sizes[
        sample_metadata$Sample_Name
      ]
    ) / 1e6,
  Zero_Count_Fraction =
    as.numeric(
      zero_count_fraction[
        sample_metadata$Sample_Name
      ]
    ),
  stringsAsFactors = FALSE
)

write.csv(
  raw_qc_metrics,
  file = file.path(
    QC_DIR,
    "GSE139028_raw_data_QC_metrics.csv"
  ),
  row.names = FALSE
)

p_raw_qc <- ggplot(
  raw_qc_metrics,
  aes(
    x = reorder(
      Sample_Name,
      Library_Size_Million
    ),
    y = Library_Size_Million,
    fill = Group
  )
) +
  geom_col(
    width = 0.75
  ) +
  coord_flip() +
  scale_fill_manual(
    values = GROUP_COLORS,
    drop = FALSE
  ) +
  labs(
    title = "GSE139028 raw library-size QC",
    subtitle = "Total raw gene counts per RNA-seq sample",
    x = "Sample",
    y = "Library size (million counts)",
    fill = "Group"
  ) +
  theme_classic(
    base_size = 12
  )

print(
  p_raw_qc
)

ggsave(
  filename = file.path(
    FIGURES_DIR,
    "01_GSE139028_raw_data_QC.png"
  ),
  plot = p_raw_qc,
  width = 8,
  height = 5.2,
  dpi = 300
)


# ---- 06. Low-expression filtering and TMM normalization / 低表达过滤与TMM标准化 ----

group <- factor(
  sample_metadata$Group,
  levels = GROUP_LEVELS
)

dge <- edgeR::DGEList(
  counts = count_matrix,
  group = group,
  genes = data.frame(
    NCBI_Gene_ID =
      rownames(count_matrix),
    stringsAsFactors = FALSE
  )
)

keep_gene <- edgeR::filterByExpr(
  dge,
  group = group
)

if (!any(keep_gene)) {
  stop(
    "No genes were retained by filterByExpr."
  )
}

dge_filtered <- dge[
  keep_gene,
  ,
  keep.lib.sizes = FALSE
]

dge_filtered <- edgeR::calcNormFactors(
  dge_filtered,
  method = "TMM"
)

normalized_expression <- edgeR::cpm(
  dge_filtered,
  log = TRUE,
  prior.count = 2
)

if (any(is.na(normalized_expression))) {
  stop(
    "The normalized expression matrix contains missing values."
  )
}

if (any(!is.finite(normalized_expression))) {
  stop(
    "The normalized expression matrix contains non-finite values."
  )
}

filtering_summary <- data.frame(
  Metric = c(
    "Input genes",
    "All-zero input genes",
    "Genes retained by filterByExpr",
    "Genes removed before testing",
    "Tested-gene fraction"
  ),
  Value = c(
    nrow(count_matrix),
    sum(
      rowSums(count_matrix) == 0
    ),
    sum(keep_gene),
    sum(!keep_gene),
    sum(keep_gene) /
      length(keep_gene)
  ),
  stringsAsFactors = FALSE
)

write.csv(
  filtering_summary,
  file = file.path(
    QC_DIR,
    "GSE139028_filtering_summary.csv"
  ),
  row.names = FALSE
)

sample_metadata$Raw_Library_Size <-
  raw_qc_metrics$Raw_Library_Size[
    match(
      sample_metadata$Sample_ID,
      raw_qc_metrics$Sample_ID
    )
  ]

sample_metadata$Filtered_Library_Size <-
  dge_filtered$samples$lib.size

sample_metadata$TMM_Normalization_Factor <-
  dge_filtered$samples$norm.factors

sample_metadata$Effective_Library_Size <- with(
  dge_filtered$samples,
  lib.size * norm.factors
)

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
    "Organism",
    "Instrument",
    "Library_Strategy",
    "Raw_Library_Size",
    "Filtered_Library_Size",
    "TMM_Normalization_Factor",
    "Effective_Library_Size"
  )
]

write.csv(
  sample_metadata_output,
  file = file.path(
    RESULTS_DIR,
    "GSE139028_sample_metadata.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

normalized_output <- data.frame(
  NCBI_Gene_ID =
    rownames(normalized_expression),
  normalized_expression,
  check.names = FALSE
)

normalized_connection <- gzfile(
  file.path(
    RESULTS_DIR,
    "GSE139028_normalized_expression.csv.gz"
  ),
  open = "wt"
)

write.csv(
  normalized_output,
  file = normalized_connection,
  row.names = FALSE
)

close(
  normalized_connection
)

normalized_long <- data.frame(
  Expression =
    as.vector(
      normalized_expression
    ),
  Sample_Name = rep(
    colnames(normalized_expression),
    each = nrow(normalized_expression)
  ),
  stringsAsFactors = FALSE
)

normalized_long$Group <- factor(
  sample_metadata$Group[
    match(
      normalized_long$Sample_Name,
      sample_metadata$Sample_Name
    )
  ],
  levels = GROUP_LEVELS
)

p_normalized <- ggplot(
  normalized_long,
  aes(
    x = Sample_Name,
    y = Expression,
    fill = Group
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.75
  ) +
  scale_fill_manual(
    values = GROUP_COLORS,
    drop = FALSE
  ) +
  labs(
    title = "GSE139028 normalized expression distribution",
    subtitle = "TMM-normalized log2 counts per million across filtered genes",
    x = "Sample",
    y = "log2 CPM",
    fill = "Group"
  ) +
  theme_classic(
    base_size = 12
  )

print(
  p_normalized
)

ggsave(
  filename = file.path(
    FIGURES_DIR,
    "02_GSE139028_normalized_expression_distribution.png"
  ),
  plot = p_normalized,
  width = 8,
  height = 5.2,
  dpi = 300
)

cat("Filtering and normalization summary:\n")

print(
  filtering_summary
)


# ---- 07. PCA and sample correlation / PCA与样本相关性 ----

gene_variance <- matrixStats::rowVars(
  normalized_expression
)

n_pca_genes <- min(
  TOP_VARIABLE_GENES_FOR_PCA,
  length(gene_variance)
)

top_variable_gene_ids <- rownames(
  normalized_expression
)[
  order(
    gene_variance,
    decreasing = TRUE
  )[seq_len(n_pca_genes)]
]

pca_result <- prcomp(
  t(
    normalized_expression[
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
  Sample_ID =
    sample_metadata$Sample_ID[
      match(
        rownames(pca_result$x),
        sample_metadata$Sample_Name
      )
    ],
  Sample_Name =
    rownames(pca_result$x),
  PC1 =
    pca_result$x[, 1],
  PC2 =
    pca_result$x[, 2],
  Group = factor(
    sample_metadata$Group[
      match(
        rownames(pca_result$x),
        sample_metadata$Sample_Name
      )
    ],
    levels = GROUP_LEVELS
  ),
  stringsAsFactors = FALSE
)

write.csv(
  pca_table,
  file = file.path(
    RESULTS_DIR,
    "GSE139028_PCA_scores.csv"
  ),
  row.names = FALSE
)

p_pca <- ggplot(
  pca_table,
  aes(
    x = PC1,
    y = PC2,
    color = Group,
    label = Sample_Name
  )
) +
  geom_point(
    size = 3
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
    title = "GSE139028 PCA",
    subtitle = paste0(
      "Top ",
      n_pca_genes,
      " variable genes"
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
    "03_GSE139028_PCA.png"
  ),
  plot = p_pca,
  width = 7.5,
  height = 5.8,
  dpi = 300
)

sample_correlation <- cor(
  normalized_expression,
  method = "pearson"
)

correlation_connection <- gzfile(
  file.path(
    QC_DIR,
    "GSE139028_sample_correlation.csv.gz"
  ),
  open = "wt"
)

write.csv(
  data.frame(
    Sample_Name =
      rownames(sample_correlation),
    sample_correlation,
    check.names = FALSE
  ),
  file = correlation_connection,
  row.names = FALSE
)

close(
  correlation_connection
)

correlation_annotation <- data.frame(
  Group = factor(
    sample_metadata$Group,
    levels = GROUP_LEVELS
  ),
  row.names =
    sample_metadata$Sample_Name
)

correlation_annotation_colors <- list(
  Group = GROUP_COLORS
)

pheatmap::pheatmap(
  sample_correlation,
  annotation_col =
    correlation_annotation,
  annotation_row =
    correlation_annotation,
  annotation_colors =
    correlation_annotation_colors,
  border_color = NA,
  color =
    EXPRESSION_HEATMAP_COLORS,
  main =
    "GSE139028 sample correlation",
  silent = FALSE
)

pheatmap::pheatmap(
  sample_correlation,
  annotation_col =
    correlation_annotation,
  annotation_row =
    correlation_annotation,
  annotation_colors =
    correlation_annotation_colors,
  border_color = NA,
  color =
    EXPRESSION_HEATMAP_COLORS,
  main =
    "GSE139028 sample correlation",
  filename = file.path(
    FIGURES_DIR,
    "04_GSE139028_sample_correlation_heatmap.png"
  ),
  width = 7,
  height = 6,
  silent = TRUE
)


# ---- 08. Differential-expression model / 差异表达模型 ----

design <- model.matrix(
  ~ 0 + group
)

colnames(design) <- levels(
  group
)

rownames(design) <-
  sample_metadata$Sample_Name

if (qr(design)$rank <
    ncol(design)) {
  stop(
    "The differential-expression design matrix is not full rank."
  )
}

cat("Design-matrix columns:\n")

print(
  colnames(design)
)

dge_filtered <- edgeR::estimateDisp(
  dge_filtered,
  design = design,
  robust = TRUE
)

ql_fit <- edgeR::glmQLFit(
  dge_filtered,
  design = design,
  robust = TRUE
)

burn_vs_normal_contrast <- c(
  Normal_skin = -1,
  Burn_eschar = 1
)

ql_test <- edgeR::glmQLFTest(
  ql_fit,
  contrast = burn_vs_normal_contrast
)

de_table <- edgeR::topTags(
  ql_test,
  n = Inf,
  sort.by = "PValue"
)$table

if (nrow(de_table) !=
    nrow(dge_filtered)) {
  stop(
    "The differential-expression result table does not contain every tested gene."
  )
}

de_table$NCBI_Gene_ID <-
  rownames(de_table)

cat("Top differential-expression results:\n")

print(
  head(
    de_table,
    10
  )
)


# ---- 09. Gene-ID mapping and complete result table / Gene ID映射与完整结果表 ----

tested_entrez_ids <- as.character(
  de_table$NCBI_Gene_ID
)

gene_symbol <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = tested_entrez_ids,
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
)

ensembl_id <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = tested_entrez_ids,
  column = "ENSEMBL",
  keytype = "ENTREZID",
  multiVals = "first"
)

gene_name <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = tested_entrez_ids,
  column = "GENENAME",
  keytype = "ENTREZID",
  multiVals = "first"
)

all_gene_results <- data.frame(
  NCBI_Gene_ID =
    tested_entrez_ids,
  Gene_Symbol = unname(
    gene_symbol[
      tested_entrez_ids
    ]
  ),
  Ensembl_ID = unname(
    ensembl_id[
      tested_entrez_ids
    ]
  ),
  Gene_Name = unname(
    gene_name[
      tested_entrez_ids
    ]
  ),
  Mapping_Status = ifelse(
    is.na(
      unname(
        gene_symbol[
          tested_entrez_ids
        ]
      )
    ),
    "Unmapped",
    "Mapped"
  ),
  GEO_ID = GEO_ID,
  Contrast_ID = CONTRAST_ID,
  Contrast_Label =
    CONTRAST_LABEL,
  Case_Group =
    "Burn eschar",
  Control_Group =
    "Normal skin",
  Case_N = 6,
  Control_N = 3,
  Case_Patient_N = 6,
  Control_Patient_N = 3,
  log2FC =
    de_table$logFC,
  Fold_Change =
    2^de_table$logFC,
  Mean_log2CPM =
    de_table$logCPM,
  Statistic =
    de_table$F,
  P_value =
    de_table$PValue,
  FDR =
    de_table$FDR,
  stringsAsFactors = FALSE
)

all_gene_results$Direction <- ifelse(
  all_gene_results$log2FC > 0,
  "Up",
  ifelse(
    all_gene_results$log2FC < 0,
    "Down",
    "No_change"
  )
)

all_gene_results$DE_Status <- ifelse(
  all_gene_results$FDR <
    FDR_CUTOFF &
    all_gene_results$log2FC >=
      LOG2FC_CUTOFF,
  "Up_significant",
  ifelse(
    all_gene_results$FDR <
      FDR_CUTOFF &
      all_gene_results$log2FC <=
        -LOG2FC_CUTOFF,
    "Down_significant",
    "Not_significant"
  )
)

all_gene_results$NegLog10_FDR <- -log10(
  pmax(
    all_gene_results$FDR,
    .Machine$double.xmin
  )
)

if (anyDuplicated(
  paste(
    all_gene_results$NCBI_Gene_ID,
    all_gene_results$Contrast_ID,
    sep = "::"
  )
) > 0) {
  stop(
    "Duplicated Gene ID and Contrast ID combinations were detected."
  )
}

all_results_connection <- gzfile(
  file.path(
    RESULTS_DIR,
    "GSE139028_all_gene_results.csv.gz"
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

cat("Differential-expression result counts:\n")

print(
  table(
    all_gene_results$DE_Status
  )
)


# ---- 10. Volcano plot / 火山图 ----

volcano_data <- all_gene_results

volcano_data$DE_Status <- factor(
  volcano_data$DE_Status,
  levels = c(
    "Down_significant",
    "Not_significant",
    "Up_significant"
  )
)

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
    size = 1.2
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
    title =
      "GSE139028: burn eschar vs normal skin",
    subtitle = paste0(
      "edgeR robust quasi-likelihood model; ",
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
      "05_GSE139028_volcano_",
      CONTRAST_ID,
      ".png"
    )
  ),
  plot = p_volcano,
  width = 8,
  height = 6,
  dpi = 300
)


# ---- 11. Top differential-gene heatmap / 主要差异基因热图 ----

ranked_gene_ids <-
  all_gene_results$NCBI_Gene_ID[
    order(
      all_gene_results$FDR,
      -abs(
        all_gene_results$log2FC
      )
    )
  ]

n_heatmap_genes <- min(
  TOP_GENES_FOR_HEATMAP,
  length(ranked_gene_ids)
)

top_heatmap_ids <- ranked_gene_ids[
  seq_len(n_heatmap_genes)
]

heatmap_matrix <- normalized_expression[
  top_heatmap_ids,
  ,
  drop = FALSE
]

heatmap_z <- t(
  scale(
    t(heatmap_matrix)
  )
)

gene_label_match <- match(
  top_heatmap_ids,
  all_gene_results$NCBI_Gene_ID
)

heatmap_labels <-
  all_gene_results$Gene_Symbol[
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
  make.unique(
    heatmap_labels
  )

heatmap_annotation <- data.frame(
  Group = factor(
    sample_metadata$Group,
    levels = GROUP_LEVELS
  ),
  row.names =
    sample_metadata$Sample_Name
)

heatmap_annotation_colors <- list(
  Group = GROUP_COLORS
)

heatmap_title <- paste0(
  "Top ",
  n_heatmap_genes,
  " genes: ",
  CONTRAST_LABEL
)

pheatmap::pheatmap(
  heatmap_z,
  annotation_col =
    heatmap_annotation,
  annotation_colors =
    heatmap_annotation_colors,
  show_colnames = TRUE,
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
  show_colnames = TRUE,
  show_rownames = TRUE,
  border_color = NA,
  color =
    EXPRESSION_HEATMAP_COLORS,
  main = heatmap_title,
  filename = file.path(
    FIGURES_DIR,
    paste0(
      "06_GSE139028_top_differential_genes_heatmap_",
      CONTRAST_ID,
      ".png"
    )
  ),
  width = 8,
  height = 10,
  silent = TRUE
)


# ---- 12. Cross-check protease genes reported in the publication / 核对论文报告的蛋白酶基因 ----

published_protease_fc <- data.frame(
  Gene_Symbol = c(
    "CTSG",
    "ELANE",
    "PRTN3",
    "CTSK",
    "CTSL",
    "CTSS",
    "CTSV",
    "MMP2",
    "MMP7",
    "MMP9",
    "MMP12",
    "MMP1",
    "MMP8",
    "MMP13"
  ),
  Published_Fold_Change = c(
    1.28,
    0.50,
    2.90,
    1.59,
    18.41,
    10.47,
    0.06,
    3.03,
    38.55,
    65.37,
    57.35,
    659.88,
    3588.45,
    245.37
  ),
  stringsAsFactors = FALSE
)

input_entrez_ids <- rownames(
  count_matrix
)

input_gene_symbols <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = input_entrez_ids,
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
)

input_symbol_set <- unique(
  unname(
    input_gene_symbols[
      !is.na(
        input_gene_symbols
      )
    ]
  )
)

tested_crosscheck <- all_gene_results[
  all_gene_results$Gene_Symbol %in%
    published_protease_fc$Gene_Symbol,
  c(
    "NCBI_Gene_ID",
    "Gene_Symbol",
    "log2FC",
    "Fold_Change",
    "P_value",
    "FDR",
    "Direction",
    "DE_Status"
  )
]

published_gene_crosscheck <- merge(
  published_protease_fc,
  tested_crosscheck,
  by = "Gene_Symbol",
  all.x = TRUE,
  sort = FALSE
)

published_gene_crosscheck$Reanalysis_Status <- ifelse(
  !is.na(
    published_gene_crosscheck$log2FC
  ),
  "Tested",
  ifelse(
    published_gene_crosscheck$Gene_Symbol %in%
      input_symbol_set,
    "Filtered_low_expression",
    "Not_found_or_unmapped_in_input"
  )
)

published_gene_crosscheck$Direction_Consistent <- ifelse(
  published_gene_crosscheck$Reanalysis_Status ==
    "Tested",
  sign(
    log2(
      published_gene_crosscheck$Published_Fold_Change
    )
  ) == sign(
    published_gene_crosscheck$log2FC
  ),
  NA
)

published_gene_crosscheck$Calculated_to_Published_FC_Ratio <- ifelse(
  published_gene_crosscheck$Reanalysis_Status ==
    "Tested",
  published_gene_crosscheck$Fold_Change /
    published_gene_crosscheck$Published_Fold_Change,
  NA_real_
)

write.csv(
  published_gene_crosscheck,
  file = file.path(
    RESULTS_DIR,
    "GSE139028_published_gene_crosscheck.csv"
  ),
  row.names = FALSE
)

cat("Published-gene cross-check:\n")

print(
  published_gene_crosscheck
)


# ---- 13. Create the BurnOmicsDB-ready result table / 创建BurnOmicsDB标准结果表 ----

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
  Organism =
    "Homo sapiens",
  Study_Population =
    "Pediatric",
  Tissue =
    "Skin",
  Sample_Context = paste0(
    "Tangentially excised burn-eschar tissue versus ",
    "normal skin collected during elective plastic surgery"
  ),
  Time_or_Stage = paste0(
    "Burn tissue excised within one week after injury; ",
    "exact AE-level day unavailable"
  ),
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
  Mean_log2CPM =
    all_gene_results$Mean_log2CPM,
  P_value =
    all_gene_results$P_value,
  FDR =
    all_gene_results$FDR,
  DE_Status =
    all_gene_results$DE_Status,
  Platform =
    "Illumina HiSeq 1000 (GPL15433)",
  Input_Data =
    "Author-provided NCBI Entrez Gene raw counts",
  Normalization =
    "TMM",
  Analysis_Method = paste0(
    "edgeR filterByExpr low-expression filtering, TMM normalization, ",
    "robust quasi-likelihood negative-binomial GLM, and unpaired ",
    "Burn_eschar minus Normal_skin contrast; positive log2FC represents ",
    "Case_Group minus Control_Group"
  ),
  Annotation_Method = paste0(
    "org.Hs.eg.db ",
    as.character(
      packageVersion(
        "org.Hs.eg.db"
      )
    ),
    "; Entrez Gene ID mapped with mapIds(multiVals='first')"
  ),
  Quality_Notes = paste0(
    "Small unpaired pediatric bulk-tissue cohort; three normal-skin and six ",
    "burn-eschar RNA-seq tissue donors are reported, but exact AE-level ",
    "patient linkage and clinical demographics are unavailable; burn type, ",
    "wound depth, TBSA, anatomical site, and exact sampling day are ",
    "heterogeneous or unavailable at the AE level; normal skin was collected ",
    "during elective plastic surgery; expression differences may reflect ",
    "intracellular regulation, leukocyte infiltration, cell-composition ",
    "changes, and major tissue-structure differences"
  ),
  Case_Patient_N =
    all_gene_results$Case_Patient_N,
  Control_Patient_N =
    all_gene_results$Control_Patient_N,
  Is_Paired_Contrast =
    FALSE,
  Is_Pooled =
    FALSE,
  Expression_Scale =
    "TMM-normalized log2 CPM with prior count 2",
  Metadata_Confidence = paste0(
    "Direct_from_GEO_SOFT_for_sample_group_and_platform; ",
    "cohort_context_from_publication"
  ),
  Quality_Grade =
    "Core_small_cohort",
  stringsAsFactors = FALSE
)

core_column_order <- c(
  "NCBI_Gene_ID",
  "Gene_Symbol",
  "Ensembl_ID",
  "Gene_Name",
  "GEO_ID",
  "Organism",
  "Study_Population",
  "Tissue",
  "Sample_Context",
  "Time_or_Stage",
  "Contrast_ID",
  "Contrast_Label",
  "Case_Group",
  "Control_Group",
  "Case_N",
  "Control_N",
  "log2FC",
  "Fold_Change",
  "Direction",
  "Mean_log2CPM",
  "P_value",
  "FDR",
  "DE_Status",
  "Platform",
  "Input_Data",
  "Normalization",
  "Analysis_Method",
  "Annotation_Method",
  "Quality_Notes"
)

if (!identical(
  colnames(database_ready)[
    seq_along(core_column_order)
  ],
  core_column_order
)) {
  stop(
    "The BurnOmicsDB core-column order is incorrect."
  )
}

database_connection <- gzfile(
  file.path(
    RESULTS_DIR,
    "GSE139028_database_ready_all_genes.csv.gz"
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
  nrow(database_ready),
  "\n",
  sep = ""
)


# ---- 14. Save analysis objects, summary, and session information / 保存分析对象、摘要与环境信息 ----

significant_count_table <- as.data.frame(
  table(
    all_gene_results$DE_Status
  ),
  stringsAsFactors = FALSE
)

colnames(significant_count_table) <- c(
  "DE_Status",
  "Count"
)

analysis_objects <- list(
  project =
    "BurnOmicsDB",
  GEO_ID =
    GEO_ID,
  counts_file =
    counts_file,
  soft_file =
    soft_file,
  paper_file =
    paper_file,
  sample_metadata =
    sample_metadata_output,
  metadata_audit =
    metadata_audit,
  input_qc =
    input_qc,
  raw_qc_metrics =
    raw_qc_metrics,
  filtering_summary =
    filtering_summary,
  raw_count_dimensions =
    dim(count_matrix),
  retained_gene_logical =
    keep_gene,
  filtered_edgeR_object =
    dge_filtered,
  normalized_expression =
    normalized_expression,
  design_matrix =
    design,
  ql_fit =
    ql_fit,
  ql_test =
    ql_test,
  all_gene_results =
    all_gene_results,
  published_gene_crosscheck =
    published_gene_crosscheck,
  database_ready =
    database_ready,
  thresholds = list(
    FDR_CUTOFF =
      FDR_CUTOFF,
    LOG2FC_CUTOFF =
      LOG2FC_CUTOFF,
    PAPER_FDR_CUTOFF =
      PAPER_FDR_CUTOFF
  )
)

saveRDS(
  analysis_objects,
  file = file.path(
    OBJECTS_DIR,
    "GSE139028_analysis_objects.rds"
  ),
  compress = "xz"
)

n_up <- sum(
  all_gene_results$DE_Status ==
    "Up_significant"
)

n_down <- sum(
  all_gene_results$DE_Status ==
    "Down_significant"
)

n_fdr_001 <- sum(
  all_gene_results$FDR <
    PAPER_FDR_CUTOFF
)

mapped_gene_n <- sum(
  all_gene_results$Mapping_Status ==
    "Mapped"
)

unmapped_gene_n <- sum(
  all_gene_results$Mapping_Status ==
    "Unmapped"
)

summary_lines <- c(
  "BurnOmicsDB - GSE139028 analysis summary",
  "",
  paste0(
    "Analysis date: ",
    Sys.Date()
  ),
  "Data type: RNA-seq raw counts summarized at NCBI Entrez Gene level",
  paste0(
    "Input genes: ",
    nrow(count_matrix)
  ),
  paste0(
    "Input samples: ",
    ncol(count_matrix)
  ),
  "RNA-seq tissue donors represented: 9",
  "Normal-skin tissue donors represented: 3",
  "Burn-eschar tissue donors represented: 6",
  "Exact AE-level patient linkage: unavailable",
  "Pooling: No",
  "Paired design: No",
  "Repeated measures: No",
  "",
  "Groups:",
  paste0(
    "  Normal skin: ",
    actual_group_counts["Normal_skin"],
    " samples"
  ),
  paste0(
    "  Burn eschar: ",
    actual_group_counts["Burn_eschar"],
    " samples"
  ),
  "",
  paste0(
    "All-zero input genes: ",
    input_qc$all_zero_gene_n
  ),
  paste0(
    "Genes retained after filterByExpr: ",
    sum(keep_gene)
  ),
  paste0(
    "Genes removed before statistical testing: ",
    sum(!keep_gene)
  ),
  paste0(
    "Mapped tested genes: ",
    mapped_gene_n
  ),
  paste0(
    "Unmapped tested genes: ",
    unmapped_gene_n
  ),
  "",
  "Normalization: TMM",
  "Normalized expression scale: log2 CPM with prior count 2",
  paste0(
    "Statistical model: edgeR robust quasi-likelihood negative-binomial GLM; ",
    "unpaired Burn_eschar minus Normal_skin contrast"
  ),
  paste0(
    "Threshold: |log2FC| >= ",
    LOG2FC_CUTOFF,
    " and FDR < ",
    FDR_CUTOFF
  ),
  "",
  "Contrast:",
  paste0(
    "  ",
    CONTRAST_ID,
    ": ",
    CONTRAST_LABEL,
    " (case samples = 6, control samples = 3)"
  ),
  "",
  "Differential-expression counts:",
  paste0(
    "  Up_significant: ",
    n_up
  ),
  paste0(
    "  Down_significant: ",
    n_down
  ),
  paste0(
    "  Not_significant: ",
    nrow(all_gene_results) -
      n_up -
      n_down
  ),
  paste0(
    "  Genes with FDR < 0.01 regardless of effect-size threshold: ",
    n_fdr_001
  ),
  "",
  paste0(
    "Positive log2FC means higher expression in burn eschar than in ",
    "normal skin."
  ),
  paste0(
    "Only genes retained by filterByExpr were statistically tested and ",
    "exported. Low-expression genes were not added back to the result tables."
  ),
  "",
  "Major limitations:",
  paste0(
    "  The cohort is small and the comparison is unpaired."
  ),
  paste0(
    "  GEO does not provide exact AE-to-patient demographic linkage."
  ),
  paste0(
    "  Burn type, wound depth, TBSA, anatomical site, and exact post-burn ",
    "sampling day are heterogeneous or unavailable at the AE level."
  ),
  paste0(
    "  Control skin was obtained during elective plastic-surgery procedures."
  ),
  paste0(
    "  Bulk-tissue expression differences can reflect intracellular ",
    "regulation, inflammatory-cell infiltration, cell composition, and ",
    "tissue-structure differences."
  ),
  "Quality grade: Core_small_cohort"
)

writeLines(
  summary_lines,
  con = file.path(
    RESULTS_DIR,
    "GSE139028_analysis_summary.txt"
  ),
  useBytes = TRUE
)

capture.output(
  sessionInfo(),
  file = file.path(
    RESULTS_DIR,
    "GSE139028_R_sessionInfo.txt"
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
  "MMP1",
  "MMP8",
  "MMP9",
  "MMP13",
  "CTSL",
  "CTSS",
  "CTSV",
  "IL6"
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
    )
  ),
]

print(
  gene_check_table
)

# End of script / 代码结束
