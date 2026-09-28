#!/usr/bin/env Rscript
# Author: Tao Xiong
# Date: 2026-09-11
# Description: 区分物种的原生分布点 (Native) 与引种分布点 (Introduced)，并分别输出到独立 CSV 文件。
#              （直接基于 WCVP Accepted Name 进行空间相交计算，跳过重复的名称标准化）

suppressPackageStartupMessages({
  library(dplyr)
  library(rWCVP)
  library(rWCVPdata)
  library(tidyverse)
  library(sf)
  library(data.table)
  library(stringr)
})

# ==============================================================================
# 1. 容错读取分布数据函数（解决 stateProvince 中未转义逗号导致提前截断的问题）
# ==============================================================================
read_distribution_csv <- function(file_path) {
  # 读取所有行以检查行数
  lines <- readLines(file_path, encoding = "UTF-8", warn = FALSE)
  total_lines <- length(lines)
  if (total_lines <= 1) stop("输入文件内容为空！")
  
  # 先尝试用 fread 快速读取
  df <- tryCatch(
    suppressWarnings(data.table::fread(file_path, data.table = FALSE, showProgress = FALSE)),
    error = function(e) NULL
  )
  
  # 如果 fread 成功读取且行数与文件一致（总行数 - 1），直接返回
  if (!is.null(df) && nrow(df) == (total_lines - 1)) {
    cat(sprintf("✓ 快速读取完成，共 %d 条记录。\n", nrow(df)))
    return(df)
  }
  
  # 否则说明因地址列包含未加双引号的逗号导致 fread 提前截断（如 line 275517），启动正则容错解析
  cat(sprintf("! 检测到地名字段包含未转义逗号导致标准读取截断（文件总行数: %d），启动正则容错解析...\n", total_lines - 1))
  
  data_lines <- lines[-1]
  data_lines <- data_lines[nzchar(trimws(data_lines))]
  
  # 格式: Species, stateProvince (可能含逗号), countryCode, Longitude, Latitude, type
  pattern <- "^([^,]+),(.*),([^,]*),([-+]?[0-9]*\\.?[0-9]+(?:[eE][-+]?[0-9]+)?),([-+]?[0-9]*\\.?[0-9]+(?:[eE][-+]?[0-9]+)?),([^,]+)$"
  m <- stringr::str_match(data_lines, pattern)
  valid_idx <- !is.na(m[, 1])
  
  parsed_df <- data.frame(
    Species       = trimws(m[valid_idx, 2]),
    stateProvince = trimws(m[valid_idx, 3]),
    countryCode   = trimws(m[valid_idx, 4]),
    Longitude     = as.numeric(m[valid_idx, 5]),
    Latitude      = as.numeric(m[valid_idx, 6]),
    type          = trimws(m[valid_idx, 7]),
    stringsAsFactors = FALSE
  )
  
  cat(sprintf("✓ 容错解析完成，成功解析并恢复 %d / %d 条完整记录。\n", nrow(parsed_df), length(data_lines)))
  return(parsed_df)
}

# ==============================================================================
# 2. 读取输入数据与准备物种名录
# ==============================================================================
args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) > 0) args[1] else "distribution_clean_distinct.csv"

if (!file.exists(input_file)) {
  stop(paste("未找到输入文件:", input_file, "\n用法: Rscript filter_native_and_introduced.R [input.csv]"), call. = FALSE)
}

cat(sprintf("==> 正在读取输入数据: %s ...\n", input_file))
data <- read_distribution_csv(input_file)

# 基础列检查
required_cols <- c("Species", "Longitude", "Latitude")
if (!all(required_cols %in% colnames(data))) {
  stop(paste("输入文件缺少必要字段:", paste(setdiff(required_cols, colnames(data)), collapse = ", ")), call. = FALSE)
}

# 格式清洗与准备：由于物种名已是 WCVP Accepted Name，直接赋值为 taxon_name，无需重复标准化
data <- data %>%
  filter(!is.na(Species) & trimws(Species) != "") %>%
  mutate(Species = trimws(Species)) %>%
  mutate(Species = str_replace(Species, "^[a-z]", toupper)) %>%
  mutate(taxon_name = Species)

cat(sprintf("==> 有效分布记录: %d 条，包含接受物种 (Accepted Species): %d 个。\n",
            nrow(data), length(unique(data$taxon_name))))

# ==============================================================================
# 3. 空间相交分析：区分 Native 与 Introduced
# ==============================================================================
dir.create("./native", showWarnings = FALSE, recursive = TRUE)
dir.create("./introduced", showWarnings = FALSE, recursive = TRUE)

unique_taxa <- sort(unique(data$taxon_name))
missing_log <- "./missing_distribution_ranges.txt"
if (file.exists(missing_log)) file.remove(missing_log)

all_native_list  <- list()
all_intro_list   <- list()
all_outside_list <- list()

# 预先获取本地 rWCVPdata 数据集指针，避免远程校验 .wcvp_fresh 失败
w_names <- rWCVPdata::wcvp_names
w_dist  <- rWCVPdata::wcvp_distributions

cat("==> 开始逐个物种进行空间分布区多边形判定...\n")

for (i in seq_along(unique_taxa)) {
  name <- unique_taxa[i]
  cat(sprintf("[%d/%d] 正在处理: %s ...\n", i, length(unique_taxa), name))
  
  sp_df <- data %>% filter(taxon_name == name)
  
  # 转为 sf 空间点对象（WGS84, EPSG:4326），remove = FALSE 保留原始经纬度字段
  occs <- tryCatch({
    st_as_sf(sp_df, coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE)
  }, error = function(e) return(NULL))
  
  if (is.null(occs) || nrow(occs) == 0) next
  
  # 下载物种在 WCVP 中的地理分布多边形（显式传入本地数据集，避免联网验证报错）
  dist_poly <- tryCatch({
    wcvp_distribution(name,
                      taxon_rank = "species",
                      introduced = TRUE,
                      extinct = FALSE,
                      location_doubtful = FALSE,
                      wcvp_names = w_names,
                      wcvp_distributions = w_dist)
  }, error = function(e) {
    write(paste("Skipped:", name, "- Error in wcvp_distribution:", e$message),
          file = missing_log, append = TRUE)
    return(NULL)
  })
  
  if (is.null(dist_poly) || nrow(dist_poly) == 0) {
    write(paste("Skipped:", name, "- WCVP 中无该物种分布区多边形"),
          file = missing_log, append = TRUE)
    next
  }
  
  # 分离原生区和引种区多边形
  native_polys <- dist_poly %>% filter(occurrence_type == "native")
  intro_polys  <- dist_poly %>% filter(occurrence_type == "introduced")
  
  # 1) 判断是否落在原生区缓冲范围内 (~1km 缓冲，0.009 度)
  is_native <- rep(FALSE, nrow(occs))
  if (nrow(native_polys) > 0) {
    buffered_native <- native_polys %>% st_union() %>% st_buffer(0.009)
    is_native <- st_intersects(occs, buffered_native, sparse = FALSE)[, 1]
  }
  
  # 2) 判断是否落在引种区缓冲范围内
  is_intro <- rep(FALSE, nrow(occs))
  if (nrow(intro_polys) > 0) {
    buffered_intro <- intro_polys %>% st_union() %>% st_buffer(0.009)
    is_intro <- st_intersects(occs, buffered_intro, sparse = FALSE)[, 1]
  }
  
  # 3) 分类判定逻辑：
  #    - 优先判定为 native
  #    - 若不在原生区但在引种区内，判定为 introduced
  #    - 既不在原生区也不在已知引种区（如近海漂移、未记录引种区等），标记为 outside_or_unresolved
  occs$distribution_status <- case_when(
    is_native ~ "native",
    !is_native & is_intro ~ "introduced",
    TRUE ~ "outside_or_unresolved"
  )
  
  # 还原为纯数据框，保留全部原始元数据字段
  sp_res <- st_drop_geometry(occs)
  
  sp_native  <- sp_res %>% filter(distribution_status == "native")
  sp_intro   <- sp_res %>% filter(distribution_status == "introduced")
  sp_outside <- sp_res %>% filter(distribution_status == "outside_or_unresolved")
  
  # 单物种分文件保存
  safe_name <- gsub(" ", "_", name)
  if (nrow(sp_native) > 0) {
    write.csv(sp_native, paste0("./native/", safe_name, "_native.csv"), row.names = FALSE)
    all_native_list[[name]] <- sp_native
  }
  if (nrow(sp_intro) > 0) {
    write.csv(sp_intro, paste0("./introduced/", safe_name, "_introduced.csv"), row.names = FALSE)
    all_intro_list[[name]] <- sp_intro
  }
  if (nrow(sp_outside) > 0) {
    all_outside_list[[name]] <- sp_outside
  }
}

# ==============================================================================
# 4. 汇总保存
# ==============================================================================
cat("\n==> 正在导出全量汇总文件...\n")

if (length(all_native_list) > 0) {
  all_native_df <- bind_rows(all_native_list)
  write.csv(all_native_df, "./native/All_Native.csv", row.names = FALSE)
  cat(sprintf("✓ 原生分布点全量汇总: ./native/All_Native.csv (共 %d 条记录)\n", nrow(all_native_df)))
}

if (length(all_intro_list) > 0) {
  all_intro_df <- bind_rows(all_intro_list)
  write.csv(all_intro_df, "./introduced/All_Introduced.csv", row.names = FALSE)
  cat(sprintf("✓ 引种分布点全量汇总: ./introduced/All_Introduced.csv (共 %d 条记录)\n", nrow(all_intro_df)))
}

if (length(all_outside_list) > 0) {
  all_outside_df <- bind_rows(all_outside_list)
  write.csv(all_outside_df, "All_Outside_Unresolved.csv", row.names = FALSE)
  cat(sprintf("ℹ 既非原生也非已知引种点: All_Outside_Unresolved.csv (共 %d 条记录)\n", nrow(all_outside_df)))
}

cat("\n==> 全部处理完成！未获取到分布区信息的物种已记录至: missing_distribution_ranges.txt\n")
