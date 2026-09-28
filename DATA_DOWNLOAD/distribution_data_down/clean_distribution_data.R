# !/usr/bin/Rscript
# Author: Tao Xiong
# Date: 2026-09-10
# 本脚本用于清理分布数据,但仍有部分要求，要求需要进行清理的数据集必须包含Species、Longitude、Latitude三个字段作为列名

suppressPackageStartupMessages({
  library(dplyr)
  library(CoordinateCleaner)
  library(data.table)
})

# 1. 读取命令行参数
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("请提供输入文件路径，例如: Rscript clean_coord.R input.csv", call. = FALSE)
}
file <- args[1]
data <- fread(file, data.table = FALSE)

# 2. 检查必须包含的列名
required_cols <- c("Species", "Longitude", "Latitude")
if (!all(required_cols %in% colnames(data))) {
  stop(paste("输入文件必须包含以下列:", paste(required_cols, collapse = ", ")), call. = FALSE)
}

# 3. 基础清洗：经纬度数值化、去空值、去除非法范围坐标
occs.coord <- data %>%
  filter(!is.na(Species) & trimws(Species) != "") %>%
  mutate(
    Longitude = suppressWarnings(as.numeric(Longitude)),
    Latitude  = suppressWarnings(as.numeric(Latitude))
  ) %>%
  filter(!is.na(Longitude) & !is.na(Latitude)) %>%
  filter(Longitude >= -180 & Longitude <= 180 & Latitude >= -90 & Latitude <= 90)

if (nrow(occs.coord) == 0) {
  stop("有效坐标记录为 0，请检查输入文件格式及经纬度列！", call. = FALSE)
}

# 4. 执行 CoordinateCleaner 清洗
# 注：如果某些物种记录过少导致 outliers 报错，可通过将 tests 里的 "outliers" 移除
# 或者使用 CoordinateCleaner 默认的行为（它会自动跳过 n < 7 的物种）
flags <- clean_coordinates(
  x = occs.coord,
  lon = "Longitude",
  lat = "Latitude",
  species = "Species",
  tests = c("capitals", "centroids", "equal", "gbif", "institutions", "outliers", "seas", "zeros"),
  value = "clean" # 直接返回清洗后的 data.frame
)

# 5. 去除物种与经纬度完全一致的重复记录，并保留其他列信息
clean_distinct <- flags %>%
  distinct(Species, Longitude, Latitude, .keep_all = TRUE)

# 6. 保存输出
output_file <- "distribution_clean_distinct.csv"
write.csv(clean_distinct, output_file, quote = TRUE, row.names = FALSE)
cat(sprintf("处理完成！原始有效记录: %d 条，清洗去重后保留: %d 条。结果已保存至: %s\n", 
            nrow(occs.coord), nrow(clean_distinct), output_file))