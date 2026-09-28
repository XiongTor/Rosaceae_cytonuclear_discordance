# !/usr/bin/Rscript
# Author: Tao Xiong
# Date: 2025.03.14
# Description: This script is used to get the distribution data of Rosaceae from GBIF, idigbio, CVH, NSII.

# 地理分布数据下载
library(dplyr)
library(sf)
library(spData)
library(tmap)
library(purrr)
library(readr)
library(magrittr) # for %T>% pipe
library(rgbif) # for occ_download
library(taxize) # for get_gbifid
library(maptools)
library(raster)
library(ggplot2)
library(ggcorrplot)
library(CoordinateCleaner)
library(zoo)
library(ridigbio)
library(readxl)  



## 下载gbif数据（提前准备物种名录）------方法一
file_url <- read.csv("rosa_species_subfam_inf.csv", header=TRUE)

#设置请求的服务器IP，即VPN挂的地址
Sys.setenv(ALL_PROXY = "127.0.0.1:7890")

#下载GBIF数据
library(rgbif)

for (i in file_url$Species) {
  tryCatch({
    print(paste("Downloading GBIF data for species:", i))
    
    # 1. 使用 occ_search 并移除已废弃的 return = "data" 参数
    item_gbif <- occ_search(
      scientificName = i,
      limit = 90000,
      hasCoordinate = TRUE
    )
    
    # 2. 检查数据是否存在
    if (!is.null(item_gbif$data) && nrow(item_gbif$data) > 0) {
      # 3. 提取实际的 data 部分
      item_gbif_data <- item_gbif$data
      
      # 4. 检查需要的列是否存在
      required_cols <- c("species", "decimalLongitude", "decimalLatitude")
      if (all(required_cols %in% colnames(item_gbif_data))) {
        # 5. 筛选有效坐标（排除 0,0 坐标）
        item_gbif_reduced_cleaned <- item_gbif_data[
          item_gbif_data$decimalLatitude != 0 & 
          item_gbif_data$decimalLongitude != 0,
          required_cols
        ]
        
        # 6. 创建结果数据框（更简洁的方式）
        combined <- data.frame(
          species_name = item_gbif_reduced_cleaned$species,
          lon = item_gbif_reduced_cleaned$decimalLongitude,
          lat = item_gbif_reduced_cleaned$decimalLatitude
        )
        
        # 7. 这里可以添加保存或合并结果的代码
        # 比如写入文件或累加到总数据框
        write.table(combined, file = paste0("./results/", gsub(" ", "_", i), ".csv", sep = ""), sep = ",", quote = FALSE, row.names = FALSE)
        
        print(paste("Successfully processed", nrow(combined), "records for", i))
      } else {
        print(paste("Missing required columns for", i))
      }
    } else {
      print(paste("No data found for", i))
    }
  }, error = function(e) {
    message(paste("Error processing", i, ":", conditionMessage(e)))
  })
}




## 下载gbif数据（提前准备物种名录）------方法二
library(dplyr)
library(purrr)
library(readr)  
library(magrittr) 
library(rgbif) 
library(taxize)


# 提供gbif的账户
user <- "xueqin" 
pwd <- "961121wxq"   
email <- "335015647@qq.com" 

# 读入物种名录
file_url <- read.csv("rosa_species_subfam_inf.csv", header=TRUE)

# 获取每个物种的taxon_key号码，用于后续下载
gbif_taxon_keys <- file_url[,3] %>%  # the first column is species names
    taxize::get_gbifid_(method="backbone") %>%  # match names to the GBIF backbone to get taxonkeys
    imap(~ .x %>% mutate(original_sciname = .y)) %>%  # add original name back into data.frame
    bind_rows() %T>% # combine all data.frames into one
    readr::write_tsv(file = "all_matches.tsv") %>% # save as side effect for you to inspect if you want
    filter(matchtype == "EXACT" & status == "ACCEPTED") %>% # get only accepted and matched names
    filter(kingdom == "Plantae") %>% # remove anything that might have matched to a non-plant
    pull(usagekey)

# 开始下载
occ_download(
    pred_in("taxonKey", gbif_taxon_keys),
    # pred_in("basisOfRecord", c('PRESERVED_SPECIMEN','HUMAN_OBSERVATION','OBSERVATION','MACHINE_OBSERVATION')),
    # pred_gt("elevation", 5000),
    # pred("country", "US"),
    pred("hasCoordinate", TRUE),
    # pred("hasGeospatialIssue", FALSE),
    # pred_gte("year", 1990),
    format = "SIMPLE_CSV",
    user=user,pwd=pwd,email=email
)

# 监听打包进度（控制台会打印当前状态，如 RUNNING，直到变成 SUCCEEDED）
# 0006806-260903145123482为提供的Download key，执行完上一步会获得这个key
occ_download_wait('0010694-260903145123482')

# 将压缩包下载到本地并读入为数据框 d
d <- occ_download_get('0010694-260903145123482', path = "./") %>%
  occ_download_import()

# 或者直接在服务器中使用aria2c下载
aria2c -s 16 -x 16 "https://api.gbif.org/v1/occurrence/download/request/0006806-260903145123482.zip"

# 查看前几行数据
head(d)

# 下载到Linux中后，过滤出特定的行和列
# https://dwc.tdwg.org/terms/#dwc:establishmentMeans 在这里查看每列的含义

# 过滤规则
# 需要同时拥有物种名和经纬度
# 只保留标本馆标本数据记录，剔除一般目击证据，化石等
# 限制了分类单元必须为SPECIES（种）、VARIETY（变种） 或 SUBSPECIES（亚种），摸棱两可的数据去掉
# 排除外来引入和人工栽培记录
# issue列排除模糊匹配与强制升阶元，即过滤掉鉴定不出种，直接归并到属或者科的物种分布信息

awk -F'\t' -v OFS='\t' '
    NR==1 {
        for (i=1; i<=NF; i++) col[$i] = i;
        # 打印 17 列新表头
        print "species", "taxonRank", "decimalLongitude", "decimalLatitude", "coordinateUncertaintyInMeters", \
              "year", "basisOfRecord", "establishmentMeans", "identifiedBy", "issue", \
              "institutionCode", "catalogNumber", "countryCode", "stateProvince", "locality", \
              "elevation", "elevationAccuracy";
        next;
    }
    {
        sp   = $col["species"];
        lon  = $col["decimalLongitude"];
        lat  = $col["decimalLatitude"];
        rec  = $col["basisOfRecord"];
        rank = $col["taxonRank"];
        est  = $col["establishmentMeans"];
        iss  = $col["issue"];

        # 严格质控过滤规则
        if (sp != "" && lon != "" && lat != "" && \
            (rank == "SPECIES" || rank == "VARIETY" || rank == "SUBSPECIES") && \
            est != "introduced" && \
            iss !~ /TAXON_MATCH_FUZZY/ && iss !~ /TAXON_MATCH_HIGHERRANK/) {
            print sp, rank, lon, lat, $col["coordinateUncertaintyInMeters"], \
                  $col["year"], rec, est, $col["identifiedBy"], iss, \
                  $col["institutionCode"], $col["catalogNumber"], $col["countryCode"], \
                  $col["stateProvince"], $col["locality"], $col["elevation"], $col["elevationAccuracy"];
        }
    }' 0006806-260903145123482.csv > final_all_data.tsv

### 查询固定的列中，不同类型的数据占比==========================
python3 -c "
	import pandas as pd
	from collections import Counter

	input_file = 'final_all_data.tsv'
	# 在这里填入你想探查统计的所有列名
	inspect_cols = ['establishmentMeans', 'basisOfRecord', 'taxonRank','issue']

	counters = {c: Counter() for c in inspect_cols}
	total_rows = 0

	# 分块流式读取（内存占用极低，绝不爆内存）
	for chunk in pd.read_csv(input_file, sep='\t', usecols=inspect_cols, chunksize=100000, low_memory=False):
	    total_rows += len(chunk)
	    for c in inspect_cols:
	        vals = chunk[c].fillna('[EMPTY/NA]').astype(str)
	        counters[c].update(vals)

	print(f'=== 数据集总记录行数: {total_rows} ===\n')
	for c in inspect_cols:
	    print(f'【 列名: {c} 】')
	    for val, count in counters[c].most_common():
	        pct = (count / total_rows) * 100
	        print(f'  - {val:<30s} : {count:8d} 条 ({pct:6.2f}%)')
	    print()
	" > summary_inf.txt


# 提取出关键信息列：
awk -F'\t' -v OFS='\t' 'NR==1{for(i=1;i<=NF;i++)a[$i]=i} {print $a["species"], $a["stateProvince"], $a["countryCode"], $a["decimalLongitude"], $a["decimalLatitude"]}' final_all_data.tsv > final_all_data_clean1.tsv

######################################################################################################################


###################################################### idigbio #######################################################
## 下载idigbio数据
library("ridigbio")
species_list <- read.csv("rosa_species_subfam_inf.csv", header=TRUE)
species_names <- species_list$accepted_name
all_results <- list()  # 初始化存储结果的列表 

already_down <- data.frame(genus=names(all_results))

remaining_sp <- species_names[!species_list$accepted_name %in% already_down$genus]


for (species in remaining_sp) {
  rq <- list(scientificname = species)  # 创建查询条件
  results <- idig_search_records(rq)    # 获取结果
  all_results[[species]] <- results     # 将结果存储在列表中
}


combined_results <- do.call(rbind, all_results)  
#过滤一下经纬度为空的数值
final_idig <- combined_results %>% 
  filter(combined_results$geopoint.lon!="NA" &combined_results$geopoint.lat!="NA")

#清理数据
data_problem <- clean_coordinates(x = final_idig,   
                                  lon = "geopoint.lon", 
                                  lat = "geopoint.lat",
                                  species = "scientificname",
                                  tests = c("capitals", "centroids", "equal", "gbif", "institutions", "outliers",
                                            "seas", "zeros")) 


summary(data_problem)
data_problematic <- data[which(data_problem$.summary== "FALSE"),]
distribution_all_clean <- data[which(data_problem$.summary== "TRUE"),]

write.csv(final_idig,file="idigbio_rosaceae_new_add.csv")

# 提取出关键信息列
awk -F',' -v OFS=',' 'NR==1{for(i=1;i<=NF;i++){gsub(/"/,"",$i); a[$i]=i}} {print $a["scientificname"], $a["country"], $a["data.dwc:eventDate"], $a["geopoint.lon"], $a["geopoint.lat"]}' idigbio_rosaceae_new_add.csv > idigbio_clean1.csv

#################################################CVH NSII###############################################
#下载CVH NSII数据
# 相关脚本调用位置在：/data/xiongtao/project/nsii_dist_test
mkdir result_nsii
while read -r line ;do
  echo y|python RetrieveHerbarChina.py -w nsii -g $line -o result -d 0
  echo y|python GeoTransfer.py -f result/${line}_nsii_list_raw.csv -d 0
  cp result/${line}_nsii.csv result_nsii
done<need_to_down_again.txt

mkdir result_cvh
while read -r line ;do
  python RetrieveHerbarChina.py -w cvh -g $line -o result -d 0 
  python RetrieveHerbarChina.py -w cvh -f result/${line}_cvh_list.txt -d 0
  cp result/${line}_cvh.csv result_cvh
done<cvh_list.txt

echo "Species,CN_Location,Count,Longitude,Latitude" >> all_nsii_distribution_data.csv
for name in `ls result_nsii/*_nsii.csv`;do 
  sed 1d $name >> all_nsii_distribution_data.csv
done