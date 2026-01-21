library(dplyr)
library(data.table)


pval_check <- function(path, key) {
  # fread로 파일 읽기
  data <- fread(path, sep = "\t", header = TRUE, check.names = FALSE, data.table = FALSE)
  
  # 첫 번째 열이 rownames일 경우 처리
  rownames(data) <- data[[1]]
  data[[1]] <- NULL
  
  number <- nrow(data)
  count <- data.frame(row.names = colnames(data))
  
  # 각 열별 계산
  for (i in colnames(data)) {
    data_ok <- sum(data[[i]] < 0.05, na.rm = TRUE)
    count[i, "True"] <- data_ok
    count[i, "False"] <- number - data_ok
    count[i, "True_rate"] <- data_ok / number * 100
    count[i, "False_rate"] <- (number - data_ok) / number * 100
    count[i, key] <- data_ok / number * 100
    count[i, "check"] <- (data_ok / number * 100) > 99.5
  }
  
  return(count)
}

# 사용 예시 ----
minfi <- pval_check("Minfi_detection-P_all_CpGs.txt", "Minfi")

# 저장 ----
# write.table(sesame, "/disk0/sm/methyl/03_pollution/total/data/pval_sesame.txt", sep = "\t", quote = FALSE)
write.table(minfi, "data/pval_minfi.txt", sep = "\t", quote = FALSE)
