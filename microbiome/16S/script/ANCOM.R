# if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
library(ANCOMBC)
library(phyloseq)
library(tidyverse)

# ==============================================================================
# 2. 설정
# ==============================================================================
# 분석할 그룹 변수 설정
args <- commandArgs(trailingOnly = TRUE)
if(length(args) > 0) {
  GROUP_COLUMNS <- args
} else {
  GROUP_COLUMNS <- c('group', 'day') # 기본값
}

# 분석할 분류 단계
TAX_LEVELS_TARGET <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# 결과 저장 폴더 생성
RESULT_DIR <- "ANCOM"
if(!dir.exists(RESULT_DIR)) dir.create(RESULT_DIR)

# ==============================================================================
# 3. 데이터 로드 및 전처리
# ==============================================================================
cat("=== 데이터 로딩 중 ===\n")

# 3-1. Feature Table
x <- read.delim("filtered/txt/denoise_noMT_table.tsv", header = TRUE, comment.char = "", row.names = 1, check.names = FALSE, skip = 1)
otu_mat <- as.matrix(x)
class(otu_mat) <- "numeric"

# 3-2. Metadata
meta_data <- read.table("/work/meta_data.tsv", sep="\t", header=TRUE, check.names=FALSE, comment.char="", row.names = 1)

# 3-3. Taxonomy
tax_df <- read.delim("taxonomy/taxonomy.tsv", header=TRUE, row.names=1, check.names=FALSE)

# Taxonomy 파싱
if (ncol(tax_df) < 6 && any(grepl(";", tax_df[,1]))) {
  cat(" -> 분류 정보 파싱 중...\n")
  tax_sep <- tax_df %>%
    separate(col = 1, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), 
             sep = ";", fill = "right", extra = "drop")
  tax_mat <- as.matrix(tax_sep)
} else {
  tax_mat <- as.matrix(tax_df)
}

# 컬럼명 강제 지정
ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
if(ncol(tax_mat) <= length(ranks)) {
  colnames(tax_mat) <- ranks[1:ncol(tax_mat)]
} else {
  colnames(tax_mat)[1:length(ranks)] <- ranks
}

tax_mat[is.na(tax_mat)] <- "Unassigned"
tax_mat[tax_mat == ""] <- "Unassigned"

# 3-4. Phyloseq 객체 생성
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)
SAMP <- sample_data(meta_data)
TAX <- tax_table(tax_mat)

ps <- phyloseq(OTU, SAMP, TAX)
ps <- prune_taxa(taxa_sums(ps) > 0, ps)

available_ranks <- rank_names(ps)
TAX_LEVELS <- intersect(TAX_LEVELS_TARGET, available_ranks)

cat(" -> Phyloseq 객체 생성 완료. \n")
cat(" -> 분석 예정 단계:", paste(TAX_LEVELS, collapse=", "), "\n")

# [이름 정리 함수]
clean_taxon_name <- function(x) {
  sapply(x, function(name) {
    if(is.na(name) || name == "" || name == "Unassigned") return(name)
    clean_name <- sub("^[kpcofgs]__", "", name)
    clean_name <- gsub("_", " ", clean_name)
    return(clean_name)
  })
}

# ==============================================================================
# 4. 분석 및 CSV 저장 루프
# ==============================================================================
set.seed(123)

for(group_col in GROUP_COLUMNS){
  cat("\n==========================================================\n")
  cat(" 처리 중인 그룹 컬럼:", group_col, "\n")
  cat("==========================================================\n")
  
  # 데이터 준비
  ps_group <- ps
  sample_data(ps_group)[[group_col]] <- as.factor(sample_data(ps_group)[[group_col]])
  ps_group <- subset_samples(ps_group, !is.na(sample_data(ps_group)[[group_col]]))
  
  all_levels_results <- list()
  
  for(tax_level in TAX_LEVELS) {
    cat(sprintf("  >> 분류 단계: %s ... ", tax_level))
    
    tryCatch({
      # 1. ANCOM-BC2 실행
      res <- ancombc2(
        data = ps_group,
        tax_level = tax_level,
        fix_formula = group_col,
        p_adj_method = "BH",
        prv_cut = 0.10, lib_cut = 1000,
        struc_zero = TRUE, neg_lb = TRUE, pseudo_sens = TRUE,
        group = group_col, global = TRUE, pairwise = TRUE,
        alpha = 0.05, verbose = FALSE 
      )
      
      out <- res$res
      
      # 2. 데이터 가공
      if(nrow(out) > 0) {
        # -------------------------------------------------------
        # (A) 통계값 추출 (LFC, P, Q)
        # -------------------------------------------------------
        stats_df <- out %>% 
          select(taxon, starts_with("lfc_"), starts_with("p_"), starts_with("q_")) %>%
          select(-contains("Intercept")) 

        stats_long <- stats_df %>%
          pivot_longer(cols = -taxon, 
                       names_to = c("metric", "group"), 
                       names_pattern = "^(lfc|p|q)_(.*)",
                       values_to = "value")
        
        # 참조 그룹 복원
        all_groups <- levels(sample_data(ps_group)[[group_col]])
        present_stat_groups <- unique(stats_long$group)
        expected_stat_groups <- paste0(group_col, all_groups)
        missing_stat_groups <- setdiff(expected_stat_groups, present_stat_groups)
        
        if(length(missing_stat_groups) > 0) {
          ref_rows <- expand.grid(taxon = unique(stats_df$taxon),
                                  metric = c("lfc", "p", "q"),
                                  group = missing_stat_groups,
                                  stringsAsFactors = FALSE) %>%
            mutate(value = case_when(
              metric == "lfc" ~ 0,
              metric == "p" ~ 1,
              metric == "q" ~ 1
            ))
          stats_long <- bind_rows(stats_long, ref_rows)
        }
        
        stats_long$group_clean <- gsub(paste0("^", group_col), "", stats_long$group)

        # Wide 포맷 변환 (통계값)
        stats_wide <- stats_long %>%
          mutate(col_name = paste(metric, group_clean, sep = "_")) %>%
          select(taxon, col_name, value) %>%
          pivot_wider(names_from = col_name, values_from = value)
        
        
        # -------------------------------------------------------
        # (B) 풍부도 데이터 준비 (Mean & Individual)
        # -------------------------------------------------------
        # 해당 레벨로 데이터 병합 (glom) 및 상대 풍부도 변환
        ps_glom <- tax_glom(ps_group, taxrank = tax_level, NArm = FALSE)
        ps_rel <- transform_sample_counts(ps_glom, function(x) x / sum(x))
        melt_df <- psmelt(ps_rel)
        
        # taxon 이름 컬럼 매칭 (phyloseq rank name 확인)
        tax_col_name <- tax_level
        if(!(tax_level %in% colnames(melt_df))) tax_col_name <- "OTU"
        
        # ANCOM 결과에 있는 Taxon만 남기기
        melt_df_filtered <- melt_df %>%
          rename(Taxon = !!sym(tax_col_name)) %>%
          filter(Taxon %in% out$taxon)

        # (B-1) 그룹별 평균 풍부도 (Mean Abundance)
        mean_abund <- melt_df_filtered %>%
          rename(Group = all_of(group_col)) %>%
          group_by(Taxon, Group) %>%
          summarise(Mean = mean(Abundance), .groups = "drop") %>%
          mutate(col_name = paste("mean", Group, sep = "_")) %>%
          select(Taxon, col_name, Mean) %>%
          pivot_wider(names_from = col_name, values_from = Mean)

        # (B-2) 개인별 풍부도 (Individual Abundance)
        # 샘플 ID가 컬럼명이 됩니다. (Sample1, Sample2, ...)
        # 에러 방지: 같은 Taxon/Sample 조합이 있을 경우 합산
        indiv_abund <- melt_df_filtered %>%
          select(Taxon, Sample, Abundance) %>%
          group_by(Taxon, Sample) %>%
          summarise(Abundance = sum(Abundance), .groups = "drop") %>%
          pivot_wider(names_from = Sample, values_from = Abundance, values_fill = 0)

        
        # -------------------------------------------------------
        # (C) 최종 데이터 병합 (Stats + Mean + Individual)
        # -------------------------------------------------------
        final_level_df <- stats_wide %>%
          left_join(mean_abund, by = c("taxon" = "Taxon")) %>%    # 평균 풍부도 결합
          left_join(indiv_abund, by = c("taxon" = "Taxon")) %>%   # 개인 풍부도 결합
          mutate(
            type = tolower(tax_level),
            taxonomy_clean = clean_taxon_name(taxon)
          ) %>%
          # 컬럼 순서 정리: type, taxonomy, taxonomy_name, 통계값..., 평균값..., 샘플값...
          select(type, taxonomy = taxon, taxonomy_name = taxonomy_clean, everything())
        
        all_levels_results[[tax_level]] <- final_level_df
        cat("완료.\n")
        
      } else {
        cat("결과 없음.\n")
      }
      
    }, error = function(e) {
      cat("실패 (", e$message, ")\n")
    })
  }
  
  # -----------------------------------------------------------
  # [최종 통합 CSV 파일 저장]
  # -----------------------------------------------------------
  output_file <- file.path(RESULT_DIR, paste0("all_results_combined_", group_col, ".csv"))
  
  if(length(all_levels_results) > 0) {
    final_df <- bind_rows(all_levels_results)
    final_df[is.na(final_df)] <- 0 
    
    write.csv(final_df, output_file, row.names = FALSE, quote = FALSE)
    cat(" [저장] 통합 파일 (Stats + Mean + Individual):", output_file, "\n")
    
  } else {
    cat("\n[경고] 저장할 결과가 없습니다. 그룹:", group_col, "\n")
  }
}

cat("\n=== 스크립트 실행 완료 ===\n")