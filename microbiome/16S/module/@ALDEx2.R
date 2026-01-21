library(ALDEx2)
library(compositions)

# 첫 줄 주석은 건너뛰고 두 번째 줄을 header로 읽기
x <- read.delim("filtered/txt/denoise_noMT_table.tsv", header = TRUE, comment.char = "", row.names = 1, check.names = FALSE, skip = 1)

conds <- as.vector(read.table("/work/meta_data.tsv", sep="\t", header=TRUE, check.names=FALSE, comment.char="")$s)


# ALDEx2 실행 (Monte Carlo 샘플 128)
x.clr <- aldex.clr(x, conds, mc.samples=128, denom="all")
x.test <- aldex.ttest(x.clr, paired=FALSE)
x.effect <- aldex.effect(x.clr, CI=TRUE, verbose=TRUE)

# 결과 확인
head(cbind(x.test, x.effect))


sig <- cbind(x.test, x.effect)
sig <- sig[sig$we.eBH < 0.1, ]
sig[order(-sig$effect), ]  # effect size 기준 내림차순





#########################################
x <- read.delim("filtered/txt/denoise_noMT_table.tsv", header = TRUE, comment.char = "", row.names = 1, check.names = FALSE, skip = 1)

conds <- as.vector(read.table("/work/meta_data.tsv", sep="\t", header=TRUE, check.names=FALSE, comment.char="")$group)


# CLR (Centered Log-Ratio) 변환 + Kruskal-Wallis 검정 + BH 보정 + 시각화
# CLR 변환
x.clr.mat <- t(apply(x, 1, function(f) clr(f + 1)))
# +1은 zero count 방지, 행(feature) 기준으로 CLR 적용

# Kruskal-Wallis 적용 (3개 이상 그룹 가능)
kw.p <- apply(x.clr.mat, 1, function(f) kruskal.test(f ~ conds)$p.value)

# BH 보정
kw.p.adj <- p.adjust(kw.p, method="BH")

# 결과 확인
head(kw.p.adj)


sig.features <- names(kw.p.adj)[kw.p.adj < 0.1]
length(sig.features)
sig.features


library(ggplot2)
df <- data.frame(
  t(clr.sig[1, , drop=FALSE]),
  group = conds
)
colnames(df)[1] <- "CLR"
ggplot(df, aes(x=group, y=CLR)) +
  geom_boxplot() +
  geom_jitter(width=0.2) +
  ggtitle(sig.features[1])
# qiime2 composition ancom