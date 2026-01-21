#!/usr/bin/env Rscript

### calculate diversity ###

# install required packages
required_pkg <- c("optparse", "ape", "rbiom", "compositions", "BiocManager")
a <- sapply(required_pkg, function(x) {
  if (!requireNamespace(x, quietly = TRUE))
    install.packages(x, repos = "http://cran.us.r-project.org")
})
if (!"microbiome" %in% installed.packages()) {
  BiocManager::install("microbiome")
}

library("optparse")

option_list = list(
  make_option(c("-f", "--file"), action="store", type="character", default=NULL, 
              help="Merged MetaPhlAn profiles. Samples as columns and species as rows required.",
              metavar="character"),
  make_option(c("-o", "--out_directory"), action="store", type="character", default="diversity_analysis",
              help="Output directory. [default = %default]"),
  make_option(c("-p", "--outfile_prefix"), action="store", type="character", default=NULL,
              help="Prefix for output files. [default = input file basename]"),
  make_option(c("-t", "--tree"), action="store", type="character", default=NULL, 
              help="Path to MetaPhlAn species Newick tree (mandatory for UniFrac).")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$file)) {
  print_help(opt_parser)
  stop('Input file must be supplied.', call.=FALSE)
}

if (!file.exists(opt$file)) {
  stop(paste0('Input file "', opt$file, '" does not exist!'), call.=FALSE)
}

if (!is.null(opt$tree) && !file.exists(opt$tree)) {
  stop(paste0('Tree file "', opt$tree, '" does not exist!'), call.=FALSE)
}

outfile_prefix <- ifelse(is.null(opt$outfile_prefix), tools::file_path_sans_ext(basename(opt$file)), opt$outfile_prefix)

current_dir <- getwd()
dir.create(file.path(current_dir, opt$out_directory), showWarnings = FALSE, recursive = TRUE)
outfile_base <- file.path(current_dir, opt$out_directory, outfile_prefix)

### Load and preprocess data ###
mpa_table <- read.table(opt$file, comment.char = '#', sep = '\t', header = TRUE, check.names=FALSE)

# Remove NCBI ID columns if present
vec <- grepl("ncbi", colnames(mpa_table), ignore.case=TRUE)
if (any(vec)) {
  mpa_table <- mpa_table[grep("t__|s__", mpa_table[,1]), !vec]
} else {
  mpa_table <- mpa_table[grep("t__|s__", mpa_table[,1]), ]
}

# Clean taxon names
taxon_sep <- ifelse(grepl("t__", mpa_table[1,1]), "t__", "s__")
if (taxon_sep == "t__") {
  mpa_table[,1] <- gsub(".+\\|t__SGB", "", mpa_table[,1])
} else {
  mpa_table[,1] <- gsub(".+\\|s__", "", mpa_table[,1])
}
mpa_table[,1] <- gsub("_group$", "", mpa_table[,1])

rownames(mpa_table) <- mpa_table[,1]
mpa_table <- mpa_table[, -1]

# Remove samples with all zeros
removed <- which(colSums(mpa_table) == 0)
if (length(removed) > 0) {
  message <- paste0("# WARNING: ", length(removed), " samples with 100% unknown species were removed.")
  write(message, stdout())
  write(paste(names(removed), collapse='\n'), stdout())
  write(message, file=paste0(outfile_base, '_samples.log'))
  write(paste(names(removed), collapse='\n'), file=paste0(outfile_base, '_samples.log'), append = TRUE)
  mpa_table <- mpa_table[, -removed]
}

# Normalize to relative abundance (sum to 1)
mpa_table <- mpa_table / 100

### ALPHA DIVERSITY ###
library(microbiome)

alpha_metrics <- c("richness_observed", "diversity_shannon", "diversity_gini_simpson", "dominance_gini")
alpha_names <- c("richness", "shannon", "simpson", "gini")

alpha_res <- sapply(alpha_metrics, function(metric) {
  alpha_val <- microbiome::alpha(mpa_table, index = metric)
  rownames(alpha_val) <- colnames(mpa_table)
  alpha_val
}, simplify=FALSE)

alpha_df <- do.call(cbind, alpha_res)
colnames(alpha_df) <- alpha_names

write.table(alpha_df, paste0(outfile_base, "_alpha_diversity.tsv"), sep="\t", quote=FALSE, col.names=NA)

### BETA DIVERSITY ###
library(rbiom)
library(compositions)
library(ape)
library(vegan)

beta_metrics <- list(
  bray_curtis = function(x) as.matrix(vegan::vegdist(t(x), method = "bray")),
  jaccard = function(x) as.matrix(vegan::vegdist(t(x), method = "jaccard")),
  weighted_unifrac = function(x, tree) {
    if (is.null(tree)) {
      warning("Tree file required for weighted_unifrac. Skipping.")
      return(NULL)
    }
    rbiom::bdiv_distmat(as.matrix(x), tree=tree, method="unifrac", weighted=TRUE)
  },
  unweighted_unifrac = function(x, tree) {
    if (is.null(tree)) {
      warning("Tree file required for unweighted_unifrac. Skipping.")
      return(NULL)
    }
    rbiom::bdiv_distmat(as.matrix(x), tree=tree, method="unifrac", weighted=FALSE)
  },
  clr = function(x) {
    clr_table <- compositions::clr(t(x))
    dist_mat <- dist(clr_table, method = "euclidean")
    as.matrix(dist_mat)
  },
  aitchison = function(x) {
    clr_table <- compositions::clr(t(x))
    dist_mat <- dist(clr_table, method = "euclidean")
    as.matrix(dist_mat)
  }
)

# Read tree if provided
mpa_tree <- NULL
if (!is.null(opt$tree)) {
  mpa_tree <- ape::read.tree(opt$tree)
  if (taxon_sep == "s__") {
    mpa_tree$tip.label <- gsub(".+\\|s__", "", mpa_tree$tip.label)
  }
  # Filter species not in tree
  removed_species <- setdiff(rownames(mpa_table), mpa_tree$tip.label)
  if (length(removed_species) > 0) {
    message <- paste0("# WARNING: ", length(removed_species), " species not in tree were removed.")
    write(message, stdout())
    write(paste(removed_species, collapse='\n'), stdout())
    write(message, file=paste0(outfile_base, '_species.log'))
    write(paste(removed_species, collapse='\n'), file=paste0(outfile_base, '_species.log'), append=TRUE)
    mpa_table <- mpa_table[setdiff(rownames(mpa_table), removed_species), ]
  }
}

beta_results <- list()
for (metric_name in names(beta_metrics)) {
  message("Calculating beta diversity: ", metric_name)
  fun <- beta_metrics[[metric_name]]
  if (metric_name %in% c("weighted_unifrac", "unweighted_unifrac")) {
    mat <- tryCatch(fun(mpa_table, mpa_tree), error = function(e) {
      warning(paste("Skipping", metric_name, "due to error:", e$message))
      return(NULL)
    })
  } else {
    mat <- tryCatch(fun(mpa_table), error = function(e) {
      warning(paste("Skipping", metric_name, "due to error:", e$message))
      return(NULL)
    })
  }
  if (!is.null(mat)) {
    beta_results[[metric_name]] <- as.matrix(mat)
  }
}

# Write beta diversity distance matrices to files
for (metric_name in names(beta_results)) {
  write.table(beta_results[[metric_name]], paste0(outfile_base, "_beta_", metric_name, ".tsv"), sep="\t", quote=FALSE, col.names=NA)
}

message("Diversity calculation finished.")
