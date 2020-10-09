suppressPackageStartupMessages({
  library(rjson)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(adegenet)
  library(phangorn)
  library(phytools)
})

dna <- read.phyDat(snakemake@input[["msa"]], format = "fasta")
tree <- phytools::force.ultrametric(read.tree(snakemake@input[["tree"]]))
mod <- modelTest(dna, unroot(tree), multicore = T, mc.cores = snakemake@params[[1]])
mod <- mod[which.min(mod$BIC), ]
fit <- optim.pml(pml(tree, dna), optRoptNni = T, optRooted = T, model = mod$model)
anc <- as.character(ancestral.pml(fit, type = "bayes", return = "phyDat"))

rownames(anc) <- c(tree$tip.label, tree$node.label)

path <- snakemake@output[["seqs"]]
write.dna(anc, path, format = "fasta")

ref <- tail(tree$tip.label, 1)

fasta2DNAbin(path, snpOnly = T) %>%
  as.character() %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  pivot_longer(., cols = tail(colnames(.), -1)) %>%
  left_join(., filter(., rowname == ref), by = "name") ->
snp

s <- split(snp, snp$rowname.x)
s <- setNames(lapply(seq_along(s), function(idx) {
  x <- s[[idx]]
  k <- names(s)[[idx]]
  muts <- with(filter(x, value.y != value.x), toupper(paste0(value.y, name, value.x)))
  list(
    muts = (if(length(muts) == 1) list(muts) else muts),
    sequence = toupper(paste0(anc[which(rownames(anc) == k), ], collapse = ""))
  )
}), names(s))
write(rjson::toJSON(
  list(
    nodes = s,
    reference = list(
      nuc = toupper(paste0(anc[which(rownames(anc) == ref), ], collapse = ""))
    )
  ),
  indent = 1
), snakemake@output[["json"]])
