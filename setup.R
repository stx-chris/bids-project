library(DESeq2)
library(openxlsx)
library(pheatmap)
library(org.Hs.eg.db)
library(ggplot2)
library(FactoMineR)
library(ggrepel)
library(apeglm)
library(ggrepel)
library(viridis)
library(tidyverse)
library(kableExtra)
library(pheatmap)
library(genefilter)
library(EnhancedVolcano)
library(UpSetR)
library(msigdbr)
library(doMC) # falls es hiermit Probleme gibt könnt ihr das Paket auch weglassen
library(foreach) # falls es hiermit Probleme gibt könnt ihr das Paket auch weglassen
library(BiocParallel) # falls es hiermit Probleme gibt könnt ihr das Paket auch weglassen
# optional -> increase speed of analysis due to multicore processing
register(MulticoreParam(4)) # -> diese Z

packages <- c(
  "doMC",
  "foreach",
  "BiocParallel",
  "UpSetR",
  "genefilter",
  "DESeq2",
  "openxlsx",
  "pheatmap",
  "org.Hs.eg.db",
  "ggplot2",
  "FactoMineR",
  "ggrepel",
  "apeglm",
  "viridis",
  "tidyverse",
  "kableExtra",
  "pheatmap",
  "EnhancedVolcano",
  "msigdbr"
)
BiocManager::install(packages)

mainDir <- "/Users/chris/Library/CloudStorage/GoogleDrive-christoph.kovacs@gmail.com/My Drive/BIDS Bioinformatik & Systembiologie 2025/KW4 Lernaufgabe/MIRACUM_BIDS_Bioinformatik_Systembiologie_RNA_Sequenzierung_Aufgabe"
analysisDir <- file.path(mainDir, "analysis")
degDIR <- file.path(analysisDir, "DEG")
gseaDIR <- file.path(analysisDir, "GSEA")
gageDIR <- file.path(analysisDir, "GSEA", "GAGE")
dir.create(degDIR, recursive = T)
dir.create(gageDIR, recursive = T)

setwd(mainDir)
getwd()

ensembl2entrez <- function(ensembl) {
  entrez <- mget(as.character(ensembl), org.Hs.egENSEMBL2EG, ifnotfound = NA)
  entrez <- lapply(entrez, function(i) return(i[1]))
  return(unlist(entrez))
}

entrez2ensembl <- function(entrez) {
  esbl <- mget(as.character(entrez), org.Hs.egENSEMBL, ifnotfound = NA)
  esbl <- lapply(esbl, function(i) return(i[1]))
  return(unlist(esbl))
}

entrez2symbol <- function(entrez) {
  symbol <- mget(as.character(entrez), org.Hs.egSYMBOL, ifnotfound = NA)
  symbol <- unlist(lapply(symbol, function(i) return(i[1])))
  return(symbol)
}

entrez2genename <- function(entrez) {
  symbol <- mget(as.character(entrez), org.Hs.egGENENAME, ifnotfound = NA)
  symbol <- unlist(lapply(symbol, function(i) return(i[1])))
  return(symbol)
}

getGeneMat <- function(ensIDs) {
  geneMat <- data.frame(ENSEMBL=ensIDs)
  geneMat$ENTREZ <- ensembl2entrez(geneMat$ENSEMBL)
  idxNA <- !is.na(geneMat$ENTREZ)
  sym <- entrez2symbol(na.omit(geneMat$ENTREZ))
  genename <- entrez2genename(na.omit(geneMat$ENTREZ))
  geneMat$Symbol <- NA
  geneMat$Symbol[idxNA] <- sym
  geneMat$Genename <- NA
  geneMat$Genename[idxNA] <- genename
  rownames(geneMat) <- geneMat$ENSEMBL
  return(geneMat)
}

library(readxl)

getwd()

# Annotation einlesen
annotation <- 
  read_excel("targets.xlsx") %>%
  mutate(group = as.factor(group)) %>%
  column_to_rownames("label") %>%
  mutate(label = rownames(.))

read_counts <- function(file) {
  read_tsv(file, skip = 4, col_names = c("GeneID", "Count1", "Count2", "Count3"), show_col_types = FALSE) %>%
    select(GeneID, Count1)
}

# Eine Liste mit Count-Daten je Datei
count_list <- map(annotation$file, read_counts)

# Count-Matrix erstellen
count_matrix <- reduce(count_list, full_join, by = "GeneID") %>%
  setNames(c("GeneID", annotation$label)) %>%
  column_to_rownames("GeneID")

# Matrix speichern
#write.table(
#  x = count_matrix,
#  file = file.path(analysisDir, "counts.txt"),
#  sep = "\t",
#  quote = FALSE
#)

# Gene (ENSEMBL-IDs) als Referenz extrahieren
gene_reference <- getGeneMat(rownames(count_matrix))

### Vorbereiten des DESeq2 Objektes
# Erstellen des DESeq2-Objekts
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData   = annotation,
  design    = ~ group
)

# Setze untreated group
dds$group <- relevel(dds$group, ref = 'DMSO') 

# Vorfilter
keep <- rowSums(counts(dds) >= 10) >= 2
#sum(keep) # keeping 18,198 records
#sum(dim(count_matrix)[1]-sum(keep)) # discarding 39,575 records
dds <- dds[keep,] 

### Normalisierung und Differentielle Expression

# Schätze size factors
dds <- estimateSizeFactors(dds)
#sizeFactors(dds)

# Schätze Dispersion
dds <- estimateDispersions(dds)
#plotDispEsts(dds)

# Führe Wald Test durch
dds <- nbinomWaldTest(dds)
#plotMA(dds)

# rlog-Transformation, behalte Gruppeninformation
rld <- rlog(dds, blind = FALSE)
rlog_matrix <- assay(rld)
#head(rlog_matrix)

### QC Plot vor und nach der Normalisierung

## Roh-Counts extrahieren
raw_counts <- counts(dds)

# Log2 transformieren
log_raw <- log2(raw_counts + 1)

# Daten für ggplot umformen
df_raw <- as.data.frame(log_raw) %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "log2count") %>%
left_join(as.data.frame(colData(dds)) %>%
rownames_to_column("sample"), by = "sample")

# Boxplot erstellen
plot_counts_trans <- df_raw %>%
  ggplot(aes(x = sample, y = log2count, fill = group)) +   # why not group?
  geom_boxplot(outlier.size = 0.3) +
  theme_bw() +
  labs(title = "Counts (trans) pro Bedingung",
       x = "Bedingung", y = "Counts (trans)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## Normalisierte Counts extrahieren
norm_counts <- counts(dds, normalized = TRUE)

# Log2 transformieren
log_norm <- log2(norm_counts + 1)

# Daten für ggplot umformen
df_norm <- as.data.frame(log_norm) %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "log2normCount") %>%
left_join(as.data.frame(colData(dds)) %>%
rownames_to_column("sample"), by = "sample")

# Boxplot erstellen
plot_counts_normtrans <- df_norm %>% 
  ggplot(aes(x = sample, y = log2normCount, fill = group)) +
  geom_boxplot(outlier.size = 0.3) +
  theme_bw() +
  labs(title = "Counts (norm, trans) pro Bedingung",
       x = "Bedingung", y = "Counts (norm, trans)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


### PCA: Principal Component Analysis (linear relationship analysis)
pca <- plotPCA(rld, intgroup = "group", returnData = TRUE)
# Variance explained
percentVar <- round(100 * attr(pca, "percentVar"))

# Plot
plot_pca <- ggplot(pca, aes(PC1, PC2, color = group, label = name)) +
  geom_point(size = 3) +
  geom_text(vjust = -1.2, size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  theme_minimal()
# Interpretation
# ==============
# PCA erklärt 87% der Gesamtvarianz in 2 Dimensionen
# DMSO bildet einen Cluster, die Samples kann als (linear) homogen angenommen werden
# Sorafenib bildet einen Cluster, die Samples können daher als (linear) homogen angenommen werden
# Trametinib bildet zwei Cluster, die Samples sind daher nicht (linear) homogen

### t-SNA: t-distributed stochastic neighbor embedding for non-linear relationship analysis
library(Rtsne)
# Extract assay data (genes x samples) and transpose to samples x genes
mat <- t(assay(rld))

# Optionally: reduce to top variable genes to speed up and denoise
top_var_genes <- head(order(apply(mat, 2, var), decreasing = TRUE), 500)
mat_top <- mat[, top_var_genes]

# Run t-SNE (you can adjust perplexity; try 5–30 depending on sample size)
set.seed(42)  # for reproducibility
tsne_out <- Rtsne(mat_top, dims = 2, perplexity = 1, verbose = TRUE)

# Combine with metadata
tsne_df <- data.frame(
  TSNE1 = tsne_out$Y[,1],
  TSNE2 = tsne_out$Y[,2],
  group = colData(rld)$group  # adjust if your metadata column is named differently
)

# Plot with ggplot2
plot_tsna <- ggplot(tsne_df, aes(x = TSNE1, y = TSNE2, color = group)) +
  geom_point(size = 3) +
  labs(title = "t-SNE of RNA-seq samples", x = "t-SNE 1", y = "t-SNE 2") +
  theme_minimal()
# Interpretation
# ==============
# Ganz deutliche Cluster zu erkennen; DMSO und Tramatinib bilden jeweils einen Cluster, zusammen einen größeren. Daher kann davon ausgegangen, dass Trametinib keinen über die Baseline hinausgehenden qualitativen Unterschied aufweist.
# Sorafenib bildet einen Cluster deutlich abseits der Kontrollgruppe und von Trametinib, das kann als Indikator für einen qualitativen Unterschied gewertet werden



if (!require("ashr", quietly = TRUE))
    install.packages("ashr")

library(ashr)

dds <- DESeq(dds)

# Sorafenib vs DMSO
res_sora <- lfcShrink(dds, contrast = c("group", "Sorafenib", "DMSO"), type = "ashr")

# Trametinib vs DMSO
res_tram <- lfcShrink(dds, contrast = c("group", "Trametinib", "DMSO"), type = "ashr")

# Hilfsfunktion für Filter & Format
filter_deg <- function(res) {
  as.data.frame(res) |>
    rownames_to_column("ENSEMBL") |>
    left_join(gene_reference[, c("ENSEMBL", "Symbol")], by = "ENSEMBL") |>
    filter(padj < 0.05 & !is.na(Symbol)) |>
    select(Symbol, log2FoldChange, padj) |>
    arrange(desc(abs(log2FoldChange)))  # nach Effektgröße sortieren
}

deg_sora <- filter_deg(res_sora)
deg_tram <- filter_deg(res_tram)

## compare up and down-regulated genes
count_regulation <- function(df) {
  df |> 
    mutate(Regulation = case_when(
      log2FoldChange > 0 ~ "Upregulated",
      log2FoldChange < 0 ~ "Downregulated"
    )) |> 
    count(Regulation)
}

sora_counts <- count_regulation(deg_sora)
tram_counts <- count_regulation(deg_tram)
 
# Visualize
sora_counts$Comparison <- "Sorafenib vs DMSO"
tram_counts$Comparison <- "Trametinib vs DMSO"
all_counts <- bind_rows(sora_counts, tram_counts)

plot_deg <- ggplot(all_counts, aes(x = Comparison, y = n, fill = Regulation)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Differenziell regulierte Gene", y = "Anzahl Gene") +
  theme_minimal()

tab_deg_sora <- kable(head(deg_sora, 10), caption = "Top DEGs: Sorafenib vs DMSO")
tab_deg_tram <- kable(head(deg_tram, 10), caption = "Top DEGs: Trametinib vs DMSO")


# **Erstellt jeweils einen Volcano Plot für die beiden Vergleiche.**
if (!requireNamespace("EnhancedVolcano", quietly = TRUE)) {
  BiocManager::install("EnhancedVolcano")
}
library(EnhancedVolcano)

plot_volcano_sora <- EnhancedVolcano(res_sora,
    lab = gene_reference$Symbol[match(rownames(res_sora), gene_reference$ENSEMBL)],
    x = 'log2FoldChange',
    y = 'padj',
    pCutoff = 0.05,
    FCcutoff = 1,
    pointSize = 2.0,
    labSize = 3.0,
    title = 'Sorafenib vs DMSO',
    subtitle = 'Differentiell exprimierte Gene',
    caption = 'log2FC > 1, FDR < 0.05',
    legendLabels = c('NS', 'Log2FC', 'FDR', 'FDR & Log2FC'),
    legendPosition = 'right')

plot_volcano_tram <- EnhancedVolcano(res_tram,
    lab = gene_reference$Symbol[match(rownames(res_tram), gene_reference$ENSEMBL)],
    x = 'log2FoldChange',
    y = 'padj',
    pCutoff = 0.05,
    FCcutoff = 1,
    pointSize = 2.0,
    labSize = 3.0,
    title = 'Trametinib vs DMSO',
    subtitle = 'Differentiell exprimierte Gene',
    caption = 'log2FC > 1, FDR < 0.05',
    legendLabels = c('NS', 'Log2FC', 'FDR', 'FDR & Log2FC'),
    legendPosition = 'right')

if (!requireNamespace("UpSetR", quietly = TRUE)) {
  install.packages("UpSetR")
}
library(UpSetR)

# function to obtain genes per set
upSetSets <- function(sets){
  list_names <- names(sets)
  attach(sets,warn.conflicts = F)
  res <- lapply(1:length(list_names),function(y){
    combinations <- combn(list_names,y)
    res<-as.list(apply(combinations,2,function(x){
      if(length(x)==1){
        p <- setdiff(get(x),unlist(sapply(setdiff(list_names,x),get)))
      }

      else if(length(x) < length(list_names)){
        p <- setdiff(Reduce(intersect,lapply(x,get)),
        Reduce(union,sapply(setdiff(list_names,x),get)))
      }

      else p <- Reduce(intersect,lapply(x,get))

      if(!identical(p,character(0))) p
      else NA
    }))

    if(y==length(list_names)) {
      res[[1]] <- unlist(res);
      res<-res[1]
    }
    names(res) <- apply(combinations,2,paste,collapse="-")
    res
  })
  result <- lapply(res, function(x) x[!is.na(x)])
  result <- unlist(result, recursive = F)
  result <- lapply(result,function(x) data.frame(ID=x))
  return(result)
  detach(sets)
}

### Gemeinsamkeiten

#### Durch beide Medikamente hoch-regulierte Gene

up_sora <- deg_sora |> filter(log2FoldChange > 0) |> pull(Symbol)
up_tram <- deg_tram |> filter(log2FoldChange > 0) |> pull(Symbol)

up_list <- list(
  Sorafenib = up_sora,
  Trametinib = up_tram
)
#### Durch beide Medikamente runter-regulierte Gene

down_sora <- deg_sora |> filter(log2FoldChange < 0) |> pull(Symbol)
down_tram <- deg_tram |> filter(log2FoldChange < 0) |> pull(Symbol)

down_list <- list(
  Sorafenib = down_sora,
  Trametinib = down_tram
)

## Plots
# Hochregulierte Gene
plot_upset_up <- upset(fromList(up_list), 
      order.by = "freq",
      mainbar.y.label = "Gemeinsame hochregulierte Gene",
      sets.x.label = "Gene pro Medikament")

# Runterregulierte Gene
plot_upset_down <- upset(fromList(down_list), 
      order.by = "freq",
      mainbar.y.label = "Gemeinsame runterregulierte Gene",
      sets.x.label = "Gene pro Medikament")

tab_upset_sora_up <- head(upSetSets(up_list)$Sorafenib)  # für hochregulierte Gene
tab_upset_tram_up <- head(upSetSets(up_list)$Trametinib)  # für hochregulierte Gene

tab_upset_sora_down <- head(upSetSets(down_list)$Sorafenib)  # für runterregulierte Gene
tab_upset_tram_down <- head(upSetSets(down_list)$Trametinib)  # für runterregulierte Gene


hyperG <- function(geneSets,DEgenes,universe, cutoff=0.1, mincount=2, parallel=T, adj.P.Val = F,
                   set.size = NULL){
  #' hyperG
  #' 
  #' @description Calculates Fisher's Exact test with the specified genes and the supplied gene-sets.
  #'
  #' @param geneSets list. Gene-Set the calculation is based on, e.g. go.bp
  #' @param DEgenes character vector. Gene IDs used for testing. Same identifiers as used for the gene-sets, e.g. ENTREZ IDs.
  #' @param universe character vector. Universe gene IDs.
  #' @param cutoff numeric. Cutoff used to identify sig. pathways. Default: 0.1.
  #' @param mincount numeric. Consider only pathways which contain at least mincount genes. Default: 2
  #' @param parallel boolean. Use parallel calculation. Default: TRUE
  #' @param adj.P.Val boolean. Use adjusted p-value for significance filtering. Is always calculated.
  #' @param set.size vector. Min and max size of allowed gene-sets. Default min:10 genes and max:500 genes.
  #'  
  #' @return the significant regualted pathways.
  #' @export
  #' @importFrom foreach, doMC
  
  require(foreach)
  require(doMC)
  if(parallel){
    registerDoMC(cores=detectCores())
    cores=detectCores()
  }else{
    cores=1
  }
  if(!is.null(set.size)){
    print('Set Size Limits')
    idx <- lapply(geneSets,function(x){length(x) <= set.size[2] & length(x) >= set.size[1]})
    geneSets <- geneSets[unlist(idx)]
  }
  l <- length(setdiff(universe,DEgenes))
  DElen <- length(DEgenes)
  results <- mclapply(1:length(geneSets), function(i){
    results <- matrix(data=NA,ncol=7,nrow = 1)
    colnames(results) <- c('Term','Count','Size','p-value','adj.P.Val','odds ratio','GeneIDs')
    geneSet <- intersect(universe, geneSets[[i]])
  e <- intersect(DEgenes,geneSet)
    a <- length(e)
    b <- DElen - a
    c <- length(geneSet) - a
    d <- l - c
    contigency.matrix <- cbind(c(a,b),c(c,d))
    res <- fisher.test(contigency.matrix,alternative = 'greater')
    results[1,'Term'] <- names(geneSets)[i]
    results[1,'Count'] <- a
    results[1,'Size'] <- length(geneSets[[i]])
    results[1,'p-value'] <- res$p.value
    results[1,'odds ratio'] <- res$estimate[[1]]
    # find genes annotated in the consensus term
    if(a > 0){
      genes <- intersect(DEgenes,geneSet)
      eid <- genes
      eid <- eid[order(eid)]
      results[1,'GeneIDs'] <- paste(eid,collapse="|")
    }
    return(results)
  }, mc.cores=cores)
    
  results <- as.data.frame(do.call(rbind, results))
  for(i in c(2, 3, 4, 5)){
    results[, i] <- as.numeric(as.character(results[, i]))
  }
  
  if(nrow(results) != 1){
    results <- results[order(results[,'p-value'],decreasing = FALSE),]
  results[,'adj.P.Val'] <- p.adjust(results[,'p-value'], 'BH')
  if(adj.P.Val){
    results <- as.data.frame(subset(results,results[,'adj.P.Val']<=cutoff))
  }else{
    results <- as.data.frame(subset(results,results[,'p-value']<=cutoff))
  }
    results <- as.data.frame(subset(results,results[,'Count']>=mincount))
  }else results <- as.data.frame(results)
  
  return(results)
}

get_geneset_ag <- function(
  species = "Homo sapiens",
  category = NULL,
  subcollection = NULL,
  format = "entrez"
) {
  require(msigdbr)
  db_df <- msigdbr(species = species ,category = category, subcollection = subcollection)
  if(format == "entrez"){
    m_list = db_df %>% split(x = as.character(.$entrez_gene), f = .$gs_name)
    for(idx in 1:length(m_list)){
      m_list[[idx]] <- unique(m_list[[idx]] )
    }
    return(m_list)
  }
  if(format == "symbol"){
    m_list = db_df %>% split(x = .$gene_symbol, f = .$gs_name)
    for(idx in 1:length(m_list)){
      m_list[[idx]] <- unique(m_list[[idx]])
    }
    return(m_list)
  }
  if(format == "df"){
    return(db_df)
  }
}

# load Hallmark gene set
hallmark <- get_geneset_ag(species = "Homo sapiens", category = "H", format = "symbol")
# Extrahiere die ENSEMBL-IDs aus dem dds-Objekt
ensembl_ids <- rownames(dds)

# Verknüpfe mit gene_reference, um die Symbolnamen zu erhalten
universe_df <- tibble(ENSEMBL = ensembl_ids) |>
  left_join(gene_reference[, c("ENSEMBL", "Symbol")], by = "ENSEMBL")

#  Erstelle das Universe
universe <- universe_df$Symbol |> unique() |> na.omit()

top_enrichment_plot <- function(ds, title) {
  top_enrichment <- ds |> 
    slice_min(order_by = `adj.P.Val`, n = 10)

  ggplot(top_enrichment, aes(x = reorder(Term, -`adj.P.Val`), y = -log10(`adj.P.Val`))) +
    geom_col(fill = "steelblue") +
    coord_flip() +
    labs(
      title = title,
      x = "Pathway",
      y = "-log10(adj. p-value)"
    ) +
    theme_minimal()
}

enrichment_up_sora <- hyperG(hallmark, up_sora, universe, cutoff = 0.05, adj.P.Val = TRUE)
plot_enrichment_up_sora <- top_enrichment_plot(enrichment_up_sora, "Top Enriched Hallmark Pathways (Sorafenib UP)")
## Interpretation (KI-unterstützt)
## ==============
## Die Behandlung mit Sorafenib führt zu einer signifikanten Hochregulation von Genen, die mit dem Zelltod (Apoptosis, p53), zellulärem Stress (zB Uv-Antwort), Immunantwort (Interferon Gamma) und Stoffwechselprozessen (Adipogenese, Häm- und Peroxisomenstoffwechsel) assoziiert sind. Dies unterstützt die bekannte cytotoxische und wachstumshemmende Wirkung von Sorafenib auf zellulärer Ebene.

enrichment_up_tram <- hyperG(hallmark, up_tram, universe, cutoff = 0.05, adj.P.Val = TRUE) # empty???
plot_enrichment_up_tram <- top_enrichment_plot(enrichment_up_tram, "Top Enriched Hallmark Pathways (Trametinib UP)")
## Interpretation
## ==============
## ?

enrichment_down_sora <- hyperG(hallmark, down_sora, universe, cutoff = 0.05, adj.P.Val = TRUE)
plot_enrichment_down_sora <- top_enrichment_plot(enrichment_down_sora, "Top Enriched Hallmark Pathways (Sorafenib DOWN)")
## Interpretation (KI-unterstützt)
## ==============
## Sorafenib führt zur signifikanten Herunterregulierung von Genen, die an Zellwachstum (MYC, mTORC1, G2M), Zellzyklusprogression, Entzündungsreaktionen (TNFα/NFκB) und Stressantworten beteiligt sind. Diese Ergebnisse sprechen für eine starke antiproliferative und antiinflammatorische Wirkung des Medikaments, ergänzt durch eine Hemmung zellulärer Überlebenssignale.

enrichment_down_tram <- hyperG(hallmark, down_tram, universe, cutoff = 0.05, adj.P.Val = TRUE)
plot_enrichment_down_tram <- top_enrichment_plot(enrichment_down_tram, "Top Enriched Hallmark Pathways (Trametinib DOWN)")
## Interpretation (KI-unterstützt)
## ==============
## Trametinib führt zu einer ausgeprägten Unterdrückung entzündungsbezogener Signalwege (TNFα/NFκB, IL-2, IL-6, STAT3), Zellüberlebenssignale sowie Komponenten der MAPK-Kaskade (KRAS). Diese Ergebnisse sind kompatibel mit dem Wirkmechanismus von Trametinib als MEK-Inhibitor und deuten auf eine kombinierte hemmende Wirkung auf Proliferation, Entzündung und Stressantwort hin.

save.image(file = "session.RData")