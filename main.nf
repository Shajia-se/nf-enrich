#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process functional_enrich {
  tag "functional_enrich"
  stageInMode 'symlink'
  stageOutMode 'move'

  publishDir "${params.project_folder}/${params.enrich_output}", mode: 'copy', overwrite: true

  input:
    path(manifest_tsv)

  output:
    path("functional_enrich_results")

  script:
  """
  set -euo pipefail
  mkdir -p tmp
  export TMPDIR=\$PWD/tmp
  export TEMP=\$PWD/tmp
  export TMP=\$PWD/tmp
  export MPLCONFIGDIR=\$PWD/tmp/mpl
  mkdir -p "\$MPLCONFIGDIR"

  cat > run.R <<'RS'
  suppressPackageStartupMessages({
    library(clusterProfiler)
    library(org.Mm.eg.db)
    library(ChIPpeakAnno)
    library(GenomicRanges)
    library(openxlsx)
    library(ggplot2)
  })

  manifest <- read.delim("${manifest_tsv}", as.is=TRUE, check.names=FALSE)
  dir.create("functional_enrich_results", showWarnings=FALSE)

  read_peak_gr <- function(path) {
    d <- read.table(path, sep="\\t", header=FALSE, as.is=TRUE, comment.char="", quote="")
    if (ncol(d) < 3) stop("Peak file must have at least 3 columns: ", path)
    GRanges(seqnames=d[[1]], ranges=IRanges(as.integer(d[[2]]) + 1L, as.integer(d[[3]])))
  }

  extract_genes <- function(path) {
    d <- read.delim(path, as.is=TRUE, check.names=FALSE)
    if ("ENTREZID" %in% colnames(d)) {
      vals <- d$ENTREZID
    } else if ("geneId" %in% colnames(d)) {
      vals <- d$geneId
    } else {
      vals <- character()
    }
    vals <- unique(unlist(strsplit(as.character(vals), "/")))
    vals <- vals[!is.na(vals) & vals != ""]
    vals
  }

  write_table_xlsx <- function(df, prefix) {
    write.table(df, paste0(prefix, ".tsv"), sep="\\t", quote=FALSE, row.names=FALSE)
    openxlsx::write.xlsx(df, paste0(prefix, ".xlsx"), rowNames=FALSE)
  }

  safe_plot <- function(expr, file, width=8, height=5) {
    pdf(file, width=width, height=height)
    ok <- TRUE
    tryCatch(force(expr), error=function(e) { plot.new(); title(conditionMessage(e)); ok <<- FALSE })
    dev.off()
    ok
  }

  run_go <- function(genes, universe, ont) {
    enrichGO(
      gene = genes,
      universe = universe,
      OrgDb = org.Mm.eg.db,
      keyType = "ENTREZID",
      ont = ont,
      pAdjustMethod = "BH",
      qvalueCutoff = ${params.qvalue_cutoff},
      readable = TRUE
    )
  }

  run_kegg <- function(genes, universe) {
    enrichKEGG(
      gene = genes,
      universe = universe,
      organism = "${params.kegg_organism}",
      pAdjustMethod = "BH",
      qvalueCutoff = ${params.qvalue_cutoff}
    )
  }

  by_source <- split(manifest, manifest$source)

  for (source in names(by_source)) {
    S <- by_source[[source]]
    source_dir <- file.path("functional_enrich_results", source)
    dir.create(source_dir, showWarnings=FALSE, recursive=TRUE)

    genes_by_sample <- lapply(S$annotated_tsv, extract_genes)
    names(genes_by_sample) <- S$sample
    universe <- sort(unique(unlist(genes_by_sample)))

    if (length(universe) == 0) next

    for (i in seq_len(nrow(S))) {
      sample <- S$sample[i]
      sample_dir <- file.path(source_dir, sample)
      dir.create(sample_dir, showWarnings=FALSE, recursive=TRUE)

      genes <- genes_by_sample[[sample]]
      if (length(genes) == 0) next

      for (ont in c("BP","MF","CC")) {
        ego <- tryCatch(run_go(genes, universe, ont), error=function(e) NULL)
        if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
          df <- as.data.frame(ego)
          write_table_xlsx(df, file.path(sample_dir, paste0("GO_", ont)))
          safe_plot(print(dotplot(ego, showCategory=${params.show_category}) + ggtitle(paste(sample, source, ont))), file.path(sample_dir, paste0("GO_", ont, ".pdf")), 9, 5)
        }
      }

      ekegg <- tryCatch(run_kegg(genes, universe), error=function(e) NULL)
      if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
        df <- as.data.frame(ekegg)
        write_table_xlsx(df, file.path(sample_dir, "KEGG"))
        safe_plot(print(dotplot(ekegg, showCategory=${params.show_category}) + ggtitle(paste(sample, source, "KEGG"))), file.path(sample_dir, "KEGG.pdf"), 9, 5)
      }
    }

    if (length(genes_by_sample) >= 2) {
      comp_dir <- file.path(source_dir, "compareCluster")
      dir.create(comp_dir, showWarnings=FALSE, recursive=TRUE)

      comp_bp <- tryCatch(compareCluster(
        geneCluster = genes_by_sample,
        fun = "enrichGO",
        universe = universe,
        OrgDb = org.Mm.eg.db,
        keyType = "ENTREZID",
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = ${params.pvalue_cutoff}
      ), error=function(e) NULL)

      if (!is.null(comp_bp) && nrow(as.data.frame(comp_bp)) > 0) {
        df <- as.data.frame(comp_bp)
        write_table_xlsx(df, file.path(comp_dir, "compareCluster_GO_BP"))
        safe_plot(print(dotplot(comp_bp, showCategory=${params.show_category}) + ggtitle(paste(source, "GO BP compareCluster"))), file.path(comp_dir, "compareCluster_GO_BP.pdf"), 10, 6)
      }

      comp_kegg <- tryCatch(compareCluster(
        geneCluster = genes_by_sample,
        fun = "enrichKEGG",
        organism = "${params.kegg_organism}",
        pAdjustMethod = "BH",
        pvalueCutoff = ${params.pvalue_cutoff}
      ), error=function(e) NULL)

      if (!is.null(comp_kegg) && nrow(as.data.frame(comp_kegg)) > 0) {
        df <- as.data.frame(comp_kegg)
        write_table_xlsx(df, file.path(comp_dir, "compareCluster_KEGG"))
        safe_plot(print(dotplot(comp_kegg, showCategory=${params.show_category}) + ggtitle(paste(source, "KEGG compareCluster"))), file.path(comp_dir, "compareCluster_KEGG.pdf"), 10, 6)
      }
    }

    if (nrow(S) == 2) {
      venn_dir <- file.path(source_dir, "overlap")
      dir.create(venn_dir, showWarnings=FALSE, recursive=TRUE)

      gr1 <- read_peak_gr(S$peak_file[1])
      gr2 <- read_peak_gr(S$peak_file[2])
      overlap <- tryCatch(findOverlapsOfPeaks(gr1, gr2, maxgap=${params.venn_maxgap}), error=function(e) NULL)
      if (!is.null(overlap)) {
        pdf(file.path(venn_dir, paste0(source, "_venn.pdf")), width=8, height=6)
        tryCatch({
          makeVennDiagram(
            overlap,
            fill = c("${params.venn_color_1}", "${params.venn_color_2}"),
            col = "black",
            cat.col = "black",
            cat.cex = 1.1,
            cex = 1.2,
            lwd = 2,
            alpha = 0.7,
            main = paste("Overlap -", source)
          )
        }, error=function(e) { plot.new(); title(conditionMessage(e)) })
        dev.off()
      }
    }
  }
  RS

  Rscript run.R
  """
}

workflow {
  def peakSources = (params.enrich_peak_sources ?: 'idr,consensus,diffbind')
    .toString()
    .split(',')
    *.trim()
    .findAll { it }
    .unique()

  def rows = []

  def resolveAnnotated = { source, sample ->
    def primary = file("${params.chipseeker_output}/${sample}/annotated_peaks.${sample}.tsv")
    if (primary.exists()) return primary

    if (source == 'idr' && sample.startsWith('idr__')) {
      def legacy = sample.replaceFirst(/^idr__/, '')
      def fallback = file("${params.chipseeker_output}/${legacy}/annotated_peaks.${legacy}.tsv")
      if (fallback.exists()) return fallback
    }
    return null
  }

  def addRows = { source, baseDir, pattern, sampleFn ->
    def dir = file(baseDir)
    assert dir.exists() : "${source} directory not found: ${baseDir}"
    def files = dir.listFiles()?.findAll { f -> f.isFile() && f.name ==~ globToRegex(pattern) }?.sort { it.name } ?: []
    files.each { f ->
      def sample = sampleFn(f)
      def annotated = resolveAnnotated(source, sample)
      if (annotated != null) {
        rows << [source, sample, annotated, file(f.toString())]
      }
    }
  }

  if (peakSources.contains('idr')) {
    addRows('idr', params.idr_output, params.idr_peak_pattern ?: "*_idr.sorted.chr.narrowPeak") { f ->
      "idr__" + f.baseName
        .replaceFirst(/_idr\.sorted\.chr$/, '')
        .replaceFirst(/_idr\.sorted$/, '')
        .replaceFirst(/\.narrowPeak$/, '')
    }
  }

  if (peakSources.contains('consensus')) {
    addRows('consensus', params.peak_consensus_output, params.consensus_peak_pattern ?: "*_consensus.bed") { f ->
      "consensus__${f.baseName}"
    }
  }

  if (peakSources.contains('diffbind')) {
    addRows('diffbind', params.diffbind_output, params.diffbind_peak_pattern ?: "*.bed") { f ->
      "diffbind__${f.baseName}"
    }
  }

  assert !rows.isEmpty() : "No enrich input rows found. Check chipseeker_output and enrich_peak_sources."

  Channel
    .fromList(['source\tsample\tannotated_tsv\tpeak_file'] + rows.collect { r ->
      "${r[0]}\t${r[1]}\t${r[2]}\t${r[3]}"
    })
    .collectFile(name: 'functional_enrich_manifest.tsv', newLine: true)
    .set { manifest_ch }

  functional_enrich(manifest_ch)
}

def globToRegex(pattern) {
  '^' + pattern
    .replace('.', '\\.')
    .replace('*', '.*')
    .replace('?', '.') + '$'
}
