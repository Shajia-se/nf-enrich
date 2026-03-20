#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process functional_enrich {
  tag "functional_enrich"
  stageInMode 'symlink'
  stageOutMode 'move'

  publishDir "${params.project_folder}/${params.enrich_output}", mode: 'copy', overwrite: true

  input:
    tuple path(manifest_tsv), path(input_files)

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

  manifest <- read.delim("${manifest_tsv}", as.is=TRUE, check.names=FALSE, quote="", comment.char="")
  # Handle possible UTF-8 BOM and surrounding whitespace in header names
  colnames(manifest) <- trimws(sub("^\\ufeff", "", colnames(manifest)))
  required_cols <- c("source", "sample", "annotated_tsv", "peak_file")
  missing_cols <- setdiff(required_cols, colnames(manifest))
  # Backward compatibility: some older manifest files were produced without header
  if (length(missing_cols) > 0 && ncol(manifest) == 4) {
    colnames(manifest) <- required_cols
    missing_cols <- setdiff(required_cols, colnames(manifest))
  }
  if (length(missing_cols) > 0) {
    stop(
      "Manifest is missing required column(s): ",
      paste(missing_cols, collapse=", "),
      "; found: ",
      paste(colnames(manifest), collapse=", ")
    )
  }
  manifest <- manifest[, required_cols, drop=FALSE]
  manifest <- manifest[complete.cases(manifest[, c("source","sample","annotated_tsv","peak_file")]), , drop=FALSE]
  manifest\$source <- trimws(as.character(manifest\$source))
  manifest\$sample <- trimws(as.character(manifest\$sample))
  manifest\$annotated_tsv <- trimws(as.character(manifest\$annotated_tsv))
  manifest\$peak_file <- trimws(as.character(manifest\$peak_file))

  # Drop accidental header-like rows carried as data
  manifest <- manifest[
    !(manifest\$source == "source" &
      manifest\$sample == "sample" &
      manifest\$annotated_tsv == "annotated_tsv" &
      manifest\$peak_file == "peak_file"),
    ,
    drop=FALSE
  ]

  # Keep only rows whose files can be resolved in work dir (absolute path or staged basename)
  can_resolve <- function(p) {
    file.exists(p) || file.exists(basename(p))
  }
  keep <- vapply(
    seq_len(nrow(manifest)),
    function(i) can_resolve(manifest\$annotated_tsv[i]) && can_resolve(manifest\$peak_file[i]),
    logical(1)
  )
  if (any(!keep)) {
    message("Dropping ", sum(!keep), " manifest rows with unresolved files")
  }
  manifest <- manifest[keep, , drop=FALSE]
  if (nrow(manifest) == 0) {
    stop("Manifest has 0 usable rows after filtering empty entries")
  }
  dir.create("functional_enrich_results", showWarnings=FALSE)
  write.table(manifest, "functional_enrich_results/manifest_used.tsv", sep="\t", quote=FALSE, row.names=FALSE)

  resolve_input_path <- function(path) {
    if (file.exists(path)) return(path)
    b <- basename(path)
    if (file.exists(b)) return(b)
    stop("Cannot open file: ", path, " (staged basename not found: ", b, ")")
  }

  read_peak_gr <- function(path) {
    path <- resolve_input_path(path)
    d <- read.table(path, sep="\\t", header=FALSE, as.is=TRUE, comment.char="", quote="")
    if (ncol(d) < 3) stop("Peak file must have at least 3 columns: ", path)
    GRanges(seqnames=d[[1]], ranges=IRanges(as.integer(d[[2]]) + 1L, as.integer(d[[3]])))
  }

  extract_genes <- function(path) {
    path <- resolve_input_path(path)
    d <- read.delim(path, as.is=TRUE, check.names=FALSE)
    if ("ENTREZID" %in% colnames(d)) {
      vals <- d\$ENTREZID
    } else if ("geneId" %in% colnames(d)) {
      vals <- d\$geneId
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

  pretty_sample <- function(x) {
    x <- sub("^(idr__|consensus__|diffbind__)", "", x)
    x
  }

  clean_term_labels <- function(x) {
    x <- as.character(x)
    x <- gsub(" - Mus musculus (house mouse)", "", x, fixed=TRUE)
    x <- gsub(" [Mus musculus]", "", x, fixed=TRUE)
    x
  }

  clean_enrich_object <- function(obj) {
    if (is.null(obj)) return(obj)
    if ("enrichResult" %in% class(obj)) {
      if (!is.null(obj@result) && "Description" %in% colnames(obj@result)) {
        obj@result[, "Description"] <- clean_term_labels(obj@result[, "Description"])
      }
    } else if ("compareClusterResult" %in% class(obj)) {
      if (!is.null(obj@compareClusterResult) && "Description" %in% colnames(obj@compareClusterResult)) {
        obj@compareClusterResult[, "Description"] <- clean_term_labels(obj@compareClusterResult[, "Description"])
      }
    }
    obj
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

  by_source <- split(manifest, manifest\$source, drop=TRUE)

  for (source in names(by_source)) {
    S <- by_source[[source]]
    source_dir <- file.path("functional_enrich_results", source)
    dir.create(source_dir, showWarnings=FALSE, recursive=TRUE)
    write.table(S, file.path(source_dir, "source_manifest.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

    genes_by_sample <- lapply(S\$annotated_tsv, extract_genes)
    names(genes_by_sample) <- S\$sample
    universe <- sort(unique(unlist(genes_by_sample)))

    if (length(universe) == 0) next

    for (i in seq_len(nrow(S))) {
      sample <- S\$sample[i]
      sample_dir <- file.path(source_dir, sample)
      dir.create(sample_dir, showWarnings=FALSE, recursive=TRUE)

      genes <- genes_by_sample[[sample]]
      if (length(genes) == 0) {
        writeLines("No input genes extracted from annotation table", file.path(sample_dir, "NO_ENRICHMENT.txt"))
        next
      }

      sample_stat <- data.frame(
        sample = sample,
        pretty_sample = pretty_sample(sample),
        n_input_genes = length(genes),
        GO_BP_terms = 0L,
        GO_MF_terms = 0L,
        GO_CC_terms = 0L,
        KEGG_terms = 0L,
        stringsAsFactors = FALSE
      )

      for (ont in c("BP","MF","CC")) {
        ego <- tryCatch(run_go(genes, universe, ont), error=function(e) NULL)
        ego <- clean_enrich_object(ego)
        if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
          df <- as.data.frame(ego)
          if ("Description" %in% colnames(df)) df[, "Description"] <- clean_term_labels(df[, "Description"])
          write_table_xlsx(df, file.path(sample_dir, paste0("GO_", ont)))
          safe_plot(print(dotplot(ego, showCategory=${params.show_category}) + ggtitle(paste(sample, source, ont))), file.path(sample_dir, paste0("GO_", ont, ".pdf")), 9, 5)
          sample_stat[[paste0("GO_", ont, "_terms")]] <- nrow(df)
        }
      }

      ekegg <- tryCatch(run_kegg(genes, universe), error=function(e) NULL)
      ekegg <- clean_enrich_object(ekegg)
      if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
        df <- as.data.frame(ekegg)
        if ("Description" %in% colnames(df)) df[, "Description"] <- clean_term_labels(df[, "Description"])
        write_table_xlsx(df, file.path(sample_dir, "KEGG"))
        safe_plot(print(dotplot(ekegg, showCategory=${params.show_category}) + ggtitle(paste(sample, source, "KEGG"))), file.path(sample_dir, "KEGG.pdf"), 9, 5)
        sample_stat\$KEGG_terms <- nrow(df)
      }

      write.table(sample_stat, file.path(sample_dir, "enrichment_status.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
      if (sum(sample_stat[1, c("GO_BP_terms","GO_MF_terms","GO_CC_terms","KEGG_terms")]) == 0) {
        writeLines("No significant enrichment terms under current cutoff", file.path(sample_dir, "NO_ENRICHMENT.txt"))
      }
    }

    if (length(genes_by_sample) >= 2) {
      comp_dir <- file.path(source_dir, "compareCluster")
      dir.create(comp_dir, showWarnings=FALSE, recursive=TRUE)

      has_comp_go <- FALSE
      for (ont in c("BP","MF","CC")) {
        comp_go <- tryCatch(compareCluster(
          geneCluster = genes_by_sample,
          fun = "enrichGO",
          universe = universe,
          OrgDb = org.Mm.eg.db,
          keyType = "ENTREZID",
          ont = ont,
          pAdjustMethod = "BH",
          pvalueCutoff = ${params.pvalue_cutoff}
        ), error=function(e) NULL)
        comp_go <- clean_enrich_object(comp_go)

        if (!is.null(comp_go) && nrow(as.data.frame(comp_go)) > 0) {
          has_comp_go <- TRUE
          df <- as.data.frame(comp_go)
          if ("Description" %in% colnames(df)) df[, "Description"] <- clean_term_labels(df[, "Description"])
          write_table_xlsx(df, file.path(comp_dir, paste0("compareCluster_GO_", ont)))
          safe_plot(print(dotplot(comp_go, showCategory=${params.show_category}) + ggtitle(paste(source, "GO", ont, "compareCluster"))), file.path(comp_dir, paste0("compareCluster_GO_", ont, ".pdf")), 10, 6)
        }
      }

      comp_kegg <- tryCatch(compareCluster(
        geneCluster = genes_by_sample,
        fun = "enrichKEGG",
        organism = "${params.kegg_organism}",
        pAdjustMethod = "BH",
        pvalueCutoff = ${params.pvalue_cutoff}
      ), error=function(e) NULL)
      comp_kegg <- clean_enrich_object(comp_kegg)

      if (!is.null(comp_kegg) && nrow(as.data.frame(comp_kegg)) > 0) {
        df <- as.data.frame(comp_kegg)
        if ("Description" %in% colnames(df)) df[, "Description"] <- clean_term_labels(df[, "Description"])
        write_table_xlsx(df, file.path(comp_dir, "compareCluster_KEGG"))
        safe_plot(print(dotplot(comp_kegg, showCategory=${params.show_category}) + ggtitle(paste(source, "KEGG compareCluster"))), file.path(comp_dir, "compareCluster_KEGG.pdf"), 10, 6)
      }

      if (!has_comp_go &&
          (is.null(comp_kegg) || nrow(as.data.frame(comp_kegg)) == 0)) {
        writeLines("No significant compareCluster enrichment under current cutoff", file.path(comp_dir, "NO_ENRICHMENT.txt"))
      }
    }

    if (nrow(S) == 2) {
      venn_dir <- file.path(source_dir, "overlap")
      dir.create(venn_dir, showWarnings=FALSE, recursive=TRUE)

      gr1 <- read_peak_gr(S\$peak_file[1])
      gr2 <- read_peak_gr(S\$peak_file[2])
      label1 <- pretty_sample(S\$sample[1])
      label2 <- pretty_sample(S\$sample[2])
      overlap <- tryCatch(findOverlapsOfPeaks(gr1, gr2, maxgap=${params.venn_maxgap}), error=function(e) NULL)
      if (!is.null(overlap)) {
        ov <- findOverlaps(gr1, gr2, maxgap=${params.venn_maxgap})
        n1 <- length(gr1)
        n2 <- length(gr2)
        n12 <- length(unique(S4Vectors::queryHits(ov)))
        n12 <- min(n12, n1, n2)
        venn_counts <- data.frame(
          source = source,
          sample_1 = label1,
          sample_2 = label2,
          n_sample_1 = n1,
          n_sample_2 = n2,
          n_overlap = n12,
          n_unique_1 = n1 - n12,
          n_unique_2 = n2 - n12,
          stringsAsFactors = FALSE
        )
        write.table(venn_counts, file.path(venn_dir, paste0(source, "_venn_counts.tsv")), sep="\t", quote=FALSE, row.names=FALSE)

        pdf(file.path(venn_dir, paste0(source, "_venn.pdf")), width=8, height=6)
        tryCatch({
          if (requireNamespace("VennDiagram", quietly=TRUE)) {
            grid::grid.newpage()
            VennDiagram::draw.pairwise.venn(
              area1 = n1,
              area2 = n2,
              cross.area = n12,
              category = c(label1, label2),
              fill = c("${params.venn_color_1}", "${params.venn_color_2}"),
              alpha = c(0.6, 0.6),
              lwd = 2,
              cex = 1.2,
              cat.cex = 1.1,
              cat.col = c("black", "black"),
              scaled = TRUE
            )
          } else {
            makeVennDiagram(
              overlap,
              NameOfPeaks = c(label1, label2),
              fill = c("${params.venn_color_1}", "${params.venn_color_2}"),
              col = "black",
              cat.col = "black",
              cat.cex = 1.1,
              cex = 1.2,
              lwd = 2,
              alpha = 0.7,
              main = paste("Overlap -", source)
            )
          }
        }, error=function(e) { plot.new(); title(conditionMessage(e)) })
        dev.off()
      }
    } else {
      venn_dir <- file.path(source_dir, "overlap")
      dir.create(venn_dir, showWarnings=FALSE, recursive=TRUE)
      writeLines(
        paste0("Overlap skipped: source=", source, " has ", nrow(S), " sample(s), requires exactly 2."),
        file.path(venn_dir, "OVERLAP_SKIPPED.txt")
      )
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
    // Try multiple naming conventions because chipseeker outputs evolved over time.
    def sampleAliases = [sample]
    if (sample.startsWith('idr__'))       sampleAliases << sample.replaceFirst(/^idr__/, '')
    if (sample.startsWith('consensus__')) sampleAliases << sample.replaceFirst(/^consensus__/, '')
    if (sample.startsWith('diffbind__'))  sampleAliases << sample.replaceFirst(/^diffbind__/, '')
    sampleAliases = sampleAliases.unique()

    for (alias in sampleAliases) {
      def p = file("${params.chipseeker_output}/${alias}/annotated_peaks.${alias}.tsv")
      if (p.exists()) return p
    }

    // Last fallback: legacy layout where directory uses original sample but filename uses trimmed alias
    for (dirAlias in sampleAliases) {
      for (fileAlias in sampleAliases) {
        def p = file("${params.chipseeker_output}/${dirAlias}/annotated_peaks.${fileAlias}.tsv")
        if (p.exists()) return p
      }
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

  def allInputFiles = rows
    .collectMany { r -> [r[2], r[3]] }
    .unique()

  Channel
    .fromList(['source\tsample\tannotated_tsv\tpeak_file'] + rows.collect { r ->
      "${r[0]}\t${r[1]}\t${r[2]}\t${r[3]}"
    })
    .collectFile(name: 'functional_enrich_manifest.tsv', newLine: true)
    .map { mf -> tuple(mf, allInputFiles) }
    .set { manifest_ch }

  functional_enrich(manifest_ch)
}

def globToRegex(pattern) {
  '^' + pattern
    .replace('.', '\\.')
    .replace('*', '.*')
    .replace('?', '.') + '$'
}
