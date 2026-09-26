# =============================================================================
# Shared helper: taxon names in figure labels are set in italics, at every rank
# (phylum down to species), as in the manuscript text. Labels that are not names
# (Other, Unknown, Mixed, ...) stay upright.
#
# Usage in a figure script:
#   source(file.path(BASE, "taxon_italics_n76.R"))
#   ... + scale_y_discrete(labels = md_taxon) +
#         theme(axis.text.y = ggtext::element_markdown(size = 7))
#   ... md_taxon_in("Enterobacteriaceae abundance")  # embedded name
# =============================================================================

suppressPackageStartupMessages(library(ggtext))

# every name GTDB assigned to the 584 MAGs, at all ranks
TAXA_VOCAB <- local({
  f <- file.path(Sys.getenv("COVID_MAG_BASE", unset = getwd()), "result2/inputfiles/new/MAG_quality_taxonomy.tsv")
  tax <- utils::read.delim(f, check.names = FALSE, quote = "", comment.char = "")
  nm  <- trimws(unlist(strsplit(as.character(tax$gtdbtk_classification), ";")))
  nm  <- sub("^[a-z]__", "", nm)
  nm  <- nm[nzchar(nm) & grepl("^[A-Z]", nm)]
  # names used in the figures that predate the GTDB release, or come from
  # MetaPhlAn rather than GTDB-Tk
  sort(unique(c(nm, "Enterobacteriaceae", "Cyanobacteria", "Melainabacteria",
                "Firmicutes", "Proteobacteria", "Actinobacteria", "Bacteroidetes",
                "Verrucomicrobia", "Euryarchaeota")))
})

NON_TAXA <- c("Other", "Others", "Unknown", "Unclassified", "Mixed", "None",
              "NA", "Bacteria", "Archaea", "Total", "All")

# a label that is itself a taxon name, e.g. an axis of phyla or species
md_taxon <- function(x) {
  x <- as.character(x)
  keep <- !is.na(x) & !(x %in% NON_TAXA) &
    (x %in% TAXA_VOCAB |
       grepl("^[A-Z][a-z]+(_[A-Z])? (sp[0-9]+|[a-z][a-z0-9_]{2,})$", x) |  # Genus species / sp900...
       grepl("^[A-Z]\\. [a-z]{3,}", x))                                    # E. coli
  ifelse(keep, paste0("*", x, "*"), x)
}

# ggtext cannot be used on an axis that ggh4x resizes (force_panelsizes drops the
# markdown grob and the asterisks print) or on a legend that patchwork merges with
# a plain element. For those, turn the markdown into a plotmath expression, which
# ggplot draws with the ordinary text element.
expr_from_md <- function(x) {
  x <- as.character(x)
  out <- lapply(x, function(s) {
    if (is.na(s) || !nzchar(s)) return(str2lang('""'))
    s <- gsub('"', "'", s, fixed = TRUE)        # a quote would close the plotmath string
    parts <- strsplit(s, "*", fixed = TRUE)[[1]]  # odd pieces upright, even italic
    if (!length(parts)) return(str2lang('""'))
    bits <- vapply(seq_along(parts), function(i)
      if (i %% 2 == 0) sprintf('italic("%s")', parts[i]) else sprintf('"%s"', parts[i]),
      character(1))
    str2lang(paste(bits[nzchar(parts)], collapse = "*"))
  })
  as.expression(out)
}

# a label that contains a taxon name inside other words,
# e.g. "Enterobacteriaceae abundance" or "E. coli biofilm formation"
md_taxon_in <- function(x, extra = character()) {
  x <- as.character(x)
  one <- TAXA_VOCAB[nchar(TAXA_VOCAB) > 5 & !grepl(" ", TAXA_VOCAB)]
  # binomials KEGG puts in pathway and module names, which are not all in the
  # catalogue (Vibrio cholerae, Bacillus anthracis, ...)
  two <- unique(c(TAXA_VOCAB[grepl(" ", TAXA_VOCAB)],
                  "Escherichia coli", "Vibrio cholerae", "Helicobacter pylori",
                  "Bacillus anthracis", "Yersinia pestis", "Staphylococcus aureus",
                  "Pseudomonas aeruginosa", "Salmonella enterica",
                  "Mycobacterium tuberculosis", "Neisseria meningitidis"))
  vapply(x, function(s) {
    if (is.na(s) || !nzchar(s)) return(s)
    s <- gsub("\\b([A-Z]\\. [a-z]{3,})", "*\\1*", s)          # abbreviated binomials
    for (w in c(extra, two)) {
      if (!grepl(w, s, fixed = TRUE)) next
      s <- gsub(paste0("(?<![*A-Za-z])", w, "(?![*A-Za-z])"), paste0("*", w, "*"),
                s, perl = TRUE)
    }
    for (w in intersect(one, unlist(strsplit(s, "[^A-Za-z_]+")))) {
      if (w %in% NON_TAXA) next
      s <- gsub(paste0("(?<![*A-Za-z])", w, "(?![*A-Za-z])"), paste0("*", w, "*"),
                s, perl = TRUE)
    }
    s
  }, character(1), USE.NAMES = FALSE)
}
