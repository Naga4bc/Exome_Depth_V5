# Exome Depth 
#Code Written By: Vyomesh J & Bhumika P
#Date: 21-Dec-24
#-----------------------------------------------------------------------------------------------------------------------#
#Code Modified By: @Nagaraj
# Date: 17-Sep-25
# Usage: Rscript Exome_Depth_CNACall_V5.R
# Changes:
#   - Removed CLI arguments.
#   - Reads sample_genelist.txt and projects.txt directly.
#   - Auto-selects BAMs, BED file, threshold = 0.0001.
#   - Writes logs to ExomeDepth_log.txt with timestamps.
#-----------------------------------------------------------------------------------------------------------------------#


# ---------------- Docker Activation ---------------- #
activate_docker <- function() {
  docker_cmd <- paste(
    "sudo docker run --rm",
    "-u $(id -u bioinfo4):$(id -g bioinfo4)",
    "-v /home/bioinfo4/Patient_Samples/Exome_Depth_Dockerization:/workspace:rw",
    "euformatics/exomedepth:v1.1",
    "Rscript /workspace/Exome_Depth_CNACall_V5.R"
  )
  system(docker_cmd)
  quit(save = "no")  # exit current R session, Docker will continue
}

# If not inside Docker, re-launch inside Docker
if (!file.exists("/.dockerenv")) {
  activate_docker()
}


library(ExomeDepth)

# ---------------- Logging Function ---------------- #
log_file <- "/workspace/ExomeDepth_log.txt"
log_message <- function(msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] %s\n", timestamp, msg))
  write(sprintf("[%s] %s", timestamp, msg), file = log_file, append = TRUE)
}

# ---------------- Function to read sample genelist ---------------- #
read_sample_genelist <- function(file_path = "/workspace/sample_genelist.txt") {
  if (!file.exists(file_path)) stop("Error: sample_genelist.txt file not found.")
  sample_genelist <- readLines(file_path)
  parsed_data <- lapply(sample_genelist, function(line) {
    strsplit(line, ":")[[1]]
  })
  return(lapply(parsed_data, function(item) {
    list(
      sample_ids = unlist(strsplit(item[1], ",")),
      genes = unlist(strsplit(item[2], ","))
    )
  }))
}

# ---------------- Function to find BAM files from projects ---------------- #
find_bam_files <- function(sample_ids, projects_file = "/workspace/projects.txt") {
  if (!file.exists(projects_file)) stop("Error: projects.txt file not found.")
  projects <- readLines(projects_file)
  base_root <- "/workspace/basespace/Projects/"
  bam_paths <- c()
  for (sid in sample_ids) {
    bam_found <- FALSE
    for (project in projects) {
      bam_path <- file.path(base_root, project, "AppResults", sid, "Files", paste0(sid, ".bam"))
      if (file.exists(bam_path)) {
        bam_paths <- c(bam_paths, bam_path)
        log_message(paste("Found BAM for sample", sid, "in project", project))
        bam_found <- TRUE
        break
      }
    }
    if (!bam_found) {
      stop(paste("Error: BAM file not found for sample ID:", sid))
    }
  }
  return(bam_paths)
}

# ---------------- Function to select BED file based on sample ID ---------------- #
select_bed_file <- function(sample_ids) {
  if (any(grepl("-CE-| -CEFu-", sample_ids))) {
    return("/workspace/Indiegene_Target_2109PD006-V1_4BaseCare_1K_DNA_GRCh37.bed")
  } else if (any(grepl("-FEV2F2both-", sample_ids))) {
    return("/workspace/TarGT_First_v2_CDS_GRCh37_13_Mar_23.bed")
  } else if (any(grepl("-CDSV36XFUS-", sample_ids))) {
    return("/workspace/T1_CDS_V3_hg19_9D24.bed")
  } else if (any(grepl("-SE8-", sample_ids))) {
    return("/workspace/SureSelect_V8_Coverde_Modified.bed")
  } else {
    stop("Error: No matching BED file found for given sample IDs.")
  }
}

# ---------------- Function to extract gene range ---------------- #
get_gene_range <- function(bed_data, gene_name) {
  gene_data <- bed_data[bed_data$gene_name == gene_name, ]
  if (nrow(gene_data) == 0) stop(paste("Gene", gene_name, "not found in the BED file."))
  if (length(unique(gene_data$chromosome)) > 1) stop(paste("Gene", gene_name, "found on multiple chromosomes."))
  return(list(
    range = paste(gene_data$chromosome[1], ":", min(gene_data$start), "-", max(gene_data$end), sep = ""),
    gene_name = gene_name,
    chromosome = gene_data$chromosome[1],
    start = min(gene_data$start),
    end = max(gene_data$end)
  ))
}

# ---------------- MAIN SCRIPT ---------------- #

# Fixed parameters
sample_genelist_path <- "/workspace/sample_genelist.txt"
threshold <- 0.0001

log_message("Starting ExomeDepth CNV pipeline")

# Load sample genelist
genelist_data <- read_sample_genelist(sample_genelist_path)
log_message(paste("Loaded sample_genelist.txt with", length(genelist_data), "entries."))

# ---------------- Process Each Sample Set ---------------- #
for (genelist_row in genelist_data) {
  sample_ids <- genelist_row$sample_ids
  genes      <- genelist_row$genes

  log_message(
    paste(
      "Processing main sample:", sample_ids[1],
      "with comparison samples:", paste(sample_ids[-1], collapse = ",")
    )
  )

  # ---------------- Find Inputs ---------------- #
  my.bam        <- find_bam_files(sample_ids)
  bed_file_path <- select_bed_file(sample_ids)

  if (!file.exists(bed_file_path)) stop("Error: BED file not found.")
  log_message(paste("Using BED file:", bed_file_path))

  bed_data <- read.table(bed_file_path, sep = "\t", header = FALSE)
  colnames(bed_data) <- c("chromosome", "start", "end", "gene_name")

  # ---------------- BAM Counts ---------------- #
  log_message("Counting reads from BAM files")

  my.counts <- getBamCounts(
    bed.frame = bed_data,
    bam.files = my.bam,
    include.chr = FALSE
  )

  my.counts$chromosome <- sub("^chr", "", my.counts$chromosome)
  my.counts$seqnames   <- paste("chr", my.counts$chromosome, sep = "")

  # ---------------- ExomeDepth Analysis ---------------- #
  log_message("Running ExomeDepth analysis")

  test <- new(
    "ExomeDepth",
    test     = my.counts[, 5],
    reference= my.counts[, 6],
    formula  = "cbind(test, reference) ~ 1",
    subset.for.speed = seq(1, nrow(my.counts), 100)
  )

  my.test        <- my.counts[, 5]
  my.ref.samples <- my.counts[, c(6, 7, 8)]

  my.choice <- select.reference.set(
    test.counts      = my.test,
    reference.counts = as.matrix(my.ref.samples),
    bin.length       = (my.counts$end - my.counts$start) / 1000,
    n.bins.reduced   = 10000
  )

  my.matrix             <- as.matrix(my.counts[, my.choice$reference.choice, drop = FALSE])
  my.reference.selected <- apply(my.matrix, 1, sum)

  all.exons <- new(
    "ExomeDepth",
    test     = my.test,
    reference= my.reference.selected,
    formula  = "cbind(test, reference) ~ 1"
  )

  all.exons <- CallCNVs(
    x     = all.exons,
    transition.probability = threshold,
    chromosome = my.counts$chromosome,
    start      = my.counts$start,
    end        = my.counts$end,
    name       = my.counts$exon
  )

  # ---------------- Output Handling ---------------- #
  setwd("/workspace")   # Ensure working directory

  # Create sample folder
  sample_folder <- file.path(getwd(), sample_ids[1])
  if (!dir.exists(sample_folder)) {
    dir.create(sample_folder, recursive = TRUE, showWarnings = FALSE)
    log_message(paste("Created folder:", sample_folder))
  }

  # Save CNV calls
  output.file <- file.path(sample_folder, paste0(sample_ids[1], "_ED_CNA.csv"))
  write.csv(file = output.file, x = all.exons@CNV.calls, row.names = FALSE)
  log_message(paste("Saved CNV calls to", output.file))

  # ---------------- Gene Plots ---------------- #
  for (gene_name in genes) {
    tryCatch({
      gene_range <- get_gene_range(bed_data, gene_name)

      # Gene-specific folder
      gene_folder <- file.path(sample_folder, gene_name)
      if (!dir.exists(gene_folder)) {
        dir.create(gene_folder, recursive = TRUE, showWarnings = FALSE)
        log_message(paste("Created gene folder:", gene_folder))
      }

      # Plot file path
      plot_file <- file.path(
        gene_folder,
        paste0(
          sample_ids[1], "_plot_", gene_name, "_",
          gene_range$start, "_", gene_range$end, "_Amplification.png"
        )
      )

      # Generate plot
      png(plot_file)
      if (length(all.exons@CNV.calls) > 0) {
        plot(
          all.exons,
          sequence = sub("^chr", "", gene_range$chromosome),
          xlim     = c(gene_range$start - 1000, gene_range$end + 1000),
          count.threshold = 50,
          main     = paste(sample_ids[1], gene_name, sep = " - "),
          cex.lab  = 0.8,
          with.gene= TRUE
        )
      }
      dev.off()
      log_message(paste("Saved plot for gene", gene_name, ":", plot_file))

    }, error = function(e) {
      log_message(paste("Error in plotting gene", gene_name, ":", e$message))
    })
  }

  log_message(paste("Completed processing for:", sample_ids[1]))
}


  
