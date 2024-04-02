#!/opt/conda/envs/is_anno/bin/Rscript

library("optparse")
library("assertthat")

# Retrieve the subcommand
subcommand <- commandArgs(trailingOnly = TRUE)[1]

###############################
## Insertion Site Annotation ##
###############################

if (subcommand == "annotate") {
    # Parse command line arguments.
    options <- list(
        make_option(
            c("-i", "--input"), type="character", default=FALSE,
            help="Input BED file."
        ),
        make_option(
            c("-o", "--output"), type="character", default=FALSE,
            help="Output TSV file."
        ),
        make_option(
            c("-g", "--genome"), type="character", default=FALSE,
            help="Genome GFF/GTF file."
        ),
        make_option(
            c("-p", "--promoter"), type="character", default="-2000,200",
            help="Promoter region. [default: %default]"
        ),
        make_option(
            c("-l", "--level"), type="character", default="transcript",
            help="Annotation level. [default: %default]"
        )
    )

    args <- parse_args(
        OptionParser(option_list=options),
        positional_arguments=1
    )

    input_file <- args$options$input
    output_file <- args$options$output
    genome_file <- args$options$genome
    promoter <- as.numeric(strsplit(args$options$promoter, ",")[[1]])
    level <- args$options$level

    # Input validation.
    assert_that(
        is.readable(input_file),
        msg="Input file is not readable."
    )
    assert_that(
        grepl(input_file, pattern="\\.bed$", ignore.case=TRUE),
        msg="Input file does not have a BED extension."
    )
    assert_that(
        is.writeable(dirname(output_file)),
        msg="Output directory is not writeable."
    )
    assert_that(
        grepl(output_file, pattern="\\.tsv$", ignore.case=TRUE),
        msg="Output file does not have a TSV extension."
    )
    assert_that(
        is.readable(genome_file),
        msg="Genome file is not readable."
    )
    assert_that(
        grepl(genome_file, pattern="\\.(gff|gtf)(\\.gz)?$", ignore.case=TRUE),
        msg="Genome file does not have a GFF/GTF extension."
    )
    assert_that(
        all(is.numeric(promoter)) & length(promoter) == 2,
        msg="Promoter region must be numeric vector of length 2."
    )
    assert_that(
        promoter[1] < promoter[2],
        msg="Promoter region start must be less than end."
    )
    assert_that(
        is.string(level) & level %in% c("gene", "transcript"),
        msg="Annotation level must be 'gene' or 'transcript'."
    )

    # Import the insertion sites using rtracklayer.
    insertion_sites <- rtracklayer::import.bed(input_file)

    # Write an output file if there are no insertion sites.
    if (length(insertion_sites) == 0) {
        readr::write_tsv(
            tibble::tibble(
                sample_id=character(),
                seqnames=character(),
                start=integer(),
                end=integer(),
                width=integer(),
                strand=character(),
                score=integer(),
                annotation=character(),
                gene_seqnames=character(),
                gene_start=integer(),
                gene_end=integer(),
                gene_width=integer(),
                gene_strand=character(),
                gene_id=character(),
                transcript_id=character(),
                tss_distance=integer(),
                gene_version=character(),
                gene_name=character(),
                gene_source=character(),
                gene_biotype=character()
            ),
            output_file
        )
        quit(save="no", status=0)
    }

    # Import the genome as a TxDb using GenomicFeatures.
    txdb <- GenomicFeatures::makeTxDbFromGFF(genome_file)

    # Annotate the insertion sites to the nearest gene.
    annotated_insertion_sites <- insertion_sites |>
        ChIPseeker::annotatePeak(
            tssRegion=promoter,
            TxDb=txdb,
            level=level,
            sameStrand=TRUE
        ) |>
        tibble::as_tibble() |>
        dplyr::select(!name) |>
        dplyr::rename(
            "gene_seqnames"="geneChr",
            "gene_start"="geneStart",
            "gene_end"="geneEnd",
            "gene_width"="geneLength",
            "gene_strand"="geneStrand",
            "gene_id"="geneId",
            "transcript_id"="transcriptId",
            "tss_distance"="distanceToTSS"
        )

    # Add additional gene info the annotated insertion sites.
    gene_info <- genome_file |>
        rtracklayer::import() |>
        tibble::as_tibble() |>
        dplyr::select(
            gene_id, gene_version, gene_name, gene_source, gene_biotype
        ) |>
        unique()

    annotated_insertion_sites <- dplyr::left_join(
        annotated_insertion_sites,
        gene_info,
        by="gene_id"
    )

    # Add a sample ID column.
    annotated_insertion_sites <- annotated_insertion_sites |>
        dplyr::mutate(
            sample_id=input_file |>
                tools::file_path_sans_ext() |>
                basename() |>
                sub(x=_, pattern="_S\\d+_L00\\d\\.insertion_sites$", replacement="")
        ) |>
        dplyr::relocate(sample_id, .before=everything())

    # Save the annotated insertion sites.
    readr::write_tsv(annotated_insertion_sites, output_file)
}

###############################
## Genomic Distribution Plot ##
###############################

process_insertion_sites <- function(insertion_site_list) {
    # Import the insertion sites.
    insertion_sites <- input_files |>
        (\(x){
            names(x) <- x |>
                basename() |>
                tools::file_path_sans_ext() |>
                sub(
                    x=_, pattern="_S\\d+_L00\\d\\.insertion_sites\\.annotated$",
                    replacement=""
                )
            return(x)
        })() |>
        lapply(readr::read_tsv, col_types="cciiicicciiiccciiccc")
    
    insertion_sites <- insertion_sites[
        stringr::str_sort(names(insertion_sites), numeric=TRUE, decreasing=!decreasing_order)
    ]

    # Discard samples that have no insertion sites.
    insertion_sites <- insertion_sites[
        vapply(insertion_sites, \(x) nrow(x) > 0, logical(1))
    ]
    
    # Clean up the annotations.
    insertion_sites <- lapply(insertion_sites, \(x) {
        x$annotation <- dplyr::case_when(
            grepl(x$annotation, pattern="^Intron") ~ "Intron",
            grepl(x$annotation, pattern="^Exon") ~ "Exon",
            grepl(x$annotation, pattern="^Promoter") ~ "Promoter",
            grepl(x$annotation, pattern="5' UTR") ~ "5'UTR",
            grepl(x$annotation, pattern="3' UTR") ~ "3'UTR",
            grepl(x$annotation, pattern="^Downstream") ~ "Downstream",
            x$annotation == "Distal Intergenic" ~ "Intergenic",
        )
        return(x)
    })

    # Return the results.
    return(insertion_sites)    
}

if (subcommand == "genomic_feature_plot") {
    library("ggplot2")

    # Parse the command line arguments.
    options <- list(
        make_option(
            c("-i", "--input"), type="character", default=FALSE,
            help="Input TSV file(s)."
        ),
        make_option(
            c("-o", "--plot_output"), type="character", default=FALSE,
            help="Output plot."
        ),
        make_option(
            c("-t", "--table_output"), type="character", default=FALSE,
            help="Output table."
        ),
        make_option(
            c("-d", "--decreasing"), action="store_true", default=FALSE,
            help="Sort the sample IDs decreasing order. [default: %default]"
        ),
        make_option(
            c("-p", "--dpi"), type="integer", default=300,
            help="DPI of the output PNG file. [default: %default]"
        )
    )

    args <- parse_args(
        OptionParser(option_list=options),
        positional_arguments=1
    )

    input_files <- strsplit(args$options$input, ",", fixed=TRUE)[[1]]
    plot_output <- args$options$plot_output
    table_output <- args$options$table_output
    decreasing_order <- args$options$decreasing
    dpi <- args$options$dpi

    # Input validation.
    assert_that(
        all(vapply(input_files, is.readable, logical(1))),
        msg="One ore more input files are not readable."
    )
    assert_that(
        all(vapply(input_files, \(x) grepl(x, pattern="\\.tsv$", ignore.case=TRUE), logical(1))),
        msg="One ore more input files do not have a TSV extension."
    )
    assert_that(
        is.writeable(dirname(plot_output)),
        msg="Plot output directory is not writeable."
    )
    assert_that(
        has_extension(plot_output, "png"),
        msg="Output plot does not have a PNG extension."
    )
    assert_that(
        is.writeable(dirname(table_output)),
        msg="Table output directory is not writeable."
    )
    assert_that(
        has_extension(table_output, "tsv"),
        msg="Output table does not have a PNG extension."
    )
    assert_that(
        is.flag(decreasing_order),
        msg="Decreasing order must be a flag."
    )
    assert_that(
        is.count(dpi) & dpi <= 300,
        msg="DPI must be a positive integer ≤ 300."
    )

    # Process the insertion sites.
    insertion_sites <- process_insertion_sites(input_files)

    # Prepare the data for plotting.
    df <- insertion_sites |>
        dplyr::bind_rows(.id="sample_id") |>
        dplyr::count(sample_id, annotation) |>
        tidyr::complete(sample_id, annotation, fill=list(n=0)) |>
        dplyr::group_by(sample_id) |>
        dplyr::mutate(frac=n/sum(n)) |>
        dplyr::ungroup() |>
        dplyr::mutate(
            sample_id=factor(
                sample_id,
                levels=stringr::str_sort(
                    unique(sample_id), numeric=TRUE, decreasing=!decreasing_order
                )
            ),
            annotation=forcats::fct_reorder(
                annotation, frac, .fun=mean, .desc=TRUE
            )
        )

    # Plot the data.    
    p <- ggplot(df, aes(x=frac, y=sample_id, fill=annotation)) +
        geom_col(width=0.85) +
        theme_classic() +
        scale_fill_viridis_d() +
        labs(
            x="Percent of Insertion Sites",
            y="Sample ID",
            fill="Genomic\nFeature",
            title="Insertion Site\nGenomic Feature Distribution"
        ) +
        scale_x_continuous(labels=scales::percent)

    plot_height <- length(insertion_sites) * 0.60    
    ggsave(
        plot_output, plot=p, height=plot_height, width=6,
        device=ragg::agg_png, dpi=300, units="in", bg="white"
    )

    # Export a table of results.
    df |>
        dplyr::select(sample_id, annotation, n) |>
        tidyr::pivot_wider(names_from=annotation, values_from=n, values_fill=0) |>
        readr::write_tsv(table_output)
}

###################
## Promoter Plot ##
###################

if (subcommand == "promoter_plot") {
    library("ggplot2")

    # Parse command line arguments
    options <- list(
        make_option(
            c("-i", "--input"), type="character", default=FALSE,
            help="Input TSV file(s)."
        ),
        make_option(
            c("-o", "--output"), type="character", default=FALSE,
            help="Output plot."
        ),
        make_option(
            c("-b", "--promoter_plot_bins"), type="integer", default=100,
            help="Number of bins for the promoter plot. [default: %default]"
        ),
        make_option(
            c("-d", "--decreasing"), action="store_true", default=FALSE,
            help="Sort the sample IDs decreasing order. [default: %default]"
        ),
        make_option(
            c("-p", "--dpi"), type="integer", default=300,
            help="DPI of the output PNG file. [default: %default]"
        )
    )

    args <- parse_args(
        OptionParser(option_list=options),
        positional_arguments=1
    )

    input_files <- strsplit(args$options$input, ",", fixed=TRUE)[[1]]
    output <- args$options$output
    n_bins <- args$options$promoter_plot_bins
    decreasing_order <- args$options$decreasing
    dpi <- args$options$dpi

    # Input validation.
    assert_that(
        all(vapply(input_files, is.readable, logical(1))),
        msg="One ore more input files are not readable."
    )
    assert_that(
        all(vapply(input_files, \(x) grepl(x, pattern="\\.tsv$", ignore.case=TRUE), logical(1))),
        msg="One ore more input files do not have a TSV extension."
    )
    assert_that(
        is.writeable(dirname(output)),
        msg="Output directory is not writeable."
    )
    assert_that(
        has_extension(output, "png"),
        msg="Output plot does not have a PNG extension."
    )
    assert_that(
        is.count(n_bins) & n_bins <= 1000,
        msg="Number of bins must be a positive integer ≤ 1000."
    )
    assert_that(
        is.flag(decreasing_order),
        msg="Decreasing order must be a flag."
    )
    assert_that(
        is.count(dpi) & dpi <= 300,
        msg="DPI must be a positive integer ≤ 300."
    )

    # Process the insertion sites.
    insertion_sites <- process_insertion_sites(input_files)

    # Create the promoter plot.
    p <- insertion_sites |>
        dplyr::bind_rows(.id="sample_id") |>
        ggplot(aes(x=tss_distance)) +
            geom_histogram(
                bins=n_bins, fill="salmon",
                aes(y=after_stat(count)/sum(after_stat(count)))
            ) +
            theme_bw() +
            facet_wrap(
                sample_id~., ncol=1, scales="free_y",
                strip.position="right"
            ) +
            scale_x_continuous(
                labels=scales::label_number(suffix=" kbp",scale=1e-3)
            ) +
            scale_y_continuous(labels=scales::percent) +
            theme(
                strip.text.y.right=element_text(angle=0),
                strip.background=element_blank(),
            ) +
            labs(
                x="Distance to Promoter",
                y="Percent of Insertion Sites",
                title="Distribution of Insertion Sites\nRelative to Promoters"
            )

    plot_height <- length(insertion_sites) * 0.75    
    ggsave(
        output, plot=p, height=plot_height, width=6,
        device=ragg::agg_png, dpi=dpi, units="in", bg="white"
    )
}

####################
## Top Sites Plot ##
####################

if (subcommand == "top_sites_plot") {
    library("ggplot2")

    # Parse the command line arguments.
    options <- list(
        make_option(
            c("-i", "--input"), type="character", default=FALSE,
            help="Input TSV file(s)."
        ),
        make_option(
            c("-o", "--plot_output"), type="character", default=FALSE,
            help="Output plot."
        ),
        make_option(
            c("-t", "--table_output"), type="character", default=FALSE,
            help="Output table."
        ),
        make_option(
            c("-n", "--top_sites"), type="integer", default=10,
            help="Number of top insertion sites to plot. [default: %default]"
        ),
        make_option(
            c("-d", "--decreasing"), action="store_true", default=FALSE,
            help="Sort the sample IDs decreasing order. [default: %default]"
        ),
        make_option(
            c("-p", "--dpi"), type="integer", default=300,
            help="DPI of the output PNG file. [default: %default]"
        )
    )

    args <- parse_args(
        OptionParser(option_list=options),
        positional_arguments=1
    )

    input_files <- strsplit(args$options$input, ",", fixed=TRUE)[[1]]
    plot_output <- args$options$plot_output
    table_output <- args$options$table_output
    top_sites <- args$options$top_sites
    decreasing_order <- args$options$decreasing
    dpi <- args$options$dpi

    # Input validation.
    assert_that(
        all(vapply(input_files, is.readable, logical(1))),
        msg="One ore more input files are not readable."
    )
    assert_that(
        all(vapply(input_files, \(x) grepl(x, pattern="\\.tsv$", ignore.case=TRUE), logical(1))),
        msg="One ore more input files do not have a TSV extension."
    )
    assert_that(
        is.writeable(dirname(plot_output)),
        msg="Plot output directory is not writeable."
    )
    assert_that(
        has_extension(plot_output, "png"),
        msg="Output plot does not have a PNG extension."
    )
    assert_that(
        is.writeable(dirname(table_output)),
        msg="Table output directory is not writeable."
    )
    assert_that(
        has_extension(table_output, "tsv"),
        msg="Output table does not have a PNG extension."
    )
    assert_that(
        is.count(top_sites) & top_sites <= 100,
        msg="Number of top sites must be a positive integer ≤ 100."
    )
    assert_that(
        is.flag(decreasing_order),
        msg="Decreasing order must be a flag."
    )
    assert_that(
        is.count(dpi) & dpi <= 300,
        msg="DPI must be a positive integer ≤ 300."
    )

    # Process the insertion sites.
    insertion_sites <- process_insertion_sites(input_files)

    # Prepare the data for plotting.
    df <- insertion_sites |>
        purrr::imap(\(x, y) {
            x  <- dplyr::mutate(x, 
                gene=dplyr::case_when(
                    annotation == "Intergenic" ~ glue::glue("chr{seqnames}:{start}-{end}:{strand}"),
                    annotation != "Intergenic" & (is.na(gene_name) | gene_name == "") ~ gene_id,
                    TRUE ~ gene_name
                ),
                sample_id=y
            )

            top_genes <- x |>
                dplyr::slice_max(score, n=top_sites, with_ties=FALSE) |>
                dplyr::arrange(desc(score)) |>
                dplyr::pull(gene)
            
            x <- x |>
                dplyr::mutate(
                    gene=ifelse(gene %in% top_genes, gene, NA)
                ) |>
                dplyr::group_by(sample_id, gene) |>
                dplyr::summarize(score=sum(score), .groups="drop") |>
                dplyr::mutate(
                    gene=forcats::fct_reorder(gene, score, .fun=unique, .desc=TRUE),
                    frac=score/sum(score)
                )
        })

    # Plot the top insertion sites.
    p <- df |>
        purrr::imap(\(x, y) {
            ggplot(x, aes(x=sample_id, y=frac, fill=gene)) +
                geom_col() +
                theme_classic() +
                scale_fill_viridis_d(na.value="lightgrey") +
                theme(axis.title.x=element_blank()) +
                labs(
                    y="Percent of Insertion Sites",
                    fill="Genomic\nLocation"
                ) +
                scale_y_continuous(labels=scales::percent)
        }) |>
        patchwork::wrap_plots(ncol=1)

    plot_height <- length(insertion_sites) * 4.5
    ggsave(
        plot_output, plot=p, height=plot_height, width=4.25, limitsize=FALSE,
        device=ragg::agg_png, dpi=300, units="in", bg="white"
    )

    # Save the table.
    df |>
        dplyr::bind_rows() |>
        dplyr::select(!frac) |>
        dplyr::mutate(
            gene=as.character(gene),
            gene=ifelse(is.na(gene), "Other", gene)
        ) |>
        tidyr::pivot_wider(names_from=gene, values_from=score) |>
        readr::write_tsv(table_output)
}

##############################
## Insertion Site Fractions ##
##############################

if (subcommand == "fractions") {
    library("ggplot2")

    # Retrieve the command line arguments.
    options <- list(
        make_option(
            c("-i", "--input"), type="character", default=FALSE,
            help="Input TSV file(s)."
        ),
        make_option(
            c("-o", "--output"), type="character", default=FALSE,
            help="Output file."
        ),
        make_option(
            c("-d", "--dpi"), type="integer", default=300,
            help="DPI of the output PNG file. [default: %default]"
        ),
        make_option(
            c("-m", "--height_multiplier"), type="numeric", default=0.85,
            help="Height multiplier for the plot. [default: %default]"
        )
    )

    args <- parse_args(
        OptionParser(option_list=options),
        positional_arguments=1
    )

    input_files <- strsplit(args$options$input, ",", fixed=TRUE)[[1]]
    output_file <- args$options$output
    dpi <- args$options$dpi
    height_multiplier <- args$options$height_multiplier

    # Input validation
    assert_that(
        all(vapply(input_files, is.readable, logical(1))),
        msg="One ore more input files are not readable."
    )
    assert_that(
        all(vapply(input_files, \(x) grepl(x, pattern="\\.tsv$", ignore.case=TRUE), logical(1))),
        msg="One ore more input files do not have a TSV extension."
    )
    assert_that(
        is.writeable(dirname(output_file)),
        msg="Output directory is not writeable."
    )
    assert_that(
        grepl(output_file, pattern="\\.png$", ignore.case=TRUE),
        msg="Output file does not have a PNG extension."
    )
    assert_that(
        is.count(dpi) & dpi <= 300,
        msg="DPI must be a positive integer ≤ 300."
    )
    assert_that(
        is.numeric(height_multiplier) & height_multiplier > 0 & height_multiplier <= 10,
        msg="Height multiplier must be a positive numeric value ≤ 10."
    )

    # Import the insertion sites.
    insertion_sites <- input_files |>
        (\(x){
            names(x) <- x |>
                basename() |>
                tools::file_path_sans_ext() |>
                sub(
                    x=_, pattern="_S\\d+_L00\\d\\.insertion_sites\\.annotated$",
                    replacement=""
                )
            return(x)
        })() |>
        lapply(readr::read_tsv, col_types="cciiicicciiiccciiccc")

    # Function to generate a random color palette.
    random_pal <- function(vals) {
        n_vals <- length(unique(vals))
        color_map <- unique(randomcoloR::randomColor(n_vals + 100))[1:n_vals]
        names(color_map) <- unique(vals)
        colors <- vapply(vals, \(x) color_map[x], character(1))
        return(colors)
    }

    # Prepare data for plotting.
    df <- insertion_sites |>
        # Summarize the fraction of reads per insertion site.
        dplyr::bind_rows() |>
        dplyr::mutate(
            sample_id=sub(
                x=sample_id, pattern="_S\\d+_L00\\d\\.insertion_sites\\.filtered$",
                replacement=""
            )
        ) |>
        dplyr::group_by(sample_id) |>
        dplyr::mutate(frac=score / sum(score)) |>
        dplyr::ungroup() |>
        # Add a column with the custom colors.
        dplyr::mutate(
            posid=paste(seqnames, start, end, strand, sep=":"),
            color=random_pal(posid),
            sample_id=factor(
                sample_id,
                levels=rev(stringr::str_sort(unique(sample_id)))
            )
        ) |>
        split(~sample_id) |>
        lapply(\(x) {
            x <- x |>
                dplyr::mutate(
                    color=forcats::fct_drop(color),
                    color=forcats::fct_reorder(color, frac, .desc=TRUE)
                )
            return(x)
        })

    # Plot the data.
    p <- 
        lapply(df, \(x) {
            ggplot(x, aes(x=frac, y=sample_id)) +
                geom_col(aes(fill=color)) +
                scale_fill_identity() +
                scale_x_continuous(labels=scales::percent) +
                theme_classic() +
                theme(
                    axis.title.x=element_blank(),
                    axis.title.y=element_blank()
                )
        }) |>
        patchwork::wrap_plots(ncol=1) |>
        patchwork::patchworkGrob() |>
        gridExtra::grid.arrange(
            left="Sample ID", bottom="Percent of Reads",
            top="Percent of Reads per Insertion Site"
        )

    plot_height= length(insertion_sites) * 0.85
    ggsave(
        output_file, plot=p, height=plot_height, width=6,
        device=ragg::agg_png, dpi=dpi, units="in", bg="white"
    )
}

######################
## Seqlogo Creation ##
######################

if (subcommand == "seqlogo") {
    library("ggseqlogo")

    # Load the command line arguments.
    # The script takes two arguments:
    #   -i/--input: Comma-separated list of input files.
    #   -o/--output: Output .png file.
    options <- list(
        make_option(
            c("-i", "--input"), type="character", default=FALSE,
            help="Comma-separated list of input files."
        ),
        make_option(
            c("-o", "--output"), type="character", default=FALSE,
            help="Output file."
        ),
        make_option(
            c("-d", "--dpi"), type="integer", default=300,
            help="DPI of the output file PNG. [default: %default]"
        ),
        make_option(
            c("-m", "--height_multiplier"), type="numeric", default=1.25,
            help="Height multiplier of the output file PNG. [default: %default]"
        )
    )

    args <- parse_args(
        OptionParser(option_list=options),
        positional_arguments=1
    )

    infiles <- strsplit(args$options$input, ",", fixed=TRUE)[[1]]
    outfile <- args$options$output

    # Validate the inputs.
    assert_that(
        is.writeable(dirname(outfile)),
        msg="Output directory is not writeable."
    )
    assert_that(
        all(vapply(infiles, \(x) grepl(x, pattern="\\.(fasta|fa)(\\.gz)?$"), logical(1))),
        msg="One or more input files are not FASTA files."
    )
    assert_that(
        all(vapply(infiles, is.readable, logical(1))),
        msg="One or more input files are not readable."
    )
    assert_that(
        has_extension(outfile, "png"),
        msg="Output file does not have a PNG extension."
    )
    assert_that(
        is.count(args$options$dpi),
        msg="DPI must be a positive integer."
    )
    assert_that(
        args$options$dpi <= 300,
        msg="DPI must be less than or equal to 300."
    )
    assert_that(
        is.numeric(args$options$height_multiplier),
        msg="Height multiplier must be a numeric value."
    )
    assert_that(
        args$options$height_multiplier > 0 & args$options$height_multiplier <= 10,
        msg="Height multiplier must be > 0 and <= 10."
    )

    # Read the input files into a list of DNAStringSets.
    # Then convert to a list of character vectors.
    seq_data <- infiles |>
        (\(x) {
            names(x) <- x |>
                basename() |>
                tools::file_path_sans_ext() |>
                sub("\\.insertion_site_seq", "", x=_)
            return(x)
        })() |>
        lapply(\(x) {
            x <- x |>
                Biostrings::readDNAStringSet() |>
                as.character() |>
                unname()
            return(x)
        })

    # Order the sequences alphanumerically.
    seq_data <- seq_data[
        stringr::str_sort(names(seq_data), numeric=TRUE, decreasing=TRUE)
    ]

    # Remove samples with no sequences.
    seq_data <- seq_data[lengths(seq_data) > 0]

    # Create the sequence logo and save as PNG.
    p <- ggseqlogo(seq_data, ncol=1, seq_type="dna") +
        ggplot2::theme(text=ggplot2::element_text(size=7)) +
        ggplot2::scale_x_continuous(breaks=seq_len(41), labels=seq(-20, 20, 1))

    fig_height <- length(infiles) * args$options$height_multiplier

    ggplot2::ggsave(
        outfile, plot=p, height=fig_height, width=8, bg="white",
        device=ragg::agg_png, dpi=args$options$dpi, units="in"
    )
}