#' @include AllClass.R
#' @importFrom dplyr %>%
NULL

#' Constructor of MetRKAT class
#'
#' @param MetRKAT
#'
#' @return
#' @export
#'
#' @examples
setMethod("initialize", signature="MetRKAT", function(.Object, ...) {
    .Object <- callNextMethod()

    # If gap value is equal to window size, set it to 0 (disabled)
    if (.Object@gap == .Object@wsize) {
        .Object@gap <- 0
    }

    if (! is.null(.Object@dbsrc)) {
        message("Reading data from the DB connection provided...", appendLF = FALSE)

        .Object@data <- dplyr::tbl(.Object@dbsrc, "data")
        .Object@annotation <- dplyr::tbl(.Object@dbsrc, "annotation")

        message("Done!")
    } else {

        message("Discarding/imputing NA values... ", appendLF = FALSE)

        # Impute NA values
        isNA_data <- rowMeans(is.na(.Object@data))

        # Remove those rows having more than the specified percentage
        keep_rows <- names(isNA_data[isNA_data < .Object@max.na])

        # Convert to numeric
        .Object@data <- as.data.frame(data.matrix(.Object@data[keep_rows, ]))

        # Impute those rows having NA values but in lower percentage
        impute_rows <- names(isNA_data[isNA_data < .Object@max.na &
                                           isNA_data > 0])

        if (length(impute_rows)) {
            # Call function "impute" from package Hmisc on every row, using the mean as
            # impute method.
            .Object@data[impute_rows, ] <- t(apply(.Object@data[impute_rows, ], 1,
                                                   Hmisc::impute, mean))
        }

        # Show message of completed NA step process
        message("Done!")

        # Show message of next step: annotation data frame
        message("Preparing annotation dataset... ", appendLF = FALSE)

        if (! nrow(.Object@annotation)) {
            # Initialize 'annotation' slot with a data frame containing annotation
            # data linked with the rows in the data provided

            # Locations: rownames, chr, pos, strand.
            # Exlude sites not present on data.
            # Convert the chr notation to numeric to proper ordering.
            chrToNum <- function(chroms) {
                chr_num <- gsub("chr", "", chroms)
                chr_num <- gsub("X", 23, chr_num)
                chr_num <- gsub("Y", 24, chr_num)

                chr_num <- as.numeric(chr_num)

                return(chr_num)
            }

            annot_loc <- data.frame(IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Locations,
                                    stringsAsFactors = FALSE) %>%
                tibble::rownames_to_column("site") %>%
                dplyr::filter(site %in% rownames(.Object@data)) %>%
                dplyr::mutate(chr = chrToNum(chr))

            # Other: contains UCSC_RefGene_Name
            # Split the column 'gene' as one site may affect to more than one gene and
            # are separeted my semicolon. In this verbose way dplyr can later group by
            # this field.
            annot_genes <- data.frame(IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Other,
                                      stringsAsFactors = FALSE) %>%
                dplyr::select(gene = UCSC_RefGene_Name) %>%
                tibble::rownames_to_column("site") %>%
                dplyr::mutate(gene = strsplit(gene, ";")) %>%
                tidyr::unnest(gene)

            # Prepare the final data frame merging all the required information, ordering
            # by chromosome and location
            .Object@annotation <- dplyr::inner_join(annot_loc,
                                                    annot_genes, by = "site")
        } else {
            # Make sure to remove those sites not present on data
            .Object@annotation <- dplyr::filter(.Object@annotation,
                                                site %in% rownames(.Object@data))
        }

        # Order by chr and position
        .Object@annotation <- dplyr::arrange(.Object@annotation, chr, pos)

        # Reorder the data rows according the the position on the annotation.
        .Object@data <- .Object@data[unique(as.character(.Object@annotation$site)), ]

        # Add a new column to annotation that links to data based on position.
        .Object@annotation$row <- match(.Object@annotation$site, rownames(.Object@data))

        # Move the site names from rownames to a proper column in the df
        # .Object@data <- dplyr::add_rownames(.Object@data, "site")

        # Show message of done process
        message("Done!")
    }

    return(.Object)
})

#' Apply the test to an already prepared method
#'
#' @param MetRKAT
#'
#' @return
#' @export
#'
#' @examples
setMethod("applyRKAT", signature="MetRKAT", function(object, overwrite = FALSE, ...) {

    # Check for existing slots on the instance instead of repeating the process.
    # This allows to "reset" SQLite queries and solve expired connections errors.
    intervals_pvals <- NULL

    if (! overwrite) {
        if (is.null(object@dbsrc)) {
            # If no database has been used, return the object directly as there is no
            # need to recreate the "query".
            if(! is.null(object@results))
                return(object@results)
        } else {
            n_results <- as.numeric(dplyr::collect(dplyr::tbl(
                object@dbsrc,
                dplyr::sql("SELECT MAX(_ROWID_) FROM results LIMIT 1")
            )))

            if (! is.na(n_results))
                intervals_pvals <- dplyr::tbl(object@dbsrc, "results")
        }
    }

    if (is.null(intervals_pvals)) {
        message("Retrieving p-value for each window...")

        .apply.rkat <- function(df_subset) {
            # Get methylation data inside that window
            # Samples need to be on rows
            window_data <- t(dplyr::slice(object@data,
                                          df_subset$first_row:df_subset$last_row))

            # Create a list containing the multiple kernels (if provided)
            window_kernels <- lapply(object@distmethod, function(dmethod) {
                # Apply the method to that window data
                window_dist <- as.matrix(vegan::vegdist(window_data,
                                                        method=dmethod))
                # Make kernel matrix
                output_kernel <- MiRKAT::D2K(window_dist)

                return(output_kernel)
            })

            # Pass the provided parameters (including y)
            # TODO: change params to be fixed instead of ellipsis?
            window_pval <- MiRKAT::MiRKAT(Ks = window_kernels, ...)

            # In case of multiple kernels, a list is returned
            if (is.list(window_pval)) {
                out_df <- data.frame(df_subset,
                                     pval = window_pval$indivP,
                                     kernel = object@distmethod,
                                     omnibus = window_pval$omnibus_p,
                                     stringsAsFactors = FALSE)

            } else {
                # Append to current row
                out_df <- data.frame(df_subset,
                                     pval = window_pval,
                                     kernel = object@distmethod,
                                     omnibus = NA,
                                     stringsAsFactors = FALSE)
            }

            return(out_df)
        }

        if (! is.null(object@dbsrc)) {
            # Retrieve the number of results in the table
            # Note: using _ROWID_ instead of COUNT to avoid the full scan table done
            # by sqlite. It should be noted that the database must always start
            # from a clean state.
            total_results <- as.numeric(dplyr::collect(
                dplyr::tbl(object@dbsrc,
                           dplyr::sql("SELECT MAX(_ROWID_) FROM results LIMIT 1"))
            ))

            # In this case we also use _ROWID_ as 'intervals' table cannot have
            # duplicates.
            total_intervals <- as.numeric(dplyr::collect(
                dplyr::tbl(object@dbsrc,
                           dplyr::sql("SELECT MAX(_ROWID_) FROM intervals LIMIT 1"))
            ))

            # Check that there are intervals on the table
            if (is.na(total_intervals))
                stop("There are no intervals on the table, create them first.")

            if (is.na(total_results)) {
                # To avoid retrieving all the rows, we manually select them in
                # chunks of 100 000 rows.
                db_limit <- 100000

                # Modify the object data to be a real data frame (supporting slice)
                object@data <- dplyr::collect(object@data, n = Inf)

                offset_list <- seq(from = 0, to = total_intervals, by = db_limit)

                for (sql_offset in offset_list) {

                    sql_query <- sprintf("SELECT * FROM intervals LIMIT %d OFFSET %d",
                                         db_limit,
                                         sql_offset)

                    message(sprintf("Analyzing chunk %d-%d from a total of %d",
                                    sql_offset,
                                    sql_offset + db_limit,
                                    total_intervals))

                    # As the limit is 100 000, the data frame is returned.
                    chunk_df <- dplyr::tbl(object@dbsrc, dplyr::sql(sql_query)) %>%
                        dplyr::collect(n = Inf) %>%
                        dplyr::rowwise() %>%
                        dplyr::do(.apply.rkat(.)) %>%
                        dplyr::ungroup()

                    # Save the 100 000 rows data frame to database
                    dplyr::db_insert_into(object@dbsrc$con, "results", chunk_df)
                }
            } else {
                # Be sure that all the results are present
                if (is.na(total_intervals) || total_intervals > total_results)
                    stop(sprintf("The number of results (%d) does not match the number of intervals (%d), manual intervention required.", total_results, total_intervals))

                message(sprintf("Using the previously stored %d results.",
                             total_results))
            }

            # Point to the results table
            intervals_pvals <- dplyr::tbl(object@dbsrc, "results")
        } else {
            # Check that there are intervals on the table
            if (is.null(object@intervals))
                stop("There are no intervals in the slot, create them first.")

            # FOR EVERY WINDOW:
            # Create the matrix distances.
            # Convert to kernel
            # Execute MiRKAT function (or adapt it before?)
            # Return a list of p-values
            intervals_pvals <- object@intervals %>%
                dplyr::rowwise() %>%
                dplyr::do(.apply.rkat(.)) %>%
                dplyr::ungroup()
        }
    }

    # Add the annotation information
    row_annotation <- object@annotation %>%
        dplyr::select(chr, pos, row) %>%
        dplyr::distinct()

    output_result <- intervals_pvals %>%
        dplyr::inner_join(row_annotation, by = c("first_row" = "row")) %>%
        dplyr::rename(start = pos) %>%
        dplyr::inner_join(row_annotation, by = c("last_row" = "row")) %>%
        dplyr::rename(end = pos, chr = chr.x) %>%
        dplyr::select(first_row, last_row, start, end, chr, pval, kernel, omnibus)

    return(output_result)
})

#' Method to create the intervals based on the configuration options
#'
#' @param MetRKAT
#'
#' @return
#' @export
#'
#' @examples
setMethod("createIntervals", signature="MetRKAT", function(object, overwrite = FALSE) {

    # Check for existing slots on the instance instead of repeating the process.
    # This allows to "reset" SQLite queries and solve expired connections errors.
    window_intervals <- NULL

    if (! overwrite) {
        if (is.null(object@dbsrc)) {
            # If no database has been used, return the object directly as there is no
            # need to recreate the "query".
            if (! is.null(object@intervals))
                return(object@intervals)
        } else {
            n_intervals <- as.numeric(dplyr::collect(dplyr::tbl(
                object@dbsrc,
                dplyr::sql("SELECT MAX(_ROWID_) FROM intervals LIMIT 1")
            )))

            if (! is.na(n_intervals))
                window_intervals <- dplyr::tbl(object@dbsrc, "intervals")
        }
    }

    if (is.null(window_intervals)) {
        message("Creating window intervals...")

        if (! is.null(object@dbsrc)) {
            # Modify the object data to be a real data frame (supporting slice)
            object@annotation <- dplyr::collect(object@annotation, n = Inf)
        }

        .overlap.generic <- function(df_subset, slide_function, group_function) {
            # Initialize empty output dataframe
            chr_df_out <- data.frame()

            # If a gap value is provided, create different overlapped windows
            slide_range <- c(0)

            if (object@gap > 0) {
                slide_range <- slide_function(df_subset)
            }

            # Some slide methods require the chromosomic position as well as
            # row position, so for those that do not return a list, include
            # dummy values.
            if (! is.list(slide_range)) {
                slide_range <- list(
                    rows = slide_range,
                    pos = rep(NA, length(slide_range))
                )
            }

            chr_df_out <- dplyr::bind_rows(mapply(function(start_row, start_pos) {
                # Select a row subset from a starting point to the end
                rows_slide <- df_subset[start_row:nrow(df_subset), ]

                df_slide <- rows_slide %>%
                        dplyr::group_by(group_function(row_subset = rows_slide, gap_pos = start_pos)) %>%
                        dplyr::summarise(first_row = min(row), last_row = max(row)) %>%
                        dplyr::select(first_row, last_row)

                # If a database connection is provided, insert the records into
                # the "intervals" table instead of returning them.
                if (! is.null(object@dbsrc)) {
                    # Use raw SQL to use the combined UNIQUE index and ignore
                    # if the combination of rows already exists on the table.
                    DBI::dbBegin(object@dbsrc$con)

                    # Note: sqlite does have a maximum number of multiple inserts
                    # per query, that cannot be changed at runtime so we need
                    # to split the rows either individually or in groups of 500.
                    slide_chunks <- split(df_slide, ceiling(1:nrow(df_slide)/450))

                    for (df_chunk in slide_chunks) {
                        sql_query <- sprintf("INSERT OR IGNORE INTO intervals (first_row, last_row) VALUES (%s);",
                                             paste(df_chunk$first_row, df_chunk$last_row, sep=",",
                                                   collapse="),("))

                        DBI::dbExecute(object@dbsrc$con, sql_query)
                    }

                    DBI::dbCommit(object@dbsrc$con)

                    # Return empty data.frame
                    df_slide <- data.frame()
                }

                return(df_slide)
            }, slide_range$rows, slide_range$pos, SIMPLIFY = FALSE))

            return(chr_df_out)
        }

        # Default interval method based on window size and overlap
        .split.default <- function(submethod = NULL) {

            .overlap.default <- function(row_subset, gap_pos = NULL) {
                return(floor(seq_along(row_subset$row)/ object@wsize))
            }

            .slide.default <- function(df_subset) {
                # Move the window using only row positions by a certain gap
                # TODO: increase "gap" + 1 (example wsize3 gap=2)
                # CHANGE FROM 0 to 1??
                slide_rows <- seq(from = 1, to = nrow(df_subset), by = object@gap)

                return(slide_rows)
            }

            # Annotation already contains a column 'row'
            default_df <- object@annotation %>%
                dplyr::select(site, chr, row) %>%
                dplyr::distinct() %>%
                dplyr::group_by(chr) %>%
                dplyr::arrange(row) %>%
                dplyr::do({
                    .overlap.generic(., .slide.default, .overlap.default)
                }) %>%
                dplyr::ungroup()

            return(default_df)
        }

        # Windows defined by genes
        .split.genes <- function(submethod = NULL) {

            # Group by gene, then select tha minimum and maximum row positions
            # of every group.
            genes_df <- object@annotation %>%
                dplyr::group_by(gene) %>%
                dplyr::summarise(first_row = min(row), last_row = max(row))

            return(genes_df)
        }

        # Windows defined by distance between positions
        .split.location <- function(submethod = NULL) {

            # Classify into groups based on cummulative sum of the diference
            # between each row (already ordered in ascending direction)
            #
            # Example:
            #   Positions: 400 1600 3600 6400
            #   Diff: 0 1200 2000 2800
            #   Cumsum: 0 1200 3200 6000
            #   Groups(Wsize=5000): 0 0 0 1
            .overlap.location <- function(row_subset, gap_pos = NULL) {
                return(floor(c(0, cumsum(diff(row_subset$pos))) / object@wsize))
            }

            .overlap.location.reset <- function(row_subset) {
                # "Temporary" implementation. Try to use primitives later
                global_diff <- 0

                # Diff between positions
                pos_diffs <- c(0, diff(row_subset))

                # Groups
                pos_groups <- c(0)

                # Current group
                current_group <- 0

                for (i in 2:length(row_subset)) {
                    # 0 250  500  600  400  500 1000
                    # 250 + 0 = 250 > 1000 --> FALSE --> current_group
                    # 500 + 250 = 750 > 1000 --> FALSE --> current group
                    # 600 + 750 = 1350 > 100 --> TRUE --> current_group + 1 & global_diff <- 0
                    # 400 + 0 = 400 > 1000 --> FALSE
                    global_diff <- global_diff + pos_diffs[i]

                    if (global_diff > object@wsize) {
                        global_diff <- 0
                        current_group <- current_group + 1
                    }

                    pos_groups[i] <- current_group
                }

                return(pos_groups)
            }

            .overlap.location.fixed <- function(row_subset, gap_pos) {
                # Consider an additional difference between the start of the window
                # and the position of the first site in the subset.
                start_dif <- row_subset$pos[1] - gap_pos

                # Append that difference to the others
                return(floor(c(start_dif, cumsum(diff(row_subset$pos) + start_dif)) / object@wsize))
            }

            .slide.location <- function(df_subset) {
                # Move the window using only row positions by a certain gap
                # Calculate the starting position point from each window, using
                # the minimum position of the chromosome as a reference.
                slide_rows_pos <- seq(from = min(df_subset$pos),
                                      to = max(df_subset$pos),
                                      by = object@gap)

                # Transform the positions to the number of rows
                slide_rows_intervals <- findInterval(df_subset$pos,
                                                      slide_rows_pos)

                # Select the first match of every interval, which will correspond
                # to the row position
                slide_rows <- match(unique(slide_rows_intervals),
                                        slide_rows_intervals)

                return(slide_rows)
            }

            .slide.location.fixed <- function(df_subset) {
                # Move the window using only row positions by a certain gap
                # Calculate the starting position point from each window, using
                # the minimum position of the chromosome as a reference.
                slide_positions <- seq(from = min(df_subset$pos),
                                       to = max(df_subset$pos),
                                       by = object@gap)

                # Find every interval to match the nearest window that will use.
                # It differs from the default location in which no gaps are discarded, even
                # when they will point to the same first site, however with this fixed
                # location option, the gap is used to distribute the groups, so the same
                # set of rows could have different window distribution.
                #
                # Add 1 to select the higher next position
                slide_rows <- findInterval(slide_positions,
                                           df_subset$pos,
                                           left.open = TRUE)

                # With fixed locations there is no skipping gap, so return all of them
                slide_return <- list(
                    rows = slide_rows,
                    pos = slide_positions
                )

                return(slide_return)
            }

            # Assign the default method and override if submethod is passed
            # and an exists
            slide_method <- get0(sprintf(".slide.location.%s", submethod),
                                 ifnotfound = .slide.location)

            overlap_method <- get0(sprintf(".overlap.location.%s", submethod),
                                  ifnotfound = .overlap.location)

            position_df <- object@annotation %>%
                dplyr::select(site, chr, row, pos) %>%
                dplyr::distinct() %>%
                dplyr::group_by(chr) %>%
                dplyr::arrange(pos) %>%
                dplyr::do({
                   .overlap.generic(., slide_method, overlap_method)
                }) %>%
                dplyr::ungroup()

            return(position_df)
        }

        # The 'wmethod' parameters defines the function to call. The underscore
        # '_' on the string is used to call 'submethods' inside it.
        window_function <- unlist(strsplit(object@wmethod, '_'))

        # Set the default submethod parameter to NULL
        if (length(window_function) < 2) {
            window_function[2] <- NA
        }

        # Call the proper function
        window_intervals <- do.call(sprintf(".split.%s", window_function[1]),
                                    list("submethod" = window_function[2]))

        # Return the table pointer if a database has been used.
        if (! is.null(object@dbsrc)) {
            window_intervals <- dplyr::tbl(object@dbsrc, "intervals")
        } else {
            # Select only "first_row" and "last_row"
            window_intervals <- dplyr::select(window_intervals, first_row, last_row)#, start_row)
        }
    }

    window_intervals <- dplyr::select(window_intervals, first_row, last_row)

    return(window_intervals)
})

#' Generates a manhattan plot
#'
#' @param MetRKAT
#'
#' @return
#' @export
#'
#' @examples
setMethod("plotManhattan", signature="MetRKAT", function(object, pvals) {
  library(dplyr)
  object@results<-as.data.frame(object@results)
    # Prepare GRanges object
    chrToNumber <- function(chroms) {
        # Remove 'chr' prefix
        chroms <- gsub("chr", "", chroms, ignore.case = TRUE)

        # Replace X, Y by numbers
        chroms[tolower(chroms) == 'x'] <- 23
        chroms[tolower(chroms) == 'y'] <- 24

        # Convert to numeric values for qqman
        chroms <- as.numeric(chroms)

        return(chroms)
    }

    # Use the omnibus p-values
    is_omnibus <- pvals == 'omnibus'

    results_filtered <- object@results

    if (! is_omnibus) {
        object@results <- dplyr::filter(object@results, kernel == pvals)
    }
   if (is_omnibus=="omnibus") {
  	P =object@results$omnibus} else{ P=object@results$pval}

manhattan_df <- object@results %>%
  dplyr::select(CHR = chr, BP = start) %>%
  dplyr::mutate(CHR = chrToNumber(CHR), SNP = NA ,P=P) %>%
  dplyr::distinct() %>%
  dplyr::collect(n = Inf)

    # Prevent 'Inf' values when doing log10
    manhattan_df$P[manhattan_df$P == 0] <- 1E-7

    qqman::manhattan(manhattan_df)
})

#' Generates a continuous window plot
#'
#' @param MetRKAT
#'
#' @return
#' @export
#'
#' @examples
setMethod("plotWindows", signature="MetRKAT", function(object, chrom,
                                                       pvals,
                                                       cutoff = 0.05,
                                                       startpos = 0,
                                                       endpos = Inf,
                                                       logp = TRUE) {
    # Search like character
    chrom <- as.character(chrom)

    # Convert to numeric
    cutoff <- as.numeric(cutoff)

    is_omnibus <- pvals == 'omnibus'

    # Select the proper kernel
    if (! is_omnibus) {
        filtered_results <- dplyr::filter(object@results, kernel == pvals)
    } else {
        # Assign the pval as the omnibus
        filtered_results <- dplyr::mutate(object@results, pval = omnibus)
    }

    # Filter results
    # TODO: limit more the results as there are too many to list.
    filtered_results <- filtered_results %>%
        dplyr::filter(chr == chrom, pval <= cutoff, start >= startpos, end <= endpos) %>%
        dplyr::distinct(chr, start, end, pval) %>%
        dplyr::collect(n = Inf)

    # Merge the start-end window position into one string for X axisss
    wplot_df <- tidyr::unite(filtered_results, "window", start, end)

    # Plot can use factors as x values but does not draw a continuous line
    wplot_df$window_num <- 1:length(wplot_df$window)

    # Disable x labels as we need to provide numeric x values to draw a
    # continuous line
    # TODO: rotate x-axis labels?
    ylab_text <- "p-value"

    if (logp) {
        wplot_df$pval[wplot_df$pval == 0] <- 1E-7
        wplot_df$pval <- -log10(wplot_df$pval)

        ylab_text <- expression(-log[10](pval))
    }

    plot(pval ~ window_num, dat=wplot_df, type="l", ylab=ylab_text, xlab="", xaxt="n" )
    axis(1, at=seq_along(wplot_df$pval), labels = FALSE)
    text(x=seq_along(wplot_df$pval), y=par()$usr[3]-0.06*(par()$usr[4]-par()$usr[3]),
         labels=wplot_df$window, srt=45, adj=1, xpd=TRUE, cex=0.5)
})

#' Generates a chromosome ideogram plot
#'
#' @param MetRKAT
#'
#' @return
#' @export
#'
#' @examples
setMethod("plotChromosome", signature="MetRKAT", function(object, chrom,
                                                          pvals,
                                                          cutoff = 0.05,
                                                          logp = TRUE,
                                                          hg_ideo = NULL,
                                                          ...) {
    # Prepare GRanges object
    chrToString <- function(chroms) {
        # Remove chr prefix in case it exists
        chroms <- gsub("chr", "", chroms, ignore.case = TRUE)

        # Make a numeric copy (if 'x/y' exists, as.numeric would turn out
        # to NA)
        chroms.num <- as.numeric(chroms)

        # Reverse numbers to X/Y if they exists
        chroms[chroms.num == 23] <- "X"
        chroms[chroms.num == 24] <- "Y"

        chroms <- paste0("chr", chroms)

        return(chroms)
    }

    # Convert it to character
    chrom <- as.character(chrom)

    # Convert to numeric
    cutoff <- as.numeric(cutoff)

    is_omnibus <- pvals == 'omnibus'

    # Select the proper kernel
    if (! is_omnibus) {
        filtered_results <- dplyr::filter(object@results, kernel == pvals)
    } else {
        # Assign the pval as the omnibus
        filtered_results <- dplyr::mutate(object@results, pval = omnibus)
    }

    filtered_results <-  filtered_results %>%
        dplyr::filter(chr == chrom, pval <= cutoff) %>%
        dplyr::distinct(chr, start, end, pval) %>%
        dplyr::collect(n = Inf)

    if (logp) {
        filtered_results$pval[filtered_results$pval == 0] <- 1E-7
        filtered_results$pval <- -log10(filtered_results$pval)
    }

    # TODO: make only the chromosome specified
    granges_result <- GenomicRanges::GRanges(seqnames = chrToString(filtered_results$chr),
                              ranges = IRanges::IRanges(start = filtered_results$start,
                                               end = filtered_results$end),
                              pval = filtered_results$pval)

    # TODO: allow multiple versions?
    if (is.null(hg_ideo)) {
        hg_ideo <- IdeoViz::getIdeo("hg19")
    }

    # Correct format of the chromosome
    if (is.numeric(chrom)) {
        chrom <- paste0("chr", chrom)
    }

    # Extra parameters
    plot_params <- list(...)

    # Plot default parameters
    plot_defaults <- list(
        ylab = if (logp) "log10(p-values)" else "p-values",
        plot_title = sprintf("%s windows", chrom),
        plotType = "lines",
        chrom = chrom,
        ideoTable = hg_ideo,
        values_GR = granges_result,
        value_cols = colnames(mcols(granges_result)),
        col = RColorBrewer::brewer.pal(n = 5, "Spectral")
    )

    # Override with the passed arguments
    plot_defaults <- modifyList(plot_defaults, plot_params)

    do.call(IdeoViz::plotOnIdeo, plot_defaults)
})

#' Transforms the object to use a SQLite database
#'
#' @param MetRKAT
#'
#' @return
#' @export
#'
#' @examples
setMethod("toSQLite", signature="MetRKAT", function(object, dbname = NULL, overwrite = FALSE) {
    # Try to restore the connection
    if (is.null(dbname)) {
        # Check that it really exists
        if (is.null(object@dbsrc))
            stop("No previous connection detected, you must provide a db name.")

        # Set the dbname to the one previously used
        dbname <- object@dbsrc$path
    }

    object@dbsrc <- dplyr::src_sqlite(dbname, create = TRUE)

    # Check for previously created tables
    db_tables <- dplyr::db_list_tables(object@dbsrc$con)

    if (overwrite || ! "data" %in% db_tables) {
        # Copy data
        message("Copying 'data' to sqlite database... ", appendLF = FALSE)

        object@data <- dplyr::copy_to(object@dbsrc, object@data,
                                      temporary = FALSE, name = "data")

        message("Done!")
    } else {
        object@data <- dplyr::tbl(object@dbsrc, "data")
    }

    if (overwrite || ! "annotation" %in% db_tables) {
        # Copy annotation
        message("Copying 'annotation' to sqlite database... ", appendLF = FALSE)

        object@annotation <- dplyr::copy_to(object@dbsrc, object@annotation, temporary = FALSE,
                                      name = "annotation",
                                      indexes = list("site", "chr", "pos", "row"))

        message("Done!")
    } else {
        object@annotation <- dplyr::tbl(object@dbsrc, "annotation")
    }

    if (overwrite || ! "intervals" %in% db_tables) {
        # Create intervals table
        message("Creating empty 'intervals' table... ", appendLF = FALSE)

        dplyr::db_create_table(object@dbsrc$con, "intervals", temporary = FALSE,
                               types = c(first_row = "number",
                                         last_row = "number"))

        dplyr::db_create_indexes(object@dbsrc$con, "intervals",
                                 indexes = list(c("first_row", "last_row")),
                                 unique = TRUE)

        message("Done!")
    } else {
        object@intervals <- createIntervals(object, overwrite = FALSE)
    }

    if (overwrite || ! "results" %in% db_tables) {
        # Create intervals table
        message("Creating empty 'results' table... ", appendLF = FALSE)

        dplyr::db_create_table(object@dbsrc$con, "results", temporary = FALSE,
                               types = c(first_row = "number",
                                         last_row = "number",
                                         pval = "real",
                                         kernel = "varchar(25)",
                                         omnibus = "real"))

        message("Done!")
    } else {
        #object@results <- applyRKAT(object, overwrite = FALSE)
    }

    return(object)
    object@results <- applyRKAT(object, overwrite = FALSE)
})

#' Main helper function to perform all the required steps for the analysis.
#'
#' @param y
#' @param data
#' @param wsize
#' @param wmethod
#' @param overlap
#' @param annotation
#'
#' @return
#' @export
#'
#' @examples
MetRKAT <- function(y, data, wsize = 1000, wmethod = "default",
                    gap = 0, distmethod = "euclidean", annotation = data.frame(),
                    max.na = 0.5, dbname = NULL, ...) {

    analysis <- new("MetRKAT",
                  data = data,
                  wsize = wsize,
                  wmethod = wmethod,
                  distmethod = distmethod,
                  gap = gap,
                  max.na = max.na,
                  annotation = annotation)

    if (! is.null(dbname)) {
        analysis <- toSQLite(analysis, dbname)
    }

    analysis@intervals <- createIntervals(analysis)

    analysis@results <- applyRKAT(analysis, y = y)

    message("Done! NOTE: the data provided probably has a different row order, so to use the row position information of the results the reordered data should be retrieved by using the slot 'data'.")

    return(analysis)
}
