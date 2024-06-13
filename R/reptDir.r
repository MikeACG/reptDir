#' @importFrom Rdpack reprompt
#' @import data.table

#' @title
#' reptDir: Replication Direction Finder
#' 
#' @description
#' Determine the replication direction of genomic ranges based on replication timing signal.
#'
#' @details
#' The method used for assigning replication direction to each range
#' is inspired by the descriptions of \insertCite{Morganella2016_repdir;textual}{reptDir}.
#' Briefly, replication time signal extrema are found by identifying the first genomic ranges switching their slope
#' from positive to negative (peaks) or visceversa (valleys) after a streak of consecutive ranges with the same slope sign. Slopes
#' are defined as the difference between the signal
#' of any range and the previous neighbor range as defined by standard genomic coordinates in the reference strand.
#' Ranges with 0 slope and whose next neighbor also has 0 slope (flat regions), as well as
#' `extrema` themselves, are assigned `NA` replication direction.
#' Ranges for which a putative direction assignment is possible are
#' therefore those between extrema or between extrema and flat regions.
#' These collections of consecutive ranges are deemed replication direction domains.
#' Additional ranges may be assigned `NA` direction by `minLen`
#' and then by `minSlope` criteria in that order. These parameters
#' correspond to filters applied by \insertCite{Morganella2016_repdir;textual}{reptDir} and
#' \insertCite{Haradhvala2016_repdir;textual}{reptDir} respectively. if a domain spans less than 
#' `minLen` nucleotides, its ranges are determined to have unreliable replication direction.
#' If a range or its next neighbor has a slope less than `minSlope`,
#' its replication direction is deemed unreliable.
#' The remaining domains with only reliable ranges are assigned replication direction
#' based on the peaks being replication origins and valleys termination zones.
#' Therefore, a domain whose ranges have a positive slope indicates replication
#' from peak to valley in the "q-to-p" direction relative to the chromosomal arms.
#' A negative slope indicates peak to valley replication in the "p-to-q" direction.
#' The leading strand is the reference strand in q-to-p replication while
#' the lagging strand is the reference in the p-to-q case.
#'
#' @param reptdt A `data.table` with at least the following columns: \cr
#'  `start`: start positions of genomic ranges (open) \cr
#'  `end`: end positions of genomic ranges (closed) \cr
#'  `signal`: replication time signal. \cr
#' The input `data.table` is assumed to have only genomic ranges
#' of a single chromosome ordered by start position.
#' The ranges are also assumed to be completely disjoint from one
#' another. The replication signal is assumed to be processed
#' already such that it meaningfully depicts the replication time profile of interest.
#' Such data is tipically found already in this format in the
#' UCSC genome browser for example. See the README at
#' \url{https://github.com/MikeACG/reptTools} for a full usage example.
#' @param minLen A scalar integer, minimum number of consecutive bases of direction domains.
#' @param minSlope A scalar numeric, minimum slope at both sides of ranges to be included in direction domains.
#'
#' @return A `data.table` with the following columns in addition to input columns: \cr
#'  `dl`: signal difference between current range and previous range to the left \cr
#'  `dr`: signal difference between next range to the right and current \cr
#'  `extrema`: whether or not the range is a replication peak or valley \cr
#'  `sdl`: sign of `dl` \cr
#'  `sdr`: sign of `dr` \cr
#'  `direction`: replication direction.
#'
#' @references \insertAllCited{}
#'
#' @export
reptDir <- function(reptdt, minLen, minSlope) {

    # copy input to new data.table
    repdDT <- data.table::copy(reptdt)

    # define slope to the left and right of each range
    d <- diff(repdDT$signal)
    repdDT[, ':=' ("dl" = c(Inf, d), "dr" = c(d, Inf))]
    repdDT[, ':=' ("sdl" = c(0L, sign(d)), "sdr" = c(sign(d), 0L))]

    # check each streak of left slopes to find extrema
    rd <- rle(repdDT$sdl)
    cs <- cumsum(rd$lengths)
    extrema <- integer(nrow(repdDT))
    for (ii in 2:(length(rd$values) - 1)) { # edges of chromosome are never extrema

        if ( (rd$values[ii] == 1L & rd$values[ii + 1L] == -1L) | (rd$values[ii] == -1L & rd$values[ii + 1L] == 1L) ) {

            extrema[cs[ii]] <- 1L # valley or peak

        }

        if ( (rd$values[ii - 1L] == 1L & rd$values[ii] == 0L & rd$values[ii + 1L] == -1L) | (rd$values[ii - 1L] == -1L & rd$values[ii] == 0L & rd$values[ii + 1L] == 1L) ) {

            extrema[cs[ii - 1L]:cs[ii]] <- 1L # flat valley or peak

        }

    }
    
    # if range has left or right non-zero slope and is not an extrema then it contributes to consecutive replication in some direction
    repdDT[, "extrema" := extrema]
    check <- ifelse((repdDT$sdl != 0 | repdDT$sdr != 0) & repdDT$extrema != 1, 1L, 0L)

    # assign an id to each consecutive group of ranges to check
    r <- rle(check)
    g <- list()
    for (ii in 1:length(r$values)) {

        if (r$values[ii] == 1L) g[[ii]] <- rep(ii, r$lengths[ii])

    }
    group <- rep(NA_integer_, nrow(repdDT))
    group[check == 1L] <- unlist(g)

    # get the width of each group
    wpg <- tapply(repdDT$end - repdDT$start, group, sum)
    # get ranges in groups with desired width and exclude the ones with insufficient slope
    keep <- group %in% as.integer(names(wpg)[wpg >= minLen])
    keep[abs(repdDT$dl) < minSlope | abs(repdDT$dr) < minSlope] <- FALSE

    # direction of replication, peaks are origins and valleys terminations
    repdDT[, "direction" := rep(NA_character_, nrow(repdDT))]
    repdDT[keep == TRUE, "direction" := ifelse(dl + dr > 0, "q2p", "p2q")]

    return(repdDT)

} 
