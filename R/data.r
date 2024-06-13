#' Replication timing data
#'
#' Chromosome 10 Wavelet-smoothed signal data from MCF7 cells as provided by the UCSC genome browser (see source)
#'
#' @format ## `mcf7Rept`
#' A `data.table` with 133,719 rows and 4 columns:
#' \describe{
#'   \item{chr:}{Chromosome name}
#'   \item{start:}{start positions of genomic ranges (open)}
#'   \item{end:}{end positions of genomic ranges (closed)}
#'   \item{signal:}{replication timing signal}
#' }
#'
#' @source <https://genome-euro.ucsc.edu/cgi-bin/hgTables?db=hg19&hgta_group=regulation&hgta_track=wgEncodeUwRepliSeq&hgta_table=wgEncodeUwRepliSeqMcf7WaveSignalRep1&hgta_doSchema=describe+table+schema>
"mcf7Rept"
