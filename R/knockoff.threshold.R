#' Threshold for the knockoff filter
#'
#' Computes the threshold for the knockoff filter.
#'
#' @param W the test statistics
#' @param fdr target false discovery rate (default: 0.1)
#' @param offset either 0 or 1 (default: 1). The offset used to compute the rejection threshold on the
#' statistics. The value 1 yields a slightly more conservative procedure ("knockoffs+") that
#' controls the FDR according to the usual definition, while an offset of 0 controls a modified FDR.
#' @return The threshold for variable selection.
#'
#' @export
#' @md
knockoff.threshold <- function(W, fdr=0.10, offset=1) {
  if(offset!=1 && offset!=0) {
    stop('Input offset must be either 0 or 1')
  }
  ts = sort(c(0, abs(W)))
  ratio = sapply(ts, function(t)
    (offset + sum(W <= -t)) / max(1, sum(W >= t)))
  ok = which(ratio <= fdr)
  ifelse(length(ok) > 0, ts[ok[1]], Inf)
}


#' Select Variables based on knockoff statistics
#'
#'
#' @param Ws the test statistics, it is either a nrep by p matrix, or list of length to be nrep, or a vector with length being p.
#' @param fdr target false discovery rate (default: 0.1)
#' @param offset either 0 or 1 (default: 1). The offset used to compute the rejection threshold on the
#' statistics. The value 1 yields a slightly more conservative procedure ("knockoffs+") that
#' controls the FDR according to the usual definition, while an offset of 0 controls a modified FDR.
#' @return An object of class "knockoff.select. This object is a list
#'  containing at least the following components:
#'  \item{W}{computed test statistics}
#'  \item{Ws}{matrix of computed test statistics}
#'  \item{thre}{computed selection threshold}
#'  \item{index}{index of selected variables}
#'
#' @export
#' @md
knockoff.select <- function(Ws, fdr=0.10, offset=1) {
  X.names = colnames(Ws)
  if(is.vector(Ws)){
    W = Ws
  }else{
    W = colMeans(Ws)
  }
  if(offset!=1 && offset!=0) {
    stop('Input offset must be either 0 or 1')
  }

  # filter
  thre <- knockoff.threshold(W,fdr,offset)
  s <- which(W >= thre)
  if (!is.null(X.names))
    names(s) = X.names[s]

  structure(list(W = W,
                 Ws = Ws,
                 t = thre,
                 s = s),
            class = 'knockoff.select')
}


#' Print results for the multiple knockoff filter
#'
#' Prints the list of variables selected by the knockoff filter and the corresponding function call.
#'
#' @param x the output of a call to knockoff.filter
#' @param ... unused
#'
#' @method print knockoff.filter
#' @export
#' @md
print.knockoff.filter<- function(x, ...) {
  cat('Call:\n')
  print(x$call)
  cat('\nSelected variables:\n')
  print(x$shat)
  if(!is.null(x$shat_list)){
    #cat('\nSelected variables for each knockoff copy:\n')
    #print(x$shat_list)
    cat('\nFrequency of selected variables from', length(x$shat_list) ,'knockoff copys:\n')
    print(colSums(x$shat_mat))
  }
}

#' Verify dependencies for chosen statistics
#'
#' @param statistic the statistic chosen by the user
#'
#' @keywords internal
verify_stat_depends <- function(statistic) {

}
