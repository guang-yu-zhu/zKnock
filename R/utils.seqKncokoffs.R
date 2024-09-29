#' Internal check of whether class of response is compatible with family
#'
#' Do not call this function on its own
#'
#' @param y response variable
#' @param type "regression", "classification", or "survival"
#'
#' @return this function doesn't return anything
#' @export
#'
#' @keywords internal
check_family <- function(y, type) {
  is_regression <- type == "regression" && is.numeric(y) && !inherits(y, "Surv")
  is_classification <- type == "classification" && is.factor(y) && length(unique(y)) == 2
  is_survival <- type == "survival" && inherits(y, "Surv")

  if (!(is_regression || is_classification || is_survival)) {
    stop("One of the following must hold for the input:\n",
         "1) type = 'regression' and class(y) = 'numeric',\n",
         "2) type = 'classification' and class(y = 'factor' with length(unique(y)) = 2, or\n",
         "3) type = 'survival' and class(y) = 'Surv' from the 'survival' package.")
  }
}

#' Internal check of whether input data frame (or tibble) is of the right format
#'
#' Do not call this function on its own
#'
#' @param X data frame or tibble
#' @param method character string, either "seq" or "mx"
#'
#' @return this function doesn't return anything
#' @export
#'
#' @keywords internal
check_design <- function(X, method="seq", check.dim=TRUE) {

  if(!("data.frame" %in% class(X))) {
    stop(paste0(deparse(substitute(X)), " should be either a data.frame or tibble"))
  }

  if(check.dim & ncol(X)<=2) {
    stop(paste0(deparse(substitute(X)), " should have ncol(X) > 2"))
  }

  if(method=="seq" & sum(!unlist(lapply(X, class)) %in% c("factor", "numeric")) > 0) {
    stop(paste0(deparse(substitute(X)), " should only contain columns of class 'numeric' or 'factor'"))
  }

  if(method=="mx" & sum(!unlist(lapply(X, class)) %in% c("numeric")) > 0) {
    stop(paste0(deparse(substitute(X)), " should only contain columns of class 'numeric'"))
  }

}

#' Internal check of normality of the numeric input covariates
#'
#' Do not call this function on its own
#'
#' @param X data frame or tibble
#'
#' @return this function doesn't return anything
#' @export
#'
#' @keywords internal
check_normality <- function(X) {

  X_numeric <- dplyr::select_if(X, is.numeric)

  is.distinct <- unlist(lapply(X_numeric, dplyr::n_distinct))==nrow(X_numeric)

  p.values <- unlist(lapply(X_numeric, function(x) suppressWarnings(ks.test(x, y="pnorm", mean=mean(x), sd=sd(x))$p.value)))

  is.normal <- p.values >= 0.05

  if (!is.null(is.normal) & sum(!is.normal) > 0 | sum(!is.distinct) > 0) {

    warning.message <- paste0("Some of the numeric input covariates may have ties and/or may not be normally distributed. This could affect the quality of corresponding knockoffs since they are sampled from a Gaussian distribution. ")

    if (!is.null(is.normal) & sum(!is.normal) > 0) {
      warning.message <- paste0(warning.message, paste(names(which(!is.normal)), collapse=", "), " had normality rejected by Kolmogorov-Smirnov test. ")
    }

    if (sum(!is.distinct & is.normal) > 0) {
      warning.message <- paste0(warning.message, paste(names(which(!is.distinct & is.normal)), collapse=", "), " had ties (but did not reject normality). ")
    }

    warning.message <- paste0(warning.message, "Please consider applying a normalizing transformation on these variables if needed.")

    warning(warning.message)

  }

}


#' Normal score transformation function
#'
#' @param y a numeric vector representing the continuous variable
#'
#' @return a vector of length(y) with the normal-score transformed variable
#' @export
#'
#' @keywords internal
ns.transform <- function(y) {

  # Normal Q-Q plot:
  yt <- qqnorm(y, plot.it=FALSE)$x

  return(yt)

}

#' Heuristic check for whether a variable can be reasonably treated as continuous
#'
#' @param X a numeric variable vector
#'
#' @return a logical TRUE or FALSE depending on whether n_distinct(x) > 30
#' @export
#'
#' @keywords internal
check_if_continuous <- function(X) {
  `%>%` <- magrittr::`%>%`
  X_numeric <- dplyr::select_if(X, is.numeric)
  is.continuous <- sum(X_numeric %>% lapply(dplyr::n_distinct) %>% unlist() <= 30) > 0
  if (is.continuous) warning("Some of the numeric columns of X have suspiciously few distinct values: n_distinct <= 30. Those columns should perhaps not be treated as continuous variables. Please review carefully and read the documentation about the gcm parameter of the knockoff.statistics function.")
}


#' Estimate adjacency matrix using graphical LASSO (glasso)
#'
#' @param X data.frame (or tibble) of covariates
#' @importFrom CVglasso CVglasso
#' @return adjacency matrix
#' @keywords internal
glasso_adjacency_matrix <- function(X){

  if (all(unlist(lapply(X, is.numeric)))){

    precision.matrix <- CVglasso::CVglasso(X, trace="none")$Omega # glasso::glasso(cov(X), rho=rho)$wi
    adjacency.matrix <- precision.matrix + t(precision.matrix)
    adjacency.matrix[adjacency.matrix !=0] = 1

  } else {

    # first, construct model matrix with dummy encoded factor variables
    X.dummy <- model.matrix(~., data=X)

    # The columns of the model matrix X correspond to these original variables (stored in vars):
    assignements <- attributes(X.dummy)$assign[-1]

    # Remove intercept of the model matrix
    X.dummy <- X.dummy[,-1]

    ## then construct its precision matrix
    precision.matrix.dummy <- CVglasso::CVglasso(X.dummy, trace="none")$Omega

    # This is done to guarantee symmetry in the adjacency matrix
    precision.matrix.dummy <- precision.matrix.dummy + t(precision.matrix.dummy)

    ## next construct adjacency matrix of original problem under the assumption that for each pair of nodes (A,B) involving at least one categorical variable, say A, there is an edge iff at least one of the dummy encoded features (of A) has an edge (to (any of the dummy encodings of ) B)
    adjacency.matrix.temp <- matrix(0, nrow=nrow(precision.matrix.dummy), ncol=ncol(X))
    ## combine entries by column
    for (i in c(1:ncol(X))){
      cols <- which(assignements==i)
      if (length(cols) == 1){
        adjacency.matrix.temp[,i] <- precision.matrix.dummy[,cols]
      } else {
        adjacency.matrix.temp[,i] = rowSums(precision.matrix.dummy[,cols])
      }
    }
    ## combine entries by row
    adjacency.matrix <- matrix(0, nrow=ncol(X), ncol=ncol(X))
    for (i in c(1:ncol(X))){
      rows <- which(assignements==i)
      if (length(rows) == 1){
        adjacency.matrix[i,] = adjacency.matrix.temp[rows,]
      } else {
        adjacency.matrix[i,] = colSums(adjacency.matrix.temp[rows,])
      }
    }

    ## convert to proper adjacency matrix (i.e. all entries either 1 or 0)
    adjacency.matrix[adjacency.matrix != 0] = 1

  }
  colnames(adjacency.matrix) <- names(X)
  return (adjacency.matrix)

}
