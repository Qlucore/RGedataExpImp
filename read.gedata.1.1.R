#' Read a .gedata file (Qlucore Omics Explorer data format) and extract the information
#' 
#' Given a .gedata file, read it and extract the data, sample annotations and
#' variable annotations. If \code{out.rds} is not NULL, save the objects in an
#' .rds file. Otherwise, return (invisibly) the data matrix and two data frames
#' containing the annotations.
#' 
#' @param in.gedata The path to a .gedata file
#' @param out.rds The name of the .rds file where the data matrix and the
#'   annotation data frames should be saved
#' @export
#' @author Charlotte Soneson, Martin Hjelmstedt
#' @return Returns invisibly a list containing three elements:
#' \itemize{
#' \item \code{data} - the data matrix (variables as rows, samples as columns)
#' \item \code{samp.annot} - a data frame with sample annotations
#' \item \code{var.annot} - a data frame with variable annotations
#' }
read.gedata <- function(in.gedata, out.rds = NULL) {
  ## Find the number of samples and variables,
  ## and the number of annotations for each of them.
  n.col <- max(count.fields(in.gedata, sep = "\t", quote = ""))
  X <- read.delim(in.gedata, sep = "\t",
                  fill = TRUE, blank.lines.skip = FALSE,
                  check.names = FALSE, 
                  as.is = TRUE, header = FALSE, 
                  col.names = 1:n.col, quote = "")

  i <- 1
  while (X[i, 2] != "samples") {
    i <- i + 1
    if (i > nrow(X)) {
      stop("Error: meta data not in gedata v 1.1 format")
    }
  }
  saminfo.rownumber <- i
  varinfo.rownumber <- saminfo.rownumber + 1
  n.metarows <- varinfo.rownumber + 1

  nbr.samples <- as.numeric(X[saminfo.rownumber, 1])
  nbr.sample.annotations <- as.numeric(X[saminfo.rownumber, 4])
  nbr.variables <- as.numeric(X[varinfo.rownumber, 1])
  nbr.variable.annotations <- as.numeric(X[varinfo.rownumber, 4])
  
  ## Extract data values
  X.data <- X[-(1:(nbr.sample.annotations + n.metarows + 1)),
              -(1:(nbr.variable.annotations + 1))]
  rownames(X.data) <- X[-(1:(nbr.sample.annotations + n.metarows + 1)), 1]

  X.data <- as.matrix(X.data)
  mode(X.data) <- 'numeric'
  
  colnames(X.data) <- X[n.metarows + 1, -(1:(nbr.variable.annotations + 1))]
  
  ## Extract variable annotations
  X.var.annot <- cbind(X[-(1:(nbr.sample.annotations + n.metarows + 1)),
                         1:(nbr.variable.annotations)])
  rownames(X.var.annot) <- X[-(1:(nbr.sample.annotations + n.metarows + 1)), 1]
  colnames(X.var.annot) <- X[nbr.sample.annotations + n.metarows + 1,
                             1:(nbr.variable.annotations)]
  X.var.annot <- as.data.frame(X.var.annot)

  ## Extract sample annotations
  X.samp.annot <- X[(n.metarows + 1):(nbr.sample.annotations + n.metarows),
                    -(1:(nbr.variable.annotations + 1))]
  rownames(X.samp.annot) <- X[(n.metarows + 1):(nbr.sample.annotations + n.metarows),
                              (nbr.variable.annotations + 1)]
  colnames(X.samp.annot) <- X[(n.metarows + 1), -(1:(nbr.variable.annotations + 1))]
  X.samp.annot <- as.data.frame(t(X.samp.annot))
  
  out.list <- list(data = X.data, var.annot = X.var.annot, samp.annot = X.samp.annot)
  
  ## Save
  if (!is.null(out.rds))
    saveRDS(out.list, file = out.rds)
  
  return(invisible(out.list))
}
