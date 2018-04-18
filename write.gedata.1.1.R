#' Write data and annotations to a .gedata v 1.1 file (Qlucore Omics Explorer data format)
#' 
#' Given a data matrix and sample and variable annotations, create a .gedata file containing this information
#' 
#' @param in.data.matrix A data matrix, with rows corresponding to variables and columns corresponding to samples
#' @param in.var.annots A data frame of variable annotations
#' @param in.sample.annots A data frame of sample annotations
#' @param outfile The output gedata file
#' @author Charlotte Soneson, Martin Hjelmstedt
#' @export
#' @return Nothing is returned, the information is written to a gedata file
write.gedata <- function(in.data.matrix, in.var.annots, in.sample.annots, outfile) {

  data.matrix <- as.matrix(as.data.frame(in.data.matrix))
  var.annots <- as.data.frame(in.var.annots)
  sample.annots <- as.data.frame(in.sample.annots)
  
  n.samples <- ncol(data.matrix)
  n.vars <- nrow(data.matrix)

  if (!is.data.frame(var.annots))
    stop("Error: var.annots must be a data.frame")
  if (!is.data.frame(sample.annots))
    stop("Error: sample.annots must be a data.frame")

  if (n.vars == 0 || n.samples == 0)
    stop('Error: Empty data matrix!')

  if (n.vars != nrow(var.annots))
    stop('Error: Dimensions of data matrix and variable annotations are not compatible!')
  if (n.samples != nrow(sample.annots))
    stop('Error: Dimensions of data matrix and sample annotations are not compatible!')

  if(!grepl(".gedata$", outfile))
    stop('Error: Output file name must be of type .gedata!')
  
  sample.annots <- t(sample.annots)
  
  fileheader2 <- paste('qlucore\tgedata\tversion 1.1\n\n\n', n.samples, 
                       '\tsamples\twith\t', nrow(sample.annots), 
                       '\tannotations\n', n.vars, '\tvariables\twith\t', 
                       ncol(var.annots), '\tannotations\n\n', sep = '')
  
  annotations <- rbind(cbind(matrix('', nrow(sample.annots), 
                                    ncol(var.annots)), 
                             rownames(sample.annots), sample.annots))
  
  thedata <- cbind(var.annots, rep('', n.vars), data.matrix)
  
  output.file <- file(outfile, 'w')
  cat(fileheader2, file = output.file)
  write.table(annotations, file = output.file, 
              quote = FALSE, row.names = FALSE, 
              col.names = FALSE, sep = '\t')
  cat(colnames(var.annots), file = output.file, sep = '\t')
  cat('\n', file = output.file)
  write.table(rbind(thedata), file = output.file, 
              quote = FALSE, row.names = FALSE, 
              col.names = FALSE, sep = '\t')
  close(output.file)
}
