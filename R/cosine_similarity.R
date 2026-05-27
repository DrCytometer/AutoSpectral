# cosine_similarity.r

#' @title Calculate Cosine Similarity
#' @description
#' Calculates the cosine similarity between columns of the input matrix (spectra).
#'
#' @param spectra The spectral matrix (or dataframe), represented as
#' fluorophores in rows and detectors in columns.
#'
#' @return The cosine similarity matrix in fluorophores x fluorophores.
#' @export

cosine.similarity <- function( spectra ) {
  spectra.t <- as.matrix( t( spectra ) )
  euclidean.norm <- sqrt( colSums( spectra.t^2 ) + 1e-9 )
  dot.product <- t( spectra.t ) %*% spectra.t
  similarity.matrix <- dot.product / outer( euclidean.norm, euclidean.norm )
  return( similarity.matrix )
}

#' @title dot Cosine Similarity Rows
#' @description
#' Cosine similarity of each row of mat against a single reference vector.
#'
#' @param mat The matrix (or dataframe), represented as events in rows and
#' detectors in columns.
#' @param ref.vec The vector (numeric) of the reference spectrum to which the
#' rows of the matrix will be compared.
#'
#' @return Returns a numeric vector of length nrow(mat).

.cosine.sim.rows <- function( mat, ref.vec ) {
  dot.prod <- mat %*% ref.vec
  mat.norm <- sqrt( rowSums( mat^2 ) )
  ref.norm <- sqrt( sum( ref.vec^2 ) )
  as.numeric( dot.prod / ( mat.norm * ref.norm + 1e-9 ) )
}
