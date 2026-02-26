#' Bilateral / anisotropic filtering of gradient field
#' 
#' Gradient fields are smoothed using bilateral filtering,
#' in which the smoothed gradient of each point is computed as
#' the weighted average of the neighbors' gradients, considering
#' both distance in space and also similarity in gradients.
#' 
#' The weight of each neighbor is computed from the product of two scores:
#' * `distance` score: Generally, closer neighbors have greater weight.
#'   * if `'euclidean'`: Gaussian transformation of the Euclidean distance
#'     of a cell from its neighbor, so that more distant neighbors have less weight.
#'   * if `'projected'`: An anisotropic filter that accounts for expected
#'     change in expression along the direction of the neighbor. The expected
#'     change in expression is calculated from the gradient field as the total
#'     derivative in the direction of the neighbor. This change in expression is
#'     then Gaussian transformed so that neighbors that are more distant along the
#'     direction of greatest change have less weight.
#'   * if `'constant'`: All neighbors have equal `distance` weights
#' * `similarity` score: Generally, neighbors with more similar gradients have
#'   greater weight
#'   * if `'euclidean'`: Gaussian transformation of the Euclidean distance
#'     between a cell's gradient field and its neighbor's gradient field.
#'   * if `'projected'`: Gaussian transformation of the cosine distance
#'     between a cell's gradient field and its neighbor's gradient field.
#'   * if `'constant'`: All neighbors have equal `similarity` weights
#' 
#' @param coords A `N` x `2` matrix of cell coordinates.
#' @param field A `2` x `D` x `N` array in column-major ordering
#'   containing the spatial gradient in expression for each of
#'   `D` latent variables at every point in space.
#' @param adj A `N` x `N` sparse adjacency matrix
#'   in dgCMatrix format.
#' @param include_self A boolean whether or not to include the
#'   each point's gradient in its own smoothed value. Defaults to TRUE.
#' @param distance Method for computing distance score in weighted average.
#'   See description for details. Defaults to `'euclidean'`.
#' @param similarity Method for computing similarity score in weighted average.
#'   See description for details. Defaults to `'euclidean'`.
#' 
#' @returns A `2` x `D` x `N` array in column-major ordering
#'   containing the smoothed spatial gradient in expression for each of
#'   `D` latent variables at every point in space.
#' 
#' @export
smooth_field = function(coords, field, adj, include_self=TRUE,
                        distance="euclidean", similarity="euclidean") {
    if (include_self) {
        diag(adj) = 1
    } else {
        diag(adj) = 0
    }

    ## method to compute distance between points 
    if (distance == "euclidean") {
        ## Euclidean distance
        dist_mode = 0
    } else if (distance == "projected") {
        ## Projected distance (anisotropic)
        dist_mode = 1
    } else if (distance == "constant") {
        dist_mode = 2
    } else {
        stop("Invalid distance method")
    }

    ## compute similarity between fields
    if (similarity == "euclidean") {
        ## Euclidean distance
        sim_mode = 0
    } else if (similarity == "projected") {
        sim_mode = 1
    } else if (similarity == "constant") {
        sim_mode = 2
    } else {
        stop("Invalid similarity method")
    }
    
    res = smooth_field_cpp(
        pvec = diff(adj@p), 
        adj_i = adj@i, ## keep it 0-indexed
        adj_p = adj@p,
        field = field,
        coords = coords,
        distance = dist_mode,
        similarity = sim_mode
    )
    return(res)
}


#' Compute a spatial gradient field at each point (cell)
#'
#' Distance between neighboring cells is normalized to unit distance
#' so that only the direction from each cell to its neighbors matters.
#' The gradient is then the average gradient in expression of each
#' embedding dimension between the index cell and its neighbors.
#'
#' @param coords A `N` x `2` matrix of cell coordinates.
#' @param embeddings A `N` x `D` matrix of cell embeddings.
#' @param adj A `N` x `N` sparse adjacency matrix
#'   in dgCMatrix format.
#'
#' @returns A `2` x `D` x `N` array in column-major ordering
#'   containing the spatial gradient in expression for each of
#'   `D` embedding dimensions at every point in space.
#'
#' @export
estimate_field = function(coords, adj, embeddings) {
    diag(adj) = 0 ## cannot include self by definition 
    res = estimate_field_cpp(
        coords = coords, 
        embeddings = embeddings, 
        adj_i = adj@i, 
        adj_p = adj@p
    )
    return(res)
}


