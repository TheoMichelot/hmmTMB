## A Matern SPDE field as an mgcv smooth. The construction is Lindgren, Rue
## and Lindstrom (2011): a Gaussian field with a Matern covariance is the
## solution of a stochastic PDE, and a finite element approximation on a
## triangulation turns that into a GMRF whose precision is sparse and known in
## closed form. 'fmesher' builds the mesh and the finite element matrices;
## nothing here needs INLA.
##
## Formally the field is just another smoother, which is why it can be written
## as an mgcv smooth and dropped into any hmmTMB formula. It differs from the
## smooths hmmTMB already handles in two ways: it carries three penalty
## matrices combined through two parameters rather than a single lambda * S,
## and its precision is proper and parameter-dependent, so its log-determinant
## is computed inside the likelihood rather than precomputed. Both are handled
## by make_matrices() and by smooth_penalty() in src/likelihood.hpp.

#' Matern SPDE smooth
#'
#' A smooth term for a spatially (or temporally) continuous Gaussian random
#' field, usable as \code{s(x, y, bs = "spde", xt = list(mesh = mesh))} in any
#' hmmTMB formula, and so on any transition probability or any parameter of any
#' observation distribution.
#'
#' @section The mesh in one dimension:
#' A one-dimensional field needs no mesh: \code{s(x, bs = "spde")} builds one
#' by default, as \code{k} quadratic B-spline basis functions on evenly spaced
#' knots, so
#'
#' \preformatted{
#' f <- ~ s(x, bs = "spde")          # and s(x, bs = "spde", k = 40) for a finer one
#' }
#'
#' is the whole of it. The knots are spread over the range of the covariate
#' widened by a fifth on each side, which is the one-dimensional version of the
#' outer extension of a two-dimensional mesh: with the knots stopping at the
#' data, the boundary condition inflates the field's variance at the two ends
#' by a third or more.
#'
#' \code{k} is a resolution, not a smoothing parameter, and this is the one
#' place where a field does not behave like the smooths around it. The range is
#' estimated and is what does the smoothing, so \code{k} is not a cap on
#' wiggliness that the penalty then pulls back from: it only has to be large
#' enough to resolve the range the data ask for, knots about a third of a range
#' apart as in two dimensions. That is why the default is \code{k = 15} rather
#' than \code{mgcv}'s 10 -- with the widening, 15 knots put the knot spacing at
#' about half the range the smoother starts from, which is a fifth of the
#' spread of the covariate.
#'
#' Too small a \code{k} does not simply oversmooth. Once the knots are further
#' apart than the range the data want, the weights are near enough independent
#' that the likelihood cannot tell one short range from another, and the range
#' collapses towards zero while the standard deviation grows to compensate. The
#' fitted curve can still look right when that happens, because the basis is
#' then doing the smoothing on its own, but the two parameters are meaningless.
#' So read \code{lambda()} before reading the curve: a range far below the knot
#' spacing, or a standard deviation orders of magnitude above the scale of the
#' linear predictor, means raise \code{k} and refit.
#'
#' Passing your own \code{fmesher::fm_mesh_1d()} overrides all of this, and is
#' worth doing when the covariate is unevenly spread and you want the knots
#' placed accordingly.
#'
#' @section Building a mesh in two dimensions:
#' In two dimensions the mesh is yours to build and pass in, because a good one
#' depends on the domain, the data and the range you expect.
#' \code{fmesher::fm_mesh_2d()} takes locations, a maximum triangle edge length
#' inside the domain and in the outer extension, and a cutoff below which
#' nearby points are merged:
#'
#' \preformatted{
#' mesh <- fmesher::fm_mesh_2d(loc = cbind(data$x, data$y),
#'                             max.edge = c(0.1, 0.3),
#'                             cutoff = 0.05, offset = c(0.1, 0.3))
#' f <- ~ s(x, y, bs = "spde", xt = list(mesh = mesh))
#' }
#'
#' A rule of thumb is \code{max.edge} no larger than a third of the range you
#' expect, and an offset about equal to the range; without the outer extension
#' the field's variance is inflated at the boundary.
#'
#' The field is evaluated at every time step, so missing covariates are
#' replaced by the last non-missing value as elsewhere in hmmTMB. If that is
#' too crude for locations, interpolate them before fitting.
#'
#' @section Parameterisation:
#' The precision is
#' \deqn{Q(\tau, \kappa) = \tau^2(\kappa^4 C + 2\kappa^2 G_1 + G_2)}
#' with \eqn{C}, \eqn{G_1}, \eqn{G_2} the finite element matrices from
#' \code{fmesher::fm_fem()}. \code{alpha = 2} throughout, so the Matern
#' smoothness is \eqn{\nu = 2 - d/2}: \eqn{\nu = 1} in two dimensions and
#' \eqn{\nu = 3/2} in one. The range is \eqn{\rho = \sqrt{8\nu}/\kappa} and
#' the marginal standard deviation is
#' \eqn{\sigma^2 = \Gamma(\nu) / (\Gamma(2)(4\pi)^{d/2}\kappa^{2\nu}\tau^2)},
#' which is \eqn{1/(\tau^2\kappa^2 4\pi)} in two dimensions and
#' \eqn{1/(4\tau^2\kappa^3)} in one. Both constants depend on \eqn{d}, so the
#' one- and two-dimensional smoothers do not share them.
#'
#' The two estimated parameters are \eqn{\log\sigma} and \eqn{\log\rho}, not
#' \eqn{\log\tau} and \eqn{\log\kappa}, and appear as \code{sd} and
#' \code{range} in \code{lambda()}. The change is exact -- a linear map of the
#' log parameters, absorbed into the constants in front of the three matrices
#' -- and worth making: the parameters mean something on their own, and the
#' likelihood's long ridge along constant marginal variance runs diagonally in
#' \eqn{(\tau, \kappa)} but along an axis in \eqn{(\sigma, \rho)}, which the
#' optimiser finds much easier.
#'
#' Both are estimated freely, with no prior. A single realisation identifies
#' them only loosely, and a field whose range approaches the size of the domain
#' is close to improper: the range then runs off and the standard deviation
#' follows it, with the fitted surface barely changing. The fit is still usable
#' if that happens, but the two parameters are not, and a finer mesh will not
#' help -- it is the data that are uninformative.
#'
#' The field is proper for any positive \eqn{\kappa}, so no identifiability
#' constraint is imposed, as in INLA. It is not centred, and its level is only
#' weakly separated from the intercept when the range is large.
#'
#' @section Bandwidth:
#' A model containing this smooth is fitted with the banded forward algorithm,
#' because the exact one makes the Hessian with respect to the field weights
#' dense and so defeats the sparsity the SPDE representation exists to provide.
#' See \code{HMM$update_bw()} and \code{HMM$check_bw()}.
#'
#' @param object,data,knots As \code{mgcv::smooth.construct()}; for
#'   \code{Predict.matrix}, \code{knots} is absent and the rest are as
#'   \code{mgcv::Predict.matrix()}.
#'
#' @return A \code{smoothCon} object of class \code{spde.smooth}, with sparse
#'   \code{X} and sparse penalties.
#'
#' @references
#' Lindgren, F., Rue, H. and Lindstrom, J. (2011). An explicit link between
#' Gaussian fields and Gaussian Markov random fields. \emph{JRSS-B} 73, 423-498.
#'
#' Fischer, J.-O. (2026). Fast and scalable inference in hidden Markov models
#' with Gaussian fields.
#'
#' @examples
#' # A one-dimensional field over a covariate, with a default mesh
#' if(requireNamespace("fmesher", quietly = TRUE)) {
#'   d <- data.frame(x = runif(100))
#'   sm <- mgcv::smoothCon(mgcv::s(x, bs = "spde"), data = d)[[1]]
#'   c(k = ncol(sm$X), sm$theta.names)
#' }
#'
#' @exportS3Method mgcv::smooth.construct
smooth.construct.spde.smooth.spec <- function(object, data, knots) {
  if(!requireNamespace("fmesher", quietly = TRUE)) {
    stop("bs = \"spde\" needs the fmesher package", call. = FALSE)
  }
  if(length(object$term) > 2 | length(object$term) < 1) {
    stop("An SPDE smooth is one- or two-dimensional, but this one has ",
         length(object$term), " terms.", call. = FALSE)
  }

  mesh <- object$xt$mesh
  if(is.null(mesh)) {
    if(length(object$term) == 2) {
      stop("A two-dimensional SPDE smooth needs a mesh: build one with ",
           "fmesher::fm_mesh_2d() and pass it as ",
           "s(x, y, bs = \"spde\", xt = list(mesh = mesh)). The right mesh ",
           "depends on the domain and on the range you expect, so there is ",
           "no useful default.", call. = FALSE)
    }
    mesh <- spde_mesh_1d(data[[object$term]], object$bs.dim)
  }
  if(!inherits(mesh, c("fm_mesh_1d", "fm_mesh_2d", "inla.mesh", "inla.mesh.1d"))) {
    stop("xt$mesh must be an fmesher mesh, from fm_mesh_2d() or fm_mesh_1d().",
         call. = FALSE)
  }
  # The Matern constants below depend on the dimension of the domain, so take
  # it from the mesh rather than from the number of terms, and make sure the
  # two agree before doing so.
  d <- ifelse(inherits(mesh, c("fm_mesh_1d", "inla.mesh.1d")), 1L, 2L)
  if(d != length(object$term)) {
    stop("A ", d, "-dimensional mesh was given for a smooth of ",
         length(object$term), " covariate", ifelse(length(object$term) > 1, "s", ""),
         ".", call. = FALSE)
  }

  object$X <- fmesher::fm_basis(mesh, spde_loc(object$term, data))
  fem <- fmesher::fm_fem(mesh)
  # The mass matrix is the lumped one, as in INLA: fm_fem already forms
  # g2 = g1 C0^-1 g1 with it, and using the consistent c1 in the kappa^4 term
  # would make the three matrices inconsistent with each other.
  # lambda = exp(L theta) with theta = (log sigma, log rho). Substituting
  # kappa = sqrt(8 nu) / rho and tau^2 = Gamma(nu) / (Gamma(2) (4 pi)^(d/2)
  # kappa^(2 nu) sigma^2) into tau^2 (kappa^4 C + 2 kappa^2 G1 + G2) leaves
  # these constants in front of the three matrices, and these powers of sigma
  # and rho. With alpha = 2 the smoothness is nu = 2 - d/2, so the constants
  # and two of the three powers of rho differ between one and two dimensions.
  if(d == 1L) {
    # nu = 3/2: kappa = 2 sqrt(3) / rho and tau^2 = 1 / (4 sigma^2 kappa^3)
    object$S <- list((sqrt(3)/2) * fem$c0, fem$g1 / (4*sqrt(3)),
                     fem$g2 / (96*sqrt(3)))
    object$L <- matrix(c(-2, -2, -2, -1, 1, 3), ncol = 2)
  } else {
    # nu = 1: kappa = 2 sqrt(2) / rho and tau^2 = 1 / (4 pi sigma^2 kappa^2)
    object$S <- list((2/pi) * fem$c0, fem$g1 / (2*pi), fem$g2 / (32*pi))
    object$L <- matrix(c(-2, -2, -2, -2, 0, 2), ncol = 2)
  }
  object$theta.names <- c("sd", "range")
  # A range of a fifth of the domain and unit marginal standard deviation:
  # smooth without being flat, on the scale of a linear predictor whose other
  # terms are of order one.
  object$theta.start <- c(0, log(spde_extent(mesh, object$X)/5))

  object$rank <- rep(ncol(object$X), length(object$S))
  object$null.space.dim <- 0    # proper for any positive kappa
  # The precision depends on both parameters through more than an overall
  # scale, so its log-determinant cannot be precomputed and has to be
  # differentiated inside the likelihood. make_matrices() reads this flag; it
  # is the only thing it needs to know about this smoother.
  object$gmrf <- TRUE
  # Two flags that keep smoothCon out of the way. 'no.rescale' suppresses its
  # penalty rescaling, which would divide the finite element matrices by a norm
  # of X and so break the meaning of tau and kappa -- and which is also the only
  # place it takes a matrix norm of X, so the design matrix can stay sparse. A
  # zero-row C says there is no centring constraint, which is right for a proper
  # field and is what INLA does.
  object$no.rescale <- TRUE
  object$C <- matrix(0, 0, ncol(object$X))
  object$df <- ncol(object$X)
  object$mesh <- mesh
  object$te.ok <- 0
  class(object) <- "spde.smooth"
  return(object)
}

#' @rdname smooth.construct.spde.smooth.spec
#' @exportS3Method mgcv::Predict.matrix
Predict.matrix.spde.smooth <- function(object, data) {
  # Dense, unlike the design matrix built above: this is only reached through
  # predict.gam(), which splices the result into a dense matrix spanning every
  # coefficient in the model and so would densify it anyway.
  as.matrix(fmesher::fm_basis(object$mesh, spde_loc(object$term, data)))
}

#' Observation locations for an SPDE smooth
#'
#' @param term Character vector of one or two covariate names
#' @param data Data frame containing those covariates
#'
#' @return Numeric vector (1d) or two-column matrix (2d) of locations
spde_loc <- function(term, data) {
  if(length(term) == 1) {
    return(as.numeric(data[[term]]))
  }
  return(cbind(as.numeric(data[[term[1]]]), as.numeric(data[[term[2]]])))
}

#' Default mesh for a one-dimensional SPDE smooth
#'
#' \code{k} quadratic B-spline basis functions on evenly spaced knots, over the
#' range of the covariate widened by a fifth on each side. The widening is the
#' one-dimensional counterpart of the outer extension of a two-dimensional
#' mesh: fmesher's default boundary condition is Neumann, under which a field
#' whose knots stop at the data has its variance inflated by a third or more at
#' the two ends. A fifth of the spread of the covariate is about one range at
#' the starting value used in \code{smooth.construct.spde.smooth.spec()}, which
#' is the same rule of thumb used for the offset in two dimensions.
#'
#' A degree-2 mesh has one basis function fewer than it has knots, so \code{k}
#' knots would not give \code{k} of them; \code{k} keeps its mgcv meaning of a
#' basis dimension here, and \code{k + 1} knots are laid down.
#'
#' The default of 15 is above mgcv's 10, because \code{k} is a resolution here
#' and not a cap on wiggliness: with the widening, 15 knots put the knot
#' spacing at about half the range the smoother starts from, in the spirit of
#' \code{max.edge} in two dimensions. Ten of them leave the knots further apart
#' than that range, and it collapses towards zero as a result.
#'
#' @param x The covariate
#' @param k Basis dimension, as passed to \code{s()}; negative when the user
#'   gave none, in which case 15 is used
#'
#' @return An \code{fm_mesh_1d} with \code{k} basis functions
spde_mesh_1d <- function(x, k) {
  if(is.null(k) || k < 0) k <- 15
  if(k < 3) {
    stop("An SPDE smooth needs k of at least 3, not ", k, ".", call. = FALSE)
  }
  r <- range(as.numeric(x), na.rm = TRUE)
  if(!all(is.finite(r)) || diff(r) <= 0) {
    stop("The covariate of an SPDE smooth has to vary and be finite, so that ",
         "a default mesh can be built over its range. Pass your own mesh as ",
         "s(x, bs = \"spde\", xt = list(mesh = mesh)) if that is not the case.",
         call. = FALSE)
  }
  pad <- diff(r) / 5
  return(fmesher::fm_mesh_1d(seq(r[1] - pad, r[2] + pad, length.out = k + 1),
                             degree = 2))
}

#' A length scale for the meshed domain, for starting values
#'
#' The larger side of the bounding box of the mesh nodes that carry basis
#' weight, which is the region the data occupy rather than the outer extension.
#' A one-dimensional mesh stores its knots in 'loc' and the midpoints of its
#' basis functions in 'mid', and only the latter are one per column of X, so
#' 'mid' is the one to use where there is a choice.
#'
#' @param mesh An fmesher mesh
#' @param X The smooth's design matrix
#'
#' @return A positive number
spde_extent <- function(mesh, X) {
  loc <- as.matrix(if(is.null(mesh$mid)) mesh$loc else mesh$mid)
  loc <- loc[, seq_len(min(2, ncol(loc))), drop = FALSE]
  used <- which(Matrix::colSums(abs(X)) > 0)
  if(length(used) > 1 & max(used) <= nrow(loc)) {
    loc <- loc[used, , drop = FALSE]
  }
  return(max(apply(loc, 2, function(z) diff(range(z))), .Machine$double.eps))
}
