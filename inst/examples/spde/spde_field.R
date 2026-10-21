### A latent Gaussian field in the transition probabilities of an HMM
###
### Simulated version of the movement-ecology example in Fischer (2026): an
### animal switches between two behaviours, and the probability of switching
### out of the resting behaviour varies smoothly over space. The spatial effect
### is modelled as a Matern field, represented on a triangulated mesh through
### the SPDE approach of Lindgren, Rue and Lindstrom (2011), which turns it
### into a Gaussian Markov random field with a sparse precision matrix.
###
### The data stream is the step length, and the track is what the steps
### produce: the animal takes a state-dependent step in a direction that
### persists from one step to the next, so its position, the data we model, and
### the covariate the field lives on are all the same process. That is the
### situation the method is for, and it is what makes the example worth
### running -- a field over positions unconnected to the modelled data would
### behave quite differently.
###
### The field enters the model as an ordinary mgcv smooth,
### s(x, y, bs = "spde"), so nothing about the modelling interface changes.
### What does change is how the likelihood is computed: hmmTMB switches to the
### banded forward algorithm, because the exact one would make the Hessian
### with respect to the several hundred field weights dense and the Laplace
### approximation unaffordable.
###
### Runs in well under a minute: the mesh follows the track rather than a
### bounding box, which keeps the number of field weights down.
###
### If you use such a field with hmmTMB, please cite Fischer (2026) as well as
### Lindgren, Rue and Lindstrom (2011).

library(hmmTMB)
library(fmesher)

set.seed(3847)


# Simulate ----------------------------------------------------------------

## The true field, added to the linear predictor of Pr(resting -> travelling).
## Its period is 40 in each direction, so the fitted range should come out at
## a few tens of units.
true_field <- function(x, y) 2 * (sin(2 * pi * x / 40) + cos(2 * pi * y / 40))

## Two states: 1 is "resting" (short steps), 2 is "travelling" (long ones).
## Step lengths are gamma, parameterised by mean and standard deviation, which
## is hmmTMB's "gamma2".
mu    <- c(0.2, 5)             # mean step length in each state
sdev  <- c(0.5, 3)             # standard deviation of step length
beta0 <- qlogis(c(0.2, 0.2))   # baseline switching probabilities

n <- 5000
kappa_pull <- 0.3              # how strongly the animal is drawn home

## Samplers for the two distributions used, written out so the simulation
## depends on nothing but base R
rgamma2 <- function(n, mean, sd) {
  scale <- sd^2 / mean
  rgamma(n, shape = mean / scale, scale = scale)
}
## von Mises by rejection, which is efficient at this concentration
rvm1 <- function(mu, kappa) {
  repeat {
    angle <- runif(1, -pi, pi)
    if (runif(1) < exp(kappa * (cos(angle - mu) - 1))) return(angle)
  }
}

## The walk. At each step the animal draws a step length from the distribution
## of its current state, and a turning angle pulled towards the origin so that
## the track stays in a bounded region rather than diffusing away. The state it
## switches to next is governed by the field at the location it has just
## reached.
state <- numeric(n)
state[1] <- 1
loc <- matrix(0, n, 2)
step <- numeric(n)
heading <- 0

for (t in 1:(n - 1)) {
  step[t] <- rgamma2(1, mu[state[t]], sdev[state[t]])
  ## Turning angle is state-independent, but biased towards c(0, 0)
  pull <- atan2(-loc[t, 2], -loc[t, 1]) - heading
  heading <- heading + rvm1(pull, kappa_pull)
  loc[t + 1, ] <- loc[t, ] + step[t] * c(cos(heading), sin(heading))

  ## Transition probabilities at the new location
  p12 <- plogis(beta0[1] + true_field(loc[t + 1, 1], loc[t + 1, 2]))
  p21 <- plogis(beta0[2])
  Gamma <- matrix(c(1 - p12, p12, p21, 1 - p21), nrow = 2, byrow = TRUE)
  state[t + 1] <- sample(1:2, size = 1, prob = Gamma[state[t], ])
}
step[n] <- rgamma2(1, mu[state[n]], sdev[state[n]])

data <- data.frame(ID = 1, step = step, x = loc[, 1], y = loc[, 2])

## The two state-dependent distributions, and the track they produced
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
curve(dgamma(x, shape = mu[1]^2 / sdev[1]^2, scale = sdev[1]^2 / mu[1]),
      xlim = c(0, 12), n = 500, bty = "n", lwd = 2, col = "#00798c",
      xlab = "Step length", ylab = "Density", main = "State-dependent steps")
curve(dgamma(x, shape = mu[2]^2 / sdev[2]^2, scale = sdev[2]^2 / mu[2]),
      add = TRUE, n = 500, lwd = 2, col = "#d1495b")
plot(data$x, data$y, type = "l", asp = 1, col = "#00798c80", bty = "n",
     xlab = "x", ylab = "y", main = "Simulated track")
points(0, 0, pch = 16)
par(mfrow = c(1, 1))


# Mesh --------------------------------------------------------------------

## The mesh is the one part of an SPDE model with no sensible default: it has
## to resolve the features you expect to see. A rule of thumb is a maximum
## edge length of about a third of the range you expect, plus an outer
## extension of about one range, without which the field's variance is
## inflated at the boundary.
##
## A real track does not fill a rectangle, so the mesh is built on a non-convex
## hull of the locations rather than on their bounding box. That keeps the
## number of field weights down, and keeps them where the data are.
locs <- cbind(data$x, data$y)
mesh <- fm_mesh_2d(loc = locs,
                   boundary = list(fm_nonconvex_hull(locs, convex = -0.03),
                                   fm_nonconvex_hull(locs, convex = -0.2)),
                   max.edge = c(7, 50),
                   cutoff = 2)
plot(mesh, asp = 1)
lines(data$x, data$y, col = "#00798c60")
mesh$n   # number of field weights, integrated out by the Laplace approximation


# Fit ---------------------------------------------------------------------

## The field goes on the 1 -> 2 transition only; the 2 -> 1 transition is
## intercept-only. Transition-specific formulas are given as a matrix, with
## "." on the reference (diagonal) entries.
form <- matrix(c(".", "~ s(x, y, bs = 'spde', xt = list(mesh = mesh))",
                 "~ 1", "."),
               nrow = 2, byrow = TRUE)

hid <- MarkovChain$new(data = data, n_states = 2, formula = form,
                       initial_state = "stationary")

obs <- Observation$new(data = data, n_states = 2,
                       dists = list(step = "gamma2"),
                       par = list(step = list(mean = c(0.2, 5),
                                              sd = c(0.5, 3))))

hmm <- HMM$new(hid = hid, obs = obs)

## The field switched the banded forward algorithm on, with the default
## bandwidth of 15
hmm$bw()

system.time(hmm$fit())


# Check the bandwidth -----------------------------------------------------

## The bandwidth controls how far back in time each log-likelihood
## contribution is allowed to look. Too small and the approximation bites;
## large enough and the log-likelihood stops moving. Profile it at the
## estimates: if the curve has flattened well before the bandwidth in use,
## the approximation is fine.
hmm$check_bw(bws = seq(5, 30, by = 5))

## If the likelihood had not stabilised, we would raise it and refit:
# hmm$update_bw(25)
# hmm$fit()


# Results -----------------------------------------------------------------

## The field's two parameters are a marginal standard deviation and a range,
## on the scale of the linear predictor and of the coordinates respectively.
## The range should be of the order of the period of the true field, 40.
## They are not variances, so sd_re() reports NA for them.
hmm$lambda()$hid

## The state-dependent step length distributions, against the mu and sdev
## used to simulate
hmm$obs()$par()[, , 1]

## And the baseline switching probabilities, against beta0
hmm$coeff_fe()$hid

## The fitted field, as a surface over the two coordinates. HMM$plot() varies
## one covariate and holds the rest at their means, which for a bivariate term
## shows a single slice through the surface; plot_2d() varies both.
hmm$plot_2d("tpm", var = "x", var2 = "y", i = 1, j = 2)

## Cells further than 10% of the plot's diagonal from the nearest observation
## are left blank by default, as in mgcv::plot.gam(): a fitted surface far from
## the data says more about the basis than about the data. The track is roughly
## circular, so this removes the corners of the bounding box. Compare with
## too_far = 0, which shows how confident the extrapolation looks.
hmm$plot_2d("tpm", var = "x", var2 = "y", i = 1, j = 2, too_far = 0)

## Where is the surface actually determined? Filling by the width of the
## confidence interval rather than by the estimate shows that: narrow near the
## middle, where the pull towards the origin made the animal spend most of its
## time, and wide out at the edges it reached only occasionally. Contours are
## off because the interval width is a Monte Carlo quantity, and contouring it
## mostly draws the simulation noise.
hmm$plot_2d("tpm", var = "x", var2 = "y", i = 1, j = 2,
            n_grid = 30, n_post = 300, show = "ci", contour = FALSE)

## Every plot_2d() call returns a ggplot, so it takes further layers and
## scales like any other. The next section does the same thing from scratch,
## for when that is not enough.


# Doing it by hand ---------------------------------------------------------

## plot_2d() is HMM$predict() on a lattice, plus a ggplot. Going through
## predict() yourself gives you the numbers, and so full control of the
## picture: your own grid, your own palette, and layers the method knows
## nothing about -- here the track itself, drawn over the surface it produced.
grid <- expand.grid(x = seq(min(data$x), max(data$x), length = 150),
                    y = seq(min(data$y), max(data$y), length = 150))
tpm <- hmm$predict("tpm", newdata = grid)
grid$p <- tpm[1, 2, ]

library(ggplot2)
ggplot(grid, aes(x, y)) +
  geom_raster(aes(fill = p)) +
  # The first 600 steps only: all 5000 cover the range and hide the surface
  geom_path(data = data[1:600, ], colour = "white", linewidth = 0.3) +
  scale_fill_viridis_c("Pr(resting -> travelling)", limits = c(0, 1)) +
  coord_equal() + theme_light() +
  ggtitle("Fitted field, with the track that produced it")

## Or without ggplot2 at all. predict() returns an array, so reshaping it to a
## matrix is all base graphics needs.
image(unique(grid$x), unique(grid$y), matrix(grid$p, 150, 150),
      col = hcl.colors(30), asp = 1, zlim = c(0, 1),
      xlab = "x", ylab = "y", main = "Estimated", bty = "n")


# Against the truth --------------------------------------------------------

## Comparing with a known truth needs the truth as well as the estimate, so
## this is a case for the manual route: both surfaces on one grid, faceted, on
## a single shared colour scale. Only the part of the grid the animal actually
## visited is worth comparing, so drop the rest.
seen <- !mgcv::exclude.too.far(grid$x, grid$y, data$x, data$y, dist = 0.1)
truth <- plogis(beta0[1] + true_field(grid$x, grid$y))

both <- rbind(
  data.frame(grid[c("x", "y")], p = ifelse(seen, truth, NA), what = "True"),
  data.frame(grid[c("x", "y")], p = ifelse(seen, grid$p, NA), what = "Estimated"))
both$what <- factor(both$what, levels = c("True", "Estimated"))

ggplot(both, aes(x, y, fill = p)) +
  geom_raster(na.rm = TRUE) +
  facet_wrap("what") +
  scale_fill_viridis_c("Pr(resting -> travelling)", limits = c(0, 1),
                       na.value = "transparent") +
  coord_equal() + theme_light()

## Correlation between the fitted and the true surface, on the linear
## predictor scale, over the region the animal visited
cor(qlogis(grid$p[seen]), qlogis(truth[seen]))

## The lattice of peaks is recovered, but the ones around the edge are
## attenuated: they sit where the animal passed only occasionally, and a
## penalised fit pulls the field towards its mean where the data are thin. The
## central peak, where the pull towards the origin kept it most of the time,
## reaches very nearly the true height -- which is why the overall range hardly
## changes even though the periphery is visibly flattened.
range(truth[seen])
range(grid$p[seen])


## Decoded states
table(hmm$viterbi(), state)
