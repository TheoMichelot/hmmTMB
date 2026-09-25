### Lions in the Central Kalahari
###
### A shortened version of the movement-ecology case study in Fischer (2026).
### Six lions in the Central Kalahari Game Reserve, Botswana, carried GPS
### collars that fixed their position roughly every hour; eight deployments
### between them span July 2009 to March 2011, and come to some fifty thousand
### hourly steps. A two-state HMM splits those steps into a resting state and
### an active one, and the question is what governs the switch between them.
###
### Two things plainly do. One is the time of day: lions are nocturnal, so the
### chance of being active runs on a 24-hour cycle. The other is where they
### are, and that is the harder one to model, because "where" is a pair of
### continuous coordinates rather than a handful of covariates. This is what a
### Gaussian field is for: the probability of settling down is allowed to vary
### smoothly over the study area, as a Matern field represented on a
### triangulated mesh through the SPDE approach of Lindgren, Rue and Lindstrom
### (2011). It enters the model as an ordinary smooth term,
### s(x, y, bs = "spde"), so nothing about the modelling interface changes;
### what changes is that hmmTMB switches to the banded forward algorithm, which
### keeps the Hessian with respect to the several hundred field weights sparse
### and so keeps the Laplace approximation affordable.
###
### The example in inst/examples/spde covers the same machinery on simulated
### data, where the true field is known. Start there if you want to see what
### the method recovers; come here for what it does to real data.
###
### Two models are fitted below, one with the field alone and one with the
### field and the daily cycle together. On a laptop the first fit takes about
### three and a half minutes and the second about ten; run start to finish,
### including the bandwidth check and the posterior simulation behind the
### confidence surfaces, the script takes some twenty-five minutes.
###
### Needs fmesher for the mesh and moveHMM for the step lengths and turning
### angles; both are in hmmTMB's Suggests.
###
### If you use such a field with hmmTMB, please cite Fischer (2026) as well as
### Lindgren, Rue and Lindstrom (2011).

library(hmmTMB)
library(fmesher)


# Data --------------------------------------------------------------------

## The data are Movebank study 3809257699:
##
##   https://www.movebank.org/cms/webapp?gwt_fragment=page=studies,path=study3809257699
##
## They are not available anonymously -- the public API returns nothing for
## this study, and the authenticated one refuses without credentials -- so this
## script cannot fetch them for you, and neither can it ship them. Create a
## (free) Movebank account, accept the study's licence terms, and download it
## as a CSV from the web interface. The lines below assume that file, under its
## Movebank name, in the working directory.
##
## move2::movebank_download_study(3809257699) will fetch it from R once
## move2::movebank_store_credentials() has your login, but it returns a move2
## object with Movebank's own column naming, so the preprocessing below would
## need adjusting.

raw <- read.csv("African lions in Central Kalahari Botswana.csv")

## Drop the fixes flagged as outliers, and keep the eight collar deployments
## used in Fischer (2026). They come from six lions -- two were collared a
## second time after the first collar came off -- and between them cover July
## 2009 to March 2011. The study holds twenty deployments in all, spread over
## 2008 to 2012; these eight are contemporaneous, which is what makes a single
## shared field worth estimating rather than eight separate ones. A different
## subset would be a different analysis.
##
## Each deployment is its own track below, even where two of them are the same
## animal: the two collars are a year apart, and nothing connects the last step
## of one to the first step of the next.
tags <- c("AL152", "AL153", "AL154", "AL155", "AL156", "AL158",
          "GSM08448", "GSM08449")
raw <- subset(raw, tag.local.identifier %in% tags &
                manually.marked.outlier != "true" &
                algorithm.marked.outlier != "true")

## The collars aimed at an hourly fix but did not hit it exactly, and
## occasionally fired twice in an hour or not at all. An HMM wants a regular
## time series, so keep the fix nearest each whole hour, discard anything more
## than fifteen minutes off it, and leave the hour empty when nothing
## qualifies. The empty hours matter: dropping them instead would close the
## gaps up silently and turn a missing night into one enormous step.
raw$time <- as.POSIXct(raw$timestamp, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
raw$hr <- round(as.numeric(raw$time) / 3600) * 3600
raw$off <- abs(as.numeric(raw$time) - raw$hr)
raw <- raw[raw$off <= 15 * 60, ]
raw <- raw[order(raw$tag.local.identifier, raw$hr, raw$off), ]
raw <- raw[!duplicated(raw[c("tag.local.identifier", "hr")]), ]

## One row per hour per tag, from that tag's first fix to its last
by_tag <- split(raw, raw$tag.local.identifier)
track <- do.call(rbind, lapply(names(by_tag), function(tag) {
  d <- by_tag[[tag]]
  hours <- data.frame(hr = seq(min(d$hr), max(d$hr), by = 3600))
  d <- merge(hours, d[c("hr", "location.long", "location.lat")], all.x = TRUE)
  data.frame(ID = tag,
             time = as.POSIXct(d$hr, origin = "1970-01-01", tz = "UTC"),
             x = d$location.long, y = d$location.lat)
}))

## Step lengths and turning angles, from great-circle distances between
## consecutive locations, so steps are in kilometres while the coordinates
## stay in degrees. Both are NA wherever a location is missing, and hmmTMB
## takes that as a missing observation and skips the likelihood contribution.
data <- moveHMM::prepData(track, type = "LL", coordNames = c("x", "y"))

## The covariates cannot be missing in the same way: the field has to be
## evaluated at every time step, and hmmTMB would otherwise carry the last
## location forward, which for a gap of several hours puts the animal
## somewhere it demonstrably was not. Linear interpolation within each track
## is a better guess and never extrapolates, since every track begins and ends
## on a real fix.
fill <- function(z) approx(seq_along(z), z, seq_along(z))$y
data$x_int <- ave(data$x, data$ID, FUN = fill)
data$y_int <- ave(data$y, data$ID, FUN = fill)

## Hour of the day, for the second model, in local time. Movebank timestamps
## are UTC and Botswana runs two hours ahead of it, with no daylight saving
## since 1943, so this is a clean two-hour shift. It is worth making: the whole
## point of the term is that lions are nocturnal, and a cycle plotted against
## UTC would put the peak of activity at 23:00 when it is really one in the
## morning. Fischer (2026) works in UTC; for that, drop the tz argument.
data$hour <- as.POSIXlt(data$time, tz = "Africa/Gaborone")$hour

nrow(data)                        # hourly steps
table(data$ID)                    # per deployment
mean(is.na(data$step))            # proportion of hours with no usable step

## Everything below is plotted in degrees, and a degree of longitude is not a
## degree of latitude: at 21.5 degrees south it covers about 7% less ground.
## asp = 1 would therefore stretch every map east-west. This is the correction,
## and it is the same one ggplot2's coord_quickmap() makes for itself.
map_asp <- 1 / cos(mean(data$y, na.rm = TRUE) * pi / 180)

## The eight tracks. They overlap heavily, which is what makes a single shared
## field worth estimating: several lions visiting the same places give it far
## more to go on than one would.
plot(data$x, data$y, asp = map_asp, type = "n", bty = "n",
     xlab = "Longitude", ylab = "Latitude")
for (id in unique(data$ID)) {
  keep <- data$ID == id
  lines(data$x[keep], data$y[keep], col = "#00798c30")
}

## Step lengths are strongly bimodal, which is the HMM's starting point: a
## great many near-zero steps, and a long tail of travelling ones.
hist(data$step, breaks = 200, xlim = c(0, 3), ylim = c(0, 2000),
     col = "grey80", border = "white", main = "", xlab = "Step length (km)")


# Mesh --------------------------------------------------------------------

## The mesh is the one part of an SPDE model with no sensible default: it has
## to resolve the features you expect to see. A rule of thumb is a maximum
## edge length of about a third of the range you expect, plus an outer
## extension of roughly one range, without which the field's variance is
## inflated at the boundary.
##
## Eight tracks do not fill a rectangle, so the mesh is built on non-convex
## hulls of the locations rather than on their bounding box. That keeps the
## number of field weights down, and keeps them where the data are.
##
## Everything spatial here is in degrees, including the range the model
## reports. The study area sits at about 21.5 degrees south, where a degree of
## longitude is some 7% shorter than a degree of latitude, so an isotropic
## field in degrees is very slightly anisotropic on the ground. At this
## latitude that is not worth projecting for; nearer the poles it would be.
loc <- cbind(data$x_int, data$y_int)
mesh <- fm_mesh_2d(loc = loc,
                   boundary = list(fm_nonconvex_hull(loc, convex = 0.04),
                                   fm_nonconvex_hull(loc, convex = 0.15)),
                   min.angle = 24,
                   max.edge = c(0.07, 1),
                   cutoff = 0.025)
mesh$n   # field weights, integrated out by the Laplace approximation

plot(mesh, asp = map_asp)
points(data$x, data$y, pch = 4, cex = 0.3, lwd = 0.5, col = "#00008b30")

## A finer mesh, max.edge = c(0.03, 1) and cutoff = 0.01, is what Fischer
## (2026) uses. It has some 2500 weights rather than 600, so it resolves more
## and costs more; this one is the compromise that keeps the script runnable.


# A field in the transition probabilities ---------------------------------

## State 1 is resting and state 2 is active, which is how the initial values
## below order them. The field goes on the 2 -> 1 transition, the probability
## of an active lion settling down, and the 1 -> 2 transition is
## intercept-only. Transition-specific formulas are given as a matrix, with
## "." on the reference (diagonal) entries.
## See ?smooth.construct.spde.smooth.spec for details on the SPDE smoother.
form <- matrix(c(".", "~ 1",
                 "~ s(x_int, y_int, bs = 'spde', xt = list(mesh = mesh))", "."),
               nrow = 2, byrow = TRUE)

hid <- MarkovChain$new(data = data, n_states = 2, formula = form,
                       initial_state = "stationary")

## Step lengths are zero-inflated gamma rather than gamma: a couple of hundred
## steps are exactly zero, because the collar returned the same coordinates
## twice, and a gamma density cannot accommodate those. Turning angles are
## wrapped Cauchy, centred near pi when resting -- consecutive short steps are
## mostly GPS error, which reverses direction -- and near 0 when active, where
## travel is directed.
obs <- Observation$new(
  data = data, n_states = 2,
  dists = list(step = "zigamma2", angle = "wrpcauchy"),
  par = list(step = list(mean = c(0.007, 0.9),
                         sd = c(0.007, 1),
                         z = c(0.002, 0.005)),
             angle = list(mu = c(3.141, 0), rho = c(0.2, 0.3))))

hmm_sp <- HMM$new(hid = hid, obs = obs)

## The field switched the banded forward algorithm on, with the default
## bandwidth of 15
hmm_sp$bw()

## Fitting the model
hmm_sp$fit()


# Check the bandwidth -----------------------------------------------------

## The bandwidth controls how far back in time each log-likelihood
## contribution is allowed to look. Too small and the approximation is poor;
## too large and computations become infeasible. Profile it at the
## estimates: if the curve has flattened before the bandwidth in use, the
## approximation is fine. The row with bw = Inf is the exact algorithm, which
## is what the banded one is being compared against.
##
## Each bandwidth means a fresh TMB object, so keep the grid short: this is
## already a minute or two, and a grid twice as fine would say the same thing.
## (The value profiled is the joint log-likelihood, at fixed parameters, so it
## is not the marginal one that llk() reports; only the differences matter.)
hmm_sp$check_bw(bws = seq(5, 25, by = 5))

## Here it goes -24668.9, -24654.6, -24653.7, and then stops: by bw = 15 the
## banded algorithm agrees with the exact one to two decimal places, so the
## default needs no changing. Had it not stabilised, we would raise it and
## refit:
# hmm_sp$update_bw(25)
# hmm_sp$fit()


# Results -----------------------------------------------------------------

## The two state-dependent distributions. The resting state has a mean step of
## about seven metres, which is GPS noise around a stationary animal rather
## than movement; the active state averages roughly nine hundred metres an
## hour. The turning-angle means come back at pi and 0 -- reversal when resting,
## persistence when active.
hmm_sp$obs()$par()[, , 1]

## The field's two parameters, a marginal standard deviation on the scale of
## the linear predictor and a range in degrees. The range comes out near 0.23
## degrees, some 23 km: comfortably inside a study area and about three times the
## mesh's inner edge length, so the rule of thumb the mesh was built on holds.
## A range that came out close to the size of the study area would be the warning
## sign: the field would then be nearly improper and its parameters meaningless,
## however good the surface looked.
hmm_sp$lambda()$hid

## The fitted surface, as Pr(active -> resting) over the two coordinates.
## HMM$plot() varies one covariate and holds the rest at their means, which for
## a bivariate term shows a single slice; plot_2d() varies both.
##
## It ranges from about 0.13 to 0.65 -- a sizeable effect. The low end is a single
## coherent region in the south-centre of the study area, where an active lion
## settles within the hour with probability around 0.15 against 0.5 or more
## over most of the rest.
##
## plot_2d() sets no coordinate system, which is right for a method that has to
## serve any pair of covariates -- a surface over temperature and wind speed has
## no business being square. When the two are coordinates it does matter, and
## since the return value is a ggplot the fix is one layer. coord_quickmap()
## is the one to reach for rather than coord_equal(): it applies the latitude
## correction above, where coord_equal() would just make the panel square.
##
## One practical note. The smooth constructor keeps the mesh on the smooth
## object it builds, but hmmTMB does not keep that object: it stores the
## assembled design and penalty matrices, and rebuilds the smooth from the
## formula whenever it needs matrices at new covariate values. The formula
## holds `mesh` by name, so every prediction below looks it up afresh. Keep the
## mesh around for as long as you want to predict -- saving the fitted model on
## its own and reloading it in a new session is not enough, and the error when
## that happens is an unhelpful "object 'mesh' not found".
hmm_sp$plot_2d("tpm", var = "x_int", var2 = "y_int", i = 2, j = 1, n_grid = 100) +
  coord_quickmap()

## Cells further than 10% of the plot's diagonal from the nearest observation
## are left blank by default, as in mgcv::plot.gam(): a fitted surface far from
## the data says more about the basis than about the data. Here that blanks the
## corners of the bounding box, which no lion visited.
##
## Where the surface is actually determined is a separate question, and filling
## by the width of the confidence interval rather than by the estimate answers
## it. The interval is under 0.1 wide over the middle, where the tracks are
## dense, and half the probability scale around the rim, where a lion passed
## once or twice -- so the region driving the interpretation above is the
## well-determined part, and the bright edges are not worth reading. Contours
## are off because the interval width is a Monte Carlo quantity, and contouring
## it mostly draws the simulation noise.
hmm_sp$plot_2d("tpm", var = "x_int", var2 = "y_int", i = 2, j = 1,
               n_grid = 30, n_post = 300, show = "ci", contour = FALSE) +
  coord_quickmap()

## Every plot_2d() call returns a ggplot, so it takes further layers like any
## other -- here the tracks themselves, over the surface they produced. (The
## lattice columns are called var and var2 inside that plot, so the layer needs
## its own aes; ggplot2 is attached already, as hmmTMB depends on it.)
hmm_sp$plot_2d("tpm", var = "x_int", var2 = "y_int", i = 2, j = 1,
               contour = FALSE) +
  geom_path(data = data, aes(x = x_int, y = y_int, group = ID),
            colour = "white", linewidth = 0.1, alpha = 0.3) +
  coord_quickmap()

## For anything plot_2d() will not do, go through predict() directly: it is
## what the method calls, on a lattice of your own making, and it returns the
## array.
grid <- expand.grid(x_int = seq(min(data$x_int), max(data$x_int), length = 200),
                    y_int = seq(min(data$y_int), max(data$y_int), length = 200))
tpm <- hmm_sp$predict("tpm", newdata = grid)
image(unique(grid$x_int), unique(grid$y_int), matrix(tpm[2, 1, ], 200, 200),
      col = hcl.colors(30), asp = map_asp, zlim = c(0, 1), bty = "n",
      xlab = "Longitude", ylab = "Latitude",
      main = expression(Pr(active %->% resting)))


# Adding the daily cycle --------------------------------------------------

## Lions are nocturnal, so both transitions should depend on the time of day.
## A cyclic cubic spline is the natural term for that, and it goes on both
## off-diagonal entries alongside the field.
##
## The knots have to be given explicitly. bs = "cc" wraps the spline between
## the outermost knots, which mgcv puts at the range of the covariate: with
## hour running 0 to 23 that would tie hour 23 to hour 0 and give the cycle a
## period of 23 hours. Setting them at 0 and 24 ties hour 24 to hour 0
## instead, which is what a day does.
form2 <- matrix(
  c(".",
    "~ s(hour, bs = 'cc')",
    "~ s(hour, bs = 'cc') + s(x_int, y_int, bs = 'spde', xt = list(mesh = mesh))",
    "."),
  nrow = 2, byrow = TRUE)

hid2 <- MarkovChain$new(data = data, n_states = 2, formula = form2,
                        initial_state = "stationary",
                        gam_args = list(knots = list(hour = c(0, 24))))

## The observation model is unchanged, but it cannot simply be handed to the
## second HMM: MarkovChain and Observation are R6 objects, held by reference
## rather than copied, so the two models would share one set of observation
## parameters and fitting the second would overwrite the first. Clone it. The
## clone carries the fitted values with it, which makes this a warm start.
hmm_td <- HMM$new(hid = hid2, obs = obs$clone())

hmm_td$fit()
## This model is pretty memory hungry which your machine might not handle.

## Both criteria prefer the larger model.
AIC(hmm_sp, hmm_td)
BIC(hmm_sp, hmm_td)

## The daily cycle, as the stationary probability of being active against the
## hour. The other covariates -- the coordinates -- are held at their means,
## so this is the cycle at the centre of the study area. It runs from about
## 0.75 through the night down to 0.10 in the early afternoon.
##
## That night-time 0.75 is higher than the 0.5 the decoded states give below,
## and the two are not in conflict: this is the distribution the chain would
## settle into if the night-time transition probabilities applied forever,
## whereas the chain actually spends the night on its way there from a day
## spent resting. A cycle shorter than the chain's mixing time never reaches
## its own stationary distribution.
hmm_td$plot("delta", var = "hour", i = 2)

## Both transition probabilities behind it. Pr(resting -> active) is the more
## interesting of the two: rather than one broad nocturnal hump it has two
## peaks, around 05:00 and again around 19:00, with a shallow dip between them
## in the middle of the night. Lions here get going at dawn and at dusk.
hmm_td$plot("tpm", var = "hour", i = 1, j = 2)
hmm_td$plot("tpm", var = "hour", i = 2, j = 1)

## The field, net of the daily cycle. The two terms are additive on the linear
## predictor, so changing the hour only shifts the surface up or down there --
## but the plot is on the probability scale, where the logistic transform then
## squashes a shifted surface differently. That makes the hour a real choice:
## the default, hour at its mean and so near midday, puts the whole surface where
## a lion is very likely to be resting anyway and flattens it against 1.
## Night, when the switch is actually in play, is the slice worth looking at.
hmm_td$plot_2d("tpm", var = "x_int", var2 = "y_int", i = 2, j = 1,
               covs = list(hour = 20)) +
  coord_quickmap()

## The field's parameters, against those of the smaller model. They hardly
## move -- sd 0.565 against 0.597, range 0.239 against 0.230 -- so the field was
## not standing in for the daily cycle in the first model: the two effects are
## picking up different things.
hmm_td$lambda()$hid
hmm_sp$lambda()$hid

## And the state-dependent distributions, which barely move either: the two
## states are pinned down by the step lengths, not by what drives the switching.
hmm_td$obs()$par()[, , 1]
