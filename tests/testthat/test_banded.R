context("Banded forward algorithm")

# An independent R implementation of both algorithms, written straight from
# the definition, to check the C++ against. The banded one mirrors
# forward_alg_banded() in src/likelihood.hpp index for index.
fwd_R <- function(prob, tpm_all, delta0, ID, bw) {
  ns <- ncol(prob)
  rho <- rep(1 / ns, ns)
  llk <- 0
  ids <- unique(ID)

  for(j in seq_along(ids)) {
    idx <- which(ID == ids[j])
    st <- min(idx)
    en <- max(idx)
    k <- if(bw < 2) en - st + 1 else bw

    # First block, from the model's own initial distribution
    phi <- delta0[j, ] * prob[st, ]
    llk <- llk + log(sum(phi))
    phi <- phi / sum(phi)
    if(st + 1 <= min(st + k - 1, en)) {
      for(r in (st + 1):min(st + k - 1, en)) {
        phi <- as.vector(phi %*% tpm_all[, , r - 1]) * prob[r, ]
        llk <- llk + log(sum(phi))
        phi <- phi / sum(phi)
      }
    }

    # Later blocks, each warmed up from rho over the preceding block
    b <- st + k
    while(b <= en) {
      phi <- rho * prob[b - k, ]
      phi <- phi / sum(phi)
      if(b - k + 1 <= b - 1) {
        for(r in (b - k + 1):(b - 1)) {
          phi <- as.vector(phi %*% tpm_all[, , r - 1]) * prob[r, ]
          phi <- phi / sum(phi)
        }
      }
      for(r in b:min(b + k - 1, en)) {
        phi <- as.vector(phi %*% tpm_all[, , r - 1]) * prob[r, ]
        llk <- llk + log(sum(phi))
        phi <- phi / sum(phi)
      }
      b <- b + k
    }
  }
  llk
}

# Overlapping states and a persistent chain, so that the observations do not
# identify the state on their own and the truncation actually bites
make_model <- function(n = 300, n_id = 1, seed = 7) {
  set.seed(seed)
  z <- numeric(n)
  s <- 1
  for(i in 1:n) {
    z[i] <- stats::rnorm(1, c(0, 1)[s], 1)
    s <- sample(c(s, 3 - s), size = 1, prob = c(0.95, 0.05))
  }
  data <- data.frame(z = z)
  if(n_id > 1) data$ID <- rep(seq_len(n_id), each = n / n_id)
  obs <- Observation$new(data = data, dists = list(z = "norm"), n_states = 2,
                         par = list(z = list(mean = c(0, 1), sd = c(1, 1))))
  hid <- MarkovChain$new(n_states = 2, data = data)
  HMM$new(obs = obs, hid = hid)
}

# Joint log-likelihood at the current parameters, at a given bandwidth
llk_at <- function(hmm, bw) {
  if(is.null(hmm$.__enclos_env__$private$tmb_args_)) hmm$setup(silent = TRUE)
  obj <- hmm$.__enclos_env__$private$make_tmb_obj(
    bw = bw, include_smooths = -1, silent = TRUE)
  -obj$fn(obj$par)
}

ingredients <- function(hmm) {
  list(prob = hmm$obs()$obs_probs(),
       tpm_all = hmm$hid()$tpm(t = "all"),
       delta0 = hmm$hid()$delta0(as_matrix = TRUE),
       ID = hmm$obs()$data()$ID)
}

test_that("the exact forward algorithm is unchanged", {
  hmm <- make_model()
  g <- ingredients(hmm)
  expect_equal(llk_at(hmm, 0), fwd_R(g$prob, g$tpm_all, g$delta0, g$ID, 0))
})

test_that("the banded forward algorithm matches its R reference", {
  hmm <- make_model()
  g <- ingredients(hmm)
  for(bw in c(2, 3, 5, 8, 13, 20)) {
    expect_equal(llk_at(hmm, bw),
                 fwd_R(g$prob, g$tpm_all, g$delta0, g$ID, bw),
                 info = paste("bw =", bw))
  }
})

test_that("the banded log-likelihood converges to the exact one", {
  hmm <- make_model()
  exact <- llk_at(hmm, 0)
  err <- sapply(c(2, 5, 10, 20, 40), function(b) abs(llk_at(hmm, b) - exact))

  # Geometric decay: nowhere near exact at bw = 2, and indistinguishable by 40
  expect_gt(err[1], 1e-2)
  expect_lt(err[5], 1e-6)
  expect_lt(err[5], err[1])
})

test_that("banding respects time series boundaries", {
  hmm <- make_model(n = 300, n_id = 3)
  g <- ingredients(hmm)
  expect_equal(length(unique(g$ID)), 3L)
  for(bw in c(3, 7, 15)) {
    expect_equal(llk_at(hmm, bw),
                 fwd_R(g$prob, g$tpm_all, g$delta0, g$ID, bw),
                 info = paste("bw =", bw))
  }
  # A bandwidth longer than each series is the exact algorithm on each
  expect_equal(llk_at(hmm, 200), llk_at(hmm, 0))
})

# The bw API --------------------------------------------------------------

test_that("update_bw validates its argument and invalidates the setup", {
  hmm <- make_model(n = 100)

  expect_equal(hmm$bw(), 0L)          # exact by default
  expect_error(hmm$update_bw(1), "should be 0")
  expect_error(hmm$update_bw(-2), "should be 0")
  expect_error(hmm$update_bw(2.5), "should be 0")
  expect_error(hmm$update_bw(c(2, 3)), "should be 0")
  expect_error(hmm$update_bw("5"), "should be 0")

  hmm$setup(silent = TRUE)
  expect_false(is.null(hmm$.__enclos_env__$private$tmb_obj_))
  hmm$update_bw(10)
  expect_equal(hmm$bw(), 10L)
  # bw is TMB data, so the taped object has to go
  expect_null(hmm$.__enclos_env__$private$tmb_obj_)

  hmm$update_bw(NULL)
  expect_equal(hmm$bw(), 0L)
})

test_that("check_bw profiles the log-likelihood", {
  hmm <- make_model(n = 200)
  prof <- hmm$check_bw(bws = c(5, 10, 20))

  expect_equal(names(prof), c("bw", "llk"))
  expect_equal(prof$bw, c(5, 10, 20, Inf))
  expect_equal(prof$llk[4], llk_at(hmm, 0))
  expect_true(all(is.finite(prof$llk)))
  expect_error(hmm$check_bw(bws = c(1, 5)), "at least 2")
})
