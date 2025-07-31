#' Create maximal antichains from interval data
#'
#' @param intvls A matrix with 2 columns: left and right endpoints of intervals.
#' @param Lopen Logical. Whether the left endpoint is open.
#' @param Ropen Logical. Whether the right endpoint is open.
#'
#' @keywords internal
#' @return A list of maximal antichains (each is a vector of row indices).
#' @export
Maclist <- function(intvls, Lopen = TRUE, Ropen = FALSE) {
  m <- nrow(intvls)
  id <- seq_len(m)
  or <- order(intvls[, 1])
  maclist <- list()
  curmac <- id[or[1]]
  minend <- intvls[curmac, 2]

  for (i in 2:m) {
    curintvl <- id[or[i]]
    if (intvls[curintvl, 1] > minend ||
        ((Lopen || Ropen) && intvls[curintvl, 1] == minend)) {
      maclist <- c(maclist, list(curmac))
      oldmac <- curmac
      curmac <- integer(0)
      for (j in oldmac) {
        if (intvls[curintvl, 1] < intvls[j, 2] ||
            (!Lopen && !Ropen && intvls[curintvl, 1] == intvls[j, 2])) {
          curmac <- c(curmac, j)
        }
      }
      curmac <- c(curmac, curintvl)
      minend <- min(intvls[curmac, 2])
    } else {
      curmac <- c(curmac, curintvl)
      minend <- min(minend, intvls[curintvl, 2])
    }
  }
  c(maclist, list(curmac))
}

#' Create clique matrix and Petrie pairs from maximal antichains
#'
#' @param ml A list of maximal antichains as returned by [Maclist()].
#'
#' @keywords internal
#' @return A list of class `petrie` containing:
#' \describe{
#'   \item{pmat}{Clique matrix}
#'   \item{ppairs}{Start and end of each interval in antichains}
#' }
#' @export
Macmat <- function(ml) {
  m <- length(ml)
  temp <- sort(unique(unlist(ml)))
  n <- length(temp)
  ppairs <- matrix(0, 2, n)
  retmat <- matrix(0, m, n)
  for (i in seq_len(m)) {
    for (j in ml[[i]]) {
      if (ppairs[1, j] == 0) ppairs[1, j] <- i
      ppairs[2, j] <- i
      retmat[i, j] <- 1
    }
  }
  dimnames(ppairs) <- list(c("Start", "End"), temp)
  dimnames(retmat) <- list(NULL, temp)
  structure(list(pmat = retmat, ppairs = ppairs), class = "petrie")
}

#' Map maximal antichains to real-valued intervals
#'
#' @param intvls Matrix of intervals (n x 2).
#' @param ml Maximal antichain list from [Maclist()]. If missing, it will be computed.
#'
#' @keywords internal
#' @return A matrix with mapped interval representations for each antichain.
#' @export
MLEintvl <- function(intvls, ml = Maclist(intvls)) {
  if (ncol(intvls) != 2 || any(intvls[, 2] < intvls[, 1])) {
    stop("only one-dimensional intervals can be handled")
  }
  m <- length(ml)
  ret <- matrix(0, m, 2)
  for (i in seq_len(m)) {
    LL <- min(intvls)
    RR <- max(intvls)
    for (j in ml[[i]]) {
      LL <- max(LL, intvls[j, 1])
      RR <- min(RR, intvls[j, 2])
    }
    ret[i, ] <- c(LL, RR)
  }
  ret
}

#' Rescale probability vector
#'
#' @param pvec A probability vector.
#' @param tiny Threshold below which values are set to zero.
#'
#' @return A rescaled probability vector.
#' @export
rescaleP <- function(pvec, tiny) {
  pvec <- ifelse(pvec < tiny, 0, pvec)
  pvec / sum(pvec)
}

#' EM Algorithm for Interval-Censored Data
#'
#' @param A Either a logical matrix or an interval matrix (n x 2).
#' @param pvec Optional initial probability vector.
#' @param maxiter Maximum number of EM iterations (default 500).
#' @param tol Convergence tolerance (default 1e-12).
#'
#' @return An object of class `"icsurv"` with elements:
#' \describe{
#'   \item{pf}{Estimated probability mass function}
#'   \item{numiter}{Number of iterations}
#'   \item{converge}{Logical, whether EM converged}
#'   \item{intmap}{Interval map (if applicable)}
#' }
#' @export
EM <- function(A, pvec, maxiter = 500, tol = 1e-12) {
  if (ncol(A) == 2 && all(A[, 2] >= A[, 1])) {
    ml <- Maclist(A)
    intmap <- t(MLEintvl(A, ml))
    A <- Macmat(ml)$pmat
  } else {
    intmap <- NULL
  }

  i <- 0
  notdone <- TRUE
  n <- ncol(A)
  Meps <- .Machine$double.eps

  if (missing(pvec)) {
    pvec <- rowSums(A) / sum(A)
  }
  pvec <- rescaleP(pvec, Meps)

  while (i < maxiter && notdone) {
    i <- i + 1
    dmat <- diag(pvec)
    t1 <- dmat %*% A
    t2 <- 1 / (t(A) %*% pvec)
    np <- rescaleP(as.vector(t1 %*% t2) / n, Meps)
    if (sum(abs(np - pvec)) < tol) notdone <- FALSE
    pvec <- np
  }

  structure(
    list(pf = pvec, numiter = i, converge = !notdone, intmap = intmap),
    class = "icsurv"
  )
}


#' @export
print.icsurv <- function(x, ...) {
  cat("Estimated probabilities:\n")
  print(x$pf)
  cat("Iterations:", x$numiter, "\n")
  cat("Converged:", x$converge, "\n")
}

#' @export
print.petrie <- function(x, ...) {
  cat("Clique matrix:\n")
  print(x$pmat)
  cat("Petrie pairs:\n")
  print(x$ppairs)
}

