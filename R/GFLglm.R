#' @title Generalized fused Lasso for generalized linear models
#' @description \code{GFLglm} This function provides generalized fused Lasso estimator for
#'     generalized linear models via coordinate descent algorithm (v0.1.0)
#'
#' @importFrom magrittr %>% add set_colnames subtract
#' @importFrom purrr map map_dbl
#' @param y a vector of a response variable
#' @param group a vector of groups
#' @param D a list of adjacent relationship among groups
#' @param family the error distribution; one of "Gaussian", "Binomial", "Poisson",
#'     "Negative.Binomial" (or "NB"), "Gamma", and "Inverse.Gaussian" (or "IG")
#' @param link the link function; the default is canonical link;
#'     this can take "log" for "Gamma", "NB", and "IG"
#' @param alpha a positive value adjusting the strength of the penalty of
#'     generalized information criterion, e.g., 2 corresponds to AIC;
#'     the default is `log(length(y))` which corresponds to BIC
#' @param offset a vector of offset
#' @param adaptive if `TRUE`, penalty weights are used
#' @param as a vector of the numbers of trials for only `family=Binomial`
#' @param nlam the number of tuning parameter candidates
#' @param maxit an upper bound for iteration
#' @param lam.cand if "uniform" is given, tuning parameter candidates are defined
#'     by uniformly dividing
#' @param dsp.nb for only "NB", this decides the estimation method of dispersion parameter;
#'     1: by Pearson estimator; 2: using Pearson estimator (default; see, Table 8.4 of
#'     Hilbe, 2011);
#'     3: Newton's method; 4: only once update by update equation of Newton's method
#' @param tol a tolerance for convergence
#' @param progress if `TRUE`, progress is displayed
#' @return a list having following elements:
#' \item{fitted.values}{fitted values}
#'
#' \item{beta}{estimates of group parameters}
#'
#' \item{eta}{linear predictors}
#'
#' \item{phi}{dispersion parameter}
#'
#' \item{cluster}{cluster of groups}
#'
#' \item{offset}{offset}
#'
#' \item{weight}{penalty weights}
#'
#' \item{results}{results for each tuning parameter candidate which have the
#'     following columns: `lambda` (tuning parameter), `MSC` (model selection criterion value),
#'     `NLL` (nagative log-likelihood function value), `DF` (degrees of freedom; the number of clusters),
#'     `dispersion` (dispersion parameter), `chisq` (Pearson's chi squared statistic),
#'     `convergence` (algorithm converged (1) or not (0))}
#'
#' \item{all.beta}{estimates of group parameters for each tuning parameter candidate}
#'
#' \item{summary}{results under the optimal tuning parameter, which has
#'     `idx` (index of the optimal tuning parameter), `lambda` (the optimal tuning parameter),
#'     `dispersion` (dispersion parameter), `df` (the number of clusters),
#'     `family`, `link`, `convergence`, `time` (total runtime)}
#' @references Hilbe, J. M. (2011). *Negative Binomial Regression 2nd Edition*.
#'   Cambridge University Press. Cambridge.
#' @references Ohishi, M. (2024).
#'   Generalized fused Lasso for grouped data in generalized linear models.
#'   *Hiroshima Statistical Research Group Technical Report*, TR-No. 24-02, Hiroshima University.
#' @export
#' @examples
#' #GFLglm(y, group, D, family)

GFLglm <- function(
    y, group, D, family, link=NULL, alpha=NULL, offset=FALSE, adaptive=TRUE,
    as=NULL, nlam=100, maxit=500, lam.cand=NULL, dsp.nb=2, tol=1e-5, progress=FALSE
){

  t1 <- proc.time()[3]

  n <- length(y)
  ns <- table(group) %>% as.numeric
  m <- length(ns)
  r <- map_dbl(D, length)

  if(family == "Binomial" & max(y) > 1)
  {
    y0 <- y
    y <- y/as
  }

  if(isFALSE(offset))
  {
    q <- numeric(n)
  } else
  {
    q <- offset
    offset <- TRUE
  }

  if(is.null(as)){as <- rep(1, n)}
  if(is.null(alpha)){alpha <- log(n)}
  isNB <- family %in% c("Negative.Binomial", "NB")

  yli <- split(y, group)
  qli <- split(q, group)
  ali <- split(as, group)

  q0s <- qli %>% map_dbl(~{
    if(length(unique(.x)) == 1)
    {
      return(.x[1])
    } else
    {
      return(NA)
    }
  })

  ##############################################################################
  ###   GLM components
  ##############################################################################

  if(isNB)
  {
    GLMcomp <- setNB(link, q0s, dsp.nb, y)
    CWMf <- setCWM.nb(GLMcomp$link)
  } else
  {
    GLMcomp <- setGLM(family, link, y, as)
    CWMf <- setCWM(GLMcomp$family, GLMcomp$link)
  }

  for(i in 1:length(GLMcomp))
  {
    assign(names(GLMcomp)[i], GLMcomp[[i]])
  }

  for(i in 1:length(CWMf))
  {
    assign(names(CWMf)[i], CWMf[[i]])
  }

  ##############################################################################
  ###   initial values (MLE)
  ##############################################################################

  if(m == n)
  {
    if((family %in% c("Poisson", "Negative.Binomial")) & (0 %in% y))
    {
      Zeta0 <- y + 0.01
    } else
    {
      Zeta0 <- y
    }

    if(isNB)
    {
      phi0 <- 1
      Beta0 <- mu.(Zeta0, phi0) - q
    } else
    {
      Beta0 <- mu.(Zeta0) - q
    }
  } else
  {
    if(isNB)
    {
      if(0 %in% y)
      {
        Zeta0 <- split(y + 0.01, group)
      } else
      {
        Zeta0 <- yli
      }

      MLE <- getMLE(
        y, q, q0s, 1, b, mu., h, mu, V, n, tol,
        m, Zeta0, qli, ns
      )

      Beta0 <- MLE$beta
      phi0 <- MLE$phi
    } else
    {
      if(family == "Gaussian")
      {
        Beta0 <- map_dbl(1:m, ~mean(yli[[.x]] - qli[[.x]]))
      } else if(family == "Poisson" | (family == "Gamma" & link == "log"))
      {
        if(0 %in% y)
        {
          Zeta0 <- split(y + 0.01, group)
        } else
        {
          Zeta0 <- yli
        }
        Beta0 <- map_dbl(1:m, ~{
          yj <- Zeta0[[.x]]; qj <- qli[[.x]]; hq <- h1(qj)
          return(sum(yj*hq)/sum(hq*exp(qj)))
        }) %>% log
      } else if(all(!is.na(q0s)))
      {
        Beta0 <- map_dbl(1:m, ~{
          aj <- ali[[.x]]; yj <- yli[[.x]]
          return(sum(aj*yj)/sum(aj))
        }) %>% mu.(.) %>% subtract(q0s)
      } else
      {
        Beta0 <- map_dbl(1:m, ~{
          if(ns[.x] == 1)
          {
            return(mu.(yli[[.x]]) - qli[[.x]])
          } else if(!is.na(q0s[.x]))
          {
            aj <- ali[[.x]]; yj <- yli[[.x]]
            return(
              mu.(sum(aj*yj)/sum(aj)) - q0s[.x]
            )
          } else
          {
            aj <- ali[[.x]]; qj <- qli[[.x]]; yj <- yli[[.x]]
            obje <- function(x){
              map_dbl(x, ~sum(
                aj * (b(h(.x+qj)) - yj*h(.x+qj))
              ))
            }

            return(optimise(obje, c(-1000, 1000))$minimum)
          }
        })
      }
    }
  }

  ##############################################################################
  ###   MLE for lam.max
  ##############################################################################

  if(isNB)
  {
    ME <- getME(y, q, Beta0, phi0, b, mu., h, mu, V, n, tol, m, ns)

    betaM <- ME$beta
    dispM <- ME$phi
  } else
  {
    if(family == "Gaussian")
    {
      betaM <- mean(y - q)
    } else if(family == "Poisson" | (family == "Gamma" & link == "log"))
    {
      betaM <- log(
        sum(y*h1(q)) / sum(h1(q)*exp(q))
      )
    } else if(all(!is.na(q0s)))
    {
      betaM <- mu.(
        sum(as*y) / sum(as)
      ) - q0s[1]
    } else
    {
      betaM <- optimise(
        function(x){
          map_dbl(x, ~sum(
            as * (b(h(.x+q)) - y*h(.x+q))
          ))
        },
        range(Beta0)
      )$minimum
    }

    ZetaM <- mu(betaM + q)
    dispM <- sum(((y - ZetaM)^2) / V(ZetaM))/(n - 1)
  }

  ##############################################################################
  ###   penalty weights
  ##############################################################################

  if(adaptive)
  {
    W <- map(1:m, ~{1/abs(Beta0[.x] - Beta0[D[[.x]]])})

    if(Inf %in% unlist(W))
    {
      idxD <- map(1:m, ~{
        wi <- W[[.x]]
        idx <- which(wi == Inf)
        if(length(idx) > 0){return(c(.x, D[[.x]][idx]))}
      }) %>% unlist %>% unique %>% setdiff(1:m, .)
    } else
    {
      idxD <- 1:m
    }
  } else
  {
    W <- map(1:m, ~rep(1, length(D[[.x]])))
    idxD <- 1:m
  }

  ##############################################################################
  ###   tuning parameter candidates
  ##############################################################################

  if(isNB)
  {
    lam.max <- max(
      abs(
        map_dbl(1:m, ~sum(
          h1(betaM+qli[[.x]], dispM)*(mu(betaM+qli[[.x]], dispM) - yli[[.x]])
        ))
      ) / (2*map_dbl(W, sum))
    )
  } else
  {
    lam.max <- max(
      abs(
        map_dbl(1:m, ~sum(
          ali[[.x]]*h1(betaM+qli[[.x]])*(mu(betaM+qli[[.x]]) - yli[[.x]])
        ))
      ) / (2*map_dbl(W, sum))
    )
  }

  if(is.null(lam.cand))
  {
    Lambda <- lam.max*( (3/4)^((nlam-1):0) )
  } else if(lam.cand == "uniform")
  {
    Lambda <- seq(0, lam.max, length=nlam+1)[-1]
  }

  lam.n <- nlam

  ##############################################################################
  ###   estimation
  ##############################################################################

  MSC <- NLL <- DF <- Chisq <- Disp <- Conv <- numeric(lam.n)
  BETA <- matrix(0, m, lam.n)

  .Beta <- Beta0
  lam.idx <- 1:lam.n

  if(isNB)
  {
    for(lam.i in lam.idx)
    {
      lambda <- Lambda[lam.i]

      res <- GFLnb.fit(
        y, D, .Beta, 1, betaM, dispM, r, maxit, tol,
        h, h1, mu, mu., b, V, lambda, W, idxD, get.sol, get.dsp, get.dsp0, PSN,
        yli, qli, q0s, n, m, ns, q
      )

      BETA[,lam.i] <- .Beta <- res$beta
      Chisq[lam.i] <- res$chisq
      Disp[lam.i] <- .phi <- res$disp

      NLL[lam.i] <- nll <- fNLL(res$eta, .phi)
      DF[lam.i] <- df <- res$df
      MSC[lam.i] <- 2*nll + alpha*df

      Conv[lam.i] <- conv <- res$convergence

      if(progress){print(paste0(lam.i, "; ", "convergence: ", conv))}
    } #end for lam.i

    dispersion <- TRUE
  } else
  {
    for(lam.i in lam.idx)
    {
      lambda <- Lambda[lam.i]

      res <- GFLglm.fit(
        y, D, .Beta, betaM, dispM, r, maxit, tol,
        h, h1, mu, mu., b, V, lambda, W, idxD, getCWS,
        yli, ali, qli, q0s, ns, get.sol, m, q, n
      )

      BETA[,lam.i] <- .Beta <- res$beta
      Chisq[lam.i] <- res$chisq
      Disp[lam.i] <- res$disp

      NLL[lam.i] <- nll <- fNLL(res$eta)
      DF[lam.i] <- df <- res$df
      MSC[lam.i] <- 2*nll + alpha*df

      Conv[lam.i] <- conv <- res$convergence

      if(progress){print(paste0(lam.i, "; ", "convergence: ", conv))}
    } #end for lam.i
  }

  ##############################################################################
  ###   tuning parameter selection
  ##############################################################################

  opt <- which.min(MSC)
  Beta.hat <- BETA[,opt]

  if(dispersion)
  {
    phi.hat <- Disp[opt]
  } else
  {
    phi.hat <- 1
  }

  t2 <- proc.time()[3]

  ##############################################################################
  ###   output
  ##############################################################################

  Eta.hat <- map(1:m, ~rep(Beta.hat[.x], ns[.x])) %>% unlist %>% add(q)

  if(isNB)
  {
    Fit <- mu(Eta.hat, phi.hat)
  } else
  {
    Fit <- mu(Eta.hat)
  }

  if(progress){message(paste0("complete", " ", Sys.time()))}

  return(list(
    fitted.values = Fit,
    beta = Beta.hat,
    eta = Eta.hat,
    phi = phi.hat,
    cluster = match(Beta.hat, unique(Beta.hat)) %>% split(1:m, .),
    offset = q,
    weight = W,
    results = data.frame(
      lambda = Lambda,
      MSC = MSC,
      NLL = NLL,
      DF = DF,
      dispersion = Disp,
      chisq = Chisq,
      convergence = Conv
    ),
    all.beta = BETA %>% set_colnames(paste0("lam", 1:lam.n)),
    summary = data.frame(
      idx = opt,
      lambda = Lambda[opt],
      dispersion = Disp[opt],
      df = DF[opt],
      family = family,
      link = link,
      convergence = Conv[opt],
      time = set_names(t2 - t1, NULL)
    )
  ))
}
