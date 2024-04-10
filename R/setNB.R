#v0.1.0

#' @importFrom magrittr %>% add
#' @importFrom purrr map map_dbl

setNB <- function(link, q, dsp, y){

  ##############################################################################
  ###   NB component
  ##############################################################################

  family <- "Negative.Binomial"
  if(is.null(dsp)){dsp <- 1}

  if(is.null(link))
  {
    link <- "canonical"

    mu <- function(eta, phi){(1/phi)*exp(eta)/(1 - exp(eta))}
    mu. <- function(zeta, phi){log(zeta) - log((1/phi) + zeta)}

    h <- function(eta, phi){eta}
    h1 <- function(eta, phi){rep(1, length(eta))}

  } else if(link == "log")
  {
    mu <- function(eta, phi){exp(eta)}
    mu. <- function(zeta, phi){log(zeta)}

    h <- function(eta, phi){eta - log((1/phi) + exp(eta))}
    h1 <- function(eta, phi){1/(1 + phi*exp(eta))}
  }

  b <- function(theta, phi){
    (1/phi)*log(
      (1/phi) + (1/phi)*exp(theta)/(1 - exp(theta))
    )
  }

  V <- function(x, phi){x + phi*(x^2)}

  fNLL <- function(Eta, phi, .y=y){
    # negative log-likelihood
    psi <- 1/phi

    return(
      - sum(
        .y*h(Eta, phi) - b(h(Eta, phi), phi) +
          lgamma(.y + psi) - lgamma(.y + 1) - lgamma(psi) +
          psi*log(psi)
      )
    )
  }

  ##############################################################################
  ###   MLE
  ##############################################################################

  getMLE <- function(
    y, q, q0s, phi0, b, mu., h, mu, V, n, tol,
    m, yli, qli, ns
  ){

    Beta1 <- map_dbl(1:m, ~{
      yj <- yli[[.x]]; qj <- qli[[.x]]
      return(sum(yj)/sum(exp(qj)))
    }) %>% log

    phi <- phi0

    run <- TRUE
    while(run)
    {
      Beta0 <- Beta1
      betaM <- abs(Beta0) %>% max
      Int <- c(-2, 2)*betaM

      phi <- get.dsp(Beta0, phi, mu, q, y, V, n, length(unique(Beta0)), tol, m, ns)

      Beta1 <- map_dbl(1:m, ~{
        if(ns[.x] == 1)
        {
          return(mu.(yli[[.x]], phi) - qli[[.x]])
        } else if(!is.na(q0s[.x]))
        {
          yj <- yli[[.x]]
          return(
            mu.(mean(yj), phi) - q0s[.x]
          )
        } else
        {
          qj <- qli[[.x]]; yj <- yli[[.x]]
          obje <- function(x){
            map_dbl(x, ~sum(
              (b(h(.x+qj, phi), phi) - yj*h(.x+qj, phi))
            ))
          }

          return(optimise(obje, Int)$minimum)
        }
      })

      run <- max((Beta1 - Beta0)^2 / Beta0^2, na.rm=T) > tol
    } # end while

    return(list(beta=Beta1, phi=phi))
  }

  if(dsp == 1)
  {
    get.dsp <- function(Beta, phi, mu, q, y, V, n, .t, tol, m, ns){

      Eta <- map(1:m, ~rep(Beta[.x], ns[.x])) %>% unlist %>% add(q)
      Zeta <- mu(Eta, phi)

      return(
        sum(((y - Zeta)^2) / V(Zeta, phi))/(n - .t)
      )
    }

    get.dsp0 <- function(Beta, phi, mu, q, y, V, n, .t, tol, m, ns){phi}
  } else if(dsp == 2)
  {
    get.dsp <- function(Beta, phi, mu, q, y, V, n, .t, tol, m, ns){

      Eta <- map(1:m, ~rep(Beta[.x], ns[.x])) %>% unlist %>% add(q)
      Zeta <- mu(Eta, phi)
      chisq <- sum(((y - Zeta)^2) / V(Zeta, phi))

      return(phi*chisq/(n - .t))
    }

    get.dsp0 <- function(Beta, phi, mu, q, y, V, n, .t, tol, m, ns){phi}
  } else if(dsp == 3)
  {
    get.dsp <- get.dsp0 <- function(Beta, phi, mu, q, y, V, n, .t, tol, m, ns){

      Eta <- map(1:m, ~rep(Beta[.x], ns[.x])) %>% unlist %>% add(q)
      Zeta <- mu(Eta, phi)

      g1 <- function(psi){
        sum(
          log((psi+Zeta)/psi) + ((y-Zeta)/(psi+Zeta)) -
            digamma(y + psi) + digamma(psi)
        )
      }

      g2 <- function(psi){
        sum(
          (-(Zeta^2)-psi*y)/(psi*((psi+Zeta)^2)) -
            trigamma(y+psi) + trigamma(psi)
        )
      }

      psi1 <- n / sum(((y-Zeta)^2)/(Zeta^2))
      run <- TRUE

      while(run)
      {
        psi0 <- psi1
        psi1 <- abs(psi0 - (g1(psi0)/g2(psi0)))

        run <- abs(psi1 - psi0) > tol
      }

      return(1/psi1)
    }
  } else if(dsp == 4)
  {
    get.dsp <- get.dsp0 <- function(Beta, phi, mu, q, y, V, n, .t, tol, m, ns){

      Eta <- map(1:m, ~rep(Beta[.x], ns[.x])) %>% unlist %>% add(q)
      Zeta <- mu(Eta, phi)

      g1 <- function(psi){
        sum(
          log((psi+Zeta)/psi) + ((y-Zeta)/(psi+Zeta)) -
            digamma(y + psi) + digamma(psi)
        )
      }

      g2 <- function(psi){
        sum(
          (-(Zeta^2)-psi*y)/(psi*((psi+Zeta)^2)) -
            trigamma(y+psi) + trigamma(psi)
        )
      }

      psi0 <- n / sum(((y-Zeta)^2)/(Zeta^2))
      psi1 <- abs(psi0 - (g1(psi0)/g2(psi0)))

      return(1/psi1)
    }
  }

  ##############################################################################
  ###   MLE for lam.max
  ##############################################################################

  if(length(unique(q)) == 1)
  {
    get.betaM <- function(y, q, Beta0, phi0, b, mu., h){
      mu.(mean(y), phi0) - q[1]
    }
  } else
  {
    get.betaM <- function(y, q, Beta0, phi0, b, mu., h){

      obje <- function(x, phi){
        map_dbl(x, ~sum(
          b(h(.x+q, phi), phi) - y*h(.x+q, phi)
        ))
      }

      return(optimise(obje, range(Beta0), phi=phi0)$minimum)
    }
  }

  getME <- function(y, q, Beta0, phi0, b, mu., h, mu, V, n, tol, m, ns){

    beta1 <- mean(Beta0)
    phi <- phi0

    run <- TRUE
    while(run)
    {
      beta0 <- beta1

      beta1 <- get.betaM(y, q, Beta0, phi, b, mu., h)
      phi <- get.dsp(rep(beta1, m), phi, mu, q, y, V, n, 1, tol, m, ns)

      run <- abs(beta1 - beta0) > tol
    } # end while

    return(list(beta=beta1, phi=phi))
  }

  return(list(
    b = b, mu = mu, mu. = mu., h = h, h1 = h1, V = V,
    family = family, link = link,
    get.betaM = get.betaM, getME = getME, get.dsp = get.dsp, get.dsp0 = get.dsp0,
    PSN = dsp %in% c(1, 2), fNLL = fNLL, getMLE = getMLE
  ))
}
