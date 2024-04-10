#v0.1.0

setGLM <- function(family, link, y, as){

  ##############################################################################
  ###   GLM component
  ##############################################################################

  if(family %in% c("Gaussian", "Binomial", "Poisson") | is.null(link))
  {
    link <- "canonical"

    h <- function(eta){eta}
    h1 <- function(eta){rep(1, length(eta))}
  }

  if(family == "Poisson" | link == "log")
  {
    mu <- function(eta){exp(eta)}
    mu. <- function(zeta){log(zeta)}
  }

  if(family %in% c("Binomial", "Poisson"))
  {
    dispersion <- FALSE
    a <- function(phi){1}
  } else
  {
    dispersion <- TRUE
  }

  if(family == "Gaussian")
  {
    mu <- function(eta){eta}
    mu. <- function(zeta){zeta}
    b <- function(theta){(1/2)*(theta^2)}
    V <- function(x){1}
    a <- function(phi){phi}
    .c <- function(phi, .y=y){-((.y^2)/(2*phi)) - (1/2)*log(2*pi*phi)}
  }

  if(family == "Binomial")
  {
    mu <- function(eta){exp(eta)/(1 + exp(eta))}
    mu. <- function(zeta){log(zeta/(1-zeta))}
    b <- function(theta){-log(1 - mu(theta))}
    V <- function(x){x*(1 - x)}
    .c <- function(phi, .a=as, .y=y){log(choose(.a, .a*.y))}

    fNLL <- function(Eta, .a=as, .y=y){
      # negative log-likelihood without constant (beta に依存しない項)
      - sum(
        .a*( .y*h(Eta) - b(h(Eta)) )
      )
    }
  }

  if(family == "Poisson")
  {
    b <- mu
    V <- function(x){x}
    .c <- function(phi, .y=y){-lgamma(.y+1)}
  }

  if(family == "Gamma")
  {
    if(link == "canonical")
    {
      mu <- function(eta){-1/eta}
      mu. <- function(zeta){-1/zeta}
    } else if(link == "log")
    {
      h <- function(eta){-exp(-eta)}
      h1 <- function(eta){1/exp(eta)}
    }

    b <- function(theta){-log(-theta)}
    V <- function(x){x^2}
    a <- function(phi){phi}
    .c <- function(phi, .y=y){
      psi <- 1/phi
      return(
        -psi*log(phi) + (psi - 1)*log(.y) - lgamma(psi)
      )
    }
  }

  if(family %in% c("Inverse.Gaussian", "IG"))
  {
    family <- "Inverse.Gaussian"

    if(link == "canonical")
    {
      mu <- function(eta){1/sqrt(-eta)}
      mu. <- function(zeta){-1/(zeta^2)}
    } else if(link == "log")
    {
      h <- function(eta){-exp(-2*eta)}
      h1 <- function(eta){2/exp(2*eta)}
    }

    b <- function(theta){-2*sqrt(-theta)}
    V <- function(x){x^3}
    a <- function(phi){2*phi}
    .c <- function(phi, .y=y){
      -(1/2)*log(2*pi*phi*(.y^3)) - (1/(2*phi*.y))
    }
  }

  if(family != "Binomial")
  {
    fNLL <- function(Eta, .y=y){
      # negative log-likelihood without constant
      - sum(
        .y*h(Eta) - b(h(Eta))
      )
    }
  }

  return(list(
    b = b, mu = mu, mu. = mu., h = h, h1 = h1, V = V,
    family = family, link = link, dispersion = dispersion,
    fNLL = fNLL
  ))
}
