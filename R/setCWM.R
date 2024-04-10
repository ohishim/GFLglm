#v0.1.0

#' @importFrom magrittr %>%
#' @importFrom purrr map map_dbl

setCWM <- function(family, link){
  # set coordinate wise minimization

  if(family == "Inverse.Gaussian" & link == "log")
  {
    get.sol <- NULL

    getCWS <- function(z, signSD, y, a, q, q0, d, r, w, lambda, b, h, mu., get.sol){
      # get coordinate wise solution

      Cand <- z[which(signSD %in% c(-1, 0, 1))]

      Rl <- c(-Inf, z); Ru <- c(z, Inf)

      u1 <- sum(y*exp(-2*q))
      u2 <- sum(exp(-q))

      Cand <- map(0:r, ~{
        if(.x == 0)
        {
          w. <- -sum(w)
        } else if(.x == r)
        {
          w. <- sum(w)
        } else
        {
          w. <- sum(w[1:.x]) - sum(w[(.x+1):r])
        }

        v <- 2*lambda*w.

        D <- (u2^2) + 2*u1*v

        if(D >= 0)
        {
          out <- log((-u2 + sqrt(D))/v)

          if(Rl[.x+1] < out & out < Ru[.x+1])
          {
            return(out)
          } else
          {
            return(NULL)
          }
        } else
        {
          return(NULL)
        }
      }) %>% unlist %>% c(Cand)

      if(length(Cand) == 1)
      {
        solution <- Cand
      } else
      {
        f <- function(x){
          map_dbl(x, ~{
            sum(
              b(h(.x+q)) - y*h(.x+q)
            ) + 2*lambda*sum(w*abs(.x - z))
          })
        }

        solution <- which.min(f(Cand)) %>% Cand[.]
      }

      return(solution)
    } # end gwtCWS
  } else
  {

    if(link == "canonical")
    {
      if(family == "Gaussian")
      {
        get.sol <- function(y, a, q, q0, lambda, w, w., d, mu., b, h, z, signSD){
          (sum(y - q) - 2*lambda*w.) / d
        }
      } else if(family == "Poisson")
      {
        get.sol <- function(y, a, q, q0, lambda, w, w., d, mu., b, h, z, signSD){
          u <- 2*lambda*w. - sum(y)
          return(log(-u) - log(sum(exp(q))))
        }
      } else
      {
        get.sol <- function(y, a, q, q0, lambda, w, w., d, mu., b, h, z, signSD){
          u <- 2*lambda*w. - sum(a*y)

          if(d == 1)
          {
            out <- mu.(-u/a) - q
          } else if(!is.na(q0))
          {
            out <- mu.(-u/sum(a)) - q0
          } else
          {
            f <- function(x){
              map_dbl(x, ~{
                sum(
                  a*( b(h(.x+q)) - y*h(.x+q) )
                ) + 2*lambda*sum(w*abs(.x - z))
              })
            }

            if(all(signSD == -2))
            {
              z1 <- max(z)
              R <- ifelse(z1 > 1, abs(z1), 1)
              z2 <- z1 + R

              l <- 1
              while(f((z2 + z1)/2) > f(z2))
              {
                l <- l + 1
                z2 <- z1 + l*R
              }

              Int <- c(z1, z2)
            } else if(all(signSD == 2))
            {
              z2 <- z[1]
              R <- ifelse(z2 > 1, abs(z2), 1)
              z1 <- z2 - R

              l <- 1
              while(f((z2 + z1)/2) > f(z1))
              {
                l <- l + 1
                z1 <- z2 - l*R
              }

              Int <- c(z1, z2)
            } else
            {
              j. <- which(signSD == -2) %>% max
              Int <- z[c(j., j.+1)]
            }

            out <- optimise(f, Int)$minimum
          }

          return(out)
        }
      }
    } else if(family == "Gamma")
    {#---   Gamma with log-link

      get.sol <- function(y, a, q, q0, lambda, w, w., d, mu., b, h, z, signSD){

        u <- sum(y/exp(q))
        v <- 2*lambda*w. + d
        return(log(u/v))
      }
    }

    getCWS <- function(z, signSD, y, a, q, q0, d, r, w, lambda, b, h, mu., get.sol){
      # get coordinate wise solution

      if(sum(signSD %in% c(-1, 0, 1)) == 1)
      {
        solution <- z[which(signSD %in% c(-1, 0, 1))]
      } else
      {
        if(all(signSD == -2))
        {
          w. <- sum(w)
        } else if(all(signSD == 2))
        {
          w. <- -sum(w)
        } else
        {
          w. <- sum(sign(-signSD)*w)
        }

        solution <- get.sol(y, a, q, q0, lambda, w, w., d, mu., b, h, z, signSD)
      }

      return(solution)
    }
  } # end ifelse

  return(list(getCWS = getCWS, get.sol = get.sol))
}
