#v0.1.0

#' @importFrom magrittr %>%
#' @importFrom purrr map_dbl

setCWM.nb <- function(link){
  # set coordinate wise minimization

  if(link == "canonical")
  {
    get.sol <- function(y, q, q0, phi, lambda, w, w., d, mu., b, h, z, signSD){

      u <- 2*lambda*w. - sum(y)

      if(d == 1)
      {
        out <- mu.(-u, phi) - q
      } else if(!is.na(q0))
      {
        out <- mu.(-u/d, phi) - q0
      } else
      {
        f <- function(x, .phi=phi){
          map_dbl(x, ~{
            sum(
              b(h(.x+q, .phi), .phi) - y*h(.x+q, .phi)
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
  } else
  {#---   NB2

    get.sol <- function(y, q, q0, phi, lambda, w, w., d, mu., b, h, z, signSD){

      v <- 2*lambda*w. - sum(y)

      if(d == 1)
      {
        out <- log(-v/phi) - log(2*lambda*w. + (1/phi)) - q
      } else if(!is.na(q0))
      {
        out <- log(-v/phi) - log(2*lambda*w. + (d/phi)) - q0
      } else
      {
        f <- function(x, .phi=phi){
          map_dbl(x, ~{
            sum(
              b(h(.x+q, .phi), .phi) - y*h(.x+q, .phi)
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

  return(list(get.sol = get.sol))
}
