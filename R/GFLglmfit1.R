#v0.1.0

#' @importFrom magrittr %>% add
#' @importFrom purrr map_dbl map_lgl


GFLglm.fit1 <- function(
    h, h1, mu, mu., b, lambda, w, z, y, a, q, q0, getCWS, get.sol,
    d
){

  r <- length(z)
  Lam <- 2*lambda*w

  SD <- cbind(-Lam, Lam) %>%
    add(
      map_dbl(z, ~{
        sum( a*h1(.x + q)*(mu(.x + q) - y) ) +
          2*lambda*sum(w*sign(.x - z))
      })
    )

  signSD <- sign(SD) %>% rowSums

  solution <- getCWS(z, signSD, y, a, q, q0, d, r, w, lambda, b, h, mu., get.sol)

  if(length(solution) > 0)
  {
    return(solution)
  } else
  {
    idx <- map_lgl(SD, ~identical(all.equal(.x, 0), TRUE)) %>% which

    if(length(idx) == 1)
    {
      SD[idx] <- 0

      signSD <- sign(SD) %>% rowSums
      solution <- getCWS(z, signSD, y, a, q, q0, d, r, w, lambda, b, h, mu., get.sol)
    } else
    {
      stop("error")
    }
  }
}
