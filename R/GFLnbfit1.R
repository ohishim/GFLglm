#v0.1.0

#' @importFrom magrittr %>% add
#' @importFrom purrr map_dbl

GFLnb.fit1 <- function(
    h, h1, mu, mu., b, lambda, w, z, y, q, q0, phi, get.sol,
    d
){

  Lam <- 2*lambda*w

  SD <- cbind(-Lam, Lam) %>%
    add(
      map_dbl(z, ~{
        sum( h1(.x + q, phi)*(mu(.x + q, phi) - y) ) +
          2*lambda*sum(w*sign(.x - z))
      })
    )

  signSD <- sign(SD) %>% rowSums

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

    solution <- get.sol(y, q, q0, phi, lambda, w, w., d, mu., b, h, z, signSD)
  }

  return(solution)
}
