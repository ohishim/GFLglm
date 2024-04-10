#v0.1.0

#' @importFrom magrittr %>% add is_greater_than set_colnames set_names
#' @importFrom purrr exec map map_dbl


GFLglm.fit <- function(
    y, D, Beta0, betaM, dispM, r, maxit, tol=1e-5,
    h, h1, mu, mu., b, V, lambda, W, idxD, getCWS,
    yli, ali, qli, q0s, ns, get.sol, m, q, n
){

  getDCP <- function(w, z, r){
    # get descent cycle penalty

    .z <- unique(z)

    if(length(.z) == r)
    {
      return(list(
        w = w[order(z)],
        z = sort(z)
      ))
    } else
    {
      lab <- match(z, .z)
      .w <- split(w, lab) %>% map_dbl(sum) %>% set_names(NULL)

      return(list(
        w = .w[order(.z)],
        z = sort(.z)
      ))
    }
  }

  getFCP <- function(Ek, D, W, Beta.aft){

    WZ <- map(Ek, ~{
      i <- .x
      Wi <- W[[i]]
      Di <- D[[i]]
      idx <- setdiff(Di, Ek)

      if(length(idx) > 0)
      {
        return(cbind(Wi[which(Di %in% idx)], Beta.aft[idx]))
      }
    }) %>% exec(rbind, !!!.)

    return(
      getDCP(WZ[,1], WZ[,2], nrow(WZ))
    )
  }

  Beta.aft <- Beta0

  run <- TRUE
  it <- 0

  while(run)
  {
    it <- it + 1
    Beta.bef <- Beta.aft

    ############################################################################
    ###   descent cycle
    ############################################################################

    for(j in idxD)
    {
      DCP <- getDCP(W[[j]], Beta.aft[D[[j]]], r[j])
      Beta.aft[j] <- GFLglm.fit1(
        h, h1, mu, mu., b, lambda, DCP$w, DCP$z, yli[[j]], ali[[j]], qli[[j]], q0s[j], getCWS, get.sol,
        d=ns[j]
      )
    } #end for i

    ############################################################################
    ###   fusion cycle
    ############################################################################

    Xi <- unique(Beta.aft)
    t <- length(Xi)

    if(t == 1)
    {
      .t <- 1

      Beta.aft[1:m] <- Xi <- betaM
    } else if(t < m)
    {
      E <- split(1:m, match(Beta.aft, Xi))
      idxF <- map_dbl(E, length) %>% is_greater_than(1) %>% which

      for(k in idxF)
      {
        Ek <- E[[k]]

        FCP <- getFCP(Ek, D, W, Beta.aft)

        Beta.aft[Ek] <- Xi[k] <- GFLglm.fit1(
          h, h1, mu, mu., b, lambda, FCP$w, FCP$z,
          yli[Ek] %>% unlist, ali[Ek] %>% unlist, qli[Ek] %>% unlist,
          ifelse(length(unique(q0s[Ek]))==1, unique(q0s[Ek]), NA),
          getCWS, get.sol,
          d=ns[Ek] %>% sum
        )
      }

      .t <- unique(Xi) %>% length
    } else
    {
      .t <- m
    }

    ############################################################################
    ###   convergence check
    ############################################################################

    if(t == .t)
    {
      dif <- max((Beta.aft - Beta.bef)^2 / Beta.bef^2, na.rm=T)
    } else
    {
      dif <- Inf
    }

    if(it == maxit)
    {
      run <- FALSE
    } else
    {
      run <- dif > tol
    }

  } #end while

  Eta <- map(1:m, ~rep(Beta.aft[.x], ns[.x])) %>% unlist %>% add(q)
  Zeta <- mu(Eta)
  chisq <- sum(((y - Zeta)^2) / V(Zeta))

  if(.t == 1)
  {
    disp <- dispM
  } else if(.t == n)
  {
    disp <- Inf
  } else
  {
    disp <- chisq / (n - .t)
  }

  return(list(
    beta = Beta.aft,
    eta = Eta,
    chisq = chisq,
    disp = disp,
    df = .t,
    convergence = dif <= tol
  ))
}
