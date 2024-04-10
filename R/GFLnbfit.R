#v0.1.0

#' @importFrom magrittr %>% add is_greater_than set_names
#' @importFrom purrr exec map map_dbl

GFLnb.fit <- function(
    y, D, Beta0, phi0, betaM, dispM, r, maxit, tol=1e-5,
    h, h1, mu, mu., b, V, lambda, W, idxD, get.sol, get.dsp, get.dsp0, PSN,
    yli, qli, q0s, n, m, ns, q
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
  phi <- phi0

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
      Beta.aft[j] <- GFLnb.fit1(
        h, h1, mu, mu., b, lambda, DCP$w, DCP$z, yli[[j]], qli[[j]], q0s[j], phi, get.sol, ns[j]
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
      phi <- dispM
    } else if(t < m)
    {
      E <- split(1:m, match(Beta.aft, Xi))
      idxF <- map_dbl(E, length) %>% is_greater_than(1) %>% which

      for(k in idxF)
      {
        Ek <- E[[k]]

        FCP <- getFCP(Ek, D, W, Beta.aft)

        Beta.aft[Ek] <- Xi[k] <- GFLnb.fit1(
          h, h1, mu, mu., b, lambda, FCP$w, FCP$z,
          yli[Ek] %>% unlist, qli[Ek] %>% unlist,
          ifelse(length(unique(q0s[Ek]))==1, unique(q0s[Ek]), NA),
          phi, get.sol, sum(ns[Ek])
        )
      }

      .t <- unique(Xi) %>% length

      phi <- get.dsp(Beta.aft, phi, mu, q, y, V, n, .t, tol, m, ns)
    } else
    {
      .t <- m

      phi <- get.dsp0(Beta.aft, phi, mu, q, y, V, n, .t, tol, m, ns)
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

  if(PSN & .t == n)
  {
    disp <- Inf
    chisq <- 0
  } else
  {
    disp <- ifelse(.t == 1, dispM, phi)

    Zeta <- mu(Eta, disp)
    chisq <- sum(((y - Zeta)^2) / V(Zeta, disp))
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
