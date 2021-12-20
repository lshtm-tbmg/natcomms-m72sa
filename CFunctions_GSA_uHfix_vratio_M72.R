FitGo <- function(cntry, Vx, Fit, InitV, TimeScale, Plot, C) {
  require(here)
  setwd(here())
  fit_env <- environment()

  FitV <- c("psz1900", "rmort", "neta", "rmortTB", "CDRscale", "CDRscaleE", "alpha")
  InitialV <- c("run", "dt", "prop")

  if (length(InitV) == 3) {
    for (i in 1:length(InitV)) {
      assign(InitialV[i], InitV[i], envir = .GlobalEnv)
    }
  } else {
    for (i in 1:2) {
      assign(InitialV[i], as.numeric(InitV[i]), envir = .GlobalEnv)
    }
    assign(InitialV[3], as.numeric(InitV[3:6]), envir = .GlobalEnv)
  }
  assign("cntry", cntry, envir = .GlobalEnv)
  assign("year1", TimeScale[1], envir = .GlobalEnv)
  assign("yearend", TimeScale[2], envir = .GlobalEnv)

  hchild <- fchild
  hadult <- fadult
  helderly <- felderly

  chiyrs <- 15
  aduyrs <- 50
  yaduyrs <- 40
  eldyrs <- (Mnage - chiyrs - aduyrs)

  voldadu <- (vadult + velderly) / 2
  noldadu <- (n + nelderly) / 2
  roldadu <- (radult + relderly) / 2
  CDRscaleO <- (CDRscale + CDRscaleE) / 2

  vHadult <- vHratio * vadult
  exit <- 0

  p <- c((rep(pchild, l = chiyrs)), (rep(padult, l = aduyrs)), (rep(pelderly, l = eldyrs)))
  f <- c((rep(fchild, l = chiyrs)), (rep(fadult, l = aduyrs)), (rep(felderly, l = eldyrs)))
  h <- c((rep(hchild, l = chiyrs)), (rep(hadult, l = aduyrs)), (rep(helderly, l = eldyrs)))
  v <- c((rep(vchild, l = chiyrs)), (rep(vadult, l = aduyrs)), (rep(velderly, l = eldyrs)))
  r <- c((rep(rchild, l = chiyrs)), (rep(radult, l = yaduyrs)), (rep(roldadu, l = (aduyrs -
    yaduyrs))), (rep(relderly, l = eldyrs)))
  n <- c((rep(n, l = chiyrs)), (rep(n, l = yaduyrs)), (rep(noldadu, l = (aduyrs -
    yaduyrs))), (rep(nelderly, l = eldyrs)))

  fH <- c((rep(fHchild, l = chiyrs)), (rep(fHadult, l = (aduyrs + eldyrs))))
  vH <- c((rep(vHchild, l = chiyrs)), (rep(vHadult, l = (aduyrs + eldyrs))))
  pH <- c((rep(pHchild, l = chiyrs)), (rep(pHadult, l = (aduyrs + eldyrs))))
  rH <- c((rep(rHchild, l = chiyrs)), (rep(rHadult, l = (aduyrs + eldyrs))))

  hH <- fH
  gH <- xH
  g <- x

  nH <- n
  rmortH <- 0

  times <- seq(year1, (yearend + (1 - dt)), dt)
  steps <- length(times)

  for (i in 1:length(Fit)) {
    assign(FitV[i], Fit[i], envir = .GlobalEnv)
  }

  if (length(Fit) > 3) {
    if (rmortTB < 0) {
      ui <- rmortTB * (ui) + ui
      uni <- rmortTB * (uni) + uni
    } else {
      ui <- rmortTB * (1 - ui) + ui
      uni <- rmortTB * (1 - uni) + uni
    }
  }

  if (length(Fit) > 3) {
    if (rmortTBH < 0) {
      uiH <- rmortTBH * (uiH) + uiH
      uniH <- rmortTBH * (uniH) + uniH
    } else {
      uiH <- rmortTBH * (1 - uiH) + uiH
      uniH <- rmortTBH * (1 - uniH) + uniH
    }
  }

  if (length(Fit) > 3) {
    if (alpha < 0) {
      pH <- alpha * pH + pH
      vH <- alpha * vH + vH
      rH <- alpha * rH + rH
    } else {
      pH <- alpha * (1 - pH) + pH
      vH <- alpha * (1 - vH) + vH
      rH <- alpha * (1 - rH) + rH
    }
  }

  uichild <- ui * uiscaleC
  uiadult <- ui * uiscaleA
  uielderly <- ui * uiscaleE
  unichild <- uni * uiscaleC
  uniadult <- uni * uiscaleA
  unielderly <- uni * uiscaleE

  uichildH <- uiH * uiscaleC
  uiadultH <- uiH * uiscaleA
  uielderlyH <- uiH * uiscaleE
  unichildH <- uniH * uiscaleC
  uniadultH <- uniH * uiscaleA
  unielderlyH <- uniH * uiscaleE

  ui <- c((rep(uichild, l = chiyrs)), (rep(uiadult, l = aduyrs)), (rep(uielderly,
    l = eldyrs
  )))
  uni <- c((rep(unichild, l = chiyrs)), (rep(uniadult, l = aduyrs)), (rep(unielderly,
    l = eldyrs
  )))

  uiH <- c((rep(uichildH, l = chiyrs)), (rep(uiadultH, l = aduyrs)), (rep(uielderlyH,
    l = eldyrs
  )))
  uniH <- c((rep(unichildH, l = chiyrs)), (rep(uniadultH, l = aduyrs)), (rep(unielderlyH,
    l = eldyrs
  )))

  is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))

  is.infinite.data.frame <- function(x) do.call(cbind, lapply(x, is.infinite))

  if (length(Fit) < 5) {
    CDRscale <- 0
  }

  if (length(Vx) > 1) {
    assign("vaccine", Vx[1], envir = .GlobalEnv)
    assign("coverage", Vx[2], envir = .GlobalEnv)
    assign("coverageM", Vx[3], envir = .GlobalEnv)
    assign("effI", Vx[4], envir = .GlobalEnv)
    assign("effD", Vx[5], envir = .GlobalEnv)
    assign("D", Vx[6], envir = .GlobalEnv)
    assign("fmass", Vx[7], envir = .GlobalEnv)
    assign("vxage", Vx[8], envir = .GlobalEnv)
    assign("covH", Vx[9], envir = .GlobalEnv)
    assign("VEH", Vx[10], envir = .GlobalEnv)

    source("#VxGenGSA_M72.R")
  } else {
    d <- matrix(0, steps, Mnage)
    dH <- matrix(0, steps, Mnage)
    thetaS <- matrix(0, steps, Mnage)
    thetaL <- matrix(0, steps, Mnage)
    thetaR <- matrix(0, steps, Mnage)
    thetaSH <- matrix(0, steps, Mnage)
    thetaLH <- matrix(0, steps, Mnage)
    thetaRH <- matrix(0, steps, Mnage)
    vaccine <- 0
    effI <- 0
    effD <- 0
    D <- 0
    coverage <- 0
    coverageM <- 0
    fmass <- 0
    vxage <- 0
    VEH <- 0
    covH <- 0
  }

  source("#InitGSA_vratio_M72.R")

  for (k in year1:(yearend)) {
    if (k <= 1950) {
      yr <- 1950
    } else {
      yr <- k
    }

    mortdrop <- 1950

    if (k <= mortdrop) {
      if (rmort < 0) {
        upop <- as.vector(rmort * (mort[1, 2:102]) + mort[1, 2:102])
      } else {
        upop <- as.vector(rmort * (1 - (mort[1, 2:102])) + mort[1, 2:102])
      }
    } else {
      if (rmort < 0) {
        upop <- as.vector(rmort * (mort[1 + yr - 1950, 2:102]) + mort[1 +
          yr - 1950, 2:102])
      } else {
        upop <- as.vector(rmort * (1 - (mort[1 + yr - 1950, 2:102])) + mort[1 +
          yr - 1950, 2:102])
      }
    }

    if (k <= 1990) {
      if (rmortH < 0) {
        uHpop <- as.vector(rmortH * (mortH[1, 2:102]) + mortH[1, 2:102])
      } else {
        uHpop <- as.vector(rmortH * (1 - (mortH[1, 2:102])) + mortH[1, 2:102])
      }
    } else {
      if (rmortH < 0) {
        uHpop <- as.vector(rmortH * (mortH[1 + yr - 1990, 2:102]) + mortH[1 +
          yr - 1990, 2:102])
      } else {
        uHpop <- as.vector(rmortH * (1 - (mortH[1 + yr - 1990, 2:102])) +
          mortH[1 + yr - 1990, 2:102])
      }
    }

    if (k <= 1990) {
      CDR_yr <- 1990
    } else {
      CDR_yr <- k
    }
    if (k <= 1994) {
      CoT_yr <- 1994
    } else {
      CoT_yr <- k
    }

    if (CDRscale < 0) {
      CDRscaled <- rep(min(((CDRscale * cdr[1 + CDR_yr - 1990, ]) + cdr[1 +
        CDR_yr - 1990, ]), 1), (chiyrs + yaduyrs))
    } else {
      CDRscaled <- rep(min((CDRscale * (1 - cdr[1 + CDR_yr - 1990, ]) + cdr[1 +
        CDR_yr - 1990, ]), 1), (chiyrs + yaduyrs))
    }

    if (CDRscaleO < 0) {
      CDRscaledO <- rep(min(((CDRscaleO * cdr[1 + CDR_yr - 1990, ]) + cdr[1 +
        CDR_yr - 1990, ]), 1), (Mnage - (eldyrs + chiyrs + yaduyrs)))
    } else {
      CDRscaledO <- rep(min((CDRscaleO * (1 - cdr[1 + CDR_yr - 1990, ]) + cdr[1 +
        CDR_yr - 1990, ]), 1), (Mnage - (eldyrs + chiyrs + yaduyrs)))
    }

    if (CDRscaleE < 0) {
      CDRscaledE <- rep(min(((CDRscaleE * cdr[1 + CDR_yr - 1990, ]) + cdr[1 +
        CDR_yr - 1990, ]), 1), eldyrs)
    } else {
      CDRscaledE <- rep(min((CDRscaleE * (1 - cdr[1 + CDR_yr - 1990, ]) + cdr[1 +
        CDR_yr - 1990, ]), 1), eldyrs)
    }

    CDR <- c(CDRscaled, CDRscaledO, CDRscaledE)

    if (k == 2010) {
      CDR2010 <- CDR
    }

    CDRH <- CDR

    CoT <- suctt[1 + CoT_yr - 1994]

    if (k < 1990) {
      ind <- 1
    } else {
      ind <- (k - 1990 + 2)
    }

    hiv <- hivall[ind, ]

    artyr <- (art[CDR_yr - 1990 + 1, 1] / 100)
    artyrc <- (art[CDR_yr - 1990 + 1, 2] / 100)

    artredu <- 0.35
    artreduC <- 0.3

    uiHa[1:15] <- ((1 - artyrc) * uiH[1:15]) + (artyrc * (((uiH[1:15] - ui[1:15]) *
      artreduC) + ui[1:15]))
    uiHa[16:Mnage] <- ((1 - artyr) * uiH[16:Mnage]) + (artyr * (((uiH[16:Mnage] -
      ui[16:Mnage]) * artredu) + ui[16:Mnage]))

    uniHa[1:15] <- ((1 - artyrc) * uniH[1:15]) + (artyrc * (((uniH[1:15] - uni[1:15]) *
      artreduC) + uni[1:15]))
    uniHa[16:Mnage] <- ((1 - artyr) * uniH[16:Mnage]) + (artyr * (((uniH[16:Mnage] -
      uni[16:Mnage]) * artredu) + uni[16:Mnage]))

    pHa[1:15] <- ((1 - artyrc) * pH[1:15]) + (artyrc * (((pH[1:15] - p[1:15]) *
      artreduC) + p[1:15]))
    pHa[16:Mnage] <- ((1 - artyr) * pH[16:Mnage]) + (artyr * (((pH[16:Mnage] -
      p[16:Mnage]) * artredu) + p[16:Mnage]))

    vHa[1:15] <- ((1 - artyrc) * vH[1:15]) + (artyrc * (((vH[1:15] - v[1:15]) *
      artreduC) + v[1:15]))
    vHa[16:Mnage] <- ((1 - artyr) * vH[16:Mnage]) + (artyr * (((vH[16:Mnage] -
      v[16:Mnage]) * artredu) + v[16:Mnage]))

    rHa[1:15] <- ((1 - artyrc) * rH[1:15]) + (artyrc * (((rH[1:15] - r[1:15]) *
      artreduC) + r[1:15]))
    rHa[16:Mnage] <- ((1 - artyr) * rH[16:Mnage]) + (artyr * (((rH[16:Mnage] -
      r[16:Mnage]) * artredu) + r[16:Mnage]))

    xHa <- ((1 - artyr) * xH) + (artyr * (((xH - x) * artredu) + x))
    gHa <- xHa

    if (k <= 1950) {
      yr <- 1950
    } else {
      yr <- k
    }

    fertdrop <- 1950

    if (k < fertdrop) {
      br <- bb[1]
      if (k == year1) {
        B <- round(br * psizeALL[1])
      } else {
        B <- round(br * psizeALL[((k - year1) * (1 / dt))])
      }
    } else {
      B <- round(bb[1 + yr - 1950] * psizeALL[((k - year1) * (1 / dt))])
    }

    if (k < 1990) {
      BH <- 0
    } else {
      brh <- H04rate[(k - 1990 + 1)]
      BH <- round(brh * psizematrix[((k - year1) * (1 / dt)), 1])
    }

    if (k > year1) {
      i <- ((1 / dt) * (k - year1) + 1)
      start <- 1

      lambda[i - 1, 1:Mnage] <- t(neta * (1 - exp(colSums(-(myneta[
        1:4,
        1:Mnage
      ]) * z * ((Imatrix[i - 1, 1:4]) / (psizematrix[i - 1, 1:4]))))))

      AIDSdeaths[i, 1] <- B * uHpop[1] * dt
      AIDSdeaths[i, 2:Mnage] <- agepopall[i - 1, 1:(Mnage - 1)] * uHpop[2:Mnage] *
        dt

      BKdeaths[i, 1] <- (B * upop[1] * dt) - AIDSdeaths[i, 1]
      BKdeaths[i, 2:Mnage] <- (agepopall[i - 1, 1:(Mnage - 1)] * upop[2:Mnage] *
        dt) - AIDSdeaths[i, 2:Mnage]

      uH[1] <- (AIDSdeaths[i, 1] / BH) / dt
      uH[2:Mnage] <- (AIDSdeaths[i, 2:Mnage] / agepopHIV[i - 1, 1:(Mnage -
        1)]) / dt

      u[1] <- (BKdeaths[i, 1] / B) / dt
      u[2:Mnage] <- (BKdeaths[i, 2:Mnage] / agepopall[i - 1, 1:(Mnage - 1)]) / dt

      uH[is.nan(uH)] <- 0
      uH[is.infinite(uH)] <- 0
      uH[is.na(uH)] <- 0

      uH[uH < 0] <- 0
      u[u < 0] <- 0

      assign("u", u, envir = .GlobalEnv)
      assign("lambda", lambda, envir = .GlobalEnv)
      assign("BKdeaths", BKdeaths, envir = .GlobalEnv)

      for (aaa in 1:Mnage) {
        aaaa <- (u[aaa] + uH[aaa] + lambda[i - 1, aaa])

        bbbb <- (u[aaa] + uH[aaa] + vHa[aaa] + (lambda[i - 1, aaa] * pHa[aaa] *
          xHa))

        cccc <- (u[aaa] + uH[aaa] + uiHa[aaa] + nH[aaa])

        dddd <- (u[aaa] + uH[aaa] + uniHa[aaa] + nH[aaa] + w)

        eeee <- (u[aaa] + uH[aaa] + rHa[aaa, 1] + (lambda[i - 1, aaa] *
          gHa))

        if (((is.nan(aaaa)) || (is.na(aaaa)) || (is.nan(bbbb)) || (is.na(bbbb)) ||
          (is.nan(cccc)) || (is.na(cccc)) || (is.nan(dddd)) || (is.na(dddd)) ||
          (is.nan(eeee)) || (is.na(eeee))) == TRUE) {
          exit <- 1
          print("nan")
        } else {
          if (((is.infinite(aaaa)) || (is.infinite(bbbb)) || (is.infinite(cccc)) ||
            (is.infinite(dddd)) || (is.infinite(eeee))) == TRUE) {
            exit <- 1
            print("inf")
          } else {
            if ((((aaaa < 0) == TRUE) || ((bbbb < 0) == TRUE) || ((cccc <
              0) == TRUE) || ((dddd < 0) == TRUE) || ((eeee < 0) == TRUE)) ==
              TRUE) {
              exit <- 1
              print("neg")
            } else {
              if ((is.nan(u[aaa])) == TRUE) {
                exit <- 1
                print("nan2")
              } else {
                if ((is.nan(exit)) == TRUE) {
                  exit <- 1
                  print("nan3")
                }
              }
            }
          }
        }

        test <- 0

        if (exit == 0) {
          ifelse(((c(aaaa, bbbb, cccc, dddd, eeee) > 2)), test <- 1, print("proceed"))
        }

        if (test == 1) {
          print("exit2")

          uHuse <- (which.max(c(aaaa, bbbb, cccc, dddd, eeee)))

          if (uHuse == 1) {
            uH[aaa] <- ((1 / dt) - (u[aaa] + lambda[i - 1, aaa]))
          }
          if (uHuse == 2) {
            uH[aaa] <- ((1 / dt) - (u[aaa] + (lambda[i - 1, aaa] * pHa[aaa] *
              xHa) + vHa[aaa, 1]))
          }
          if (uHuse == 3) {
            uH[aaa] <- ((1 / dt) - (u[aaa] + uiHa[aaa] + nH[aaa]))
          }
          if (uHuse == 4) {
            uH[aaa] <- ((1 / dt) - (u[aaa] + uniHa[aaa] + nH[aaa] + w))
          }
          if (uHuse == 5) {
            uH[aaa] <- ((1 / dt) - (u[aaa] + (lambda[i - 1, aaa] * gHa) +
              rHa[aaa, 1]))
          }

          aaaa <- (u[aaa] + uH[aaa] + lambda[i - 1, aaa])

          bbbb <- (u[aaa] + uH[aaa] + vHa[aaa] + (lambda[i - 1, aaa] *
            pHa[aaa] * xHa))

          cccc <- (u[aaa] + uH[aaa] + uiHa[aaa] + nH[aaa])

          dddd <- (u[aaa] + uH[aaa] + uniHa[aaa] + nH[aaa] + w)

          eeee <- (u[aaa] + uH[aaa] + rHa[aaa, 1] + (lambda[i - 1, aaa] *
            gHa))

          if ((((aaaa < 0) == TRUE) || ((bbbb < 0) == TRUE) || ((cccc <
            0) == TRUE) || ((dddd < 0) == TRUE) || ((eeee < 0) == TRUE)) ==
            TRUE) {
            exit <- 1
          }
        }
      }

      upopstore[i, ] <- upop[1:Mnage]
      uHpopstore[i, ] <- uHpop[1:Mnage]
      ustore[i, ] <- u[1:Mnage]
      uHstore[i, ] <- uH[1:Mnage]

      j <- 1
      S[i, j] <- B - BH
      SH[i, j] <- BH

      S[i, 2:Mnage] <- S[i - 1, 1:(Mnage - 1)] - (u[1:(Mnage - 1)] + lambda[i -
        1, 1:(Mnage - 1)]) * S[i - 1, 1:(Mnage - 1)] * dt - hiv[1:(Mnage -
        1)] * S[i - 1, 1:(Mnage - 1)]
      L[i, 2:Mnage] <- L[i - 1, 1:(Mnage - 1)] + lambda[i - 1, 1:(Mnage -
        1)] * (1 - p[1:(Mnage - 1)]) * (S[i - 1, 1:(Mnage - 1)] + g * R[i -
        1, 1:(Mnage - 1)]) * dt - (v[1:(Mnage - 1)] + lambda[i - 1, 1:(Mnage -
        1)] * p[1:(Mnage - 1)] * x + u[1:(Mnage - 1)]) * L[i - 1, 1:(Mnage -
        1)] * dt - hiv[1:(Mnage - 1)] * L[i - 1, 1:(Mnage - 1)]

      new_infect[i, 2:Mnage] <- lambda[i - 1, 1:(Mnage - 1)] * S[
        i - 1,
        1:(Mnage - 1)
      ] * dt + lambda[i - 1, 1:(Mnage - 1)] * (x * L[i -
        1, 1:(Mnage - 1)] + g * R[i - 1, 1:(Mnage - 1)]) * dt

      new_I_react[i, 2:Mnage] <- v[1:(Mnage - 1)] * f[1:(Mnage - 1)] * (L[i -
        1, 1:(Mnage - 1)]) * dt + r[1:(Mnage - 1)] * h[1:(Mnage - 1)] *
        R[i - 1, 1:(Mnage - 1)] * dt
      new_NI_react[i, 2:Mnage] <- v[1:(Mnage - 1)] * (1 - f[1:(Mnage - 1)]) *
        L[i - 1, 1:(Mnage - 1)] * dt + r[1:(Mnage - 1)] * (1 - h[1:(Mnage -
          1)]) * R[i - 1, 1:(Mnage - 1)] * dt
      new_actv_react[i, 2:Mnage] <- v[1:(Mnage - 1)] * (L[i - 1, 1:(Mnage -
        1)]) * dt + r[1:(Mnage - 1)] * R[i - 1, 1:(Mnage - 1)] * dt
      new_actv_inf[i, 2:Mnage] <- lambda[i - 1, 1:(Mnage - 1)] * p[1:(Mnage -
        1)] * S[i - 1, 1:(Mnage - 1)] * dt + lambda[i - 1, 1:(Mnage - 1)] *
        p[1:(Mnage - 1)] * x * (L[i - 1, 1:(Mnage - 1)]) * dt + lambda[i -
        1, 1:(Mnage - 1)] * p[1:(Mnage - 1)] * g * R[i - 1, 1:(Mnage -
        1)] * dt
      new_I[i, 2:Mnage] <- lambda[i - 1, 1:(Mnage - 1)] * p[1:(Mnage - 1)] *
        f[1:(Mnage - 1)] * (S[i - 1, 1:(Mnage - 1)] + g * R[i - 1, 1:(Mnage -
          1)]) * dt + (v[1:(Mnage - 1)] + lambda[i - 1, 1:(Mnage - 1)] *
          p[1:(Mnage - 1)] * x) * f[1:(Mnage - 1)] * L[i - 1, 1:(Mnage -
          1)] * dt + r[1:(Mnage - 1)] * h[1:(Mnage - 1)] * R[i - 1, 1:(Mnage -
          1)] * dt + w * NI[i - 1, 1:(Mnage - 1)] * dt
      new_I_noconv[i, 2:Mnage] <- lambda[i - 1, 1:(Mnage - 1)] * p[1:(Mnage -
        1)] * f[1:(Mnage - 1)] * (S[i - 1, 1:(Mnage - 1)] + g * R[
        i - 1,
        1:(Mnage - 1)
      ]) * dt + (v[1:(Mnage - 1)] + lambda[i - 1, 1:(Mnage -
        1)] * p[1:(Mnage - 1)] * x) * f[1:(Mnage - 1)] * L[i - 1, 1:(Mnage -
        1)] * dt + r[1:(Mnage - 1)] * h[1:(Mnage - 1)] * R[i - 1, 1:(Mnage -
        1)] * dt
      new_NI[i, 2:Mnage] <- lambda[i - 1, 1:(Mnage - 1)] * p[1:(Mnage -
        1)] * (1 - f[1:(Mnage - 1)]) * (S[i - 1, 1:(Mnage - 1)] + g * R[i -
        1, 1:(Mnage - 1)]) * dt + (v[1:(Mnage - 1)] + lambda[i - 1, 1:(Mnage -
        1)] * p[1:(Mnage - 1)] * x) * (1 - f[1:(Mnage - 1)]) * L[
        i - 1,
        1:(Mnage - 1)
      ] * dt + r[1:(Mnage - 1)] * (1 - h[1:(Mnage - 1)]) *
        R[i - 1, 1:(Mnage - 1)] * dt
      new_actv[i, 2:Mnage] <- lambda[i - 1, 1:(Mnage - 1)] * p[1:(Mnage -
        1)] * S[i - 1, 1:(Mnage - 1)] * dt + (v[1:(Mnage - 1)] + lambda[i -
        1, 1:(Mnage - 1)] * p[1:(Mnage - 1)] * x) * L[i - 1, 1:(Mnage -
        1)] * dt + (r[1:(Mnage - 1)] + lambda[i - 1, 1:(Mnage - 1)] * p[1:(Mnage -
        1)] * g) * R[i - 1, 1:(Mnage - 1)] * dt
      new_actv_chk[i, 2:Mnage] <- new_actv_react[i, 2:Mnage] + new_actv_inf[
        i,
        2:Mnage
      ]

      new_notif[i, 2:Mnage] <- CDR[2:(Mnage)] * (new_I[i, 2:Mnage] + e *
        new_NI[i, 2:Mnage])

      R[i, 2:Mnage] <- R[i - 1, 1:(Mnage - 1)] + n[1:(Mnage - 1)] * (I[i -
        1, 1:(Mnage - 1)] + NI[i - 1, 1:(Mnage - 1)]) * dt + CDR[2:(Mnage)] *
        CoT * (new_I[i, 2:Mnage] + e * new_NI[i, 2:Mnage]) - (r[1:(Mnage -
        1)] + g * lambda[i - 1, 1:(Mnage - 1)] + u[1:(Mnage - 1)]) * R[i -
        1, 1:(Mnage - 1)] * dt - hiv[1:(Mnage - 1)] * R[i - 1, 1:(Mnage -
        1)]
      I[i, 2:Mnage] <- I[i - 1, 1:(Mnage - 1)] + (1 - CDR[2:(Mnage)] * CoT) *
        (new_I[i, 2:Mnage]) - (n[1:(Mnage - 1)] + u[1:(Mnage - 1)] + ui[1:(Mnage -
        1)]) * I[i - 1, 1:(Mnage - 1)] * dt - hiv[1:(Mnage - 1)] * I[i -
        1, 1:(Mnage - 1)]
      NI[i, 2:Mnage] <- NI[i - 1, 1:(Mnage - 1)] + ((1 - CDR[2:(Mnage)] *
        CoT * e) * new_NI[i, 2:Mnage]) - (n[1:(Mnage - 1)] + u[1:(Mnage -
        1)] + uni[1:(Mnage - 1)] + w) * NI[i - 1, 1:(Mnage - 1)] * dt -
        hiv[1:(Mnage - 1)] * NI[i - 1, 1:(Mnage - 1)]

      SH[i, 2:Mnage] <- SH[i - 1, 1:(Mnage - 1)] - (u[1:(Mnage - 1)] + uH[1:(Mnage -
        1)] + lambda[i - 1, 1:(Mnage - 1)]) * SH[i - 1, 1:(Mnage - 1)] *
        dt + hiv[1:(Mnage - 1)] * S[i - 1, 1:(Mnage - 1)]

      LH[i, 2:Mnage] <- LH[i - 1, 1:(Mnage - 1)] + lambda[i - 1, 1:(Mnage -
        1)] * (1 - pHa[1:(Mnage - 1)]) * (SH[i - 1, 1:(Mnage - 1)] + gHa *
        RH[i - 1, 1:(Mnage - 1)]) * dt - (vHa[1:(Mnage - 1)] + lambda[i -
        1, 1:(Mnage - 1)] * pHa[1:(Mnage - 1)] * xHa + u[1:(Mnage - 1)] +
        uH[1:(Mnage - 1)]) * LH[i - 1, 1:(Mnage - 1)] * dt + hiv[1:(Mnage -
        1)] * L[i - 1, 1:(Mnage - 1)]

      new_IH_react[i, 2:Mnage] <- vHa[1:(Mnage - 1)] * fH[1:(Mnage - 1)] *
        (LH[i - 1, 1:(Mnage - 1)]) * dt + rHa[1:(Mnage - 1)] * hH[1:(Mnage -
          1)] * RH[i - 1, 1:(Mnage - 1)] * dt

      new_NIH_react[i, 2:Mnage] <- vHa[1:(Mnage - 1)] * (1 - fH[1:(Mnage -
        1)]) * LH[i - 1, 1:(Mnage - 1)] * dt + rHa[1:(Mnage - 1)] * (1 -
        hH[1:(Mnage - 1)]) * RH[i - 1, 1:(Mnage - 1)] * dt

      new_actvH_react[i, 2:Mnage] <- vHa[1:(Mnage - 1)] * (LH[i - 1, 1:(Mnage -
        1)]) * dt + rHa[1:(Mnage - 1)] * RH[i - 1, 1:(Mnage - 1)] * dt

      new_IH[i, 2:Mnage] <- lambda[i - 1, 1:(Mnage - 1)] * pHa[1:(Mnage -
        1)] * fH[1:(Mnage - 1)] * (SH[i - 1, 1:(Mnage - 1)] + gHa * RH[i -
        1, 1:(Mnage - 1)]) * dt + (vHa[1:(Mnage - 1)] + lambda[i - 1, 1:(Mnage -
        1)] * pHa[1:(Mnage - 1)] * xHa) * fH[1:(Mnage - 1)] * LH[
        i - 1,
        1:(Mnage - 1)
      ] * dt + rHa[1:(Mnage - 1)] * hH[1:(Mnage - 1)] *
        RH[i - 1, 1:(Mnage - 1)] * dt + w * NIH[i - 1, 1:(Mnage - 1)] *
        dt

      new_IH_noconv[i, 2:Mnage] <- lambda[i - 1, 1:(Mnage - 1)] * pHa[1:(Mnage -
        1)] * fH[1:(Mnage - 1)] * (SH[i - 1, 1:(Mnage - 1)] + gHa * RH[i -
        1, 1:(Mnage - 1)]) * dt + (vHa[1:(Mnage - 1)] + lambda[i - 1, 1:(Mnage -
        1)] * pHa[1:(Mnage - 1)] * xHa) * fH[1:(Mnage - 1)] * LH[
        i - 1,
        1:(Mnage - 1)
      ] * dt + rHa[1:(Mnage - 1)] * hH[1:(Mnage - 1)] *
        RH[i - 1, 1:(Mnage - 1)] * dt

      new_NIH[i, 2:Mnage] <- lambda[i - 1, 1:(Mnage - 1)] * pHa[1:(Mnage -
        1)] * (1 - fH[1:(Mnage - 1)]) * (SH[i - 1, 1:(Mnage - 1)] + gHa *
        RH[i - 1, 1:(Mnage - 1)]) * dt + (vHa[1:(Mnage - 1)] + lambda[i -
        1, 1:(Mnage - 1)] * pHa[1:(Mnage - 1)] * xHa) * (1 - fH[1:(Mnage -
        1)]) * LH[i - 1, 1:(Mnage - 1)] * dt + rHa[1:(Mnage - 1)] * (1 -
        hH[1:(Mnage - 1)]) * RH[i - 1, 1:(Mnage - 1)] * dt

      new_actvH[i, 2:Mnage] <- lambda[i - 1, 1:(Mnage - 1)] * pHa[1:(Mnage -
        1)] * SH[i - 1, 1:(Mnage - 1)] * dt + (vHa[1:(Mnage - 1)] + lambda[i -
        1, 1:(Mnage - 1)] * pHa[1:(Mnage - 1)] * xHa) * LH[i - 1, 1:(Mnage -
        1)] * dt + (rHa[1:(Mnage - 1)] + lambda[i - 1, 1:(Mnage - 1)] *
        pHa[1:(Mnage - 1)] * gHa) * RH[i - 1, 1:(Mnage - 1)] * dt

      new_notifH[i, 2:Mnage] <- CDRH[2:(Mnage)] * (new_IH[i, 2:Mnage] +
        e * new_NIH[i, 2:Mnage])

      RH[i, 2:Mnage] <- RH[i - 1, 1:(Mnage - 1)] + nH[1:(Mnage - 1)] * (IH[i -
        1, 1:(Mnage - 1)] + NIH[i - 1, 1:(Mnage - 1)]) * dt + CDRH[2:(Mnage)] *
        CoT * (new_IH[i, 2:Mnage] + e * new_NIH[i, 2:Mnage]) - (rHa[1:(Mnage -
        1)] + gHa * lambda[i - 1, 1:(Mnage - 1)] + u[1:(Mnage - 1)] + uH[1:(Mnage -
        1)]) * RH[i - 1, 1:(Mnage - 1)] * dt + hiv[1:(Mnage - 1)] * R[i -
        1, 1:(Mnage - 1)]
      IH[i, 2:Mnage] <- IH[i - 1, 1:(Mnage - 1)] + (1 - CDRH[2:(Mnage)] *
        CoT) * (new_IH[i, 2:Mnage]) - (nH[1:(Mnage - 1)] + u[1:(Mnage -
        1)] + uH[1:(Mnage - 1)] + uiHa[1:(Mnage - 1)]) * IH[i - 1, 1:(Mnage -
        1)] * dt + hiv[1:(Mnage - 1)] * I[i - 1, 1:(Mnage - 1)]
      NIH[i, 2:Mnage] <- NIH[i - 1, 1:(Mnage - 1)] + ((1 - CDRH[2:(Mnage)] *
        CoT * e) * new_NIH[i, 2:Mnage]) - (nH[1:(Mnage - 1)] + u[1:(Mnage -
        1)] + uH[1:(Mnage - 1)] + uniHa[1:(Mnage - 1)] + w) * NIH[
        i - 1,
        1:(Mnage - 1)
      ] * dt + hiv[1:(Mnage - 1)] * NI[i - 1, 1:(Mnage -
        1)]

      Sv[i, 1] <- 0
      Lv[i, 1] <- 0
      Rv[i, 1] <- 0
      Iv[i, 1] <- 0
      NIv[i, 1] <- 0
      new_Iv[i, 1] <- 0
      new_NIv[i, 1] <- 0
      new_Iv_noconv[i, 1] <- 0
      new_notifv[i, 1] <- 0
      SvH[i, 1] <- 0
      LvH[i, 1] <- 0
      RvH[i, 1] <- 0
      IvH[i, 1] <- 0
      NIvH[i, 1] <- 0
      new_IvH[i, 1] <- 0
      new_NIvH[i, 1] <- 0
      new_IvH_noconv[i, 1] <- 0
      new_notifvH[i, 1] <- 0

      Sv[i, 2:Mnage] <- Sv[i - 1, 1:(Mnage - 1)] - (u[1:(Mnage - 1)] + ((1 -
        effI) * lambda[i - 1, 1:(Mnage - 1)])) * Sv[i - 1, 1:(Mnage - 1)] *
        dt - hiv[1:(Mnage - 1)] * Sv[i - 1, 1:(Mnage - 1)]

      Lv[i, 2:Mnage] <- Lv[i - 1, 1:(Mnage - 1)] + ((1 - effI) * lambda[i -
        1, 1:(Mnage - 1)]) * (1 - ((1 - effD) * p[1:(Mnage - 1)])) * (Sv[i -
        1, 1:(Mnage - 1)] + g * Rv[i - 1, 1:(Mnage - 1)]) * dt - ((1 -
        effD) * v[1:(Mnage - 1)] + ((1 - effI) * lambda[i - 1, 1:(Mnage -
        1)]) * ((1 - effD) * p[1:(Mnage - 1)]) * x + u[1:(Mnage - 1)]) *
        Lv[i - 1, 1:(Mnage - 1)] * dt - hiv[1:(Mnage - 1)] * Lv[
        i - 1,
        1:(Mnage - 1)
      ]

      new_Iv[i, 2:Mnage] <- ((1 - effI) * lambda[i - 1, 1:(Mnage - 1)]) *
        ((1 - effD) * p[1:(Mnage - 1)]) * f[1:(Mnage - 1)] * (Sv[
          i - 1,
          1:(Mnage - 1)
        ] + g * Rv[i - 1, 1:((Mnage - 1))]) * dt + ((1 - effD) *
          v[1:(Mnage - 1)] + ((1 - effI) * lambda[i - 1, 1:(Mnage - 1)]) *
            ((1 - effD) * p[1:(Mnage - 1)]) * x) * f[1:(Mnage - 1)] * Lv[i -
          1, 1:(Mnage - 1)] * dt + (1 - effD) * r[1:(Mnage - 1)] * h[1:(Mnage -
          1)] * Rv[i - 1, 1:(Mnage - 1)] * dt + w * NIv[i - 1, 1:(Mnage -
          1)] * dt
      new_Iv_noconv[i, 2:Mnage] <- ((1 - effI) * lambda[i - 1, 1:(Mnage -
        1)]) * ((1 - effD) * p[1:(Mnage - 1)]) * f[1:(Mnage - 1)] * (Sv[i -
        1, 1:(Mnage - 1)] + g * Rv[i - 1, 1:(Mnage - 1)]) * dt + ((1 -
        effD) * v[1:(Mnage - 1)] + ((1 - effI) * lambda[i - 1, 1:(Mnage -
        1)]) * ((1 - effD) * p[1:(Mnage - 1)]) * x) * f[1:(Mnage - 1)] *
        Lv[i - 1, 1:(Mnage - 1)] * dt + (1 - effD) * r[1:(Mnage - 1)] *
        h[1:(Mnage - 1)] * Rv[i - 1, 1:(Mnage - 1)] * dt
      new_NIv[i, 2:Mnage] <- ((1 - effI) * lambda[i - 1, 1:(Mnage - 1)]) *
        ((1 - effD) * p[1:(Mnage - 1)]) * (1 - f[1:(Mnage - 1)]) * (Sv[i -
          1, 1:(Mnage - 1)] + g * Rv[i - 1, 1:(Mnage - 1)]) * dt + ((1 -
          effD) * v[1:(Mnage - 1)] + ((1 - effI) * lambda[i - 1, 1:(Mnage -
          1)]) * ((1 - effD) * p[1:(Mnage - 1)]) * x) * (1 - f[1:(Mnage -
          1)]) * Lv[i - 1, 1:(Mnage - 1)] * dt + (1 - effD) * r[1:(Mnage -
          1)] * (1 - h[1:(Mnage - 1)]) * Rv[i - 1, 1:(Mnage - 1)] * dt
      new_notifv[i, 2:Mnage] <- CDR[2:(Mnage)] * (new_Iv[i, 2:Mnage] + e *
        new_NIv[i, 2:Mnage])
      new_actvv[i, 2:Mnage] <- ((1 - effI) * lambda[i - 1, 1:(Mnage - 1)]) *
        ((1 - effD) * p[1:(Mnage - 1)]) * Sv[i - 1, 1:(Mnage - 1)] * dt +
        (((1 - effD) * v[1:(Mnage - 1)]) + ((1 - effI) * lambda[
          i - 1,
          1:(Mnage - 1)
        ]) * ((1 - effD) * p[1:(Mnage - 1)]) * x) * Lv[i -
          1, 1:(Mnage - 1)] * dt + (((1 - effD) * r[1:(Mnage - 1)]) + ((1 -
          effI) * lambda[i - 1, 1:(Mnage - 1)]) * ((1 - effD) * p[1:(Mnage -
          1)]) * g) * Rv[i - 1, 1:(Mnage - 1)] * dt

      Rv[i, 2:Mnage] <- Rv[i - 1, 1:(Mnage - 1)] + n[1:(Mnage - 1)] * (Iv[i -
        1, 1:(Mnage - 1)] + NIv[i - 1, 1:(Mnage - 1)]) * dt + CDR[2:(Mnage)] *
        CoT * (new_Iv[i, 2:Mnage] + e * new_NIv[i, 2:Mnage]) - ((1 - effD) *
        r[1:(Mnage - 1)] + g * (1 - effI) * lambda[i - 1, 1:(Mnage - 1)] +
        u[1:(Mnage - 1)]) * Rv[i - 1, 1:(Mnage - 1)] * dt - hiv[1:(Mnage -
        1)] * Rv[i - 1, 1:(Mnage - 1)]

      Iv[i, 2:Mnage] <- Iv[i - 1, 1:(Mnage - 1)] + (1 - CDR[2:(Mnage)] *
        CoT) * new_Iv[i, 2:Mnage] - (n[1:(Mnage - 1)] + u[1:((Mnage - 1))] +
        ui[1:(Mnage - 1)]) * Iv[i - 1, 1:(Mnage - 1)] * dt - hiv[1:(Mnage -
        1)] * Iv[i - 1, 1:(Mnage - 1)]

      NIv[i, 2:Mnage] <- NIv[i - 1, 1:(Mnage - 1)] + (1 - CDR[2:(Mnage)] *
        CoT * e) * new_NIv[i, 2:Mnage] - (n[1:(Mnage - 1)] + u[1:(Mnage -
        1)] + uni[1:(Mnage - 1)] + w) * NIv[i - 1, 1:(Mnage - 1)] * dt -
        hiv[1:(Mnage - 1)] * NIv[i - 1, 1:(Mnage - 1)]

      SvH[i, 2:Mnage] <- SvH[i - 1, 1:(Mnage - 1)] - (u[1:(Mnage - 1)] +
        uH[1:(Mnage - 1)] + ((1 - (effI * VEH)) * lambda[i - 1, 1:(Mnage -
          1)])) * SvH[i - 1, 1:(Mnage - 1)] * dt + hiv[1:(Mnage - 1)] * Sv[i -
        1, 1:(Mnage - 1)]

      LvH[i, 2:Mnage] <- LvH[i - 1, 1:(Mnage - 1)] + ((1 - (effI * VEH)) *
        lambda[i - 1, 1:(Mnage - 1)]) * (1 - ((1 - (effD * VEH)) * pHa[1:(Mnage -
        1)])) * (SvH[i - 1, 1:(Mnage - 1)] + gHa * RvH[i - 1, 1:(Mnage -
        1)]) * dt - ((1 - (effD * VEH)) * vHa[1:(Mnage - 1)] + ((1 - (effI *
        VEH)) * lambda[i - 1, 1:(Mnage - 1)]) * ((1 - (effD * VEH)) * pHa[1:(Mnage -
        1)]) * xHa + u[1:(Mnage - 1)] + uH[1:(Mnage - 1)]) * LvH[
        i - 1,
        1:(Mnage - 1)
      ] * dt + hiv[1:(Mnage - 1)] * Lv[i - 1, 1:(Mnage -
        1)]

      new_IvH[i, 2:Mnage] <- ((1 - (effI * VEH)) * lambda[i - 1, 1:(Mnage -
        1)]) * ((1 - (effD * VEH)) * pHa[1:(Mnage - 1)]) * fH[1:(Mnage -
        1)] * (SvH[i - 1, 1:(Mnage - 1)] + gHa * RvH[i - 1, 1:((Mnage -
        1))]) * dt + ((1 - (effD * VEH)) * vHa[1:(Mnage - 1)] + ((1 - (effI *
        VEH)) * lambda[i - 1, 1:(Mnage - 1)]) * ((1 - (effD * VEH)) * pHa[1:(Mnage -
        1)]) * xHa) * fH[1:(Mnage - 1)] * LvH[i - 1, 1:(Mnage - 1)] * dt +
        (1 - (effD * VEH)) * rHa[1:(Mnage - 1)] * hH[1:(Mnage - 1)] * RvH[i -
          1, 1:(Mnage - 1)] * dt + w * NIvH[i - 1, 1:(Mnage - 1)] * dt
      new_IvH_noconv[i, 2:Mnage] <- ((1 - (effI * VEH)) * lambda[
        i - 1,
        1:(Mnage - 1)
      ]) * ((1 - (effD * VEH)) * pHa[1:(Mnage - 1)]) * fH[1:(Mnage -
        1)] * (SvH[i - 1, 1:(Mnage - 1)] + gHa * RvH[i - 1, 1:(Mnage -
        1)]) * dt + ((1 - (effD * VEH)) * vHa[1:(Mnage - 1)] + ((1 - (effI *
        VEH)) * lambda[i - 1, 1:(Mnage - 1)]) * ((1 - (effD * VEH)) * pHa[1:(Mnage -
        1)]) * xHa) * fH[1:(Mnage - 1)] * LvH[i - 1, 1:(Mnage - 1)] * dt +
        (1 - (effD * VEH)) * rHa[1:(Mnage - 1)] * hH[1:(Mnage - 1)] * RvH[i -
          1, 1:(Mnage - 1)] * dt
      new_NIvH[i, 2:Mnage] <- ((1 - (effI * VEH)) * lambda[i - 1, 1:(Mnage -
        1)]) * ((1 - (effD * VEH)) * pHa[1:(Mnage - 1)]) * (1 - fH[1:(Mnage -
        1)]) * (SvH[i - 1, 1:(Mnage - 1)] + gHa * RvH[i - 1, 1:(Mnage -
        1)]) * dt + ((1 - (effD * VEH)) * vHa[1:(Mnage - 1)] + ((1 - (effI *
        VEH)) * lambda[i - 1, 1:(Mnage - 1)]) * ((1 - (effD * VEH)) * pHa[1:(Mnage -
        1)]) * xHa) * (1 - fH[1:(Mnage - 1)]) * LvH[i - 1, 1:(Mnage - 1)] *
        dt + (1 - (effD * VEH)) * rHa[1:(Mnage - 1)] * (1 - hH[1:(Mnage -
        1)]) * RvH[i - 1, 1:(Mnage - 1)] * dt
      new_notifvH[i, 2:Mnage] <- CDRH[2:(Mnage)] * (new_IvH[i, 2:Mnage] +
        e * new_NIvH[i, 2:Mnage])
      new_actvvH[i, 2:Mnage] <- ((1 - (effI * VEH)) * lambda[i - 1, 1:(Mnage -
        1)]) * ((1 - (effD * VEH)) * pHa[1:(Mnage - 1)]) * SvH[i - 1, 1:(Mnage -
        1)] * dt + (((1 - (effD * VEH)) * vHa[1:(Mnage - 1)]) + ((1 - (effI *
        VEH)) * lambda[i - 1, 1:(Mnage - 1)]) * ((1 - (effD * VEH)) * pHa[1:(Mnage -
        1)]) * xHa) * LvH[i - 1, 1:(Mnage - 1)] * dt + (((1 - (effD * VEH)) *
        rHa[1:(Mnage - 1)]) + ((1 - (effI * VEH)) * lambda[i - 1, 1:(Mnage -
        1)]) * ((1 - (effD * VEH)) * pHa[1:(Mnage - 1)]) * gHa) * RvH[i -
        1, 1:(Mnage - 1)] * dt

      RvH[i, 2:Mnage] <- RvH[i - 1, 1:(Mnage - 1)] + nH[1:(Mnage - 1)] *
        (IvH[i - 1, 1:(Mnage - 1)] + NIvH[i - 1, 1:(Mnage - 1)]) * dt +
        CDRH[2:(Mnage)] * CoT * (new_IvH[i, 2:Mnage] + e * new_NIvH[
          i,
          2:Mnage
        ]) - ((1 - (effD * VEH)) * rHa[1:(Mnage - 1)] + gHa *
          (1 - (effI * VEH)) * lambda[i - 1, 1:(Mnage - 1)] + u[1:(Mnage -
          1)] + uH[1:(Mnage - 1)]) * RvH[i - 1, 1:(Mnage - 1)] * dt + hiv[1:(Mnage -
          1)] * Rv[i - 1, 1:(Mnage - 1)]

      IvH[i, 2:Mnage] <- IvH[i - 1, 1:(Mnage - 1)] + (1 - CDRH[2:(Mnage)] *
        CoT) * new_IvH[i, 2:Mnage] - (nH[1:(Mnage - 1)] + u[1:(Mnage -
        1)] + uH[1:(Mnage - 1)] + uiHa[1:(Mnage - 1)]) * IvH[i - 1, 1:(Mnage -
        1)] * dt + hiv[1:(Mnage - 1)] * Iv[i - 1, 1:(Mnage - 1)]

      NIvH[i, 2:Mnage] <- NIvH[i - 1, 1:(Mnage - 1)] + (1 - CDRH[2:(Mnage)] *
        CoT * e) * new_NIvH[i, 2:Mnage] - (nH[1:(Mnage - 1)] + u[1:(Mnage -
        1)] + uH[1:(Mnage - 1)] + uniHa[1:(Mnage - 1)] + w) * NIvH[i -
        1, 1:(Mnage - 1)] * dt + hiv[1:(Mnage - 1)] * NIv[i - 1, 1:(Mnage -
        1)]

      S2 <- S[i, ] + Sv[i, ] * (d[i, ] * (1 - thetaS[i, ])) - thetaS[i, ] * S[i, ]
      L2 <- L[i, ] + Lv[i, ] * (d[i, ] * (1 - thetaL[i, ])) - thetaL[i, ] * L[i, ]
      R2 <- R[i, ] + Rv[i, ] * (d[i, ] * (1 - thetaR[i, ])) - thetaR[i, ] * R[i, ]
      I2 <- I[i, ] + Iv[i, ] * d[i, ]
      NI2 <- NI[i, ] + NIv[i, ] * d[i, ]

      Sv2 <- Sv[i, ] - Sv[i, ] * (d[i, ] * (1 - thetaS[i, ])) + thetaS[i, ] * S[i, ]
      Lv2 <- Lv[i, ] - Lv[i, ] * (d[i, ] * (1 - thetaL[i, ])) + thetaL[i, ] * L[i, ]
      Rv2 <- Rv[i, ] - Rv[i, ] * (d[i, ] * (1 - thetaR[i, ])) + thetaR[i, ] * R[i, ]
      Iv2 <- Iv[i, ] - Iv[i, ] * (d[i, ])
      NIv2 <- NIv[i, ] - NIv[i, ] * (d[i, ])

      SH2 <- SH[i, ] + SvH[i, ] * (dH[i, ] * (1 - thetaSH[i, ])) - thetaSH[i, ] * SH[i, ]
      LH2 <- LH[i, ] + LvH[i, ] * (dH[i, ] * (1 - thetaLH[i, ])) - thetaLH[i, ] * LH[i, ]
      RH2 <- RH[i, ] + RvH[i, ] * (dH[i, ] * (1 - thetaRH[i, ])) - thetaRH[i, ] * RH[i, ]
      IH2 <- IH[i, ] + IvH[i, ] * dH[i, ]
      NIH2 <- NIH[i, ] + NIvH[i, ] * dH[i, ]

      SvH2 <- SvH[i, ] - SvH[i, ] * (dH[i, ] * (1 - thetaSH[i, ])) + thetaSH[i, ] * SH[i, ]
      LvH2 <- LvH[i, ] - LvH[i, ] * (dH[i, ] * (1 - thetaLH[i, ])) + thetaLH[i, ] * LH[i, ]
      RvH2 <- RvH[i, ] - RvH[i, ] * (dH[i, ] * (1 - thetaRH[i, ])) + thetaRH[i, ] * RH[i, ]
      IvH2 <- IvH[i, ] - IvH[i, ] * (dH[i, ])
      NIvH2 <- NIvH[i, ] - NIvH[i, ] * (dH[i, ])

      if (vaccine == 1) {
        num_vac[i, ] <- (S[i, ] + L[i, ] + R[i, ] + Sv[i, ] + Lv[i, ] +
          Rv[i, ] + SH[i, ] + LH[i, ] + RH[i, ] + SvH[i, ] + LvH[i, ] +
          RvH[i, ]) * thetaV2[i, ]
      }
      if (vaccine == 2) {
        num_vac[i, ] <- (S[i, ] + L[i, ] + R[i, ] + Sv[i, ] + Lv[i, ] +
          Rv[i, ] + SH[i, ] + LH[i, ] + RH[i, ] + SvH[i, ] + LvH[i, ] +
          RvH[i, ]) * thetaV2[i, ]
      }
      if (vaccine == 3) {
        num_vac[i, ] <- (S[i, ] + L[i, ] + R[i, ] + Sv[i, ] + Lv[i, ] +
          Rv[i, ] + SH[i, ] + LH[i, ] + RH[i, ] + SvH[i, ] + LvH[i, ] +
          RvH[i, ]) * thetaV2[i, ]
      }

      S[i, ] <- S2
      L[i, ] <- L2
      R[i, ] <- R2
      I[i, ] <- I2
      NI[i, ] <- NI2
      Sv[i, ] <- Sv2
      Lv[i, ] <- Lv2
      Rv[i, ] <- Rv2
      Iv[i, ] <- Iv2
      NIv[i, ] <- NIv2

      SH[i, ] <- SH2
      LH[i, ] <- LH2
      RH[i, ] <- RH2
      IH[i, ] <- IH2
      NIH[i, ] <- NIH2
      SvH[i, ] <- SvH2
      LvH[i, ] <- LvH2
      RvH[i, ] <- RvH2
      IvH[i, ] <- IvH2
      NIvH[i, ] <- NIvH2

      psize[i] <- sum(
        S[i, ], L[i, ], R[i, ], I[i, ], NI[i, ], Sv[i, ],
        Lv[i, ], Rv[i, ], Iv[i, ], NIv[i, ]
      )
      psizeH[i] <- sum(SH[i, ], LH[i, ], RH[i, ], IH[i, ], NIH[i, ], SvH[i, ], LvH[i, ], RvH[i, ], IvH[i, ], NIvH[i, ])
      psizeALL[i] <- sum(S[i, ], L[i, ], R[i, ], I[i, ], NI[i, ], Sv[i, ], Lv[i, ], Rv[i, ], Iv[i, ], NIv[i, ], SH[i, ], LH[i, ], RH[i, ], IH[i, ], NIH[i, ], SvH[i, ], LvH[i, ], RvH[i, ], IvH[i, ], NIvH[i, ])

      psize014[i] <- sum(
        S[i, 1:15], L[i, 1:15], R[i, 1:15], I[i, 1:15],
        NI[i, 1:15], Sv[i, 1:15], Lv[i, 1:15], Iv[i, 1:15], NIv[i, 1:15],
        Rv[i, 1:15]
      )
      psize1554[i] <- sum(
        S[i, 16:55], L[i, 16:55], R[i, 16:55], I[i, 16:55],
        NI[i, 16:55], Sv[i, 16:55], Lv[i, 16:55], Iv[i, 16:55], NIv[
          i,
          16:55
        ], Rv[i, 16:55]
      )
      psize5564[i] <- sum(
        S[i, 56:65], L[i, 56:65], R[i, 56:65], I[i, 56:65],
        NI[i, 56:65], Sv[i, 56:65], Lv[i, 56:65], Iv[i, 56:65], NIv[
          i,
          56:65
        ], Rv[i, 56:65]
      )
      psize65plus[i] <- sum(
        S[i, 66:Mnage], L[i, 66:Mnage], R[i, 66:Mnage],
        I[i, 66:Mnage], NI[i, 66:Mnage], Sv[i, 66:Mnage], Lv[i, 66:Mnage],
        Rv[i, 66:Mnage], Iv[i, 66:Mnage], NIv[i, 66:Mnage]
      )
      psize55plus[i] <- sum(
        S[i, 56:Mnage], L[i, 56:Mnage], R[i, 56:Mnage],
        I[i, 56:Mnage], NI[i, 56:Mnage], Sv[i, 56:Mnage], Lv[i, 56:Mnage],
        Rv[i, 56:Mnage], Iv[i, 56:Mnage], NIv[i, 56:Mnage]
      )
      psize55minus[i] <- sum(
        S[i, 1:55], L[i, 1:55], R[i, 1:55], I[i, 1:55],
        NI[i, 1:55], Sv[i, 1:55], Lv[i, 1:55], Rv[i, 1:55], Iv[i, 1:55],
        NIv[i, 1:55]
      )

      psizeVX[i] <- sum(Sv[i, ], Lv[i, ], Rv[i, ], Iv[i, ], NIv[i, ])
      psizeVX9[i] <- sum(Sv[i, 10], Lv[i, 10], Iv[i, 10], NIv[i, 10], Rv[
        i,
        10
      ])
      psizeVX9100[i] <- sum(
        Sv[i, 10:101], Lv[i, 10:101], Iv[i, 10:101],
        NIv[i, 10:101], Rv[i, 10:101]
      )
      psizeNOVX[i] <- sum(S[i, ], L[i, ], R[i, ], I[i, ], NI[i, ])
      psizeNOVX9[i] <- sum(S[i, 10], L[i, 10], I[i, 10], NI[i, 10], R[
        i,
        10
      ])
      psizeNOVX9100[i] <- sum(
        S[i, 10:101], L[i, 10:101], I[i, 10:101],
        NI[i, 10:101], R[i, 10:101]
      )

      psize1559[i] <- sum(
        S[i, 16:60], L[i, 16:60], R[i, 16:60], I[i, 16:60],
        NI[i, 16:60], Sv[i, 16:60], Lv[i, 16:60], Rv[i, 16:60], Iv[i, 16:60],
        NIv[i, 16:60]
      )
      psize1529[i] <- sum(
        S[i, 16:30], L[i, 16:30], R[i, 16:30], I[i, 16:30],
        NI[i, 16:30], Sv[i, 16:30], Lv[i, 16:30], Rv[i, 16:30], Iv[i, 16:30],
        NIv[i, 16:30]
      )
      psize3044[i] <- sum(
        S[i, 31:45], L[i, 31:45], R[i, 31:45], I[i, 31:45],
        NI[i, 31:45], Sv[i, 31:45], Lv[i, 31:45], Rv[i, 31:45], Iv[i, 31:45],
        NIv[i, 31:45]
      )
      psize4559[i] <- sum(
        S[i, 46:60], L[i, 46:60], R[i, 46:60], I[i, 46:60],
        NI[i, 46:60], Sv[i, 46:60], Lv[i, 46:60], Rv[i, 46:60], Iv[i, 46:60],
        NIv[i, 46:60]
      )
      psize60plus[i] <- sum(
        S[i, 61:Mnage], L[i, 61:Mnage], R[i, 61:Mnage],
        I[i, 61:Mnage], NI[i, 61:Mnage], Sv[i, 61:Mnage], Lv[i, 61:Mnage],
        Rv[i, 61:Mnage], Iv[i, 61:Mnage], NIv[i, 61:Mnage]
      )

      psize0509[i] <- sum(
        S[i, 6:10], L[i, 6:10], R[i, 6:10], I[i, 6:10],
        NI[i, 6:10], Sv[i, 6:10], Lv[i, 6:10], Rv[i, 6:10], Iv[i, 6:10],
        NIv[i, 6:10]
      )
      psize1019[i] <- sum(
        S[i, 11:20], L[i, 11:20], R[i, 11:20], I[i, 11:20],
        NI[i, 11:20], Sv[i, 11:20], Lv[i, 11:20], Rv[i, 11:20], Iv[i, 11:20],
        NIv[i, 11:20]
      )
      psize2029[i] <- sum(
        S[i, 21:30], L[i, 21:30], R[i, 21:30], I[i, 21:30],
        NI[i, 21:30], Sv[i, 21:30], Lv[i, 21:30], Rv[i, 21:30], Iv[i, 21:30],
        NIv[i, 21:30]
      )
      psize3039[i] <- sum(
        S[i, 31:40], L[i, 31:40], R[i, 31:40], I[i, 31:40],
        NI[i, 31:40], Sv[i, 31:40], Lv[i, 31:40], Rv[i, 31:40], Iv[i, 31:40],
        NIv[i, 31:40]
      )
      psize4049[i] <- sum(
        S[i, 41:50], L[i, 41:50], R[i, 41:50], I[i, 41:50],
        NI[i, 41:50], Sv[i, 41:50], Lv[i, 41:50], Rv[i, 41:50], Iv[i, 41:50],
        NIv[i, 41:50]
      )
      psize5059[i] <- sum(
        S[i, 51:60], L[i, 51:60], R[i, 51:60], I[i, 51:60],
        NI[i, 51:60], Sv[i, 51:60], Lv[i, 51:60], Rv[i, 51:60], Iv[i, 51:60],
        NIv[i, 51:60]
      )
      psize6069[i] <- sum(
        S[i, 61:70], L[i, 61:70], R[i, 61:70], I[i, 61:70],
        NI[i, 61:70], Sv[i, 61:70], Lv[i, 61:70], Rv[i, 61:70], Iv[i, 61:70],
        NIv[i, 61:70]
      )
      psize70plus[i] <- sum(
        S[i, 71:Mnage], L[i, 71:Mnage], R[i, 71:Mnage],
        I[i, 71:Mnage], NI[i, 71:Mnage], Sv[i, 71:Mnage], Lv[i, 71:Mnage],
        Rv[i, 71:Mnage], Iv[i, 71:Mnage], NIv[i, 71:Mnage]
      )
      psize5574[i] <- sum(
        S[i, 56:75], L[i, 56:75], R[i, 56:75], I[i, 56:75],
        NI[i, 56:75], Sv[i, 56:75], Lv[i, 56:75], Rv[i, 56:75], Iv[i, 56:75],
        NIv[i, 56:75]
      )
      psize75plus[i] <- sum(
        S[i, 76:Mnage], L[i, 76:Mnage], R[i, 76:Mnage],
        I[i, 76:Mnage], NI[i, 76:Mnage], Sv[i, 76:Mnage], Lv[i, 76:Mnage],
        Rv[i, 76:Mnage], Iv[i, 76:Mnage], NIv[i, 76:Mnage]
      )
      psize1524[i] <- sum(
        S[i, 16:25], L[i, 16:25], R[i, 16:25], I[i, 16:25],
        NI[i, 16:25], Sv[i, 16:25], Lv[i, 16:25], Rv[i, 16:25], Iv[i, 16:25],
        NIv[i, 16:25]
      )
      psize2554[i] <- sum(
        S[i, 26:55], L[i, 26:55], R[i, 26:55], I[i, 26:55],
        NI[i, 26:55], Sv[i, 26:55], Lv[i, 26:55], Rv[i, 26:55], Iv[i, 26:55],
        NIv[i, 26:55]
      )
      psize15plus[i] <- sum(
        S[i, 16:Mnage], L[i, 16:Mnage], R[i, 16:Mnage],
        I[i, 16:Mnage], NI[i, 16:Mnage], Sv[i, 16:Mnage], Lv[i, 16:Mnage],
        Rv[i, 16:Mnage], Iv[i, 16:Mnage], NIv[i, 16:Mnage]
      )

      psize1549[i] <- psize1559[i] - psize5059[i]

      agepopall[i, 1:Mnage] <- S[i, 1:Mnage] + L[i, 1:Mnage] + R[i, 1:Mnage] +
        I[i, 1:Mnage] + NI[i, 1:Mnage] + Sv[i, 1:Mnage] + Lv[i, 1:Mnage] +
        Rv[i, 1:Mnage] + Iv[i, 1:Mnage] + NIv[i, 1:Mnage] + SH[i, 1:Mnage] +
        LH[i, 1:Mnage] + RH[i, 1:Mnage] + IH[i, 1:Mnage] + NIH[i, 1:Mnage] +
        SvH[i, 1:Mnage] + LvH[i, 1:Mnage] + RvH[i, 1:Mnage] + IvH[i, 1:Mnage] +
        NIvH[i, 1:Mnage]

      psize04H[i] <- sum(
        SH[i, 1:5], LH[i, 1:5], RH[i, 1:5], IH[i, 1:5],
        NIH[i, 1:5], SvH[i, 1:5], LvH[i, 1:5], RvH[i, 1:5], IvH[i, 1:5],
        NIvH[i, 1:5]
      )
      agepopHIV[i, 1:Mnage] <- SH[i, 1:Mnage] + LH[i, 1:Mnage] + RH[
        i,
        1:Mnage
      ] + IH[i, 1:Mnage] + NIH[i, 1:Mnage] + SvH[i, 1:Mnage] +
        LvH[i, 1:Mnage] + RvH[i, 1:Mnage] + IvH[i, 1:Mnage] + NIvH[i, 1:Mnage]

      psize014H[i] <- sum(
        SH[i, 1:15], LH[i, 1:15], RH[i, 1:15], IH[
          i,
          1:15
        ], NIH[i, 1:15], SvH[i, 1:15], LvH[i, 1:15], IvH[i, 1:15],
        NIvH[i, 1:15], RvH[i, 1:15]
      )
      psize1554H[i] <- sum(
        SH[i, 16:55], LH[i, 16:55], RH[i, 16:55], IH[
          i,
          16:55
        ], NIH[i, 16:55], SvH[i, 16:55], LvH[i, 16:55], IvH[i, 16:55],
        NIvH[i, 16:55], RvH[i, 16:55]
      )
      psize5564H[i] <- sum(
        SH[i, 56:65], LH[i, 56:65], RH[i, 56:65], IH[
          i,
          56:65
        ], NIH[i, 56:65], SvH[i, 56:65], LvH[i, 56:65], IvH[i, 56:65],
        NIvH[i, 56:65], RvH[i, 56:65]
      )
      psize65plusH[i] <- sum(
        SH[i, 66:Mnage], LH[i, 66:Mnage], RH[i, 66:Mnage],
        IH[i, 66:Mnage], NIH[i, 66:Mnage], SvH[i, 66:Mnage], LvH[i, 66:Mnage],
        RvH[i, 66:Mnage], IvH[i, 66:Mnage], NIvH[i, 66:Mnage]
      )

      psizeVXH[i] <- sum(SvH[i, ], LvH[i, ], RvH[i, ], IvH[i, ], NIvH[i, ])
      psizeVX9H[i] <- sum(
        SvH[i, 10], LvH[i, 10], IvH[i, 10], NIvH[i, 10],
        RvH[i, 10]
      )
      psizeVX9100H[i] <- sum(
        SvH[i, 10:101], LvH[i, 10:101], IvH[i, 10:101],
        NIvH[i, 10:101], RvH[i, 10:101]
      )

      psize1559H[i] <- sum(
        SH[i, 16:60], LH[i, 16:60], RH[i, 16:60], IH[
          i,
          16:60
        ], NIH[i, 16:60], SvH[i, 16:60], LvH[i, 16:60], RvH[i, 16:60],
        IvH[i, 16:60], NIvH[i, 16:60]
      )
      psize60plusH[i] <- sum(
        SH[i, 61:Mnage], LH[i, 61:Mnage], RH[i, 61:Mnage],
        IH[i, 61:Mnage], NIH[i, 61:Mnage], SvH[i, 61:Mnage], LvH[i, 61:Mnage],
        RvH[i, 61:Mnage], IvH[i, 61:Mnage], NIvH[i, 61:Mnage]
      )
      psize15plusH[i] <- sum(
        SH[i, 16:Mnage], LH[i, 16:Mnage], RH[i, 16:Mnage],
        IH[i, 16:Mnage], NIH[i, 16:Mnage], SvH[i, 16:Mnage], LvH[i, 16:Mnage],
        RvH[i, 16:Mnage], IvH[i, 16:Mnage], NIvH[i, 16:Mnage]
      )

      psize0509H[i] <- sum(
        SH[i, 6:10], LH[i, 6:10], RH[i, 6:10], IH[
          i,
          6:10
        ], NIH[i, 6:10], SvH[i, 6:10], LvH[i, 6:10], RvH[i, 6:10],
        IvH[i, 6:10], NIvH[i, 6:10]
      )
      psize1019H[i] <- sum(
        SH[i, 11:20], LH[i, 11:20], RH[i, 11:20], IH[
          i,
          11:20
        ], NIH[i, 11:20], SvH[i, 11:20], LvH[i, 11:20], RvH[i, 11:20],
        IvH[i, 11:20], NIvH[i, 11:20]
      )
      psize2029H[i] <- sum(
        SH[i, 21:30], LH[i, 21:30], RH[i, 21:30], IH[
          i,
          21:30
        ], NIH[i, 21:30], SvH[i, 21:30], LvH[i, 21:30], RvH[i, 21:30],
        IvH[i, 21:30], NIvH[i, 21:30]
      )
      psize3039H[i] <- sum(
        SH[i, 31:40], LH[i, 31:40], RH[i, 31:40], IH[
          i,
          31:40
        ], NIH[i, 31:40], SvH[i, 31:40], LvH[i, 31:40], RvH[i, 31:40],
        IvH[i, 31:40], NIvH[i, 31:40]
      )
      psize4049H[i] <- sum(
        SH[i, 41:50], LH[i, 41:50], RH[i, 41:50], IH[
          i,
          41:50
        ], NIH[i, 41:50], SvH[i, 41:50], LvH[i, 41:50], RvH[i, 41:50],
        IvH[i, 41:50], NIvH[i, 41:50]
      )
      psize5059H[i] <- sum(
        SH[i, 51:60], LH[i, 51:60], RH[i, 51:60], IH[
          i,
          51:60
        ], NIH[i, 51:60], SvH[i, 51:60], LvH[i, 51:60], RvH[i, 51:60],
        IvH[i, 51:60], NIvH[i, 51:60]
      )
      psize6069H[i] <- sum(
        SH[i, 61:70], LH[i, 61:70], RH[i, 61:70], IH[
          i,
          61:70
        ], NIH[i, 61:70], SvH[i, 61:70], LvH[i, 61:70], RvH[i, 61:70],
        IvH[i, 61:70], NIvH[i, 61:70]
      )
      psize70plusH[i] <- sum(
        SH[i, 71:Mnage], LH[i, 71:Mnage], RH[i, 71:Mnage],
        IH[i, 71:Mnage], NIH[i, 71:Mnage], SvH[i, 71:Mnage], LvH[i, 71:Mnage],
        RvH[i, 71:Mnage], IvH[i, 71:Mnage], NIvH[i, 71:Mnage]
      )

      psize1549H[i] <- psize1559H[i] - psize5059H[i]

      artsize[i, 1] <- psize014H[i] * artyrc
      artsize[i, 2] <- psize15plusH[i] * artyr

      psizematrix[i, 1] <- sum(
        S[i, 1:5], L[i, 1:5], R[i, 1:5], I[i, 1:5],
        NI[i, 1:5], Sv[i, 1:5], Lv[i, 1:5], Rv[i, 1:5], Iv[i, 1:5], NIv[
          i,
          1:5
        ], SH[i, 1:5], LH[i, 1:5], RH[i, 1:5], IH[i, 1:5], NIH[
          i,
          1:5
        ], SvH[i, 1:5], LvH[i, 1:5], RvH[i, 1:5], IvH[i, 1:5], NIvH[
          i,
          1:5
        ]
      )
      psizematrix[i, 2] <- sum(
        S[i, 6:10], L[i, 6:10], R[i, 6:10], I[
          i,
          6:10
        ], NI[i, 6:10], Sv[i, 6:10], Lv[i, 6:10], Rv[i, 6:10], Iv[
          i,
          6:10
        ], NIv[i, 6:10], SH[i, 6:10], LH[i, 6:10], RH[i, 6:10], IH[
          i,
          6:10
        ], NIH[i, 6:10], SvH[i, 6:10], LvH[i, 6:10], RvH[i, 6:10],
        IvH[i, 6:10], NIvH[i, 6:10]
      )
      psizematrix[i, 3] <- sum(
        S[i, 11:15], L[i, 11:15], R[i, 11:15], I[
          i,
          11:15
        ], NI[i, 11:15], Sv[i, 11:15], Lv[i, 11:15], Rv[i, 11:15],
        Iv[i, 11:15], NIv[i, 11:15], SH[i, 11:15], LH[i, 11:15], RH[
          i,
          11:15
        ], IH[i, 11:15], NIH[i, 11:15], SvH[i, 11:15], LvH[i, 11:15],
        RvH[i, 11:15], IvH[i, 11:15], NIvH[i, 11:15]
      )
      psizematrix[i, 4] <- sum(
        S[i, 16:Mnage], L[i, 16:Mnage], R[i, 16:Mnage],
        I[i, 16:Mnage], NI[i, 16:Mnage], Sv[i, 16:Mnage], Lv[i, 16:Mnage],
        Rv[i, 16:Mnage], Iv[i, 16:Mnage], NIv[i, 16:Mnage], SH[i, 16:Mnage],
        LH[i, 16:Mnage], RH[i, 16:Mnage], IH[i, 16:Mnage], NIH[i, 16:Mnage],
        SvH[i, 16:Mnage], LvH[i, 16:Mnage], RvH[i, 16:Mnage], IvH[i, 16:Mnage],
        NIvH[i, 16:Mnage]
      )

      Imatrix[i, 1] <- sum(I[i, 1:5], Iv[i, 1:5], IH[i, 1:5], IvH[i, 1:5])
      Imatrix[i, 2] <- sum(I[i, 6:10], Iv[i, 6:10], IH[i, 6:10], IvH[
        i,
        6:10
      ])
      Imatrix[i, 3] <- sum(I[i, 11:15], Iv[i, 11:15], IH[i, 11:15], IvH[
        i,
        11:15
      ])
      Imatrix[i, 4] <- sum(
        I[i, 16:Mnage], Iv[i, 16:Mnage], IH[i, 16:Mnage],
        IvH[i, 16:Mnage]
      )

      TBDeaths[i, 1:Mnage] <- dt * ((ui[1:Mnage] * I[i - 1, 1:Mnage]) +
        (uni[1:Mnage] * NI[i - 1, 1:Mnage]) + (ui[1:Mnage] * Iv[
          i - 1,
          1:Mnage
        ]) + (uni[1:Mnage] * NIv[i - 1, 1:Mnage]) + (uiHa[1:Mnage] *
          IH[i - 1, 1:Mnage]) + (uniHa[1:Mnage] * NIH[i - 1, 1:Mnage]) +
        (uiHa[1:Mnage] * IvH[i - 1, 1:Mnage]) + (uniHa[1:Mnage] * NIvH[i -
          1, 1:Mnage]))
      TBDeathsN[i, 1:Mnage] <- dt * ((ui[1:Mnage] * I[i - 1, 1:Mnage]) +
        (uni[1:Mnage] * NI[i - 1, 1:Mnage]) + (ui[1:Mnage] * Iv[
          i - 1,
          1:Mnage
        ]) + (uni[1:Mnage] * NIv[i - 1, 1:Mnage]))
      TBDeathsH[i, 1:Mnage] <- dt * ((uiHa[1:Mnage] * IH[i - 1, 1:Mnage]) +
        (uniHa[1:Mnage] * NIH[i - 1, 1:Mnage]) + (uiHa[1:Mnage] * IvH[i -
          1, 1:Mnage]) + (uniHa[1:Mnage] * NIvH[i - 1, 1:Mnage]))

      ADeaths[i, ] <- dt * (u * S[i - 1, ] + u * L[i - 1, ] + (u + ui) *
        I[i - 1, ] + (u + uni) * NI[i - 1, ] + u * R[i - 1, ] + u * Sv[i -
        1, ] + u * Lv[i - 1, ] + u * Rv[i - 1, ] + (u + ui) * Iv[i - 1, ] + (u + uni) * NIv[i - 1, ] + (u + uH) * SH[i - 1, ] + (u + uH) *
        LH[i - 1, ] + (u + uH + uiHa) * IH[i - 1, ] + (u + uH + uniHa) *
        NIH[i - 1, ] + (u + uH) * RH[i - 1, ] + (u + uH) * SvH[i - 1, ] +
        (u + uH) * LvH[i - 1, ] + (u + uH) * RvH[i - 1, ] + (u + uH + uiHa) *
          IvH[i - 1, ] + (u + uH + uniHa) * NIvH[i - 1, ])
      ADeathsN[i, ] <- dt * (u * S[i - 1, ] + u * L[i - 1, ] + (u + ui) *
        I[i - 1, ] + (u + uni) * NI[i - 1, ] + u * R[i - 1, ] + u * Sv[i -
        1, ] + u * Lv[i - 1, ] + u * Rv[i - 1, ] + (u + ui) * Iv[i - 1, ] + (u + uni) * NIv[i - 1, ])
      ADeathsH[i, ] <- dt * ((u + uH) * SH[i - 1, ] + (u + uH) * LH[i -
        1, ] + (u + uH + uiHa) * IH[i - 1, ] + (u + uH + uniHa) * NIH[i -
        1, ] + (u + uH) * RH[i - 1, ] + (u + uH) * SvH[i - 1, ] + (u +
        uH) * LvH[i - 1, ] + (u + uH) * RvH[i - 1, ] + (u + uH + uiHa) *
        IvH[i - 1, ] + (u + uH + uniHa) * NIvH[i - 1, ])
    }

    for (i in (2 + (1 / dt) * (k - year1)):((1 / dt) * (k - year1) + 1 / dt)) {
      start <- 0

      lambda[i - 1, 1:Mnage] <- t(neta * (1 - exp(colSums(-(myneta[1:4, 1:Mnage]) *
        z * ((Imatrix[i - 1, 1:4]) / (psizematrix[i - 1, 1:4]))))))

      AIDSdeaths[i, 1:Mnage] <- agepopall[i - 1, 1:Mnage] * uHpop[1:Mnage] *
        dt

      BKdeaths[i, 1:Mnage] <- (agepopall[i - 1, 1:Mnage] * upop[1:Mnage] *
        dt) - AIDSdeaths[i, 1:Mnage]

      uH[1:Mnage] <- (AIDSdeaths[i, 1:Mnage] / agepopHIV[i - 1, 1:Mnage]) / dt

      u[1:Mnage] <- (BKdeaths[i, 1:Mnage] / agepopall[i - 1, 1:Mnage]) / dt

      uH[is.nan(uH)] <- 0
      uH[is.na(uH)] <- 0
      uH[is.infinite(uH)] <- 0

      uH[uH < 0] <- 0
      u[u < 0] <- 0

      assign("u", u, envir = .GlobalEnv)
      assign("lambda", lambda, envir = .GlobalEnv)
      assign("BKdeaths", BKdeaths, envir = .GlobalEnv)

      for (aaa in 1:Mnage) {
        aaaa <- (u[aaa] + uH[aaa] + lambda[i - 1, aaa])

        bbbb <- (u[aaa] + uH[aaa] + vHa[aaa] + (lambda[i - 1, aaa] * pHa[aaa] *
          xHa))

        cccc <- (u[aaa] + uH[aaa] + uiHa[aaa] + nH[aaa])

        dddd <- (u[aaa] + uH[aaa] + uniHa[aaa] + nH[aaa] + w)

        eeee <- (u[aaa] + uH[aaa] + rHa[aaa, 1] + (lambda[i - 1, aaa] * gHa))

        if (((is.nan(aaaa)) || (is.na(aaaa)) || (is.nan(bbbb)) || (is.na(bbbb)) ||
          (is.nan(cccc)) || (is.na(cccc)) || (is.nan(dddd)) || (is.na(dddd)) ||
          (is.nan(eeee)) || (is.na(eeee))) == TRUE) {
          exit <- 1
          print("nan")
        } else {
          if (((is.infinite(aaaa)) || (is.infinite(bbbb)) || (is.infinite(cccc)) ||
            (is.infinite(dddd)) || (is.infinite(eeee))) == TRUE) {
            exit <- 1
            print("inf")
          } else {
            if ((((aaaa < 0) == TRUE) || ((bbbb < 0) == TRUE) || ((cccc <
              0) == TRUE) || ((dddd < 0) == TRUE) || ((eeee < 0) == TRUE)) ==
              TRUE) {
              exit <- 1
              print("neg")
            } else {
              if ((is.nan(u[aaa])) == TRUE) {
                exit <- 1
                print("nan2")
              } else {
                if ((is.nan(exit)) == TRUE) {
                  exit <- 1
                  print("nan3")
                }
              }
            }
          }
        }

        test <- 0

        if (exit == 0) {
          ifelse(((c(aaaa, bbbb, cccc, dddd, eeee) > 2)), test <- 1, print("proceed"))
        }

        if (test == 1) {
          print("exit2")

          uHuse <- (which.max(c(aaaa, bbbb, cccc, dddd, eeee)))

          if (uHuse == 1) {
            uH[aaa] <- ((1 / dt) - (u[aaa] + lambda[i - 1, aaa]))
          }
          if (uHuse == 2) {
            uH[aaa] <- ((1 / dt) - (u[aaa] + (lambda[i - 1, aaa] * pHa[aaa] *
              xHa) + vHa[aaa, 1]))
          }
          if (uHuse == 3) {
            uH[aaa] <- ((1 / dt) - (u[aaa] + uiHa[aaa] + nH[aaa]))
          }
          if (uHuse == 4) {
            uH[aaa] <- ((1 / dt) - (u[aaa] + uniHa[aaa] + nH[aaa] + w))
          }
          if (uHuse == 5) {
            uH[aaa] <- ((1 / dt) - (u[aaa] + (lambda[i - 1, aaa] * gHa) +
              rHa[aaa, 1]))
          }

          aaaa <- (u[aaa] + uH[aaa] + lambda[i - 1, aaa])

          bbbb <- (u[aaa] + uH[aaa] + vHa[aaa] + (lambda[i - 1, aaa] * pHa[aaa] *
            xHa))

          cccc <- (u[aaa] + uH[aaa] + uiHa[aaa] + nH[aaa])

          dddd <- (u[aaa] + uH[aaa] + uniHa[aaa] + nH[aaa] + w)

          eeee <- (u[aaa] + uH[aaa] + rHa[aaa, 1] + (lambda[i - 1, aaa] *
            gHa))

          if ((((aaaa < 0) == TRUE) || ((bbbb < 0) == TRUE) || ((cccc < 0) ==
            TRUE) || ((dddd < 0) == TRUE) || ((eeee < 0) == TRUE)) == TRUE) {
            exit <- 1
          }
        }
      }

      upopstore[i, ] <- upop[1:Mnage]
      uHpopstore[i, ] <- uHpop[1:Mnage]
      ustore[i, ] <- u[1:Mnage]
      uHstore[i, ] <- uH[1:Mnage]

      if (k == 1950) {
        demscale <- psize[i - 1, 1] / psz1900
        S[i - 1, ] <- S[i - 1, ] / demscale
        L[i - 1, ] <- L[i - 1, ] / demscale
        I[i - 1, ] <- I[i - 1, ] / demscale
        NI[i - 1, ] <- NI[i - 1, ] / demscale
        R[i - 1, ] <- R[i - 1, ] / demscale
      }

      S[i, 1:Mnage] <- S[i - 1, 1:Mnage] - (u[1:Mnage] + lambda[i - 1, 1:Mnage]) *
        S[i - 1, 1:Mnage] * dt - hiv[1:Mnage] * S[i - 1, 1:Mnage]

      L[i, 1:Mnage] <- L[i - 1, 1:Mnage] + lambda[i - 1, 1:Mnage] * (1 - p[1:Mnage]) *
        (S[i - 1, 1:Mnage] + g * R[i - 1, 1:Mnage]) * dt - (v[1:Mnage] +
        lambda[i - 1, 1:Mnage] * p[1:Mnage] * x + u[1:Mnage]) * L[
        i - 1,
        1:Mnage
      ] * dt - hiv[1:Mnage] * L[i - 1, 1:Mnage]
      new_infect[i, 1:Mnage] <- lambda[i - 1, 1:(Mnage)] * S[i - 1, 1:(Mnage)] *
        dt + lambda[i - 1, 1:(Mnage)] * (x * L[i - 1, 1:(Mnage)] + g * R[i -
          1, 1:(Mnage)]) * dt

      new_I_react[i, 1:Mnage] <- v[1:Mnage] * f[1:(Mnage)] * (L[i - 1, 1:(Mnage)]) *
        dt + r[1:Mnage] * h[1:(Mnage)] * R[i - 1, 1:(Mnage)] * dt
      new_NI_react[i, 1:Mnage] <- v[1:Mnage] * (1 - f[1:(Mnage)]) * L[
        i - 1,
        1:(Mnage)
      ] * dt + r[1:Mnage] * (1 - h[1:(Mnage)]) * R[i - 1, 1:(Mnage)] *
        dt

      new_actv_react[i, 1:Mnage] <- v[1:Mnage] * (L[i - 1, 1:Mnage]) * dt +
        r[1:Mnage] * R[i - 1, 1:Mnage] * dt
      new_actv_inf[i, 1:Mnage] <- lambda[i - 1, 1:Mnage] * p[1:Mnage] * S[i -
        1, 1:Mnage] * dt + lambda[i - 1, 1:Mnage] * p[1:Mnage] * x * (L[i -
        1, 1:Mnage]) * dt + lambda[i - 1, 1:Mnage] * p[1:Mnage] * g * R[i -
        1, 1:Mnage] * dt

      new_I[i, 1:Mnage] <- lambda[i - 1, 1:Mnage] * p[1:Mnage] * f[1:Mnage] *
        (S[i - 1, 1:Mnage] + g * R[i - 1, 1:(Mnage)]) * dt + (v[1:Mnage] +
          lambda[i - 1, 1:Mnage] * p[1:(Mnage)] * x) * f[1:(Mnage)] * L[i -
          1, 1:(Mnage)] * dt + r[1:Mnage] * h[1:(Mnage)] * R[i - 1, 1:(Mnage)] *
          dt + w * NI[i - 1, 1:(Mnage)] * dt
      new_I_noconv[i, 1:Mnage] <- lambda[i - 1, 1:Mnage] * p[1:Mnage] * f[1:Mnage] *
        (S[i - 1, 1:Mnage] + g * R[i - 1, 1:Mnage]) * dt + (v[1:Mnage] +
          lambda[i - 1, 1:Mnage] * p[1:Mnage] * x) * f[1:Mnage] * L[
          i - 1,
          1:Mnage
        ] * dt + r[1:Mnage] * h[1:Mnage] * R[i - 1, 1:Mnage] * dt

      new_NI[i, 1:Mnage] <- lambda[i - 1, 1:Mnage] * p[1:(Mnage)] * (1 - f[1:(Mnage)]) *
        (S[i - 1, 1:(Mnage)] + g * R[i - 1, 1:(Mnage)]) * dt + (v[1:Mnage] +
          lambda[i - 1, 1:Mnage] * p[1:(Mnage)] * x) * (1 - f[1:(Mnage)]) *
          L[i - 1, 1:(Mnage)] * dt + r[1:Mnage] * (1 - h[1:(Mnage)]) * R[i -
          1, 1:(Mnage)] * dt

      new_actv[i, 1:Mnage] <- lambda[i - 1, 1:Mnage] * p[1:Mnage] * S[
        i - 1,
        1:Mnage
      ] * dt + (v[1:Mnage] + lambda[i - 1, 1:Mnage] * p[1:Mnage] *
        x) * L[i - 1, 1:Mnage] * dt + (r[1:Mnage] + lambda[i - 1, 1:Mnage] *
        p[1:Mnage] * g) * R[i - 1, 1:Mnage] * dt

      new_actv_chk[i, 1:Mnage] <- new_actv_react[i, 1:Mnage] + new_actv_inf[
        i,
        1:Mnage
      ]

      new_notif[i, 1:Mnage] <- CDR[1:(Mnage)] * (new_I[i, 1:Mnage] + e * new_NI[
        i,
        1:Mnage
      ])
      R[i, 1:Mnage] <- R[i - 1, 1:Mnage] + n[1:Mnage] * (I[i - 1, 1:(Mnage)] +
        NI[i - 1, 1:(Mnage)]) * dt + CDR[1:Mnage] * CoT * (new_I[i, 1:Mnage] +
        e * new_NI[i, 1:Mnage]) - (r[1:Mnage] + g * lambda[i - 1, 1:Mnage] +
        u[1:(Mnage)]) * R[i - 1, 1:(Mnage)] * dt - hiv[1:Mnage] * R[
        i - 1,
        1:Mnage
      ]
      I[i, 1:Mnage] <- I[i - 1, 1:Mnage] + (1 - CDR[1:Mnage] * CoT) * (new_I[
        i,
        1:Mnage
      ]) - (n[1:Mnage] + u[1:(Mnage)] + ui[1:Mnage]) * I[
        i - 1,
        1:Mnage
      ] * dt - hiv[1:Mnage] * I[i - 1, 1:Mnage]
      NI[i, 1:Mnage] <- NI[i - 1, 1:Mnage] + ((1 - CDR[1:Mnage] * CoT * e) *
        new_NI[i, 1:Mnage]) - (n[1:Mnage] + u[1:Mnage] + uni[1:Mnage] + w) *
        NI[i - 1, 1:(Mnage)] * dt - hiv[1:Mnage] * NI[i - 1, 1:Mnage]

      SH[i, 1:Mnage] <- SH[i - 1, 1:Mnage] - (u[1:Mnage] + uH[1:Mnage] + lambda[i -
        1, 1:Mnage]) * SH[i - 1, 1:Mnage] * dt + hiv[1:Mnage] * S[
        i - 1,
        1:Mnage
      ]

      LH[i, 1:Mnage] <- LH[i - 1, 1:Mnage] + lambda[i - 1, 1:Mnage] * (1 - pHa[1:Mnage]) *
        (SH[i - 1, 1:Mnage] + gHa * RH[i - 1, 1:Mnage]) * dt - (vHa[1:Mnage] +
        lambda[i - 1, 1:Mnage] * pHa[1:Mnage] * xHa + u[1:Mnage] + uH[1:Mnage]) *
        LH[i - 1, 1:Mnage] * dt + hiv[1:Mnage] * L[i - 1, 1:Mnage]

      new_IH_react[i, 1:Mnage] <- vHa[1:Mnage] * fH[1:Mnage] * (LH[i - 1, 1:Mnage]) *
        dt + rHa[1:Mnage] * hH[1:Mnage] * RH[i - 1, 1:Mnage] * dt
      new_NIH_react[i, 1:Mnage] <- vHa[1:Mnage] * (1 - fH[1:Mnage]) * LH[i -
        1, 1:Mnage] * dt + rHa[1:Mnage] * (1 - hH[1:Mnage]) * RH[i - 1, 1:Mnage] *
        dt
      new_actvH_react[i, 1:Mnage] <- vHa[1:Mnage] * (LH[i - 1, 1:Mnage]) * dt +
        rHa[1:Mnage] * RH[i - 1, 1:Mnage] * dt
      new_IH[i, 1:Mnage] <- lambda[i - 1, 1:Mnage] * pHa[1:Mnage] * fH[1:Mnage] *
        (SH[i - 1, 1:Mnage] + gHa * RH[i - 1, 1:Mnage]) * dt + (vHa[1:Mnage] +
          lambda[i - 1, 1:Mnage] * pHa[1:Mnage] * xHa) * fH[1:Mnage] * LH[i -
          1, 1:Mnage] * dt + rHa[1:Mnage] * hH[1:Mnage] * RH[i - 1, 1:Mnage] *
          dt + w * NIH[i - 1, 1:Mnage] * dt
      new_IH_noconv[i, 1:Mnage] <- lambda[i - 1, 1:Mnage] * pHa[1:Mnage] * fH[1:Mnage] *
        (SH[i - 1, 1:Mnage] + gHa * RH[i - 1, 1:Mnage]) * dt + (vHa[1:Mnage] +
          lambda[i - 1, 1:Mnage] * pHa[1:Mnage] * xHa) * fH[1:Mnage] * LH[i -
          1, 1:Mnage] * dt + rHa[1:Mnage] * hH[1:Mnage] * RH[i - 1, 1:Mnage] *
          dt
      new_NIH[i, 1:Mnage] <- lambda[i - 1, 1:Mnage] * pHa[1:Mnage] * (1 - fH[1:Mnage]) *
        (SH[i - 1, 1:Mnage] + gHa * RH[i - 1, 1:Mnage]) * dt + (vHa[1:Mnage] +
          lambda[i - 1, 1:Mnage] * pHa[1:Mnage] * xHa) * (1 - fH[1:Mnage]) *
          LH[i - 1, 1:Mnage] * dt + rHa[1:Mnage] * (1 - hH[1:Mnage]) * RH[i -
          1, 1:Mnage] * dt
      new_actvH[i, 1:Mnage] <- lambda[i - 1, 1:Mnage] * pHa[1:Mnage] * SH[i -
        1, 1:Mnage] * dt + (vHa[1:Mnage] + lambda[i - 1, 1:Mnage] * pHa[1:Mnage] *
        xHa) * LH[i - 1, 1:Mnage] * dt + (rHa[1:Mnage] + lambda[i - 1, 1:Mnage] *
        pHa[1:Mnage] * gHa) * RH[i - 1, 1:Mnage] * dt

      new_notifH[i, 1:Mnage] <- CDRH[1:(Mnage)] * (new_IH[i, 1:Mnage] + e *
        new_NIH[i, 1:Mnage])

      RH[i, 1:Mnage] <- RH[i - 1, 1:Mnage] + nH[1:Mnage] * (IH[i - 1, 1:Mnage] +
        NIH[i - 1, 1:Mnage]) * dt + CDRH[1:(Mnage)] * CoT * (new_IH[i, 1:Mnage] +
        e * new_NIH[i, 1:Mnage]) - (rHa[1:Mnage] + gHa * lambda[i - 1, 1:Mnage] +
        u[1:Mnage] + uH[1:Mnage]) * RH[i - 1, 1:Mnage] * dt + hiv[1:Mnage] *
        R[i - 1, 1:Mnage]
      IH[i, 1:Mnage] <- IH[i - 1, 1:Mnage] + (1 - CDRH[1:(Mnage)] * CoT) * (new_IH[
        i,
        1:Mnage
      ]) - (nH[1:Mnage] + u[1:Mnage] + uH[1:Mnage] + uiHa[1:Mnage]) *
        IH[i - 1, 1:Mnage] * dt + hiv[1:Mnage] * I[i - 1, 1:Mnage]
      NIH[i, 1:Mnage] <- NIH[i - 1, 1:Mnage] + ((1 - CDRH[1:(Mnage)] * CoT *
        e) * new_NIH[i, 1:Mnage]) - (nH[1:Mnage] + u[1:Mnage] + uH[1:Mnage] +
        uniHa[1:Mnage] + w) * NIH[i - 1, 1:Mnage] * dt + hiv[1:Mnage] * NI[i -
        1, 1:Mnage]

      Sv[i, 1:Mnage] <- Sv[i - 1, 1:Mnage] - (u[1:Mnage] + ((1 - effI) * lambda[i -
        1, 1:Mnage])) * Sv[i - 1, 1:Mnage] * dt - hiv[1:Mnage] * Sv[
        i - 1,
        1:Mnage
      ]

      Lv[i, 1:Mnage] <- Lv[i - 1, 1:Mnage] + ((1 - effI) * lambda[i - 1, 1:Mnage]) *
        (1 - ((1 - effD) * p[1:Mnage])) * (Sv[i - 1, 1:Mnage] + g * Rv[i -
          1, 1:Mnage]) * dt - ((1 - effD) * v[1:Mnage] + ((1 - effI) * lambda[i -
        1, 1:Mnage]) * ((1 - effD) * p[1:Mnage]) * x + u[1:Mnage]) * Lv[i -
        1, 1:Mnage] * dt - hiv[1:Mnage] * Lv[i - 1, 1:Mnage]

      new_Iv[i, 1:Mnage] <- ((1 - effI) * lambda[i - 1, 1:Mnage]) * ((1 - effD) *
        p[1:Mnage]) * f[1:Mnage] * (Sv[i - 1, 1:Mnage] + g * Rv[i - 1, 1:(Mnage)]) *
        dt + ((1 - effD) * v[1:Mnage] + ((1 - effI) * lambda[i - 1, 1:Mnage]) *
          ((1 - effD) * p[1:Mnage]) * x) * f[1:(Mnage)] * Lv[i - 1, 1:Mnage] *
          dt + (1 - effD) * r[1:Mnage] * h[1:Mnage] * Rv[i - 1, 1:Mnage] *
          dt + w * NIv[i - 1, 1:Mnage] * dt
      new_Iv_noconv[i, 1:Mnage] <- ((1 - effI) * lambda[i - 1, 1:Mnage]) * ((1 -
        effD) * p[1:Mnage]) * f[1:Mnage] * (Sv[i - 1, 1:Mnage] + g * Rv[i -
        1, 1:Mnage]) * dt + ((1 - effD) * v[1:Mnage] + ((1 - effI) * lambda[i -
        1, 1:Mnage]) * ((1 - effD) * p[1:Mnage]) * x) * f[1:Mnage] * Lv[i -
        1, 1:Mnage] * dt + (1 - effD) * r[1:Mnage] * h[1:Mnage] * Rv[i -
        1, 1:Mnage] * dt
      new_NIv[i, 1:Mnage] <- ((1 - effI) * lambda[i - 1, 1:Mnage]) * ((1 - effD) *
        p[1:Mnage]) * (1 - f[1:Mnage]) * (Sv[i - 1, 1:Mnage] + g * Rv[i -
        1, 1:Mnage]) * dt + ((1 - effD) * v[1:Mnage] + ((1 - effI) * lambda[i -
        1, 1:Mnage]) * ((1 - effD) * p[1:Mnage]) * x) * (1 - f[1:Mnage]) *
        Lv[i - 1, 1:Mnage] * dt + (1 - effD) * r[1:Mnage] * (1 - h[1:Mnage]) *
        Rv[i - 1, 1:Mnage] * dt
      new_notifv[i, 1:Mnage] <- CDR[1:(Mnage)] * (new_Iv[i, 1:Mnage] + e * new_NIv[
        i,
        1:Mnage
      ])
      new_actvv[i, 1:Mnage] <- ((1 - effI) * lambda[i - 1, 1:Mnage]) * ((1 -
        effD) * p[1:Mnage]) * Sv[i - 1, 1:Mnage] * dt + (((1 - effD) * v[1:Mnage]) +
        ((1 - effI) * lambda[i - 1, 1:Mnage]) * ((1 - effD) * p[1:Mnage]) *
          x) * Lv[i - 1, 1:Mnage] * dt + (((1 - effD) * r[1:Mnage]) + ((1 -
        effI) * lambda[i - 1, 1:Mnage]) * ((1 - effD) * p[1:Mnage]) * g) *
        Rv[i - 1, 1:Mnage] * dt

      Rv[i, 1:Mnage] <- Rv[i - 1, 1:Mnage] + n[1:Mnage] * (Iv[i - 1, 1:(Mnage)] +
        NIv[i - 1, 1:(Mnage)]) * dt + CDR[1:Mnage] * CoT * (new_Iv[i, 1:Mnage] +
        e * new_NIv[i, 1:Mnage]) - ((1 - effD) * r[1:Mnage] + g * (1 - effI) *
        lambda[i - 1, 1:Mnage] + u[1:Mnage]) * Rv[i - 1, 1:(Mnage)] * dt -
        hiv[1:Mnage] * Rv[i - 1, 1:Mnage]

      Iv[i, 1:Mnage] <- Iv[i - 1, 1:Mnage] + (1 - CDR[1:Mnage] * CoT) * new_Iv[
        i,
        1:Mnage
      ] - (n[1:Mnage] + u[1:(Mnage)] + ui[1:Mnage]) * Iv[
        i - 1,
        1:Mnage
      ] * dt - hiv[1:Mnage] * Iv[i - 1, 1:Mnage]

      NIv[i, 1:Mnage] <- NIv[i - 1, 1:Mnage] + (1 - CDR[1:Mnage] * CoT * e) *
        new_NIv[i, 1:Mnage] - (n[1:Mnage] + u[1:Mnage] + uni[1:Mnage] + w) *
        NIv[i - 1, 1:(Mnage)] * dt - hiv[1:Mnage] * NIv[i - 1, 1:Mnage]

      SvH[i, 1:Mnage] <- SvH[i - 1, 1:Mnage] - (u[1:Mnage] + uH[1:Mnage] + ((1 -
        (effI * VEH)) * lambda[i - 1, 1:Mnage])) * SvH[i - 1, 1:Mnage] *
        dt + hiv[1:Mnage] * Sv[i - 1, 1:Mnage]

      LvH[i, 1:Mnage] <- LvH[i - 1, 1:Mnage] + ((1 - (effI * VEH)) * lambda[i -
        1, 1:Mnage]) * (1 - ((1 - (effD * VEH)) * pHa[1:Mnage])) * (SvH[i -
        1, 1:Mnage] + gHa * RvH[i - 1, 1:Mnage]) * dt - ((1 - (effD * VEH)) *
        vHa[1:Mnage] + ((1 - (effI * VEH)) * lambda[i - 1, 1:Mnage]) * ((1 -
          (effD * VEH)) * pHa[1:Mnage]) * xHa + u[1:Mnage] + uH[1:Mnage]) *
        LvH[i - 1, 1:Mnage] * dt + hiv[1:Mnage] * Lv[i - 1, 1:Mnage]

      new_IvH[i, 1:Mnage] <- ((1 - (effI * VEH)) * lambda[i - 1, 1:Mnage]) *
        ((1 - (effD * VEH)) * pHa[1:Mnage]) * fH[1:Mnage] * (SvH[i - 1, 1:Mnage] +
          gHa * RvH[i - 1, 1:(Mnage)]) * dt + ((1 - (effD * VEH)) * vHa[1:Mnage] +
          ((1 - (effI * VEH)) * lambda[i - 1, 1:Mnage]) * ((1 - (effD * VEH)) *
            pHa[1:Mnage]) * xHa) * fH[1:Mnage] * LvH[i - 1, 1:Mnage] * dt +
        (1 - (effD * VEH)) * rHa[1:Mnage] * hH[1:Mnage] * RvH[i - 1, 1:Mnage] *
          dt + w * NIvH[i - 1, 1:Mnage] * dt
      new_IvH_noconv[i, 1:Mnage] <- ((1 - (effI * VEH)) * lambda[i - 1, 1:Mnage]) *
        ((1 - (effD * VEH)) * pHa[1:Mnage]) * fH[1:Mnage] * (SvH[i - 1, 1:Mnage] +
          gHa * RvH[i - 1, 1:Mnage]) * dt + ((1 - (effD * VEH)) * vHa[1:Mnage] +
          ((1 - (effI * VEH)) * lambda[i - 1, 1:Mnage]) * ((1 - (effD * VEH)) *
            pHa[1:Mnage]) * xHa) * fH[1:Mnage] * LvH[i - 1, 1:Mnage] * dt +
        (1 - (effD * VEH)) * rHa[1:Mnage] * hH[1:Mnage] * RvH[i - 1, 1:Mnage] *
          dt
      new_NIvH[i, 1:Mnage] <- ((1 - (effI * VEH)) * lambda[i - 1, 1:Mnage]) *
        ((1 - (effD * VEH)) * pHa[1:Mnage]) * (1 - fH[1:Mnage]) * (SvH[i -
          1, 1:Mnage] + gHa * RvH[i - 1, 1:Mnage]) * dt + ((1 - (effD * VEH)) *
          vHa[1:Mnage] + ((1 - (effI * VEH)) * lambda[i - 1, 1:Mnage]) * ((1 -
            (effD * VEH)) * pHa[1:Mnage]) * xHa) * (1 - fH[1:Mnage]) * LvH[i -
          1, 1:Mnage] * dt + (1 - (effD * VEH)) * rHa[1:Mnage] * (1 - hH[1:Mnage]) *
          RvH[i - 1, 1:Mnage] * dt
      new_notifvH[i, 1:Mnage] <- CDRH[1:(Mnage)] * (new_IvH[i, 1:Mnage] + e *
        new_NIvH[i, 1:Mnage])
      new_actvvH[i, 1:Mnage] <- ((1 - (effI * VEH)) * lambda[i - 1, 1:Mnage]) *
        ((1 - (effD * VEH)) * pHa[1:Mnage]) * SvH[i - 1, 1:Mnage] * dt +
        (((1 - (effD * VEH)) * vHa[1:Mnage]) + ((1 - (effI * VEH)) * lambda[i -
          1, 1:Mnage]) * ((1 - (effD * VEH)) * pHa[1:Mnage]) * xHa) * LvH[i -
          1, 1:Mnage] * dt + (((1 - (effD * VEH)) * rHa[1:Mnage]) + ((1 -
          (effI * VEH)) * lambda[i - 1, 1:Mnage]) * ((1 - (effD * VEH)) * pHa[1:Mnage]) *
          gHa) * RvH[i - 1, 1:Mnage] * dt

      RvH[i, 1:Mnage] <- RvH[i - 1, 1:Mnage] + nH[1:Mnage] * (IvH[i - 1, 1:Mnage] +
        NIvH[i - 1, 1:Mnage]) * dt + CDRH[1:(Mnage)] * CoT * (new_IvH[
        i,
        1:Mnage
      ] + e * new_NIvH[i, 1:Mnage]) - ((1 - (effD * VEH)) * rHa[1:Mnage] +
        gHa * (1 - (effI * VEH)) * lambda[i - 1, 1:Mnage] + u[1:Mnage] +
        uH[1:Mnage]) * RvH[i - 1, 1:Mnage] * dt + hiv[1:Mnage] * Rv[
        i - 1,
        1:Mnage
      ]

      IvH[i, 1:Mnage] <- IvH[i - 1, 1:Mnage] + (1 - CDRH[1:(Mnage)] * CoT) *
        new_IvH[i, 1:Mnage] - (nH[1:Mnage] + u[1:Mnage] + uH[1:Mnage] + uiHa[1:Mnage]) *
        IvH[i - 1, 1:Mnage] * dt + hiv[1:Mnage] * Iv[i - 1, 1:Mnage]

      NIvH[i, 1:Mnage] <- NIvH[i - 1, 1:Mnage] + (1 - CDRH[1:(Mnage)] * CoT *
        e) * new_NIvH[i, 1:Mnage] - (nH[1:Mnage] + u[1:Mnage] + uH[1:Mnage] +
        uniHa[1:Mnage] + w) * NIvH[i - 1, 1:Mnage] * dt + hiv[1:Mnage] *
        NIv[i - 1, 1:Mnage]

      if (vaccine == 0) {
        VX[i, 1:3] <- 0
      }

      if (vaccine == 2) {
        VX[i, 1] <- (sum(thetaV2a[i, ] * (S[i, ] + L[i, ] + R[i, ])) + sum(thetaV2m[i, ] * (S[i, ] + L[i, ] + R[i, ])))
        VX[i, 2] <- sum(thetaV2a[i, ] * (S[i, ] + L[i, ] + R[i, ]))
        VX[i, 3] <- sum(thetaV2m[i, ] * (S[i, ] + L[i, ] + R[i, ]))
      }
      if (vaccine == 3) {
        VX[i, 1] <- (sum(thetaV2a[i, ] * (S[i, ] + L[i, ] + R[i, ])) + sum(thetaV2m[i, ] * (S[i, ] + L[i, ] + R[i, ])))
        VX[i, 2] <- sum(thetaV2a[i, ] * (S[i, ] + L[i, ] + R[i, ]))
        VX[i, 3] <- sum(thetaV2m[i, ] * (S[i, ] + L[i, ] + R[i, ]))
      }
      if (vaccine == 4) {
        VX[i, 1] <- (sum(thetaV2a[i, ] * (S[i, ] + L[i, ] + R[i, ])) + sum(thetaV2m[i, ] * (S[i, ] + L[i, ] + R[i, ])))
        VX[i, 2] <- sum(thetaV2a[i, ] * (S[i, ] + L[i, ] + R[i, ]))
        VX[i, 3] <- sum(thetaV2m[i, ] * (S[i, ] + L[i, ] + R[i, ]))
      }
      if (vaccine == 5) {
        VX[i, 1] <- (sum(thetaV4a[i, ] * (S[i, ] + L[i, ] + R[i, ])) + sum(thetaV4m[i, ] * (S[i, ] + L[i, ] + R[i, ])))
        VX[i, 2] <- sum(thetaV4a[i, ] * (S[i, ] + L[i, ] + R[i, ]))
        VX[i, 3] <- sum(thetaV4m[i, ] * (S[i, ] + L[i, ] + R[i, ]))
      }
      if (vaccine == 6) {
        VX[i, 1] <- (sum(thetaV4a[i, ] * (S[i, ] + L[i, ] + R[i, ])) + sum(thetaV4m[i, ] * (S[i, ] + L[i, ] + R[i, ])))
        VX[i, 2] <- sum(thetaV4a[i, ] * (S[i, ] + L[i, ] + R[i, ]))
        VX[i, 3] <- sum(thetaV4m[i, ] * (S[i, ] + L[i, ] + R[i, ]))
      }
      if (vaccine == 7) {
        VX[i, 1] <- (sum(thetaV4a[i, ] * (S[i, ] + L[i, ] + R[i, ])) + sum(thetaV4m[i, ] * (S[i, ] + L[i, ] + R[i, ])))
        VX[i, 2] <- sum(thetaV4a[i, ] * (S[i, ] + L[i, ] + R[i, ]))
        VX[i, 3] <- sum(thetaV4m[i, ] * (S[i, ] + L[i, ] + R[i, ]))
      }
      if (vaccine == 8) {
        VX[i, 1] <- (sum(thetaV2a[i, ] * (S[i, ] + L[i, ] + R[i, ])) + sum(thetaV2m[i, ] * (S[i, ] + L[i, ] + R[i, ])))
        VX[i, 2] <- sum(thetaV2a[i, ] * (S[i, ] + L[i, ] + R[i, ]))
        VX[i, 3] <- sum(thetaV2m[i, ] * (S[i, ] + L[i, ] + R[i, ]))
      }
      if (vaccine == 9) {
        VX[i, 1] <- (sum(thetaV4a[i, ] * (S[i, ] + L[i, ] + R[i, ])) + sum(thetaV4m[i, ] * (S[i, ] + L[i, ] + R[i, ])))
        VX[i, 2] <- sum(thetaV4a[i, ] * (S[i, ] + L[i, ] + R[i, ]))
        VX[i, 3] <- sum(thetaV4m[i, ] * (S[i, ] + L[i, ] + R[i, ]))
      }

      S2 <- S[i, ] + Sv[i, ] * (d[i, ] * (1 - thetaS[i, ])) - thetaS[i, ] *
        S[i, ]
      L2 <- L[i, ] + Lv[i, ] * (d[i, ] * (1 - thetaL[i, ])) - thetaL[i, ] *
        L[i, ]
      R2 <- R[i, ] + Rv[i, ] * (d[i, ] * (1 - thetaR[i, ])) - thetaR[i, ] *
        R[i, ]
      I2 <- I[i, ] + Iv[i, ] * d[i, ]
      NI2 <- NI[i, ] + NIv[i, ] * d[i, ]

      Sv2 <- Sv[i, ] - Sv[i, ] * (d[i, ] * (1 - thetaS[i, ])) + thetaS[i, ] *
        S[i, ]
      Lv2 <- Lv[i, ] - Lv[i, ] * (d[i, ] * (1 - thetaL[i, ])) + thetaL[i, ] *
        L[i, ]
      Rv2 <- Rv[i, ] - Rv[i, ] * (d[i, ] * (1 - thetaR[i, ])) + thetaR[i, ] *
        R[i, ]
      Iv2 <- Iv[i, ] - Iv[i, ] * (d[i, ])
      NIv2 <- NIv[i, ] - NIv[i, ] * (d[i, ])

      SH2 <- SH[i, ] + SvH[i, ] * (dH[i, ] * (1 - thetaSH[i, ])) - thetaSH[i, ] * SH[i, ]
      LH2 <- LH[i, ] + LvH[i, ] * (dH[i, ] * (1 - thetaLH[i, ])) - thetaLH[i, ] * LH[i, ]
      RH2 <- RH[i, ] + RvH[i, ] * (dH[i, ] * (1 - thetaRH[i, ])) - thetaRH[i, ] * RH[i, ]
      IH2 <- IH[i, ] + IvH[i, ] * dH[i, ]
      NIH2 <- NIH[i, ] + NIvH[i, ] * dH[i, ]

      SvH2 <- SvH[i, ] - SvH[i, ] * (dH[i, ] * (1 - thetaSH[i, ])) + thetaSH[i, ] * SH[i, ]
      LvH2 <- LvH[i, ] - LvH[i, ] * (dH[i, ] * (1 - thetaLH[i, ])) + thetaLH[i, ] * LH[i, ]
      RvH2 <- RvH[i, ] - RvH[i, ] * (dH[i, ] * (1 - thetaRH[i, ])) + thetaRH[i, ] * RH[i, ]
      IvH2 <- IvH[i, ] - IvH[i, ] * (dH[i, ])
      NIvH2 <- NIvH[i, ] - NIvH[i, ] * (dH[i, ])

      if (vaccine == 1) {
        num_vac[i, ] <- (S[i, ] + L[i, ] + R[i, ] + Sv[i, ] + Lv[i, ] + Rv[i, ] + SH[i, ] + LH[i, ] + RH[i, ] + SvH[i, ] + LvH[i, ] + RvH[i, ]) * thetaV2[i, ]
      }
      if (vaccine == 2) {
        num_vac[i, ] <- (S[i, ] + L[i, ] + R[i, ] + Sv[i, ] + Lv[i, ] + Rv[i, ] + SH[i, ] + LH[i, ] + RH[i, ] + SvH[i, ] + LvH[i, ] + RvH[i, ]) * thetaV2[i, ]
      }
      if (vaccine == 3) {
        num_vac[i, ] <- (S[i, ] + L[i, ] + R[i, ] + Sv[i, ] + Lv[i, ] + Rv[i, ] + SH[i, ] + LH[i, ] + RH[i, ] + SvH[i, ] + LvH[i, ] + RvH[i, ]) * thetaV2[i, ]
      }

      S[i, ] <- S2
      L[i, ] <- L2
      R[i, ] <- R2
      I[i, ] <- I2
      NI[i, ] <- NI2
      Sv[i, ] <- Sv2
      Lv[i, ] <- Lv2
      Rv[i, ] <- Rv2
      Iv[i, ] <- Iv2
      NIv[i, ] <- NIv2

      SH[i, ] <- SH2
      LH[i, ] <- LH2
      RH[i, ] <- RH2
      IH[i, ] <- IH2
      NIH[i, ] <- NIH2
      SvH[i, ] <- SvH2
      LvH[i, ] <- LvH2
      RvH[i, ] <- RvH2
      IvH[i, ] <- IvH2
      NIvH[i, ] <- NIvH2

      psize[i] <- sum(S[i, ], L[i, ], R[i, ], I[i, ], NI[i, ], Sv[i, ], Lv[i, ], Rv[i, ], Iv[i, ], NIv[i, ])
      psizeH[i] <- sum(SH[i, ], LH[i, ], RH[i, ], IH[i, ], NIH[i, ], SvH[i, ], LvH[i, ], RvH[i, ], IvH[i, ], NIvH[i, ])
      psizeALL[i] <- sum(
        S[i, ], L[i, ], R[i, ], I[i, ], NI[i, ], Sv[i, ],
        Lv[i, ], Rv[i, ], Iv[i, ], NIv[i, ], SH[i, ], LH[i, ], RH[i, ], IH[i, ], NIH[i, ], SvH[i, ], LvH[i, ], RvH[i, ], IvH[i, ], NIvH[i, ]
      )

      psize014[i] <- sum(S[i, 1:15], L[i, 1:15], R[i, 1:15], I[i, 1:15], NI[
        i,
        1:15
      ], Sv[i, 1:15], Lv[i, 1:15], Iv[i, 1:15], NIv[i, 1:15], Rv[
        i,
        1:15
      ])
      psize1554[i] <- sum(
        S[i, 16:55], L[i, 16:55], R[i, 16:55], I[i, 16:55],
        NI[i, 16:55], Sv[i, 16:55], Lv[i, 16:55], Iv[i, 16:55], NIv[i, 16:55],
        Rv[i, 16:55]
      )
      psize5564[i] <- sum(
        S[i, 56:65], L[i, 56:65], R[i, 56:65], I[i, 56:65],
        NI[i, 56:65], Sv[i, 56:65], Lv[i, 56:65], Iv[i, 56:65], NIv[i, 56:65],
        Rv[i, 56:65]
      )
      psize65plus[i] <- sum(
        S[i, 66:Mnage], L[i, 66:Mnage], R[i, 66:Mnage],
        I[i, 66:Mnage], NI[i, 66:Mnage], Sv[i, 66:Mnage], Lv[i, 66:Mnage],
        Rv[i, 66:Mnage], Iv[i, 66:Mnage], NIv[i, 66:Mnage]
      )
      psize55plus[i] <- sum(
        S[i, 56:Mnage], L[i, 56:Mnage], R[i, 56:Mnage],
        I[i, 56:Mnage], NI[i, 56:Mnage], Sv[i, 56:Mnage], Lv[i, 56:Mnage],
        Rv[i, 56:Mnage], Iv[i, 56:Mnage], NIv[i, 56:Mnage]
      )
      psize55minus[i] <- sum(
        S[i, 1:55], L[i, 1:55], R[i, 1:55], I[i, 1:55],
        NI[i, 1:55], Sv[i, 1:55], Lv[i, 1:55], Rv[i, 1:55], Iv[i, 1:55],
        NIv[i, 1:55]
      )

      psizeVX[i] <- sum(Sv[i, ], Lv[i, ], Rv[i, ], Iv[i, ], NIv[i, ])
      psizeVX9[i] <- sum(Sv[i, 10], Lv[i, 10], Iv[i, 10], NIv[i, 10], Rv[
        i,
        10
      ])
      psizeVX9100[i] <- sum(Sv[i, 10:101], Lv[i, 10:101], Iv[i, 10:101], NIv[
        i,
        10:101
      ], Rv[i, 10:101])
      psizeNOVX[i] <- sum(S[i, ], L[i, ], R[i, ], I[i, ], NI[i, ])
      psizeNOVX9[i] <- sum(S[i, 10], L[i, 10], I[i, 10], NI[i, 10], R[i, 10])
      psizeNOVX9100[i] <- sum(S[i, 10:101], L[i, 10:101], I[i, 10:101], NI[
        i,
        10:101
      ], R[i, 10:101])

      psize1559[i] <- sum(
        S[i, 16:60], L[i, 16:60], R[i, 16:60], I[i, 16:60],
        NI[i, 16:60], Sv[i, 16:60], Lv[i, 16:60], Rv[i, 16:60], Iv[i, 16:60],
        NIv[i, 16:60]
      )
      psize1529[i] <- sum(
        S[i, 16:30], L[i, 16:30], R[i, 16:30], I[i, 16:30],
        NI[i, 16:30], Sv[i, 16:30], Lv[i, 16:30], Rv[i, 16:30], Iv[i, 16:30],
        NIv[i, 16:30]
      )
      psize3044[i] <- sum(
        S[i, 31:45], L[i, 31:45], R[i, 31:45], I[i, 31:45],
        NI[i, 31:45], Sv[i, 31:45], Lv[i, 31:45], Rv[i, 31:45], Iv[i, 31:45],
        NIv[i, 31:45]
      )
      psize4559[i] <- sum(
        S[i, 46:60], L[i, 46:60], R[i, 46:60], I[i, 46:60],
        NI[i, 46:60], Sv[i, 46:60], Lv[i, 46:60], Rv[i, 46:60], Iv[i, 46:60],
        NIv[i, 46:60]
      )
      psize60plus[i] <- sum(
        S[i, 61:Mnage], L[i, 61:Mnage], R[i, 61:Mnage],
        I[i, 61:Mnage], NI[i, 61:Mnage], Sv[i, 61:Mnage], Lv[i, 61:Mnage],
        Rv[i, 61:Mnage], Iv[i, 61:Mnage], NIv[i, 61:Mnage]
      )
      psize5574[i] <- sum(
        S[i, 56:75], L[i, 56:75], R[i, 56:75], I[i, 56:75],
        NI[i, 56:75], Sv[i, 56:75], Lv[i, 56:75], Rv[i, 56:75], Iv[i, 56:75],
        NIv[i, 56:75]
      )
      psize75plus[i] <- sum(
        S[i, 76:Mnage], L[i, 76:Mnage], R[i, 76:Mnage],
        I[i, 76:Mnage], NI[i, 76:Mnage], Sv[i, 76:Mnage], Lv[i, 76:Mnage],
        Rv[i, 76:Mnage], Iv[i, 76:Mnage], NIv[i, 76:Mnage]
      )
      psize1524[i] <- sum(
        S[i, 16:25], L[i, 16:25], R[i, 16:25], I[i, 16:25],
        NI[i, 16:25], Sv[i, 16:25], Lv[i, 16:25], Rv[i, 16:25], Iv[i, 16:25],
        NIv[i, 16:25]
      )
      psize2554[i] <- sum(
        S[i, 26:55], L[i, 26:55], R[i, 26:55], I[i, 26:55],
        NI[i, 26:55], Sv[i, 26:55], Lv[i, 26:55], Rv[i, 26:55], Iv[i, 26:55],
        NIv[i, 26:55]
      )
      psize15plus[i] <- sum(
        S[i, 16:Mnage], L[i, 16:Mnage], R[i, 16:Mnage],
        I[i, 16:Mnage], NI[i, 16:Mnage], Sv[i, 16:Mnage], Lv[i, 16:Mnage],
        Rv[i, 16:Mnage], Iv[i, 16:Mnage], NIv[i, 16:Mnage]
      )

      psize0509[i] <- sum(S[i, 6:10], L[i, 6:10], R[i, 6:10], I[i, 6:10], NI[
        i,
        6:10
      ], Sv[i, 6:10], Lv[i, 6:10], Rv[i, 6:10], Iv[i, 6:10], NIv[
        i,
        6:10
      ])
      psize1019[i] <- sum(
        S[i, 11:20], L[i, 11:20], R[i, 11:20], I[i, 11:20],
        NI[i, 11:20], Sv[i, 11:20], Lv[i, 11:20], Rv[i, 11:20], Iv[i, 11:20],
        NIv[i, 11:20]
      )
      psize2029[i] <- sum(
        S[i, 21:30], L[i, 21:30], R[i, 21:30], I[i, 21:30],
        NI[i, 21:30], Sv[i, 21:30], Lv[i, 21:30], Rv[i, 21:30], Iv[i, 21:30],
        NIv[i, 21:30]
      )
      psize3039[i] <- sum(
        S[i, 31:40], L[i, 31:40], R[i, 31:40], I[i, 31:40],
        NI[i, 31:40], Sv[i, 31:40], Lv[i, 31:40], Rv[i, 31:40], Iv[i, 31:40],
        NIv[i, 31:40]
      )
      psize4049[i] <- sum(
        S[i, 41:50], L[i, 41:50], R[i, 41:50], I[i, 41:50],
        NI[i, 41:50], Sv[i, 41:50], Lv[i, 41:50], Rv[i, 41:50], Iv[i, 41:50],
        NIv[i, 41:50]
      )
      psize5059[i] <- sum(
        S[i, 51:60], L[i, 51:60], R[i, 51:60], I[i, 51:60],
        NI[i, 51:60], Sv[i, 51:60], Lv[i, 51:60], Rv[i, 51:60], Iv[i, 51:60],
        NIv[i, 51:60]
      )
      psize6069[i] <- sum(
        S[i, 61:70], L[i, 61:70], R[i, 61:70], I[i, 61:70],
        NI[i, 61:70], Sv[i, 61:70], Lv[i, 61:70], Rv[i, 61:70], Iv[i, 61:70],
        NIv[i, 61:70]
      )
      psize70plus[i] <- sum(
        S[i, 71:Mnage], L[i, 71:Mnage], R[i, 71:Mnage],
        I[i, 71:Mnage], NI[i, 71:Mnage], Sv[i, 71:Mnage], Lv[i, 71:Mnage],
        Rv[i, 71:Mnage], Iv[i, 71:Mnage], NIv[i, 71:Mnage]
      )

      psize1549[i] <- psize1559[i] - psize5059[i]

      agepopall[i, 1:Mnage] <- S[i, 1:Mnage] + L[i, 1:Mnage] + R[i, 1:Mnage] +
        I[i, 1:Mnage] + NI[i, 1:Mnage] + Sv[i, 1:Mnage] + Lv[i, 1:Mnage] +
        Rv[i, 1:Mnage] + Iv[i, 1:Mnage] + NIv[i, 1:Mnage] + SH[i, 1:Mnage] +
        LH[i, 1:Mnage] + RH[i, 1:Mnage] + IH[i, 1:Mnage] + NIH[i, 1:Mnage] +
        SvH[i, 1:Mnage] + LvH[i, 1:Mnage] + RvH[i, 1:Mnage] + IvH[i, 1:Mnage] +
        NIvH[i, 1:Mnage]

      psize04H[i] <- sum(SH[i, 1:5], LH[i, 1:5], RH[i, 1:5], IH[i, 1:5], NIH[
        i,
        1:5
      ], SvH[i, 1:5], LvH[i, 1:5], RvH[i, 1:5], IvH[i, 1:5], NIvH[
        i,
        1:5
      ])
      agepopHIV[i, 1:Mnage] <- SH[i, 1:Mnage] + LH[i, 1:Mnage] + RH[i, 1:Mnage] +
        IH[i, 1:Mnage] + NIH[i, 1:Mnage] + SvH[i, 1:Mnage] + LvH[i, 1:Mnage] +
        RvH[i, 1:Mnage] + IvH[i, 1:Mnage] + NIvH[i, 1:Mnage]

      psize014H[i] <- sum(
        SH[i, 1:15], LH[i, 1:15], RH[i, 1:15], IH[i, 1:15],
        NIH[i, 1:15], SvH[i, 1:15], LvH[i, 1:15], IvH[i, 1:15], NIvH[i, 1:15],
        RvH[i, 1:15]
      )
      psize1554H[i] <- sum(
        SH[i, 16:55], LH[i, 16:55], RH[i, 16:55], IH[
          i,
          16:55
        ], NIH[i, 16:55], SvH[i, 16:55], LvH[i, 16:55], IvH[i, 16:55],
        NIvH[i, 16:55], RvH[i, 16:55]
      )
      psize5564H[i] <- sum(
        SH[i, 56:65], LH[i, 56:65], RH[i, 56:65], IH[
          i,
          56:65
        ], NIH[i, 56:65], SvH[i, 56:65], LvH[i, 56:65], IvH[i, 56:65],
        NIvH[i, 56:65], RvH[i, 56:65]
      )
      psize65plusH[i] <- sum(
        SH[i, 66:Mnage], LH[i, 66:Mnage], RH[i, 66:Mnage],
        IH[i, 66:Mnage], NIH[i, 66:Mnage], SvH[i, 66:Mnage], LvH[i, 66:Mnage],
        RvH[i, 66:Mnage], IvH[i, 66:Mnage], NIvH[i, 66:Mnage]
      )

      psizeVXH[i] <- sum(SvH[i, ], LvH[i, ], RvH[i, ], IvH[i, ], NIvH[i, ])
      psizeVX9H[i] <- sum(
        SvH[i, 10], LvH[i, 10], IvH[i, 10], NIvH[i, 10],
        RvH[i, 10]
      )
      psizeVX9100H[i] <- sum(
        SvH[i, 10:101], LvH[i, 10:101], IvH[i, 10:101],
        NIvH[i, 10:101], RvH[i, 10:101]
      )

      psize1559H[i] <- sum(
        SH[i, 16:60], LH[i, 16:60], RH[i, 16:60], IH[
          i,
          16:60
        ], NIH[i, 16:60], SvH[i, 16:60], LvH[i, 16:60], RvH[i, 16:60],
        IvH[i, 16:60], NIvH[i, 16:60]
      )
      psize60plusH[i] <- sum(
        SH[i, 61:Mnage], LH[i, 61:Mnage], RH[i, 61:Mnage],
        IH[i, 61:Mnage], NIH[i, 61:Mnage], SvH[i, 61:Mnage], LvH[i, 61:Mnage],
        RvH[i, 61:Mnage], IvH[i, 61:Mnage], NIvH[i, 61:Mnage]
      )
      psize15plusH[i] <- sum(
        SH[i, 16:Mnage], LH[i, 16:Mnage], RH[i, 16:Mnage],
        IH[i, 16:Mnage], NIH[i, 16:Mnage], SvH[i, 16:Mnage], LvH[i, 16:Mnage],
        RvH[i, 16:Mnage], IvH[i, 16:Mnage], NIvH[i, 16:Mnage]
      )

      psize0509H[i] <- sum(
        SH[i, 6:10], LH[i, 6:10], RH[i, 6:10], IH[i, 6:10],
        NIH[i, 6:10], SvH[i, 6:10], LvH[i, 6:10], RvH[i, 6:10], IvH[i, 6:10],
        NIvH[i, 6:10]
      )
      psize1019H[i] <- sum(
        SH[i, 11:20], LH[i, 11:20], RH[i, 11:20], IH[
          i,
          11:20
        ], NIH[i, 11:20], SvH[i, 11:20], LvH[i, 11:20], RvH[i, 11:20],
        IvH[i, 11:20], NIvH[i, 11:20]
      )
      psize2029H[i] <- sum(
        SH[i, 21:30], LH[i, 21:30], RH[i, 21:30], IH[
          i,
          21:30
        ], NIH[i, 21:30], SvH[i, 21:30], LvH[i, 21:30], RvH[i, 21:30],
        IvH[i, 21:30], NIvH[i, 21:30]
      )
      psize3039H[i] <- sum(
        SH[i, 31:40], LH[i, 31:40], RH[i, 31:40], IH[
          i,
          31:40
        ], NIH[i, 31:40], SvH[i, 31:40], LvH[i, 31:40], RvH[i, 31:40],
        IvH[i, 31:40], NIvH[i, 31:40]
      )
      psize4049H[i] <- sum(
        SH[i, 41:50], LH[i, 41:50], RH[i, 41:50], IH[
          i,
          41:50
        ], NIH[i, 41:50], SvH[i, 41:50], LvH[i, 41:50], RvH[i, 41:50],
        IvH[i, 41:50], NIvH[i, 41:50]
      )
      psize5059H[i] <- sum(
        SH[i, 51:60], LH[i, 51:60], RH[i, 51:60], IH[
          i,
          51:60
        ], NIH[i, 51:60], SvH[i, 51:60], LvH[i, 51:60], RvH[i, 51:60],
        IvH[i, 51:60], NIvH[i, 51:60]
      )
      psize6069H[i] <- sum(
        SH[i, 61:70], LH[i, 61:70], RH[i, 61:70], IH[
          i,
          61:70
        ], NIH[i, 61:70], SvH[i, 61:70], LvH[i, 61:70], RvH[i, 61:70],
        IvH[i, 61:70], NIvH[i, 61:70]
      )
      psize70plusH[i] <- sum(
        SH[i, 71:Mnage], LH[i, 71:Mnage], RH[i, 71:Mnage],
        IH[i, 71:Mnage], NIH[i, 71:Mnage], SvH[i, 71:Mnage], LvH[i, 71:Mnage],
        RvH[i, 71:Mnage], IvH[i, 71:Mnage], NIvH[i, 71:Mnage]
      )

      psize1549H[i] <- psize1559H[i] - psize5059H[i]

      artsize[i, 1] <- psize014H[i] * artyrc
      artsize[i, 2] <- psize15plusH[i] * artyr

      psizematrix[i, 1] <- sum(
        S[i, 1:5], L[i, 1:5], R[i, 1:5], I[i, 1:5],
        NI[i, 1:5], Sv[i, 1:5], Lv[i, 1:5], Rv[i, 1:5], Iv[i, 1:5], NIv[
          i,
          1:5
        ], SH[i, 1:5], LH[i, 1:5], RH[i, 1:5], IH[i, 1:5], NIH[i, 1:5],
        SvH[i, 1:5], LvH[i, 1:5], RvH[i, 1:5], IvH[i, 1:5], NIvH[i, 1:5]
      )
      psizematrix[i, 2] <- sum(
        S[i, 6:10], L[i, 6:10], R[i, 6:10], I[i, 6:10],
        NI[i, 6:10], Sv[i, 6:10], Lv[i, 6:10], Rv[i, 6:10], Iv[i, 6:10],
        NIv[i, 6:10], SH[i, 6:10], LH[i, 6:10], RH[i, 6:10], IH[i, 6:10],
        NIH[i, 6:10], SvH[i, 6:10], LvH[i, 6:10], RvH[i, 6:10], IvH[i, 6:10],
        NIvH[i, 6:10]
      )
      psizematrix[i, 3] <- sum(
        S[i, 11:15], L[i, 11:15], R[i, 11:15], I[
          i,
          11:15
        ], NI[i, 11:15], Sv[i, 11:15], Lv[i, 11:15], Rv[i, 11:15], Iv[
          i,
          11:15
        ], NIv[i, 11:15], SH[i, 11:15], LH[i, 11:15], RH[i, 11:15],
        IH[i, 11:15], NIH[i, 11:15], SvH[i, 11:15], LvH[i, 11:15], RvH[
          i,
          11:15
        ], IvH[i, 11:15], NIvH[i, 11:15]
      )
      psizematrix[i, 4] <- sum(
        S[i, 16:Mnage], L[i, 16:Mnage], R[i, 16:Mnage],
        I[i, 16:Mnage], NI[i, 16:Mnage], Sv[i, 16:Mnage], Lv[i, 16:Mnage],
        Rv[i, 16:Mnage], Iv[i, 16:Mnage], NIv[i, 16:Mnage], SH[i, 16:Mnage],
        LH[i, 16:Mnage], RH[i, 16:Mnage], IH[i, 16:Mnage], NIH[i, 16:Mnage],
        SvH[i, 16:Mnage], LvH[i, 16:Mnage], RvH[i, 16:Mnage], IvH[i, 16:Mnage],
        NIvH[i, 16:Mnage]
      )

      Imatrix[i, 1] <- sum(I[i, 1:5], Iv[i, 1:5], IH[i, 1:5], IvH[i, 1:5])
      Imatrix[i, 2] <- sum(I[i, 6:10], Iv[i, 6:10], IH[i, 6:10], IvH[i, 6:10])
      Imatrix[i, 3] <- sum(I[i, 11:15], Iv[i, 11:15], IH[i, 11:15], IvH[
        i,
        11:15
      ])
      Imatrix[i, 4] <- sum(
        I[i, 16:Mnage], Iv[i, 16:Mnage], IH[i, 16:Mnage],
        IvH[i, 16:Mnage]
      )

      TBDeaths[i, 1:Mnage] <- dt * ((ui[1:Mnage] * I[i - 1, 1:Mnage]) + (uni[1:Mnage] *
        NI[i - 1, 1:Mnage]) + (ui[1:Mnage] * Iv[i - 1, 1:Mnage]) + (uni[1:Mnage] *
        NIv[i - 1, 1:Mnage]) + (uiHa[1:Mnage] * IH[i - 1, 1:Mnage]) + (uniHa[1:Mnage] *
        NIH[i - 1, 1:Mnage]) + (uiHa[1:Mnage] * IvH[i - 1, 1:Mnage]) + (uniHa[1:Mnage] *
        NIvH[i - 1, 1:Mnage]))
      TBDeathsN[i, 1:Mnage] <- dt * ((ui[1:Mnage] * I[i - 1, 1:Mnage]) + (uni[1:Mnage] *
        NI[i - 1, 1:Mnage]) + (ui[1:Mnage] * Iv[i - 1, 1:Mnage]) + (uni[1:Mnage] *
        NIv[i - 1, 1:Mnage]))
      TBDeathsH[i, 1:Mnage] <- dt * ((uiHa[1:Mnage] * IH[i - 1, 1:Mnage]) +
        (uniHa[1:Mnage] * NIH[i - 1, 1:Mnage]) + (uiHa[1:Mnage] * IvH[i -
          1, 1:Mnage]) + (uniHa[1:Mnage] * NIvH[i - 1, 1:Mnage]))

      ADeaths[i, ] <- dt * (u * S[i - 1, ] + u * L[i - 1, ] + (u + ui) * I[i -
        1, ] + (u + uni) * NI[i - 1, ] + u * R[i - 1, ] + u * Sv[i - 1, ] +
        u * Lv[i - 1, ] + u * Rv[i - 1, ] + (u + ui) * Iv[i - 1, ] + (u +
          uni) * NIv[i - 1, ] + (u + uH) * SH[i - 1, ] + (u + uH) * LH[i -
          1, ] + (u + uH + uiHa) * IH[i - 1, ] + (u + uH + uniHa) * NIH[i -
          1, ] + (u + uH) * RH[i - 1, ] + (u + uH) * SvH[i - 1, ] + (u + uH) *
          LvH[i - 1, ] + (u + uH) * RvH[i - 1, ] + (u + uH + uiHa) * IvH[i -
          1, ] + (u + uH + uniHa) * NIvH[i - 1, ])
      ADeathsN[i, ] <- dt * (u * S[i - 1, ] + u * L[i - 1, ] + (u + ui) * I[i -
        1, ] + (u + uni) * NI[i - 1, ] + u * R[i - 1, ] + u * Sv[i - 1, ] +
        u * Lv[i - 1, ] + u * Rv[i - 1, ] + (u + ui) * Iv[i - 1, ] + (u +
          uni) * NIv[i - 1, ])
      ADeathsH[i, ] <- dt * ((u + uH) * SH[i - 1, ] + (u + uH) * LH[i - 1, ] +
        (u + uH + uiHa) * IH[i - 1, ] + (u + uH + uniHa) * NIH[i - 1, ] +
        (u + uH) * RH[i - 1, ] + (u + uH) * SvH[i - 1, ] + (u + uH) * LvH[i -
          1, ] + (u + uH) * RvH[i - 1, ] + (u + uH + uiHa) * IvH[i - 1, ] +
        (u + uH + uniHa) * NIvH[i - 1, ])

      if (i == ((1 / dt) * (k - year1) + 1 / dt)) {
        i1 <- ((1 / dt) * (k - year1) + 1)
        i2 <- ((1 / dt) * (k - year1) + 1 / dt)
        i3 <- ((1 / dt) * (k - year1) + 2)

        PSIZEy[(k - year1 + 1), 1] <- mean(psize[i1:i2])
        PSIZEy[(k - year1 + 1), 2] <- mean(psize014[i1:i2])
        PSIZEy[(k - year1 + 1), 3] <- mean(psize1554[i1:i2])
        PSIZEy[(k - year1 + 1), 4] <- mean(psize5564[i1:i2])
        PSIZEy[(k - year1 + 1), 5] <- mean(psize65plus[i1:i2])
        PSIZEy[(k - year1 + 1), 6] <- mean(psize1559[i1:i2])
        PSIZEy[(k - year1 + 1), 7] <- mean(psize1529[i1:i2])
        PSIZEy[(k - year1 + 1), 8] <- mean(psize3044[i1:i2])
        PSIZEy[(k - year1 + 1), 9] <- mean(psize4559[i1:i2])
        PSIZEy[(k - year1 + 1), 10] <- mean(psize60plus[i1:i2])
        PSIZEy[(k - year1 + 1), 11] <- mean(psize55plus[i1:i2])
        PSIZEy[(k - year1 + 1), 12] <- mean(psize0509[i1:i2])
        PSIZEy[(k - year1 + 1), 13] <- mean(psize1019[i1:i2])
        PSIZEy[(k - year1 + 1), 14] <- mean(psize2029[i1:i2])
        PSIZEy[(k - year1 + 1), 15] <- mean(psize3039[i1:i2])
        PSIZEy[(k - year1 + 1), 16] <- mean(psize4049[i1:i2])
        PSIZEy[(k - year1 + 1), 17] <- mean(psize5059[i1:i2])
        PSIZEy[(k - year1 + 1), 18] <- mean(psize6069[i1:i2])
        PSIZEy[(k - year1 + 1), 19] <- mean(psize70plus[i1:i2])
        PSIZEy[(k - year1 + 1), 20] <- mean(psize5574[i1:i2])
        PSIZEy[(k - year1 + 1), 21] <- mean(psize75plus[i1:i2])
        PSIZEy[(k - year1 + 1), 22] <- mean(psize1524[i1:i2])
        PSIZEy[(k - year1 + 1), 23] <- mean(psize2554[i1:i2])
        PSIZEy[(k - year1 + 1), 24] <- mean(psize15plus[i1:i2])
        PSIZEy[(k - year1 + 1), 25] <- mean(psize1549[i1:i2])

        PSIZEyH[(k - year1 + 1), 1] <- mean(psizeALL[i1:i2])
        PSIZEyH[(k - year1 + 1), 2] <- mean(psize014H[i1:i2])
        PSIZEyH[(k - year1 + 1), 3] <- mean(psize1554H[i1:i2])
        PSIZEyH[(k - year1 + 1), 4] <- mean(psize5564H[i1:i2])
        PSIZEyH[(k - year1 + 1), 5] <- mean(psize65plusH[i1:i2])
        PSIZEyH[(k - year1 + 1), 6] <- mean(psize1559H[i1:i2])
        PSIZEyH[(k - year1 + 1), 7] <- mean(psize60plusH[i1:i2])
        PSIZEyH[(k - year1 + 1), 8] <- mean(psize15plusH[i1:i2])
        PSIZEyH[(k - year1 + 1), 9] <- mean(psize1549H[i1:i2])
        PSIZEyH[(k - year1 + 1), 10] <- mean(psizeH[i1:i2])
        PSIZEyH[(k - year1 + 1), 11] <- mean(psizeVXH[i1:i2])
        PSIZEyH[(k - year1 + 1), 12] <- mean(psize0509H[i1:i2])
        PSIZEyH[(k - year1 + 1), 13] <- mean(psize1019H[i1:i2])
        PSIZEyH[(k - year1 + 1), 14] <- mean(psize2029H[i1:i2])
        PSIZEyH[(k - year1 + 1), 15] <- mean(psize3039H[i1:i2])
        PSIZEyH[(k - year1 + 1), 16] <- mean(psize4049H[i1:i2])
        PSIZEyH[(k - year1 + 1), 17] <- mean(psize5059H[i1:i2])
        PSIZEyH[(k - year1 + 1), 18] <- mean(psize6069H[i1:i2])
        PSIZEyH[(k - year1 + 1), 19] <- mean(psize70plusH[i1:i2])

        colnames(PSIZEyH) <- c(
          "AllagesALL", "0-14H", "15-54H", "55-64H",
          "65+H", "15-59H", "60+H", "15+H", "1549H", "AllagesH", "VAXH",
          "5-9H", "10-19H", "20-29H", "30-39H", "40-49H", "50-59H", "60-69H",
          "70+H"
        )

        HIVP[(k - year1 + 1), 1] <- 100 * mean(psizeH[i1:i2]) / mean(psizeALL[i1:i2])
        HIVP[(k - year1 + 1), 2] <- 100 * mean(psize014H[i1:i2]) / mean(psize014[i1:i2] +
          psize014H[i1:i2])
        HIVP[(k - year1 + 1), 3] <- 100 * mean(psize15plusH[i1:i2]) / mean(psize15plusH[i1:i2] +
          psize15plus[i1:i2])
        HIVP[(k - year1 + 1), 4] <- 100 * mean(psize1549H[i1:i2]) / mean(psize1549H[i1:i2] +
          psize1549[i1:i2])

        TBI[(k - year1 + 1), 1] <- 1e+05 * sum(new_I_noconv[i1:i2, ], new_NI[i1:i2, ], new_Iv_noconv[i1:i2, ], new_NIv[i1:i2, ]) / mean(psize[i1:i2])
        TBI[(k - year1 + 1), 2] <- 1e+05 * sum(
          new_I_noconv[i1:i2, 1:15],
          new_NI[i1:i2, 1:15], new_Iv_noconv[i1:i2, 1:15], new_NIv[
            i1:i2,
            1:15
          ]
        ) / mean(psize014[i1:i2])
        TBI[(k - year1 + 1), 3] <- 1e+05 * sum(
          new_I_noconv[i1:i2, 16:55],
          new_NI[i1:i2, 16:55], new_Iv_noconv[i1:i2, 16:55], new_NIv[
            i1:i2,
            16:55
          ]
        ) / mean(psize1554[i1:i2])
        TBI[(k - year1 + 1), 4] <- 1e+05 * sum(
          new_I_noconv[i1:i2, 56:65],
          new_NI[i1:i2, 56:65], new_Iv_noconv[i1:i2, 56:65], new_NIv[
            i1:i2,
            56:65
          ]
        ) / mean(psize5564[i1:i2])
        TBI[(k - year1 + 1), 5] <- 1e+05 * sum(
          new_I_noconv[i1:i2, 66:Mnage],
          new_NI[i1:i2, 66:Mnage], new_Iv_noconv[i1:i2, 66:Mnage], new_NIv[
            i1:i2,
            66:Mnage
          ]
        ) / mean(psize65plus[i1:i2])
        TBI[(k - year1 + 1), 6] <- 1e+05 * sum(
          new_I_noconv[i1:i2, 56:Mnage],
          new_NI[i1:i2, 56:Mnage], new_Iv_noconv[i1:i2, 56:Mnage], new_NIv[
            i1:i2,
            56:Mnage
          ]
        ) / mean(psize55plus[i1:i2])
        TBI[(k - year1 + 1), 7] <- 1e+05 * sum(
          new_I_noconv[i1:i2, 1:55],
          new_NI[i1:i2, 1:55], new_Iv_noconv[i1:i2, 1:55], new_NIv[
            i1:i2,
            1:55
          ]
        ) / mean(psize55minus[i1:i2])
        TBI[(k - year1 + 1), 8] <- 1e+05 * sum(
          new_I_noconv[i1:i2, 16:25],
          new_NI[i1:i2, 16:25], new_Iv_noconv[i1:i2, 16:25], new_NIv[
            i1:i2,
            16:25
          ]
        ) / mean(psize1524[i1:i2])
        TBI[(k - year1 + 1), 9] <- 1e+05 * sum(
          new_I_noconv[i1:i2, 26:55],
          new_NI[i1:i2, 26:55], new_Iv_noconv[i1:i2, 26:55], new_NIv[
            i1:i2,
            26:55
          ]
        ) / mean(psize2554[i1:i2])
        TBI[(k - year1 + 1), 10] <- 1e+05 * sum(
          new_I_noconv[i1:i2, 16:Mnage],
          new_NI[i1:i2, 16:Mnage], new_Iv_noconv[i1:i2, 16:Mnage], new_NIv[
            i1:i2,
            16:Mnage
          ]
        ) / mean(psize15plus[i1:i2])

        TBIH[(k - year1 + 1), 1] <- 1e+05 * sum(new_IH_noconv[i1:i2, ], new_NIH[i1:i2, ], new_IvH_noconv[i1:i2, ], new_NIvH[i1:i2, ]) / mean(psizeH[i1:i2])
        TBIH[(k - year1 + 1), 2] <- 1e+05 * sum(
          new_IH_noconv[i1:i2, 1:15],
          new_NIH[i1:i2, 1:15], new_IvH_noconv[i1:i2, 1:15], new_NIvH[
            i1:i2,
            1:15
          ]
        ) / mean(psize014H[i1:i2])
        TBIH[(k - year1 + 1), 3] <- 1e+05 * sum(
          new_IH_noconv[i1:i2, 16:Mnage],
          new_NIH[i1:i2, 16:Mnage], new_IvH_noconv[i1:i2, 16:Mnage], new_NIvH[
            i1:i2,
            16:Mnage
          ]
        ) / mean(psize15plusH[i1:i2])

        TBIHda[(k - year1 + 1), 1] <- 1e+05 * sum(
          new_IH_noconv[i1:i2, ],
          new_NIH[i1:i2, ], new_IvH_noconv[i1:i2, ], new_NIvH[i1:i2, ]
        ) / mean(psizeALL[i1:i2])
        TBIHda[(k - year1 + 1), 2] <- 1e+05 * sum(
          new_IH_noconv[i1:i2, 1:15],
          new_NIH[i1:i2, 1:15], new_IvH_noconv[i1:i2, 1:15], new_NIvH[
            i1:i2,
            1:15
          ]
        ) / mean(psize014H[i1:i2] + psize014[i1:i2])
        TBIHda[(k - year1 + 1), 3] <- 1e+05 * sum(
          new_IH_noconv[i1:i2, 16:Mnage],
          new_NIH[i1:i2, 16:Mnage], new_IvH_noconv[i1:i2, 16:Mnage], new_NIvH[
            i1:i2,
            16:Mnage
          ]
        ) / mean(psize15plusH[i1:i2] + psize15plus[i1:i2])

        TBIall[(k - year1 + 1), 1] <- 1e+05 * sum(
          new_IH_noconv[i1:i2, ],
          new_NIH[i1:i2, ], new_IvH_noconv[i1:i2, ], new_NIvH[i1:i2, ], new_I_noconv[i1:i2, ], new_NI[i1:i2, ], new_Iv_noconv[i1:i2, ], new_NIv[i1:i2, ]
        ) / mean(psizeALL[i1:i2])
        TBIall[(k - year1 + 1), 2] <- 1e+05 * sum(
          new_IH_noconv[i1:i2, 1:15],
          new_NIH[i1:i2, 1:15], new_IvH_noconv[i1:i2, 1:15], new_NIvH[
            i1:i2,
            1:15
          ], new_I_noconv[i1:i2, 1:15], new_NI[i1:i2, 1:15], new_Iv_noconv[
            i1:i2,
            1:15
          ], new_NIv[i1:i2, 1:15]
        ) / mean(psize014H[i1:i2] + psize014[i1:i2])
        TBIall[(k - year1 + 1), 3] <- 1e+05 * sum(
          new_IH_noconv[i1:i2, 16:Mnage],
          new_NIH[i1:i2, 16:Mnage], new_IvH_noconv[i1:i2, 16:Mnage], new_NIvH[
            i1:i2,
            16:Mnage
          ], new_I_noconv[i1:i2, 16:Mnage], new_NI[i1:i2, 16:Mnage],
          new_Iv_noconv[i1:i2, 16:Mnage], new_NIv[i1:i2, 16:Mnage]
        ) / mean(psize15plusH[i1:i2] +
          psize15plus[i1:i2])

        TBIpcHIV[(k - year1 + 1), 1] <- 100 * (sum(
          new_IH_noconv[i1:i2, ],
          new_NIH[i1:i2, ], new_IvH_noconv[i1:i2, ], new_NIvH[i1:i2, ]
        ) / sum(new_I_noconv[i1:i2, ], new_NI[i1:i2, ], new_Iv_noconv[i1:i2, ], new_NIv[i1:i2, ], new_IH_noconv[i1:i2, ], new_NIH[i1:i2, ], new_IvH_noconv[i1:i2, ], new_NIvH[i1:i2, ]))

        newinc[(k - year1 + 1)] <- 1000 * sum(new_I_noconv[i1:i2, ], new_NI[i1:i2, ], new_IH_noconv[i1:i2, ], new_NIH[i1:i2, ], new_Iv_noconv[i1:i2, ], new_NIv[i1:i2, ], new_IvH_noconv[i1:i2, ], new_NIvH[i1:i2, ])

        artnum[(k - year1 + 1), 1] <- 1000 * (sum(artsize[i1:i2, 1]) * dt)
        artnum[(k - year1 + 1), 2] <- 1000 * (sum(artsize[i1:i2, 2]) * dt)

        TBN[(k - year1 + 1), 1] <- 1e+05 * sum(new_notif[i1:i2, ], new_notifv[i1:i2, ]) / mean(psize[i1:i2])
        TBN[(k - year1 + 1), 2] <- 1e+05 * sum(new_notif[i1:i2, 1:15], new_notifv[
          i1:i2,
          1:15
        ]) / mean(psize014[i1:i2])
        TBN[(k - year1 + 1), 3] <- 1e+05 * sum(new_notif[i1:i2, 16:55], new_notifv[
          i1:i2,
          16:55
        ]) / mean(psize1554[i1:i2])
        TBN[(k - year1 + 1), 4] <- 1e+05 * sum(new_notif[i1:i2, 56:65], new_notifv[
          i1:i2,
          56:65
        ]) / mean(psize5564[i1:i2])
        TBN[(k - year1 + 1), 5] <- 1e+05 * sum(
          new_notif[i1:i2, 66:Mnage],
          new_notifv[i1:i2, 66:Mnage]
        ) / mean(psize65plus[i1:i2])
        TBN[(k - year1 + 1), 6] <- 1e+05 * sum(
          new_notif[i1:i2, 56:Mnage],
          new_notifv[i1:i2, 56:Mnage]
        ) / mean(psize55plus[i1:i2])
        TBN[(k - year1 + 1), 7] <- 1e+05 * sum(new_notif[i1:i2, 1:55], new_notifv[
          i1:i2,
          1:55
        ]) / mean(psize55minus[i1:i2])
        TBN[(k - year1 + 1), 8] <- 1e+05 * sum(
          new_notif[i1:i2, 16:Mnage],
          new_notifv[i1:i2, 16:Mnage]
        ) / mean(psize15plus[i1:i2])

        TBNHda[(k - year1 + 1), 1] <- 1e+05 * sum(new_notifH[i1:i2, ], new_notifvH[i1:i2, ]) / mean(psizeALL[i1:i2])
        TBNHda[(k - year1 + 1), 2] <- 1e+05 * sum(
          new_notifH[i1:i2, 1:15],
          new_notifvH[i1:i2, 1:15]
        ) / mean(psize014H[i1:i2] + psize014[i1:i2])
        TBNHda[(k - year1 + 1), 3] <- 1e+05 * sum(
          new_notifH[i1:i2, 16:Mnage],
          new_notifvH[i1:i2, 16:Mnage]
        ) / mean(psize15plusH[i1:i2] + psize15plus[i1:i2])

        TBNall[(k - year1 + 1), 1] <- 1e+05 * sum(new_notifH[i1:i2, ], new_notifvH[i1:i2, ], new_notif[i1:i2, ], new_notifv[i1:i2, ]) / mean(psizeALL[i1:i2])
        TBNall[(k - year1 + 1), 2] <- 1e+05 * sum(
          new_notifH[i1:i2, 1:15],
          new_notifvH[i1:i2, 1:15], new_notif[i1:i2, 1:15], new_notifv[
            i1:i2,
            1:15
          ]
        ) / mean(psize014H[i1:i2] + psize014[i1:i2])
        TBNall[(k - year1 + 1), 3] <- 1e+05 * sum(
          new_notifH[i1:i2, 16:Mnage],
          new_notifvH[i1:i2, 16:Mnage], new_notif[i1:i2, 16:Mnage], new_notifv[
            i1:i2,
            16:Mnage
          ]
        ) / mean(psize15plusH[i1:i2] + psize15plus[i1:i2])

        newtreat[(k - year1 + 1), 1] <- 1000 * sum(new_notif[i1:i2, ], new_notifv[i1:i2, ], new_notifH[i1:i2, ], new_notifvH[i1:i2, ])

        TBP[(k - year1 + 1), 1] <- 1e+05 * (sum(
          I[i1:i2, ], NI[i1:i2, ],
          Iv[i1:i2, ], NIv[i1:i2, ]
        ) / (1 / dt)) / mean(psize[i1:i2])
        TBP[(k - year1 + 1), 2] <- 1e+05 * (sum(I[i1:i2, 1:15], NI[
          i1:i2,
          1:15
        ], Iv[i1:i2, 1:15], NIv[i1:i2, 1:15]) / (1 / dt)) / mean(psize014[i1:i2])
        TBP[(k - year1 + 1), 3] <- 1e+05 * (sum(I[i1:i2, 16:30], NI[
          i1:i2,
          16:30
        ], Iv[i1:i2, 16:30], NIv[i1:i2, 16:30]) / (1 / dt)) / mean(psize1529[i1:i2])
        TBP[(k - year1 + 1), 4] <- 1e+05 * (sum(I[i1:i2, 31:45], NI[
          i1:i2,
          31:45
        ], Iv[i1:i2, 31:45], NIv[i1:i2, 31:45]) / (1 / dt)) / mean(psize3044[i1:i2])
        TBP[(k - year1 + 1), 5] <- 1e+05 * (sum(I[i1:i2, 46:60], NI[
          i1:i2,
          46:60
        ], Iv[i1:i2, 46:60], NIv[i1:i2, 46:60]) / (1 / dt)) / mean(psize4559[i1:i2])
        TBP[(k - year1 + 1), 6] <- 1e+05 * (sum(I[i1:i2, 61:Mnage], NI[
          i1:i2,
          61:Mnage
        ], Iv[i1:i2, 61:Mnage], NIv[i1:i2, 61:Mnage]) / (1 / dt)) / mean(psize60plus[i1:i2])
        TBP[(k - year1 + 1), 7] <- 1e+05 * (sum(I[i1:i2, 56:Mnage], NI[
          i1:i2,
          56:Mnage
        ], Iv[i1:i2, 56:Mnage], NIv[i1:i2, 56:Mnage]) / (1 / dt)) / mean(psize55plus[i1:i2])

        TBPHda[(k - year1 + 1), 1] <- 1e+05 * (sum(IH[i1:i2, ], NIH[i1:i2, ], IvH[i1:i2, ], NIvH[i1:i2, ]) / (1 / dt)) / mean(psizeALL[i1:i2])
        TBPHda[(k - year1 + 1), 2] <- 1e+05 * (sum(IH[i1:i2, 1:15], NIH[
          i1:i2,
          1:15
        ], IvH[i1:i2, 1:15], NIvH[i1:i2, 1:15]) / (1 / dt)) / mean(psize014H[i1:i2] +
          psize014[i1:i2])
        TBPHda[(k - year1 + 1), 3] <- 1e+05 * (sum(IH[i1:i2, 16:Mnage], NIH[
          i1:i2,
          16:Mnage
        ], IvH[i1:i2, 16:Mnage], NIvH[i1:i2, 16:Mnage]) / (1 / dt)) / mean(psize15plusH[i1:i2] +
          psize15plus[i1:i2])

        TBPall[(k - year1 + 1), 1] <- 1e+05 * (sum(I[i1:i2, ], NI[i1:i2, ], Iv[i1:i2, ], NIv[i1:i2, ], IH[i1:i2, ], NIH[i1:i2, ], IvH[i1:i2, ], NIvH[i1:i2, ]) / (1 / dt)) / mean(psizeALL[i1:i2])

        TBPb[(k - year1 + 1), 1] <- 1e+05 * ((sum(I[i1:i2, ], Iv[i1:i2, ])) / (1 / dt)) / mean(psize[i1:i2])
        TBPb[(k - year1 + 1), 2] <- 1e+05 * ((sum(I[i1:i2, 1:15], Iv[
          i1:i2,
          1:15
        ])) / (1 / dt)) / mean(psize014[i1:i2])
        TBPb[(k - year1 + 1), 3] <- 1e+05 * ((sum(I[i1:i2, 16:30], Iv[
          i1:i2,
          16:30
        ])) / (1 / dt)) / mean(psize1529[i1:i2])
        TBPb[(k - year1 + 1), 4] <- 1e+05 * ((sum(I[i1:i2, 31:45], Iv[
          i1:i2,
          31:45
        ])) / (1 / dt)) / mean(psize3044[i1:i2])
        TBPb[(k - year1 + 1), 5] <- 1e+05 * ((sum(I[i1:i2, 46:60], Iv[
          i1:i2,
          46:60
        ])) / (1 / dt)) / mean(psize4559[i1:i2])
        TBPb[(k - year1 + 1), 6] <- 1e+05 * ((sum(I[i1:i2, 61:Mnage], Iv[
          i1:i2,
          61:Mnage
        ])) / (1 / dt)) / mean(psize60plus[i1:i2])
        TBPb[(k - year1 + 1), 7] <- 1e+05 * ((sum(I[i1:i2, 56:Mnage], Iv[
          i1:i2,
          56:Mnage
        ])) / (1 / dt)) / mean(psize55plus[i1:i2])
        TBPb[(k - year1 + 1), 8] <- 1e+05 * ((sum(I[i1:i2, 16:Mnage], Iv[
          i1:i2,
          16:Mnage
        ])) / (1 / dt)) / mean(psize15plus[i1:i2])
        TBPb[(k - year1 + 1), 9] <- 1e+05 * ((sum(I[i1:i2, 31:60], Iv[
          i1:i2,
          31:60
        ])) / (1 / dt)) / ((mean(psize3044[i1:i2])) + (mean(psize4559[i1:i2])))

        TBM[(k - year1 + 1), 1] <- 1e+05 * sum(TBDeaths[i1:i2, ]) / mean(psizeALL[i1:i2])
        TBM[(k - year1 + 1), 2] <- 1e+05 * sum(TBDeaths[i1:i2, 1:15]) / mean(psize014H[i1:i2] +
          psize014[i1:i2])
        TBM[(k - year1 + 1), 3] <- 1e+05 * sum(TBDeaths[i1:i2, 16:55]) / mean(psize1554H[i1:i2] +
          psize1554[i1:i2])
        TBM[(k - year1 + 1), 4] <- 1e+05 * sum(TBDeaths[i1:i2, 56:65]) / mean(psize5564H[i1:i2] +
          psize5564[i1:i2])
        TBM[(k - year1 + 1), 5] <- 1e+05 * sum(TBDeaths[i1:i2, 66:Mnage]) / mean(psize65plusH[i1:i2] +
          psize65plus[i1:i2])
        TBM[(k - year1 + 1), 6] <- 1e+05 * sum(TBDeaths[i1:i2, 16:60]) / mean(psize1559H[i1:i2] +
          psize1559[i1:i2])
        TBM[(k - year1 + 1), 7] <- 1e+05 * sum(TBDeaths[i1:i2, 61:Mnage]) / mean(psize60plus[i1:i2])
        TBM[(k - year1 + 1), 8] <- 1e+05 * sum(TBDeaths[i1:i2, 56:Mnage]) / mean(psize55plus[i1:i2])
        TBM[(k - year1 + 1), 9] <- 1e+05 * sum(TBDeaths[i1:i2, 1:55]) / mean(psize55minus[i1:i2])
        TBM[(k - year1 + 1), 10] <- 1e+05 * sum(TBDeaths[i1:i2, 16:25]) / mean(psize1524[i1:i2])
        TBM[(k - year1 + 1), 11] <- 1e+05 * sum(TBDeaths[i1:i2, 26:55]) / mean(psize2554[i1:i2])

        TBMHda[(k - year1 + 1), 1] <- 1e+05 * sum(TBDeathsH[i1:i2, ]) / mean(psizeALL[i1:i2])

        TBPI[(k - year1 + 1), 1] <- 100 * (((sum(L[i1, ]) / psize[i1]) + (sum(L[i2, ]) / psize[i2])) / 2)

        TBPI[(k - year1 + 1), 2] <- 100 * (((sum(L[i1, 1:15]) / psize014[i1]) +
          (sum(L[i2, 1:15]) / psize014[i2])) / 2)
        TBPI[(k - year1 + 1), 3] <- 100 * (((sum(L[i1, 16:55]) / psize1554[i1]) +
          (sum(L[i2, 16:55]) / psize1554[i2])) / 2)
        TBPI[(k - year1 + 1), 4] <- 100 * (((sum(L[i1, 56:65]) / psize5564[i1]) +
          (sum(L[i2, 56:65]) / psize5564[i2])) / 2)
        TBPI[(k - year1 + 1), 5] <- 100 * (((sum(L[i1, 66:Mnage]) / psize65plus[i1]) +
          (sum(L[i2, 66:Mnage]) / psize65plus[i2])) / 2)
        TBPI[(k - year1 + 1), 6] <- 100 * (((sum(L[i1, 56:Mnage]) / psize55plus[i1]) +
          (sum(L[i2, 56:Mnage]) / psize55plus[i2])) / 2)
        TBPI[(k - year1 + 1), 7] <- 100 * (((sum(L[i1, 6:10]) / psize0509[i1]) +
          (sum(L[i2, 6:10]) / psize0509[i2])) / 2)
        TBPI[(k - year1 + 1), 8] <- 100 * (((sum(L[i1, 11:20]) / psize1019[i1]) +
          (sum(L[i2, 11:20]) / psize1019[i2])) / 2)
        TBPI[(k - year1 + 1), 9] <- 100 * (((sum(L[i1, 21:30]) / psize2029[i1]) +
          (sum(L[i2, 21:30]) / psize2029[i2])) / 2)
        TBPI[(k - year1 + 1), 10] <- 100 * (((sum(L[i1, 31:40]) / psize3039[i1]) +
          (sum(L[i2, 31:40]) / psize3039[i2])) / 2)
        TBPI[(k - year1 + 1), 11] <- 100 * (((sum(L[i1, 41:50]) / psize4049[i1]) +
          (sum(L[i2, 41:50]) / psize4049[i2])) / 2)
        TBPI[(k - year1 + 1), 12] <- 100 * (((sum(L[i1, 51:60]) / psize5059[i1]) +
          (sum(L[i2, 51:60]) / psize5059[i2])) / 2)
        TBPI[(k - year1 + 1), 13] <- 100 * (((sum(L[i1, 61:70]) / psize6069[i1]) +
          (sum(L[i2, 61:70]) / psize6069[i2])) / 2)
        TBPI[(k - year1 + 1), 14] <- 100 * (((sum(L[i1, 71:Mnage]) / psize70plus[i1]) +
          (sum(L[i2, 71:Mnage]) / psize70plus[i2])) / 2)
        TBPI[(k - year1 + 1), 15] <- 100 * (((sum(L[i1, 66:75]) / (psize5574[i1] -
          psize5564[i1])) + (sum(L[i2, 66:75]) / (psize5574[i1] - psize5564[i1]))) / 2)
        TBPI[(k - year1 + 1), 16] <- 100 * (((sum(L[i1, 76:Mnage]) / psize75plus[i1]) +
          (sum(L[i2, 76:Mnage]) / psize75plus[i2])) / 2)
        TBPI[(k - year1 + 1), 17] <- 100 * (((sum(L[i1, 16:25]) / psize1524[i1]) +
          (sum(L[i2, 16:25]) / psize1524[i2])) / 2)

        TBPIage[(k - year1 + 1), ] <- 100 * (((L[i1, ] / (L[i1, ] + R[i1, ] +
          S[i1, ] + I[i1, ] + NI[i1, ])) + (L[i2, ] / (L[i2, ] + R[i2, ] +
          S[i2, ] + I[i2, ] + NI[i2, ]))) / 2)

        TBRa[(k - year1 + 1), 1] <- 100 * ((sum(new_actv_react[i1:i2, ] +
          new_actvH_react[i1:i2, ])) / (sum(new_actv[i1:i2, ] + new_actvH[i1:i2, ])))
        TBRa[(k - year1 + 1), 2] <- 100 * ((sum(new_actv_react[i1:i2, 1:15] +
          new_actvH_react[i1:i2, 1:15])) / (sum(new_actv[i1:i2, 1:15] + new_actvH[
          i1:i2,
          1:15
        ])))
        TBRa[(k - year1 + 1), 3] <- 100 * ((sum(new_actv_react[i1:i2, 16:65] +
          new_actvH_react[i1:i2, 16:65])) / (sum(new_actv[i1:i2, 16:65] + new_actvH[
          i1:i2,
          16:65
        ])))
        TBRa[(k - year1 + 1), 4] <- 100 * ((sum(new_actv_react[i1:i2, 56:65] +
          new_actvH_react[i1:i2, 56:65])) / (sum(new_actv[i1:i2, 56:65] + new_actvH[
          i1:i2,
          56:65
        ])))
        TBRa[(k - year1 + 1), 5] <- 100 * ((sum(new_actv_react[i1:i2, 66:Mnage] +
          new_actvH_react[i1:i2, 66:Mnage])) / (sum(new_actv[i1:i2, 66:Mnage] +
          new_actvH[i1:i2, 66:Mnage])))
        TBRa[(k - year1 + 1), 6] <- 100 * ((sum(new_actv_react[i1:i2, 56:Mnage] +
          new_actvH_react[i1:i2, 56:Mnage])) / (sum(new_actv[i1:i2, 56:Mnage] +
          new_actvH[i1:i2, 56:Mnage])))
        TBRa[(k - year1 + 1), 7] <- 100 * ((sum(new_actv_react[i1:i2, ])) / (sum(new_actv[i1:i2, ])))
        TBRa[(k - year1 + 1), 8] <- 100 * ((sum(new_actvH_react[i1:i2, ])) / (sum(new_actvH[i1:i2, ])))

        TBRi[(k - year1 + 1), 1] <- 100 - (TBRa[(k - year1 + 1), 1])
        TBRi[(k - year1 + 1), 2] <- 100 - (TBRa[(k - year1 + 1), 2])
        TBRi[(k - year1 + 1), 3] <- 100 - (TBRa[(k - year1 + 1), 3])
        TBRi[(k - year1 + 1), 4] <- 100 - (TBRa[(k - year1 + 1), 4])
        TBRi[(k - year1 + 1), 5] <- 100 - (TBRa[(k - year1 + 1), 5])
        TBRi[(k - year1 + 1), 6] <- 100 - (TBRa[(k - year1 + 1), 6])
        TBRi[(k - year1 + 1), 7] <- 100 - (TBRa[(k - year1 + 1), 7])
        TBRi[(k - year1 + 1), 8] <- 100 - (TBRa[(k - year1 + 1), 8])

        TBInew[(k - year1 + 1), 1] <- sum(new_infect[i1:i2, ] + new_infectH[i1:i2, ])
        TBInew[(k - year1 + 1), 2] <- sum(new_infect[i1:i2, 1:15] + new_infectH[
          i1:i2,
          1:15
        ])
        TBInew[(k - year1 + 1), 3] <- sum(new_infect[i1:i2, 16:55] + new_infectH[
          i1:i2,
          16:55
        ])
        TBInew[(k - year1 + 1), 4] <- sum(new_infect[i1:i2, 56:65] + new_infectH[
          i1:i2,
          56:65
        ])
        TBInew[(k - year1 + 1), 5] <- sum(new_infect[i1:i2, 66:Mnage] + new_infectH[
          i1:i2,
          66:Mnage
        ])
        TBInew[(k - year1 + 1), 6] <- sum(new_infect[i1:i2, 56:Mnage] + new_infect[
          i1:i2,
          56:Mnage
        ])

        ARI[(k - year1 + 1), 1] <- (sum(new_infect[i1:i2, ])) / mean(psize[i1:i2]) *
          100
        ARI[(k - year1 + 1), 2] <- (sum(new_infect[i1:i2, 1:15])) / mean(psize014[i1:i2]) *
          100
        ARI[(k - year1 + 1), 3] <- (sum(new_infect[i1:i2, 16:55])) / mean(psize1554[i1:i2]) *
          100
        ARI[(k - year1 + 1), 4] <- (sum(new_infect[i1:i2, 56:Mnage])) / mean(psize55plus[i1:i2]) *
          100

        TBProp[(k - year1 + 1), 1] <- sum(I[i1:i2, 1:15], Iv[i1:i2, 1:15]) / sum(
          I[i1:i2],
          Iv[i1:i2]
        )
        TBProp[(k - year1 + 1), 2] <- sum(I[i1:i2, 16:55], Iv[i1:i2, 16:55]) / sum(
          I[i1:i2],
          Iv[i1:i2]
        )
        TBProp[(k - year1 + 1), 3] <- sum(I[i1:i2, 56:65], Iv[i1:i2, 56:65]) / sum(
          I[i1:i2],
          Iv[i1:i2]
        )
        TBProp[(k - year1 + 1), 4] <- sum(I[i1:i2, 66:Mnage], Iv[i1:i2, 66:Mnage]) / sum(
          I[i1:i2],
          Iv[i1:i2]
        )

        I2050[1, ] <- sum(new_infect[((2000 - year1) * (1 / dt) + 1):((2000 +
          1 - year1) * (1 / dt)), ])

        TBAc[(k - year1 + 1), 1] <- sum(new_actv[i1:i2, ], new_actvv[i1:i2, ], new_actvH[i1:i2, ], new_actvvH[i1:i2, ])
        TBAc_age[(k - year1 + 1), 1] <- sum(new_actv[i1:i2, 1:15], new_actvv[
          i1:i2,
          1:15
        ], new_actvH[i1:i2, 1:15], new_actvvH[i1:i2, 1:15])
        TBAc_age[(k - year1 + 1), 2] <- sum(new_actv[i1:i2, 16:25], new_actvv[
          i1:i2,
          16:25
        ], new_actvH[i1:i2, 16:25], new_actvvH[i1:i2, 16:25])
        TBAc_age[(k - year1 + 1), 3] <- sum(new_actv[i1:i2, 26:55], new_actvv[
          i1:i2,
          26:55
        ], new_actvH[i1:i2, 26:55], new_actvvH[i1:i2, 26:55])
        TBAc_age[(k - year1 + 1), 4] <- sum(new_actv[i1:i2, 56:65], new_actvv[
          i1:i2,
          56:65
        ], new_actvH[i1:i2, 56:65], new_actvvH[i1:i2, 56:65])
        TBAc_age[(k - year1 + 1), 5] <- sum(new_actv[i1:i2, 66:Mnage], new_actvv[
          i1:i2,
          66:Mnage
        ], new_actvH[i1:i2, 66:Mnage], new_actvvH[i1:i2, 66:Mnage])

        TBMo[(k - year1 + 1), 1] <- sum(TBDeaths[i1:i2, ])

        NV[(k - year1 + 1), 1] <- sum(num_vac[i1:i2, ])

        COVtot[(k - year1 + 1), 1] <- 100 * mean(psizeVX[i1:i2]) / mean(psize[i1:i2])
        COVtottest[(k - year1 + 1), 1] <- 100 * mean(psizeVX[i1:i2]) / (sum(
          mean(psizeVX[i1:i2]),
          mean(psizeNOVX[i1:i2])
        ))
        COV9[(k - year1 + 1), 1] <- 100 * mean(psizeVX9[i1:i2]) / (sum(
          mean(psizeVX9[i1:i2]),
          mean(psizeNOVX9[i1:i2])
        ))
        COV9100[(k - year1 + 1), 1] <- 100 * mean(psizeVX9100[i1:i2]) / (sum(
          mean(psizeVX9100[i1:i2]),
          mean(psizeNOVX9100[i1:i2])
        ))

        TBP_age[(k - year1 + 1), 1:Mnage] <- 1000 * (((I[i1, 1:Mnage] + I[
          i2,
          1:Mnage
        ]) / (1 / dt)) + ((Iv[i1, 1:Mnage] + Iv[i2, 1:Mnage]) / (1 / dt)) +
          ((NI[i1, 1:Mnage] + NI[i2, 1:Mnage]) / (1 / dt)) + ((NIv[i1, 1:Mnage] +
            NIv[i2, 1:Mnage]) / (1 / dt)))
        TBPH_age[(k - year1 + 1), 1:Mnage] <- 1000 * (((IH[i1, 1:Mnage] +
          IH[i2, 1:Mnage]) / (1 / dt)) + ((IvH[i1, 1:Mnage] + IvH[i2, 1:Mnage]) / (1 / dt)) +
          ((NIH[i1, 1:Mnage] + NIH[i2, 1:Mnage]) / (1 / dt)) + ((NIvH[i1, 1:Mnage] +
            NIvH[i2, 1:Mnage]) / (1 / dt)))
        TBMo_age[(k - year1 + 1), 1:Mnage] <- 1000 * ((TBDeaths[i1, 1:Mnage] +
          TBDeaths[i2, 1:Mnage]) / (1 / dt))

        COVage[(k - year1 + 1), ] <- 100 * (mapply(
          sum, Sv[i, ], Lv[i, ],
          Rv[i, ], Iv[i, ], NIv[i, ], SvH[i, ], LvH[i, ], RvH[i, ], IvH[i, ], NIvH[i, ]
        ) / mapply(sum, S[i, ], L[i, ], R[i, ], I[i, ], NI[i, ], Sv[i, ], Lv[i, ], Rv[i, ], Iv[i, ], NIv[i, ], SH[i, ], LH[i, ], RH[i, ], IH[i, ], NIH[i, ], SvH[i, ], LvH[i, ], RvH[i, ], IvH[i, ], NIvH[i, ]))

        if (i == (length(seq(year1, (yearend + (1 - dt)), dt)))) {
          vaccgive <- colSums(VX)
        }
      }
    }
  }

  assign("S", S, envir = .GlobalEnv)
  assign("L", L, envir = .GlobalEnv)
  assign("I", I, envir = .GlobalEnv)
  assign("NI", NI, envir = .GlobalEnv)
  assign("R", R, envir = .GlobalEnv)
  assign("new_I", new_I, envir = .GlobalEnv)
  assign("new_I_noconv", new_I_noconv, envir = .GlobalEnv)
  assign("new_NI", new_NI, envir = .GlobalEnv)
  assign("new_notif", new_notif, envir = .GlobalEnv)
  assign("NBirths", BIRTHS, envir = .GlobalEnv)
  assign("NBirthsH", BIRTHSH, envir = .GlobalEnv)
  assign("brate", brate, envir = .GlobalEnv)
  assign("Sv", Sv, envir = .GlobalEnv)
  assign("Lv", Lv, envir = .GlobalEnv)
  assign("Rv", Rv, envir = .GlobalEnv)
  assign("Iv", Iv, envir = .GlobalEnv)
  assign("NIv", NIv, envir = .GlobalEnv)
  assign("new_Iv", new_Iv, envir = .GlobalEnv)
  assign("new_Iv_noconv", new_Iv_noconv, envir = .GlobalEnv)
  assign("new_NIv", new_NIv, envir = .GlobalEnv)
  assign("new_notifv", new_notifv, envir = .GlobalEnv)
  assign("lambda", lambda, envir = .GlobalEnv)
  assign("thetaS", thetaS, envir = .GlobalEnv)
  assign("thetaL", thetaL, envir = .GlobalEnv)
  assign("thetaR", thetaR, envir = .GlobalEnv)
  assign("d", d, envir = .GlobalEnv)
  assign("SH", SH, envir = .GlobalEnv)
  assign("LH", LH, envir = .GlobalEnv)
  assign("IH", IH, envir = .GlobalEnv)
  assign("NIH", NIH, envir = .GlobalEnv)
  assign("RH", RH, envir = .GlobalEnv)
  assign("new_IH", new_IH, envir = .GlobalEnv)
  assign("new_IH_noconv", new_IH_noconv, envir = .GlobalEnv)
  assign("new_NIH", new_NIH, envir = .GlobalEnv)
  assign("new_notifH", new_notifH, envir = .GlobalEnv)
  assign("NBirthsH", BIRTHSH, envir = .GlobalEnv)
  assign("SvH", SvH, envir = .GlobalEnv)
  assign("LvH", LvH, envir = .GlobalEnv)
  assign("RvH", RvH, envir = .GlobalEnv)
  assign("IvH", IvH, envir = .GlobalEnv)
  assign("NIvH", NIvH, envir = .GlobalEnv)
  assign("new_IvH", new_IvH, envir = .GlobalEnv)
  assign("new_IvH_noconv", new_IvH_noconv, envir = .GlobalEnv)
  assign("new_NIvH", new_NIvH, envir = .GlobalEnv)
  assign("new_notifvH", new_notifvH, envir = .GlobalEnv)
  assign("thetaSH", thetaSH, envir = .GlobalEnv)
  assign("thetaLH", thetaLH, envir = .GlobalEnv)
  assign("thetaRH", thetaRH, envir = .GlobalEnv)
  assign("dH", dH, envir = .GlobalEnv)

  assign("TBDeaths", TBDeaths, envir = .GlobalEnv)
  assign("TBDeathsH", TBDeathsH, envir = .GlobalEnv)
  assign("ADeaths", ADeaths, envir = .GlobalEnv)
  assign("ADeathsH", ADeathsH, envir = .GlobalEnv)
  assign("ADeathsN", ADeathsN, envir = .GlobalEnv)
  assign("u", u, envir = .GlobalEnv)
  assign("hiv", hiv, envir = .GlobalEnv)
  assign("psize", psize, envir = .GlobalEnv)
  assign("psizeH", psizeH, envir = .GlobalEnv)
  assign("psizeall", psizeALL, envir = .GlobalEnv)
  assign("bv", bv, envir = .GlobalEnv)
  assign("psize014", psize014, envir = .GlobalEnv)

  assign("psize1529", psize1529, envir = .GlobalEnv)
  assign("psize1554", psize1554, envir = .GlobalEnv)
  assign("psize1559", psize1559, envir = .GlobalEnv)
  assign("psize3044", psize3044, envir = .GlobalEnv)
  assign("psize4559", psize4559, envir = .GlobalEnv)
  assign("psize5564", psize5564, envir = .GlobalEnv)
  assign("psize60plus", psize60plus, envir = .GlobalEnv)
  assign("psize65plus", psize65plus, envir = .GlobalEnv)
  assign("psize55plus", psize55plus, envir = .GlobalEnv)
  assign("psize55minus", psize55minus, envir = .GlobalEnv)
  assign("psize0509", psize0509, envir = .GlobalEnv)
  assign("psize1019", psize1019, envir = .GlobalEnv)
  assign("psize2029", psize2029, envir = .GlobalEnv)
  assign("psize3039", psize3039, envir = .GlobalEnv)
  assign("psize4049", psize4049, envir = .GlobalEnv)
  assign("psize5059", psize5059, envir = .GlobalEnv)
  assign("psize6069", psize6069, envir = .GlobalEnv)
  assign("psize70plus", psize70plus, envir = .GlobalEnv)
  assign("psize5574", psize5574, envir = .GlobalEnv)
  assign("psize75plus", psize75plus, envir = .GlobalEnv)
  assign("psize1524", psize1524, envir = .GlobalEnv)
  assign("psize2554", psize2554, envir = .GlobalEnv)
  assign("psize15plus", psize15plus, envir = .GlobalEnv)
  assign("COVtot", COVtot, envir = .GlobalEnv)
  assign("COVtottest", COVtot, envir = .GlobalEnv)
  assign("COV9", COV9, envir = .GlobalEnv)
  assign("COV9100", COV9100, envir = .GlobalEnv)
  assign("psizeVX", psizeVX, envir = .GlobalEnv)
  assign("psizeVX9100", psizeVX, envir = .GlobalEnv)
  assign("psizeNOVX9100", psizeNOVX, envir = .GlobalEnv)
  assign("psizeVX9", psizeVX, envir = .GlobalEnv)
  assign("psize1549", psizeVX, envir = .GlobalEnv)
  assign("psize1549H", psizeVX, envir = .GlobalEnv)

  assign("HIVP", HIVP, envir = .GlobalEnv)
  assign("TBI", TBI, envir = .GlobalEnv)
  assign("TBIH", TBIH, envir = .GlobalEnv)
  assign("TBIHda", TBIHda, envir = .GlobalEnv)
  assign("TBIall", TBIall, envir = .GlobalEnv)
  assign("TBIpcHIV", TBIpcHIV, envir = .GlobalEnv)
  assign("TBN", TBN, envir = .GlobalEnv)
  assign("TBNHda", TBNHda, envir = .GlobalEnv)
  assign("TBM", TBM, envir = .GlobalEnv)
  assign("TBMHda", TBMHda, envir = .GlobalEnv)
  assign("TBRx", TBRx, envir = .GlobalEnv)
  assign("VX", VX, envir = .GlobalEnv)
  assign("TBP", TBP, envir = .GlobalEnv)
  assign("TBPb", TBPb, envir = .GlobalEnv)
  assign("TBPI", TBPI, envir = .GlobalEnv)
  assign("PSIZEy", PSIZEy, envir = .GlobalEnv)
  assign("PSIZEyH", PSIZEyH, envir = .GlobalEnv)
  assign("TBRa", TBRa, envir = .GlobalEnv)
  assign("TBRi", TBRi, envir = .GlobalEnv)
  assign("TBRa2", TBRa2, envir = .GlobalEnv)
  assign("TBRi2", TBRi2, envir = .GlobalEnv)
  assign("TBInew", TBInew, envir = .GlobalEnv)
  assign("ARI", ARI, envir = .GlobalEnv)
  assign("TBAc", TBAc, envir = .GlobalEnv)
  assign("TBAc_age", TBAc_age, envir = .GlobalEnv)
  assign("TBMo", TBMo, envir = .GlobalEnv)
  assign("num_vac", num_vac, envir = .GlobalEnv)
  assign("NV", NV, envir = .GlobalEnv)
  assign("CDR", CDR, envir = .GlobalEnv)
  assign("CDR2010", CDR2010, envir = .GlobalEnv)
  assign("TBProp", TBProp, envir = .GlobalEnv)
  assign("new_actv", new_actv, envir = .GlobalEnv)
  assign("new_actvv", new_actvv, envir = .GlobalEnv)
  assign("I2050", I2050, envir = .GlobalEnv)
  assign("CDR_av", CDR_av, envir = .GlobalEnv)

  assign("Imatrix", Imatrix, envir = .GlobalEnv)
  assign("psizematrix", psizematrix, envir = .GlobalEnv)

  assign("vaccgive", vaccgive, envir = .GlobalEnv)
  assign("vaccgiveyr", vaccgiveyr, envir = .GlobalEnv)

  assign("AIDSdeaths", AIDSdeaths, envir = .GlobalEnv)
  assign("alldeath", alldeath, envir = .GlobalEnv)
  assign("u", u, envir = .GlobalEnv)
  assign("upop", upop, envir = .GlobalEnv)
  assign("BKdeaths", BKdeaths, envir = .GlobalEnv)

  assign("ustore", ustore, envir = .GlobalEnv)
  assign("upopstore", upopstore, envir = .GlobalEnv)

  assign("uHstore", uHstore, envir = .GlobalEnv)
  assign("uHpopstore", uHpopstore, envir = .GlobalEnv)

  assign("TBP_age", TBP_age, envir = .GlobalEnv)
  assign("TBPH_age", TBPH_age, envir = .GlobalEnv)
  assign("TBMo_age", TBMo_age, envir = .GlobalEnv)
  assign("COVage", COVage, envir = .GlobalEnv)
  assign("newinc", newinc, envir = .GlobalEnv)
  assign("newtreat", newtreat, envir = .GlobalEnv)
  assign("artnum", artnum, envir = .GlobalEnv)
  assign("TBPIage", TBPIage, envir = .GlobalEnv)

  totmort[, 1] <- sum(TBDeaths[((2025 - year1) * (1 / dt) + 1):((2050 + 1 - year1) *
    (1 / dt)), ])
  totmort[, 2] <- sum(TBDeaths[((2025 - year1) * (1 / dt) + 1):((2050 + 1 - year1) *
    (1 / dt)), 66:Mnage])
  totmort[, 3] <- sum(TBDeaths[((2025 - year1) * (1 / dt) + 1):((2050 + 1 - year1) *
    (1 / dt)), 56:Mnage])
  totmort[, 4] <- sum(TBDeaths[((2025 - year1) * (1 / dt) + 1):((2050 + 1 - year1) *
    (1 / dt)), 1:55])

  totcase[, 1] <- sum(new_I_noconv[((2025 - year1) * (1 / dt) + 1):((2050 + 1 - year1) *
    (1 / dt)), ], new_NI[((2025 - year1) * (1 / dt) + 1):((2050 + 1 - year1) * (1 / dt)), ])
  totcase[, 2] <- sum(new_I_noconv[((2025 - year1) * (1 / dt) + 1):((2050 + 1 - year1) *
    (1 / dt)), 66:Mnage], new_NI[((2025 - year1) * (1 / dt) + 1):((2050 + 1 - year1) *
    (1 / dt)), 66:Mnage])
  totcase[, 3] <- sum(new_I_noconv[((2025 - year1) * (1 / dt) + 1):((2050 + 1 - year1) *
    (1 / dt)), 56:Mnage], new_NI[((2025 - year1) * (1 / dt) + 1):((2050 + 1 - year1) *
    (1 / dt)), 56:Mnage])
  totcase[, 4] <- sum(new_I_noconv[((2025 - year1) * (1 / dt) + 1):((2050 + 1 - year1) *
    (1 / dt)), 1:55], new_NI[((2025 - year1) * (1 / dt) + 1):((2050 + 1 - year1) *
    (1 / dt)), 1:55])

  cumulout <- cbind(totmort, totcase)
  colnames(cumulout) <- c(
    "All agesM", "65+M", "55+M", "<55M", "All agesI", "65+I",
    "55+I", "<55I"
  )
  assign("cumulout", cumulout, envir = .GlobalEnv)
  assign("totmort", totmort, envir = .GlobalEnv)
  assign("totcase", totcase, envir = .GlobalEnv)
  assign("cumuloutyr", cumuloutyr, envir = .GlobalEnv)
  assign("CFR", CFR, envir = .GlobalEnv)

  LTBIc <- sum(L[231, 1:15], LH[231, 1:15]) / (psize014[231] + psize014H[231]) *
    100
  LTBIa <- sum(L[231, 16:Mnage], LH[231, 16:Mnage]) / (psize15plus[231] + psize15plusH[231]) *
    100
  assign("LTBIc", LTBIc, envir = .GlobalEnv)
  assign("LTBIa", LTBIa, envir = .GlobalEnv)
  assign("uH", uH, envir = .GlobalEnv)

  colnames(TBI) <- c(
    "All ages", "0-14", "15-54", "55-64", "65+", "55+", "<55",
    "15-24", "25-54", "15+"
  )
  colnames(HIVP) <- c("All ages", "0-14", "15+", "1549")
  colnames(TBN) <- c(
    "All ages", "0-14", "15-54", "55-64", "65+", "55+", "<55",
    "15+"
  )
  colnames(TBM) <- c(
    "All ages", "0-14", "15-54", "55-64", "65+", "15-59", "60+",
    "55+", "<55", "15-24", "25-54"
  )
  colnames(TBP) <- c("All ages", "0-14", "15-29", "30-44", "45-59", "60+", "55+")
  colnames(TBPb) <- c(
    "All ages", "0-14", "15-29", "30-44", "45-59", "60+", "55+",
    "15+", "30-59"
  )
  colnames(TBPI) <- c(
    "All ages", "0-14", "15-54", "55-64", "65+", "55+", "5-9",
    "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70+", "6574", "75+",
    "1524"
  )
  colnames(PSIZEy) <- c(
    "All ages", "0-14", "15-54", "55-64", "65+", "15-59", "15-29",
    "30-44", "45-59", "60+", "55+", "5-9", "10-19", "20-29", "30-39", "40-49",
    "50-59", "60-69", "70+", "55-74", "75+", "15-24", "25-54", "15+", "1549"
  )
  colnames(PSIZEyH) <- c(
    "AllagesALL", "0-14H", "15-54H", "55-64H", "65+H", "15-59H",
    "60+H", "15+H", "1549", "AllagesH", "VAXH", "5-9H", "10-19H", "20-29H", "30-39H",
    "40-49H", "50-59H", "60-69H", "70+H"
  )

  X <- cbind(
    psize, rowSums(S), BIRTHS, rowSums(I), rowSums(NI), rowSums(L), rowSums(R),
    rowSums(new_I), rowSums(new_NI), rowSums(new_I_react), rowSums(new_NI_react),
    rowSums(Sv), rowSums(Lv), rowSums(Rv), TBI[, 1], TBI[, 2], TBI[, 3], TBI[
      ,
      4
    ], TBI[, 5], TBI[, 6], TBI[, 7], TBN[, 1], TBN[, 2], TBN[, 3], TBN[
      ,
      4
    ], TBN[, 5], TBN[, 6], TBN[, 7], TBM[, 1], TBM[, 2], TBM[, 3], TBM[
      ,
      4
    ], TBM[, 5], TBM[, 6], TBM[, 7], TBM[, 8], TBM[, 9], TBP[, 1], TBP[
      ,
      2
    ], TBP[, 3], TBP[, 4], TBP[, 5], TBP[, 6], TBP[, 7], TBPb[, 1], TBPb[
      ,
      2
    ], TBPb[, 3], TBPb[, 4], TBPb[, 5], TBPb[, 6], TBPb[, 7], TBPI[, 1],
    TBPI[, 2], TBPI[, 3], TBPI[, 4], TBPI[, 5], TBPI[, 6], TBPI[, 7], TBPI[
      ,
      8
    ], TBPI[, 9], TBPI[, 10], TBPI[, 11], TBPI[, 12], TBPI[, 13], TBPI[
      ,
      14
    ], PSIZEy[, 1], PSIZEy[, 2], PSIZEy[, 3], PSIZEy[, 4], PSIZEy[, 5],
    PSIZEy[, 6], PSIZEy[, 7], PSIZEy[, 8], PSIZEy[, 9], PSIZEy[, 10], PSIZEy[
      ,
      11
    ], PSIZEy[, 12], PSIZEy[, 13], PSIZEy[, 14], PSIZEy[, 15], PSIZEy[
      ,
      16
    ], PSIZEy[, 17], PSIZEy[, 18], PSIZEy[, 19], PSIZEy[, 20], PSIZEy[
      ,
      21
    ], PSIZEy[, 22], PSIZEy[, 23], TBI[, 8], TBI[, 9], TBM[, 10], TBM[
      ,
      11
    ], TBPb[, 8], PSIZEy[, 24], TBRa[, 1], TBRa[, 2], TBRa[, 3], TBRa[
      ,
      4
    ], TBRa[, 5], TBRi[, 1], TBRi[, 2], TBRi[, 3], TBRi[, 4], TBRi[, 5],
    TBPb[, 9], rowSums(Iv), rowSums(NIv), rowSums(new_Iv), rowSums(new_NIv),
    rowSums(new_Iv_noconv), rowSums(new_notifv), TBI[, 10], TBN[, 8], BIRTHSH,
    HIVP[, 1], HIVP[, 2], HIVP[, 3], HIVP[, 4], rowSums(SH), rowSums(IH), rowSums(NIH),
    rowSums(LH), rowSums(RH), rowSums(new_IH), rowSums(new_NIH), rowSums(new_IH_react),
    rowSums(new_NIH_react), rowSums(SvH), rowSums(IvH), rowSums(NIvH), rowSums(LvH),
    rowSums(RvH), PSIZEyH[, 1], PSIZEyH[, 2], PSIZEyH[, 3], PSIZEyH[, 4], PSIZEyH[
      ,
      5
    ], PSIZEyH[, 6], PSIZEyH[, 7], PSIZEyH[, 8], PSIZEyH[, 9], PSIZEyH[
      ,
      10
    ], PSIZEyH[, 11], PSIZEyH[, 12], PSIZEyH[, 13], PSIZEyH[, 14], PSIZEyH[
      ,
      15
    ], PSIZEyH[, 16], PSIZEyH[, 17], PSIZEyH[, 18], PSIZEyH[, 19], TBIH[
      ,
      1
    ], TBIH[, 2], TBIH[, 3], TBIHda[, 1], TBIHda[, 2], TBIHda[, 3], TBIall[
      ,
      1
    ], TBIall[, 2], TBIall[, 3], TBIpcHIV[, 1], TBNHda[, 1], TBNHda[, 2],
    TBNHda[, 3], TBNall[, 1], TBNall[, 2], TBNall[, 3], TBMHda[, 1], TBPall[
      ,
      1
    ], TBRa[, 7], TBRa[, 8], TBRi[, 7], TBRi[, 8]
  )
  colnames(X) <- c(
    "PSIZE", "S", "Births", "I", "NI", "L", "R", "new_I", "new_NI",
    "new_I_react", "new_NI_react", "Sv", "Lv", "Rv", "TBItot", "TBI0-14", "TBI15-54",
    "TBI55-64", "TBI65+", "TBI55+", "TBI<55", "TBNtot", "TBN0-14", "TBN15-54",
    "TBN55-64", "TBN65+", "TBN55+", "TBN<55", "TBMtot", "TBM0-14", "TBM15-54",
    "TBM55-64", "TBM65+", "TBM15-59", "TBM60+", "TBM55+", "TBM<55", "TBPtot",
    "TBP0-14", "TBP15-29", "TBP30-44", "TBP45-59", "TBP60+", "TBP55+", "TBPbtot",
    "TBPb0-14", "TBPb15-29", "TBPb30-44", "TBPb45-59", "TBPb60+", "TBPb55+",
    "TBPItot", "TBPI0-14", "TBPI15-54", "TBPI55-64", "TBPI65+", "TBPI55+", "TBPI5-9",
    "TBPI10-19", "TBPI20-29", "TBPI30-39", "TBPI40-49", "TBPI50-59", "TBPI60-69",
    "TBPI70+", "YearPsizetot", "YearPsize0-14", "YearPsize15-54", "YearPsize55-64",
    "YearPsize65+", "YearPsize15-59", "YearPsize15-29", "YearPsize30-44", "YearPsize45-59",
    "YearPsize60+", "YearPsize55+", "YearPsize5-9", "YearPsize10-19", "YearPsize20-29",
    "YearPsize30-39", "YearPsize40-49", "YearPsize50-59", "YearPsize60-69", "YearPsize70+",
    "YearPsize55-74", "YearPsize75+", "YearPsize15-24", "YearPsize25-54", "TBI15-24",
    "TBI25-54", "TBM15-24", "TBM25-54", "TBPb15+", "YearPsize15plus", "TBRatot",
    "TBRa0-14", "TBRa15-64", "TBRa55-64", "TBRa65+", "TBRitot", "TBRi0-14", "TBRi15-64",
    "TBRi55-64", "TBRi65+", "TBPb30-59", "Iv", "NIv", "new_Iv", "new_NIv", "new_Iv_noconv",
    "new_notifv", "TBI15+", "TBN15+", "BirthsH", "HIVPtot", "HIVP014", "HIVP15plus",
    "HIVP1549", "SH", "IH", "NIH", "LH", "RH", "new_IH", "new_NIH", "new_I_reactH",
    "new_NI_reactH", "SvH", "IvH", "NIvH", "LvH", "RvH", "PSIZEALL", "YearPsize0-14H",
    "YearPsize15-54H", "YearPsize55-64H", "YearPsize65+H", "YearPsize15-59H",
    "YearPsize60+H", "YearPsize15+H", "YearPsize15-49H", "PSIZEH", "VAXH", "YearPsize5-9H",
    "YearPsize10-19H", "YearPsize20-29H", "YearPsize30-39H", "YearPsize40-49H",
    "YearPsize50-59H", "YearPsize60-69H", "YearPsize70+H", "TBIHtot", "TBIH0-14",
    "TBIH15p", "TBIHdatot", "TBIHda0-14", "TBIHda15p", "TBIalltot", "TBIall0-14",
    "TBIall15p", "TBpcHIV", "TBNHdatot", "TBNHda0-14", "TBNHda15p", "TBNalltot",
    "TBNall0-14", "TBNall15p", "TBMHdatot", "TBPalltot", "TBRatotN", "TBRatotH",
    "TBRitotN", "TBRitotH"
  )

  return(X)
}
