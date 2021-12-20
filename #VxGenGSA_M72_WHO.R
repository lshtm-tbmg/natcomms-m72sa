yrintro <- 2028

times <- seq(year1, (yearend + (1 - dt)), dt)
steps <- length(times)

dV2 <- matrix(0, steps, Mnage)
thetaV2 <- matrix(0, steps, Mnage)
thetaV2a <- matrix(0, steps, Mnage)
thetaV2m <- matrix(0, steps, Mnage)
thetablank <- matrix(0, steps, Mnage)
if (yrintro < yearend) {
  startyr <- seq((yrintro - year1) * (1 / dt) + 1, (yearend - year1) * (1 / dt) + 1, (1 / dt))

  if (D < fmass) {
    spacing <- fmass
  } else {
    spacing <- D
  }

  massvxyr <- seq((yrintro - year1) * (1 / dt) + 1, (yearend - year1) * (1 / dt) + 1, spacing * (1 / dt))

  thetaV2[startyr, (vxage + 1)] <- coverage
  thetaV2a[startyr, (vxage + 1)] <- coverage

  if (Mmax < Mnage) {
    thetaV2[massvxyr, (vxage + 2):(Mmax + 1)] <- coverageM
    thetaV2m[massvxyr, (vxage + 2):(Mmax + 1)] <- coverageM
  } else {
    thetaV2[massvxyr, (vxage + 2):Mnage] <- coverageM
    thetaV2m[massvxyr, (vxage + 2):Mnage] <- coverageM
  }

  if (D <= (yearend - yrintro)) {
    dV2[startyr[-1:-D], (vxage + 1 + D)] <- 1

    if (coverageM > 0) {
      newvaxprop <- ((1 - coverage) * coverageM) + coverage

      routinewaneprop <- (newvaxprop - coverageM) / newvaxprop

      dV2[startyr[-1:-D], (vxage + 1 + D)] <- routinewaneprop

      massexityr <- massvxyr + (D * (1 / dt))
      massexityr <- massexityr[massexityr <= (yearend - year1 + 1) * (1 / dt)]

      if (Mmax < (Mnage - D)) {
        dV2[massexityr, (vxage + 1 + D):(Mmax + 1 + D)] <- 1
      } else {
        dV2[massexityr, (vxage + 1 + D):Mnage] <- 1
      }
    }
  }
}

if (vaccine == 1) {
  thetaS <- thetaV2
  thetaL <- thetaV2
  thetaR <- thetaV2
  d <- dV2
} else if (vaccine == 2) {
  thetaS <- thetaV2
  thetaL <- thetablank
  thetaR <- thetablank
  d <- dV2
} else if (vaccine == 3) {
  thetaS <- thetablank
  thetaL <- thetaV2
  thetaR <- thetaV2
  d <- dV2
}

thetaSH <- thetaS * covH
thetaLH <- thetaL * covH
thetaRH <- thetaR * covH

dH <- d

assign("thetaV2a", thetaV2a, envir = .GlobalEnv)
assign("thetaV2m", thetaV2m, envir = .GlobalEnv)
assign("thetaV2", thetaV2, envir = .GlobalEnv)
assign("d", d, envir = .GlobalEnv)
assign("dH", dH, envir = .GlobalEnv)

assign("thetaS", thetaS, envir = .GlobalEnv)
assign("thetaSH", thetaSH, envir = .GlobalEnv)
assign("thetaL", thetaL, envir = .GlobalEnv)
assign("thetaLH", thetaLH, envir = .GlobalEnv)
assign("thetaR", thetaR, envir = .GlobalEnv)
assign("thetaRH", thetaRH, envir = .GlobalEnv)
