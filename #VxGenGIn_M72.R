yrintro <- 2025

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

  thetaV2[massvxyr, (vxage + 2):Mmax] <- coverageM
  thetaV2m[massvxyr, (vxage + 2):Mmax] <- coverageM

  if (D <= (yearend - yrintro)) {
    dV2[startyr, (vxage + 1 + D)] <- 1

    if (coverageM > 0) {
      massexityr <- massvxyr + (D * (1 / dt))
      massexityr <- massexityr[massexityr <= (yearend - year1 + 1) * (1 / dt)]

      dV2[massexityr, (vxage + 2 + D):Mnage] <- 1

      dropyrs <- c()
      for (i in 1:length(massvxyr)) {
        dropyrs <- c(dropyrs, ((massvxyr[i] + 1):min(steps, (massvxyr[i] + ((D - 1) / dt)))))
      }

      dV2[dropyrs, (vxage + 1 + D)] <- 0
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

assign("thetaV2a", thetaV2a, envir = .GlobalEnv)
assign("thetaV2m", thetaV2m, envir = .GlobalEnv)
assign("thetaV2", thetaV2, envir = .GlobalEnv)
assign("d", d, envir = .GlobalEnv)
