## ----echo=FALSE---------------------------------------------------------------
## To allow nlmixr to reload runs without large run times
## To run the actual models on your system, take the save options off.
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message = FALSE,
  warning = FALSE,
  out.width = "100%"
  )
options(huxtable.knit_print_df = FALSE)

## -----------------------------------------------------------------------------
## Load libraries
library(ggplot2)
library(nlmixr2)
str(theo_sd)

ggplot(theo_sd, aes(TIME, DV)) + geom_line(aes(group=ID), col="red") +
  scale_x_continuous("Time (h)") + scale_y_continuous("Concentration") +
  labs(title="Theophylline single-dose", subtitle="Concentration vs. time by individual")


## -----------------------------------------------------------------------------
one.cmt <- function() {
  ini({
    ## You may label each parameter with a comment
    tka <- 0.45 # Ka
    tcl <- log(c(0, 2.7, 100)) # Log Cl
    ## This works with interactive models
    ## You may also label the preceding line with label("label text")
    tv <- 3.45; label("log V")
    ## the label("Label name") works with all models
    eta.ka ~ 0.6
    eta.cl ~ 0.3
    eta.v ~ 0.1
    add.sd <- 0.7
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv + eta.v)
    linCmt() ~ add(add.sd)
  })
}

f <- nlmixr(one.cmt)

## -----------------------------------------------------------------------------
fit <- nlmixr(one.cmt, theo_sd, est="focei",
              control=list(print=0))

print(fit)

## -----------------------------------------------------------------------------
one.compartment <- function() {
  ini({
    tka <- 0.45 # Log Ka
    tcl <- 1 # Log Cl
    tv <- 3.45    # Log V
    eta.ka ~ 0.6
    eta.cl ~ 0.3
    eta.v ~ 0.1
    add.sd <- 0.7
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv + eta.v)
    d/dt(depot) = -ka * depot
    d/dt(center) = ka * depot - cl / v * center
    cp = center / v
    cp ~ add(add.sd)
  })
}

## -----------------------------------------------------------------------------
fit2 <- nlmixr(one.compartment, theo_sd,  est="saem",
               control=list(print=0))
print(fit2)

## -----------------------------------------------------------------------------
fitN <- nlmixr(one.compartment, theo_sd, list(pnlsTol=0.5), est="nlme")
print(fitN)

## -----------------------------------------------------------------------------
plot(fit)

## -----------------------------------------------------------------------------
print(fit)

## -----------------------------------------------------------------------------
fit$eta

## -----------------------------------------------------------------------------
traceplot(fit)

## -----------------------------------------------------------------------------
f <- function(){
  ini({
    lCl <- 1.6      #log Cl (L/hr)
    lVc <- log(90)  #log Vc (L)
    lKA <- 0.1      #log Ka (1/hr)
    prop.err <- c(0, 0.2, 1)
    eta.Cl ~ 0.1   # BSV Cl
    eta.Vc ~ 0.1   # BSV Vc
    eta.KA ~ 0.1   # BSV Ka
  })
  model({
    Cl <- exp(lCl + eta.Cl)
    Vc = exp(lVc + eta.Vc)
    KA <- exp(lKA + eta.KA)
    ## Instead of specifying the ODEs, you can use
    ## the linCmt() function to use the solved system.
    ##
    ## This function determines the type of PK solved system
    ## to use by the parameters that are defined.  In this case
    ## it knows that this is a one-compartment model with first-order
    ## absorption.
    linCmt() ~ prop(prop.err)
  })
}

## -----------------------------------------------------------------------------
nlmixr(f)

