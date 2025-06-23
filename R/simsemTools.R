#Functions for book simsem features

#create a function that will allow us to populate k=k datasets, save those datasets in a list, 
#which we can then call into further functions. ##PROBLEM. THESE DATASETS DON'T ACCOUNT FOR CENTRAL LIMIT THEOREM, REDUCING VARIANCES AND OVERESTIMATING VALUES.

# simulateDataSave_old <- function(model, nSamples = 1, group1.nobs, group2.nobs, miss = NULL, seed = seed) {
#   sampdat <- NULL
#   seed <- 12345
#   for (i in 1:nSamples){
#     sampdat[[i]] <- simulateData(model, sampdat.nobs = c(group1.nobs, group2.nobs), seed = seed) #to do: update to include 1 to n groups
#     sampdat[[i]] <- as.data.frame(sampdat[[i]])
#     sampdat[[i]]$iteration <- i
#     sampdat[[i]] <- if(is.null(miss)) {
#       sampdat[[i]]
#     } else {
#       impose(miss, sampdat[[i]])
#     }
#   }
#   #sampdat <- do.call("rbind", sampdat) #can be used to combine datasets into a dataframe
#   return(sampdat)
# }


##### my functions #####

#Function that creates population dataset for j groups
#Draws k samples of (Nj)-size dataset for group 1:j
#Impose missingness on each sample dataset
#Returns the list of sample datasets which we can then call into further functions
simulateDataSets <- function(model, nIterations = 1, nObservations, miss = NULL, seed = seed, estimator = "ML", meanstructure = FALSE) {
  seed <- 12345
  fixed.x <- FALSE
  sampdat <- rep(list(list(nObservations)), nIterations)
  nGroups = as.numeric(length(nObservations)) #will create the number of groups based on the length of the sampdat.nobs vector; e.g., c(20, 20, 20) will lead to 3 groups
  popObs <- list(rep(10000, nGroups))  #will create a placeholder population dataset for each group, based on model specifications
  popObs <- paste0("c(", do.call(paste0, c(popObs, collapse = ", ")), ")")
  population <- eval(parse(text = paste0("simulateData(model, sample.nobs = ", paste0(popObs), ", seed = ", paste0(seed), ", meanstructure = ", paste0(meanstructure) , ", estimator = \"", paste0(estimator), "\")")))
  if(nGroups > 1) {
    for (i in 1:nIterations){
      for (j in 1:nGroups){
        sampdat[[i]][[j]] <- population[population$group == j,][sample(nrow(population[population$group == j,]), nObservations[j]), ]
        sampdat[[i]][[j]]$iteration <- i
      }
      sampdat[[i]] <- do.call("rbind", sampdat[[i]])
      sampdat[[i]] <- if(is.null(miss)) {
        sampdat[[i]]
      } else {
        impose(miss, sampdat[[i]])
      }
    }
    #sampdat <- do.call("rbind", sampdat) #can be used to combine datasets into a dataframe
    return(sampdat)
  } else {
    for (i in 1:nIterations){
      sampdat[[i]] <- population[sample(nrow(population), nObservations), ]
      sampdat[[i]]$iteration <- i
      #sampdat[[i]] <- do.call("rbind", sampdat[[i]]) #not needed for the one group option
      sampdat[[i]] <- if(is.null(miss)) {
        sampdat[[i]]
      } else {
        impose(miss, sampdat[[i]])
      }
    }
    #sampdat <- do.call("rbind", sampdat) #can be used to combine datasets into a dataframe
    return(sampdat)
  }
}

# simsemPower <- function(model, 
#                         nIterations = 1, 
#                         nObservations, 
#                         miss = NULL, 
#                         seed = seed, 
#                         estimator = "ML")
#   #population model using lavaan or simsem syntax
#   #starting N, n groups, group proportions, parameter to power on,
#   #mini-sim (n = 20) for starting N
#   #assess power estimate for parameter of interest
#   #incremental N increase, conditional on parameter power estimate(s) 
#     #for minimum power estimate, (>.70 n = 20, >.75 n = 50, > .79 n = 100, <.81 n = 1000)
#     #option for any power estimate (only one of the effects need to be power .80 or above)
#   #option to change power level; in this case, the lower threshold will be power-.1, and increase accordingly

positive_definite <- function(PopulationModel){
  if(det(fitted(sem(PopulationModel))$cov) > 0) {
    text = paste0("The population model is positive definite. The determinant of the covariance matrix is ", 
                  paste0(det(fitted(sem(PopulationModel))$cov)), ".")
  } else {
    text = paste0("The population model is not-positive definite. The determinant of the covariance matrix is ", 
                  paste0(det(fitted(sem(PopulationModel))$cov)), ". Check your model before proceeding.")
  }
  print(text)
}




#Functions that aggregate additional summary outputs to better examine data coverage (missingness)
cov_coverage <- function(fit) {
  coverage <- lavInspect(fit, "coverage")
}

groupavg_cov_coverage <- function(output, group) {
  round(tril(apply(simplify2array(simplify2array(getExtraOutput(output))[group,]), 1:2, mean)), 3)
}

groupsd_cov_coverage <- function(output, group) {
  round(tril(apply(simplify2array(simplify2array(getExtraOutput(output))[group,]), 1:2, sd)), 3)
}

avg_cov_coverage <- function(output) {
  round(tril(apply(simplify2array(getExtraOutput(output)), 1:2, mean)), 3)
}
sd_cov_coverage <- function(output) {
  round(tril(apply(simplify2array(getExtraOutput(output)), 1:2, sd)), 3)
}

# library(simstudy)
# genCorGen(n = 10, nvars = 2, params1 = c(.5, 1-.5), dist = "binary", rho = .5, corstr = "cs", method = "copula", wide = T)
# 
# simbinary <- function(R, p) {
#   x = genCorGen(n = 100000, nvars = 2, params1 = c(p, 1-p), dist = "binary", rho = R, corstr = "cs", method = "ep", wide = T)
#   corX = cor(x$V1, x$V2)
#   print(corX)
# }
# simbinary(.3, .5)
#
# V <- dTemp[, stats::qbinom(Unew, 1, param1)]
# dbinom(.5, 100000, .3) / 100000

#adapted from StupidWolf response to https://stackoverflow.com/questions/65628556/generate-a-binary-variable-with-a-predefined-correlation-to-an-already-existing
simcorcovars = function(rho, x, type){
  n = length(x)
  x_bar = mean(x)
  xy_bar = rho * x_bar + (1-rho)*x_bar^2
  toflip = sum(x == 1) - round(n * xy_bar)
  if(type == "bothbinary") {
    y = x
    y[sample(which(x==0),toflip)] = 1
    y[sample(which(x==1),toflip)] = 0
    return(y)
  } else if (type == "binarycont") {
    y = rnorm(n, rho*x, sqrt(1-rho^2))
    return(y)
  }
}
# 
# simpointbiserial = function(rho, x){
# 
#   n = length(x)
#   x_bar = mean(x)
#   xy_bar = rho * x_bar + (1-rho)*x_bar^2
#   toflip = sum(x == 1) - round(n * xy_bar)
#   y = x
#   y[sample(which(x==0),toflip)] = 1
#   y[sample(which(x==1),toflip)] = 0
#   return(y)
# }
# 
simbinomial = function(rho, p, n, binvarname, contvarname) {
  sigma = matrix(c(1,rho,rho,1), ncol=2)
  s = chol(sigma)
  z = s%*%matrix(rnorm(n*2), nrow=2)
  u = pnorm(z)
  contvar = qgamma(u[1,], 15, .5)
  binvar = ifelse(u[2,] > p, 1, 0)
  y = cbind(contvar, binvar)
  y = as.data.frame(y)
  colnames(y) <- c(contvarname, binvarname)
  return(y)
}
# 
# 
# simbinomial = function(rho, data, binvarname, contvarname) {
#   sigma = matrix(c(1,rho,rho,1), ncol=2)
#   s = chol(sigma)
#   z = s%*%matrix(rnorm(n*2), nrow=2)
#   u = pnorm(z)
#   contvar = qgamma(u[1,], 1)
#   binvar = ifelse(u[2,] > p, 1, 0)
#   y = cbind(contvar, binvar)
#   y = as.data.frame(y)
#   colnames(y) <- c(contvarname, binvarname)
#   return(y)
# }

##### simulateData2: adapted very closely from simulateData from lavaan, but change the tol in the mvrnorm from 1e-6 to 1e-10 #####

# ALL CODE, WITH THE EXCEPTION OF THE MINOR CHANGES IN THE TOL IN MVRNORM ARE ORIGINALLY FROM LAVANN, NOT simsemTools.

# simulate data starting from a user-specified model
#
# initial version: YR 24 jan 2011
# revision for 0.4-11: YR 21 okt 2011
simulateData2 <- function(
                         # user-specified model
                         model           = NULL,
                         model.type      = "sem",

                         # model modifiers
                         meanstructure   = FALSE,
                         int.ov.free     = TRUE,
                         int.lv.free     = FALSE,
                         conditional.x   = FALSE,
                         fixed.x         = FALSE,
                         orthogonal      = FALSE,
                         std.lv          = TRUE,

                         auto.fix.first  = FALSE,
                         auto.fix.single = FALSE,
                         auto.var        = TRUE,
                         auto.cov.lv.x   = TRUE,
                         auto.cov.y      = TRUE,
                         ...,

                         # data properties
                         sample.nobs     = 500L,
                         ov.var          = NULL,
                         group.label     = paste("G", 1:ngroups, sep=""),
                         skewness        = NULL,
                         kurtosis        = NULL,

                         # control
                         seed            = NULL,
                         empirical       = FALSE,

                         return.type     = "data.frame",
                         return.fit      = FALSE,
                         debug           = FALSE,
                         standardized    = FALSE
                        )
{
    if(!is.null(seed)) set.seed(seed)
    #if(!exists(".Random.seed", envir = .GlobalEnv))
    #    runif(1)               # initialize the RNG if necessary
    #RNGstate <- .Random.seed

    # lavaanify
    if(is.list(model)) {
        # two possibilities: either model is already lavaanified
        # or it is something else...
        if(!is.null(model$lhs) && !is.null(model$op)  &&
           !is.null(model$rhs) && !is.null(model$free)) {
            lav <- model

            # until 0.6-5, we only used the 'ustart' column
            # but what if 'lav' is a fitted lavaan object -> use 'est'
            if(!is.null(lav$est)) {
                lav$ustart <- lav$est
                lav$se <- NULL
                lav$est <- NULL
                lav$start <- NULL
            }
        } else if(is.character(model[[1]])) {
            stop("lavaan ERROR: model is a list, but not a parameterTable?")
        }
    } else {
        lav <- lavaanify(model = model,
                         meanstructure=meanstructure,
                         int.ov.free=int.ov.free,
                         int.lv.free=int.lv.free,
                         conditional.x=conditional.x,
                         fixed.x=fixed.x,
                         orthogonal=orthogonal,
                         std.lv=std.lv,
                         auto.fix.first=auto.fix.first,
                         auto.fix.single=auto.fix.single,
                         auto.var=auto.var,
                         auto.cov.lv.x=auto.cov.lv.x,
                         auto.cov.y=auto.cov.y,
                         ngroups=length(sample.nobs))
    }

    group.values <- lav_partable_group_values(lav)
    if(debug) {
        cat("initial lav\n")
        print(as.data.frame(lav))
    }

    # fill in any remaining NA values (needed for unstandardize)
    # 1 for variances and (unstandardized) factor loadings, 0 otherwise
    idx <- which(lav$op == "=~" & is.na(lav$ustart))
    if(length(idx) > 0L) {
        if(standardized) {
             lav$ustart[idx] <- 0.7
        } else {
             lav$ustart[idx] <- 1.0
        }
    }

    idx <- which(lav$op == "~~" & is.na(lav$ustart) & lav$lhs == lav$rhs)
    if(length(idx) > 0L) lav$ustart[idx] <- 1.0

    idx <- which(lav$op == "~" & is.na(lav$ustart))
    if(length(idx) > 0L) {
        warning("lavaan WARNING: some regression coefficients are unspecified and will be set to zero")
    }

    idx <- which(is.na(lav$ustart))
    if(length(idx) > 0L) lav$ustart[idx] <- 0.0

    if(debug) {
        cat("lav + default values\n")
        print(as.data.frame(lav))
    }

    # set residual variances to enforce a standardized solution
    # but only if no *residual* variances have been specified in the syntax

    if(standardized) {
        # check if factor loadings are smaller than 1.0
        lambda.idx <- which(lav$op == "=~")
        if(any(lav$ustart[lambda.idx] >= 1.0)) {
            warning("lavaan WARNING: standardized=TRUE but factor loadings are >= 1.0")
        }

        # check if regression coefficients are smaller than 1.0
        reg.idx <- which(lav$op == "~")
        if(any(lav$ustart[reg.idx] >= 1.0)) {
            warning("lavaan WARNING: standardized=TRUE but regression coefficients are >= 1.0")
        }

        # for ordered observed variables, we will get '0.0', but that is ok
        # so there is no need to make a distinction between numeric/ordered
        # here??
        ngroups <- lav_partable_ngroups(lav)
        ov.names <- vnames(lav, "ov")
        ov.nox   <- vnames(lav, "ov.nox")
        lv.names <- vnames(lav, "lv")
        lv.y     <- vnames(lav, "lv.y")
        lv.nox   <- vnames(lav, "lv.nox")
        ov.var.idx <- which(lav$op == "~~" & lav$lhs %in% ov.nox &
                            lav$rhs == lav$lhs)
        lv.var.idx <- which(lav$op == "~~" & lav$lhs %in% lv.nox &
                            lav$rhs == lav$lhs)
        if(any(lav$user[c(ov.var.idx, lv.var.idx)] > 0L)) {
            warning("lavaan WARNING: if residual variances are specified, please use standardized=FALSE")
        }
        lav$ustart[c(ov.var.idx,lv.var.idx)] <- 0.0
        fit <- lavaan(model=lav, sample.nobs=sample.nobs, ...)
        Sigma.hat <- computeSigmaHat(lavmodel = fit@Model)
        ETA <- computeVETA(lavmodel = fit@Model)

        if(debug) {
            cat("Sigma.hat:\n"); print(Sigma.hat)
            cat("Eta:\n"); print(ETA)
        }

        # stage 1: standardize LV
        if(length(lv.nox) > 0L) {
            for(g in 1:ngroups) {
                var.group <- which(lav$op == "~~" & lav$lhs %in% lv.nox &
                                   lav$rhs == lav$lhs &
                                   lav$group == group.values[g])
                eta.idx <- match(lv.nox, lv.names)
                lav$ustart[var.group] <- 1 - diag(ETA[[g]])[eta.idx]
            }
        }
        # refit
        fit <- lavaan(model=lav, sample.nobs=sample.nobs, ...)
        Sigma.hat <- computeSigmaHat(lavmodel = fit@Model)

        if(debug) {
            cat("after stage 1:\n")
            cat("Sigma.hat:\n"); print(Sigma.hat)
        }

        # stage 2: standardize OV
        for(g in 1:ngroups) {
            var.group <- which(lav$op == "~~" & lav$lhs %in% ov.nox &
                               lav$rhs == lav$lhs &
                               lav$group == group.values[g])
            ov.idx <- match(ov.nox, ov.names)
            lav$ustart[var.group] <- 1 - diag(Sigma.hat[[g]])[ov.idx]
        }

        if(debug) {
            cat("after standardisation lav\n")
            print(as.data.frame(lav))
        }
    }


    # unstandardize
    if(!is.null(ov.var)) {
        # FIXME: if ov.var is named, check the order of the elements

        # 1. unstandardize observed variables
        lav$ustart <- lav_unstandardize_ov(partable = lav, ov.var = ov.var)

        # 2. unstandardized latent variables

        if(debug) {
            cat("after unstandardisation lav\n")
            print(as.data.frame(lav))
        }
    }

    # fit the model without data
    fit <- lavaan(model=lav, sample.nobs=sample.nobs,  ...)

    # the model-implied moments for the population
    Sigma.hat <- computeSigmaHat(lavmodel = fit@Model)
       Mu.hat <- computeMuHat(lavmodel = fit@Model)
    if(fit@Model@categorical) {
       TH <- computeTH(lavmodel = fit@Model)
    }

    if(debug) {
        cat("\nModel-implied moments (before Vale-Maurelli):\n")
        print(Sigma.hat)
        print(Mu.hat)
        if(exists("TH")) print(TH)
    }

    # ngroups
    ngroups <- length(sample.nobs)

    # prepare
    X <- vector("list", length=ngroups)
    out <- vector("list", length=ngroups)

    for(g in 1:ngroups) {
        COV <- Sigma.hat[[g]]

        # if empirical = TRUE, rescale by N/(N-1), so that estimator=ML
        # returns exact results
        if(empirical) {
            COV <- COV * sample.nobs[g] / (sample.nobs[g] - 1)
        }

        # FIXME: change to rmvnorm once we include the library?
        if(is.null(skewness) && is.null(kurtosis) && is.null(bin_x)) {
            X[[g]] <- MASS::mvrnorm(n = sample.nobs[g],
                                    mu = Mu.hat[[g]],
                                    Sigma = COV,
                                    empirical = empirical,
                                    tol = 1e-10)
        } else if (!is.null(skewness) && (!is.null(kurtosis)) && is.null(bin_x)) {
            # first generate Z
            Z <- ValeMaurelli1983(n        = sample.nobs[g],
                                  COR      = cov2cor(COV),
                                  skewness = skewness,  # FIXME: per group?
                                  kurtosis = kurtosis,
                                  debug    = debug)
            # rescale
            # Note: 'scale()' will first center, and then scale
            # but we need to first scale, and then center...
            # this was reported by Jordan Brace (9 may 2014)
            #X[[g]] <- scale(Z, center = -Mu.hat[[g]],
            #                   scale  = 1/sqrt(diag(COV)))

            # first, we scale
            TMP <- scale(Z, center = FALSE,
                         scale = 1/sqrt(diag(COV)))[,,drop=FALSE]

            # then, we center
            X[[g]] <- sweep(TMP, MARGIN=2, STATS=Mu.hat[[g]], FUN="+")
        } else {                                                                    ###THIS IS WHERE I'VE TRIED UPDATING. NOT WORKING
            # include a way to account for binary predictor variables
            # here, we can update the TMP file to 
            TMP <- scale(Z, center = FALSE,
                         scale = 1/sqrt(diag(COV)))[,,drop=FALSE]

            # then, we center
            X[[g]] <- sweep(TMP, MARGIN=2, STATS=Mu.hat[[g]], FUN="+")
        }                                                                           ###END OF UPDATES

        # any categorical variables?
        ov.ord <- vnames(lav, type="ov.ord", group = group.values[g])
        if(length(ov.ord) > 0L) {
            ov.names <- vnames(lav, type="ov", group = group.values[g])
            # use thresholds to cut
            for(o in ov.ord) {
                o.idx <- which(o == ov.names)
                th.idx <- which(lav$op == "|" & lav$lhs == o &
                                lav$group == group.values[g])
                th.val <- c(-Inf,sort(lav$ustart[th.idx]),+Inf)
                X[[g]][,o.idx] <- as.integer(cut(X[[g]][,o.idx], th.val))
            }
        }

        if(return.type == "data.frame") X[[g]] <- as.data.frame(X[[g]])
    }

    if(return.type == "matrix") {
        if(ngroups == 1L) {
            return(X[[1L]])
        } else {
            return(X)
        }

    } else if (return.type == "data.frame") {
        Data <- X[[1L]]

        # if multiple groups, add group column
        if(ngroups > 1L) {
            for(g in 2:ngroups) {
                Data <- rbind(Data, X[[g]])
            }
            Data$group <- rep(1:ngroups, times=sample.nobs)
        }
        var.names <- vnames(fit@ParTable, type="ov", group=1L)
        if(ngroups > 1L) var.names <- c(var.names, "group")
        names(Data) <- var.names
        if(return.fit) {
            attr(Data, "fit") <- fit
        }
        return(Data)

    } else if (return.type == "cov") {
        if(ngroups == 1L) {
            return(cov(X[[1L]]))
        } else {
            cov.list <- lapply(X, cov)
            return(cov.list)
        }
    }
}

Skewness <- function(x., N1=TRUE) {
    x <- x.; x <- x[!is.na(x)]; N <- length(x)
    mean.x <- mean(x); xc <- x - mean.x; var.x <- var(x)
    if(!N1) var.x <- var.x * (N-1)/N
    sd.x <- sqrt(var.x)
    sk <- sum(xc*xc*xc)/(sd.x*sd.x*sd.x)
    skewness <- N*sk/((N-1)*(N-2))
    skewness
}

Kurtosis <- function(x., N1=TRUE) {
    x <- x.; x <- x[!is.na(x)]; N <- length(x)
    mean.x <- mean(x); xc <- x - mean.x; var.x <- var(x)
    if(!N1) var.x <- var.x * (N-1)/N
    k <- sum(xc*xc*xc*xc)/(var.x*var.x)
    kurtosis <- N*(N+1)*k/((N-1)*(N-2)*(N-3))-3*(N-1)*(N-1)/((N-2)*(N-3))
    kurtosis
}

# NOTE: as pointed out in Fleishman (1978), a real solution does not
# always exist (for a/b/c/d) for all values of skew/kurtosis
#
# for example: skew = 3, only valid if kurtosis > 14 (approximately)
#
# fleishman eq 21 suggests: skew^2 < 0.0629576*kurtosis + 0.0717247
# see figure 1 page 527
#
# note also that the a/b/c/d solution is not unique, although this seems
# not to matter for generating the data

# Fleishman (1978) cubic transformation method
lav_fleishman1978 <- function(n=100, skewness=0, kurtosis=0, verbose=FALSE) {

    system.function <- function(x, skewness, kurtosis) {
        b=x[1L]; c=x[2L]; d=x[3L]
        eq1 <- b*b + 6*b*d + 2*c*c + 15*d*d - 1
        eq2 <- 2*c*(b*b + 24*b*d + 105*d*d + 2) - skewness
        eq3 <- 24*(b*d + c*c*(1 + b*b + 28*b*d) +
                   d*d*(12 + 48*b*d + 141*c*c + 225*d*d)) - kurtosis
        eq <- c(eq1,eq2,eq3)
        sum(eq*eq) ## SS
    }

    out <- nlminb(start=c(1,0,0), objective=system.function,
                  scale=10,
                  control=list(trace=ifelse(verbose,1,0), rel.tol=1e-10),
                  skewness=skewness, kurtosis=kurtosis)
    if(out$convergence != 0 || out$objective > 1e-5) warning("no convergence")
    b <- out$par[1L]; c <- out$par[2L]; d <- out$par[3L]; a <- -c

    Z <- rnorm(n=n)
    Y <- a + b*Z + c*Z*Z + d*Z*Z*Z
    Y
}

ValeMaurelli1983 <- function(n=100L, COR, skewness, kurtosis, debug = FALSE) {

    fleishman1978_abcd <- function(skewness, kurtosis) {
        system.function <- function(x, skewness, kurtosis) {
            b.=x[1L]; c.=x[2L]; d.=x[3L]
            eq1 <- b.*b. + 6*b.*d. + 2*c.*c. + 15*d.*d. - 1
            eq2 <- 2*c.*(b.*b. + 24*b.*d. + 105*d.*d. + 2) - skewness
            eq3 <- 24*(b.*d. + c.*c.*(1 + b.*b. + 28*b.*d.) +
                       d.*d.*(12 + 48*b.*d. + 141*c.*c. + 225*d.*d.)) - kurtosis
            eq <- c(eq1,eq2,eq3)
            sum(eq*eq) ## SS
        }

        out <- nlminb(start=c(1,0,0), objective=system.function,
                      scale=10,
                      control=list(trace=0),
                      skewness=skewness, kurtosis=kurtosis)
        if(out$convergence != 0 || out$objective > 1e-5) {
            warning("lavaan WARNING: ValeMaurelli1983 method did not convergence, or it did not find the roots")
        }
        b. <- out$par[1L]; c. <- out$par[2L]; d. <- out$par[3L]; a. <- -c.
        c(a.,b.,c.,d.)
    }

    getICOV <- function(b1, c1, d1, b2, c2, d2, R) {
        objectiveFunction <- function(x, b1, c1, d1, b2, c2, d2, R) {
            rho=x[1L]
            eq <- rho*(b1*b2 + 3*b1*d2 + 3*d1*b2 + 9*d1*d2) +
                  rho*rho*(2*c1*c2) + rho*rho*rho*(6*d1*d2) - R
            eq*eq
        }

        #gradientFunction <- function(x, bcd1, bcd2, R) {
        #
        #}

        out <- nlminb(start=R, objective=objectiveFunction,
                      scale=10, control=list(trace=0),
                      b1=b1, c1=c1, d1=d1, b2=b2, c2=c2, d2=d2, R=R)
        if(out$convergence != 0 || out$objective > 1e-5) warning("no convergence")
        rho <- out$par[1L]
        rho
    }

    # number of variables
    nvar <- ncol(COR)
    # check skewness
    if(is.null(skewness)) {
        SK <- rep(0, nvar)
    } else if(length(skewness) == nvar) {
        SK <- skewness
    } else if(length(skewness) == 1L) {
        SK <- rep(skewness, nvar)
    } else {
        stop("skewness has wrong length")
    }

    if(is.null(kurtosis)) {
        KU <- rep(0, nvar)
    } else if(length(kurtosis) == nvar) {
        KU <- kurtosis
    } else if(length(kurtosis) == 1L) {
        KU <- rep(kurtosis, nvar)
    } else {
        stop("kurtosis has wrong length")
    }

    # create Fleishman table
    FTable <- matrix(0, nvar, 4L)
    for(i in 1:nvar) {
        FTable[i,] <- fleishman1978_abcd(skewness=SK[i], kurtosis=KU[i])
    }

    # compute intermediate correlations between all pairs
    ICOR <- diag(nvar)
    for(j in 1:(nvar-1L)) {
        for(i in (j+1):nvar) {
            if(COR[i,j] == 0) next
            ICOR[i,j] <- ICOR[j,i] <-
                getICOV(FTable[i,2], FTable[i,3], FTable[i,4],
                        FTable[j,2], FTable[j,3], FTable[j,4], R=COR[i,j])
        }
    }

    if(debug) {
         cat("\nOriginal correlations (for Vale-Maurelli):\n")
         print(COR)
         cat("\nIntermediate correlations (for Vale-Maurelli):\n")
         print(ICOR)
         cat("\nEigen values ICOR:\n")
         print( eigen(ICOR)$values )
    }

    # generate Z ## FIXME: replace by rmvnorm once we use that package
    X <- Z <- MASS::mvrnorm(n=n, mu=rep(0,nvar), Sigma=ICOR, tol = 1e-10)

    # transform Z using Fleishman constants
    for(i in 1:nvar) {
        X[,i] <- FTable[i,1L] + FTable[i,2L]*Z[,i] + FTable[i,3L]*Z[,i]*Z[,i] +
                 FTable[i,4L]*Z[,i]*Z[,i]*Z[,i]
    }

    X
}



