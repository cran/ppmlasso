\name{ppmlasso}
\alias{ppmlasso}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ppmlasso(formula, sp.xy, env.grid, sp.scale, coord = c("X", "Y"), data = ppmdat(sp.xy = sp.xy, sp.scale = sp.scale, back.xy = env.grid, coord = c("X", "Y"), sp.file = NA, quad.file = NA, datfilename = "PPMDat", writefile = writefile), lamb = NA, n.fits = 200, ob.wt = NA, criterion = "bic", alpha = 1, family = "poisson", tol = 1e-09, gamma = 0, init.coef = NA, mu.min = 1e-16, mu.max = 1/mu.min, r = NA, interactions = NA, availability = NA, max.it = 25, min.lamb = -10, standardise = TRUE, n.blocks = NA, block.size = sp.scale * 100, seed = 1, writefile = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
%%     ~~Describe \code{formula} here~~
}
  \item{sp.xy}{
%%     ~~Describe \code{sp.xy} here~~
}
  \item{env.grid}{
%%     ~~Describe \code{env.grid} here~~
}
  \item{sp.scale}{
%%     ~~Describe \code{sp.scale} here~~
}
  \item{coord}{
%%     ~~Describe \code{coord} here~~
}
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{lamb}{
%%     ~~Describe \code{lamb} here~~
}
  \item{n.fits}{
%%     ~~Describe \code{n.fits} here~~
}
  \item{ob.wt}{
%%     ~~Describe \code{ob.wt} here~~
}
  \item{criterion}{
%%     ~~Describe \code{criterion} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{family}{
%%     ~~Describe \code{family} here~~
}
  \item{tol}{
%%     ~~Describe \code{tol} here~~
}
  \item{gamma}{
%%     ~~Describe \code{gamma} here~~
}
  \item{init.coef}{
%%     ~~Describe \code{init.coef} here~~
}
  \item{mu.min}{
%%     ~~Describe \code{mu.min} here~~
}
  \item{mu.max}{
%%     ~~Describe \code{mu.max} here~~
}
  \item{r}{
%%     ~~Describe \code{r} here~~
}
  \item{interactions}{
%%     ~~Describe \code{interactions} here~~
}
  \item{availability}{
%%     ~~Describe \code{availability} here~~
}
  \item{max.it}{
%%     ~~Describe \code{max.it} here~~
}
  \item{min.lamb}{
%%     ~~Describe \code{min.lamb} here~~
}
  \item{standardise}{
%%     ~~Describe \code{standardise} here~~
}
  \item{n.blocks}{
%%     ~~Describe \code{n.blocks} here~~
}
  \item{block.size}{
%%     ~~Describe \code{block.size} here~~
}
  \item{seed}{
%%     ~~Describe \code{seed} here~~
}
  \item{writefile}{
%%     ~~Describe \code{writefile} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (formula, sp.xy, env.grid, sp.scale, coord = c("X", 
    "Y"), data = ppmdat(sp.xy = sp.xy, sp.scale = sp.scale, back.xy = env.grid, 
    coord = c("X", "Y"), sp.file = NA, quad.file = NA, datfilename = "PPMDat", 
    writefile = writefile), lamb = NA, n.fits = 200, ob.wt = NA, 
    criterion = "bic", alpha = 1, family = "poisson", tol = 1e-09, 
    gamma = 0, init.coef = NA, mu.min = 1e-16, mu.max = 1/mu.min, 
    r = NA, interactions = NA, availability = NA, max.it = 25, 
    min.lamb = -10, standardise = TRUE, n.blocks = NA, block.size = sp.scale * 
        100, seed = 1, writefile = TRUE) 
{
    error.flag = FALSE
    formula.out = formula
    if (is(data, "list")) {
        use.form = as.character(formula)[2]
        cat.names = setdiff(unique(data$cat.names), NA)
        for (i in 1:length(cat.names)) {
            use.form = gsub(cat.names[i], paste(names(data$X)[which(data$cat.names == 
                cat.names[i])], collapse = " + "), use.form)
        }
        formula = as.formula(paste("~", use.form))
        data = data$X
    }
    wt.calc = FALSE
    if (is.na(ob.wt) == TRUE) {
        wt.calc = TRUE
        ob.wt = data$wt
    }
    call = match.call()
    mf = model.frame(formula, data = data)
    mt = attr(mf, "terms")
    y = data$Pres/data$wt
    X = if (!is.empty.model(mt)) 
        model.matrix(mt, mf)
    else matrix(, NROW(y), 0L)
    if (any(X[, 1] != 1)) {
        X = as.matrix(cbind(1, X))
    }
    cut.var = which(apply(X, 2, max) == 0)
    if (length(cut.var) > 0) {
        X = X[, -cut.var]
    }
    Xnames = polyNames(X)
    if (is.null(Xnames) == FALSE) {
        dimnames(X)[[2]][Xnames$ids] = Xnames$names
    }
    dimnames(X)[[2]][1] = "Intercept"
    area.int = FALSE
    raw.int = NA
    if (family == "area.inter") {
        family = "poisson"
        area.int = TRUE
        if (is.na(interactions) == TRUE) {
            interactions = pointInteractions(data, r, availability)
        }
        raw.int = interactions
    }
    s.means = NULL
    s.sds = NULL
    if (standardise == TRUE) {
        stand.X = standardiseX(X[, -1])
        X = as.matrix(cbind(1, stand.X$X))
        dimnames(X)[[2]][1] = "Intercept"
        s.means = stand.X$dat.means
        s.sds = stand.X$dat.sds
        if (area.int == TRUE) {
            interactions = standardiseX(interactions)$X
        }
    }
    if (family == "poisson") {
        family = get(family, mode = "function", envir = parent.frame())
    }
    if (is.function(family)) {
        family = family()
    }
    if (family$family == "binomial") {
        mu.max = 1 - mu.min
    }
    score.0 = scoreIntercept(y, X, ob.wt = ob.wt, area.int = area.int, 
        int = interactions, family)
    if (is.na(init.coef)) {
        gamma = 0
    }
    if (gamma == 0) {
        adapt.weights = rep(1, dim(X)[2])
    }
    if (gamma != 0) {
        adapt.weights = 1/abs(init.coef)^gamma
    }
    cv = rep(0, dim(data)[1])
    if (criterion == "blockCV") {
        cv = blocks(n.blocks, block.size, data, seed = seed)
        pred.mu = matrix(NA, dim(data)[1], n.fits)
    }
    data.all = data
    y.all = y
    X.all = X
    wt.all = ob.wt
    if (area.int == TRUE) {
        interactions.all = interactions
    }
    for (cv.i in 1:length(unique(cv))) {
        dat.test = data[cv == cv.i, ]
        data = data[cv != cv.i, ]
        data$wt = weights(data[data$Pres == 1, ], data[data$Pres == 
            0, ], coord)
        y = data$Pres/data$wt
        X = X[cv != cv.i, ]
        ob.wt = ob.wt[cv != cv.i]
        if (wt.calc == TRUE) {
            ob.wt = data$wt
        }
        if (area.int == TRUE) {
            interactions = interactions[cv != cv.i]
        }
        if (is.na(lamb) == TRUE) {
            new.score = abs(score.0/adapt.weights)
            sub.score = new.score[-1]
            max.lambda = max(sub.score[is.infinite(sub.score) == 
                FALSE])
            if (is.na(max.lambda) == FALSE) {
                lambs = sort(exp(seq(min.lamb, log(max.lambda + 
                  1e-05), length.out = n.fits)), decreasing = TRUE)
            }
        }
        if (is.na(lamb) == FALSE) {
            lambs = lamb
        }
        if (criterion == "blockCV") {
            lambda.max = max(abs(scoreIntercept(data$Pres/data$wt, 
                X, ob.wt = data$wt, area.int = area.int, int = interactions, 
                family = poisson())[-1]))
            lambs = exp(seq(0, -12, length.out = n.fits)) * lambda.max
        }
        n.pres = sum(y > 1e-08)
        n.var = dim(X)[2] - 1
        if (criterion == "msi") {
            lambs = max.lambda/sqrt(n.pres)
        }
        if (area.int == TRUE) {
            X.0 = as.matrix(cbind(X, interactions))
        }
        if (area.int != TRUE) {
            X.0 = X
        }
        mod.0 = glm(y ~ X.0[, -1], family = family, weights = ob.wt)
        coefs = matrix(NA, (dim(X)[2]), (length(lambs) + 1), 
            dimnames = list(dimnames(X)[[2]]))
        if (area.int == 1) {
            coefs = matrix(NA, (dim(X)[2] + area.int), (length(lambs) + 
                1), dimnames = list(c(dimnames(X)[[2]], "Interaction")))
        }
        num.param = rep(NA, (length(lambs) + 1))
        gcvs = rep(NA, (length(lambs) + 1))
        aics = rep(NA, (length(lambs) + 1))
        hqcs = rep(NA, (length(lambs) + 1))
        bics = rep(NA, (length(lambs) + 1))
        devs = rep(NA, (length(lambs) + 1))
        ll = rep(NA, (length(lambs) + 1))
        pll = rep(NA, (length(lambs) + 1))
        nlgcvs = rep(NA, (length(lambs) + 1))
        offset = c(log(mean(y)), rep(0, n.var), rep(1, area.int))
        if (any(is.na(adapt.weights) == FALSE)) {
            if (sum(is.infinite(adapt.weights[-1])) != (length(adapt.weights) - 
                1 - area.int)) {
                it.max = 100
                for (reg.path in 1:length(lambs)) {
                  mod = try(singleLasso(y, X, lamb = lambs[reg.path], 
                    ob.wt = ob.wt, alpha = alpha, b.init = offset, 
                    family = family, tol = 1e-09, gamma = gamma, 
                    init.coef = init.coef, area.int = area.int, 
                    interactions = interactions, max.it = it.max, 
                    standardise = FALSE), TRUE)
                  if (is(mod, "try-error")) {
                    break
                  }
                  if (any(is.na(mod$b))) {
                    break
                  }
                  coefs[, reg.path] = mod$b
                  gcvs[reg.path] = mod$GCV
                  aics[reg.path] = mod$AIC
                  hqcs[reg.path] = mod$HQC
                  bics[reg.path] = mod$BIC
                  devs[reg.path] = mod$dev
                  ll[reg.path] = mod$like
                  pll[reg.path] = mod$pen
                  num.param[reg.path] = sum(abs(mod$b) > 1e-07) - 
                    1 - area.int
                  offset = mod$b
                  it.max = max.it
                  cat(paste("Fitting Models:", reg.path, "of", 
                    length(lambs), "\n"))
                  flush.console()
                }
                num.done = length(lambs) + 1 - sum(is.na(aics))
                s.denom = sum(abs(mod.0$coef[-c(1, (n.var + 2))]))
                if (n.var > n.pres) {
                  s.denom = sum(abs(coefs[-c(1, (n.var + 2)), 
                    num.done]))
                }
                if (any(is.na(mod.0$coef)) == TRUE) {
                  s.denom = sum(abs(coefs[-c(1, (n.var + 2)), 
                    num.done]))
                }
                ss = rep(NA, (length(lambs) + 1))
                nlgcvs = rep(NA, (length(lambs) + 1))
                for (i.nlgcv in 1:num.done) {
                  s.num = sum(abs(coefs[-c(1, (n.var + 2)), i.nlgcv]))
                  s = s.num/s.denom
                  ss[i.nlgcv] = s
                  nlgcv = devs[i.nlgcv]/(n.pres * (1 - (area.int + 
                    n.var * s)/n.pres)^2)
                  if (area.int + n.var * s > n.pres) {
                    nlgcv = NA
                  }
                  nlgcvs[i.nlgcv] = nlgcv
                }
            }
            if (sum(is.infinite(adapt.weights[-1])) == (length(adapt.weights) - 
                1 - area.int)) {
                coefs = matrix(rep(init.coef, (length(lambs) + 
                  1)), length(adapt.weights), (length(lambs) + 
                  1))
                gcvs = rep(1, (length(lambs) + 1))
                aics = rep(1, (length(lambs) + 1))
                bics = rep(1, (length(lambs) + 1))
                hqcs = rep(1, (length(lambs) + 1))
                nlgcvs = rep(1, (length(lambs) + 1))
                devs = rep(1, (length(lambs) + 1))
                ll = rep(1, (length(lambs) + 1))
                num.param = rep(0, (length(lambs) + 1))
                num.done = length(lambs) + 1 - sum(is.na(aics))
            }
            criterion.matrix = data.frame(aics, bics, hqcs, gcvs, 
                nlgcvs)
            names(criterion.matrix) = c("AIC", "BIC", "HQC", 
                "GCV", "NLGCV")
            lambs[(length(lambs) + 1)] = 0
            coefs[, length(lambs)] = mod.0$coef
            if (gamma != 0) {
                coefs[, length(lambs)] = init.coef
            }
            num.param[(length(lambs) + 1)] = length(adapt.weights) - 
                sum(is.infinite(adapt.weights[-c(1, (n.var + 
                  2))])) - 1 - area.int
            meth.id = paste(criterion, "s", sep = "")
            if (criterion == "msi" | criterion == "blockCV") {
                choice.id = 1
            }
            if (criterion != "msi" & criterion != "blockCV") {
                choice.id = max(which.min(get(meth.id)))
            }
            lambda.hat = lambs[choice.id]
            beta.hat = coefs[, choice.id]
            eta.hat = X.0 \%*\% beta.hat
            mu.hat = family$linkinv(eta.hat)
            like.hat = ll[choice.id]
            assign(paste("coefs.", cv.i, sep = ""), coefs)
            if (criterion == "blockCV") {
                for (i in 1:length(lambs) - 1) {
                  fam.fit = poisson()
                  if (area.int == TRUE) {
                    fam.fit = "area.inter"
                  }
                  cv.fit = list(beta = coefs[, i], s.means = s.means, 
                    s.sds = s.sds, family = fam.fit, pt.interactions = interactions, 
                    formula = formula, mu = rep(0, dim(data)[1]))
                  if (area.int == TRUE) {
                    pred.mu[cv == cv.i, i] = predict.ppmlasso(cv.fit, 
                      newdata = dat.test, interactions = interactions.all[cv == 
                        cv.i])
                  }
                  if (area.int != TRUE) {
                    pred.mu[cv == cv.i, i] = predict.ppmlasso(cv.fit, 
                      newdata = dat.test)
                  }
                }
            }
        }
        data = data.all
        y = y.all
        X = X.all
        ob.wt = wt.all
        if (area.int == TRUE) {
            interactions = interactions.all
        }
    }
    if (criterion == "blockCV") {
        l.vec = apply(pred.mu, 2, calculateUnpenalisedLikelihood, 
            ob.wt = data$wt, y = data$Pres/data$wt)
        coef.mat = matrix(NA, n.blocks, dim(coefs)[1])
        for (i in 1:n.blocks) {
            coefs.i = get(paste("coefs.", i, sep = ""))
            coef.mat[i, ] = coefs.i[, which.max(l.vec)]
        }
        coef.av = apply(coef.mat, 2, mean, na.rm = TRUE)
        lambda.mult = exp(seq(0, -12, length.out = n.fits))[which.max(l.vec)]
        if (area.int != TRUE) {
            lambda.max = max(abs(scoreIntercept(data$Pres/data$wt, 
                X, ob.wt = data$wt, family = poisson())[-1]))
        }
        if (area.int == TRUE) {
            lambda.max = max(abs(scoreIntercept(data$Pres/data$wt, 
                X, ob.wt = data$wt, family = poisson(), area.int = TRUE, 
                int = scale(interactions))[-1]))
        }
        final.fit = try(singleLasso(y.all, X.all, lamb = lambda.mult * 
            lambda.max, ob.wt = wt.all, alpha = alpha, b.init = coef.av, 
            family = family, tol = 1e-09, gamma = gamma, init.coef = init.coef, 
            area.int = area.int, interactions = interactions, 
            max.it = it.max, standardise = FALSE), TRUE)
        coefs = final.fit$b
        lambs = lambda.mult * lambda.max
        gcvs = NA
        aics = NA
        hqcs = NA
        bics = NA
        nlgcvs = NA
        devs = final.fit$dev
        ll = final.fit$like
        pll = final.fit$pen
        num.param = sum(abs(final.fit$b) > 1e-07) - 1 - area.int
        if (area.int == TRUE) {
            X.0 = as.matrix(cbind(X, interactions))
        }
        if (area.int != TRUE) {
            X.0 = X
        }
        lambda.hat = lambs
        beta.hat = coefs
        eta.hat = X.0 \%*\% beta.hat
        mu.hat = family$linkinv(eta.hat)
        like.hat = ll
        criterion.matrix = data.frame(aics, bics, hqcs, gcvs, 
            nlgcvs)
    }
    family.out = family$family
    if (area.int == TRUE) {
        family.out = "area.inter"
    }
    output = list(betas = coefs, lambdas = lambs, likelihoods = ll, 
        pen.likelihoods = pll, lambda = lambda.hat, beta = beta.hat, 
        mu = mu.hat, likelihood = like.hat, criterion = criterion, 
        family = family.out, gamma = gamma, alpha = alpha, init.coef = init.coef, 
        criterion.matrix = criterion.matrix, data = X.0, pt.interactions = raw.int, 
        wt = ob.wt, pres = data$Pres, x = data$X, y = data$Y, 
        r = r, call = call, formula = formula.out, s.means = s.means, 
        s.sds = s.sds, cv.group = cv, n.blocks = n.blocks)
    class(output) = c("ppmlasso", "list")
    output
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
