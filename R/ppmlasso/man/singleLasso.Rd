\name{singleLasso}
\alias{singleLasso}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
singleLasso(y, X, lamb, ob.wt = rep(1, length(y)), alpha = 1, b.init = NA, intercept = NA, family = "gaussian", tol = 1e-09, gamma = 0, init.coef = NA, mu.min = 1e-16, mu.max = 1/mu.min, area.int = FALSE, interactions, max.it = 25, standardise = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{lamb}{
%%     ~~Describe \code{lamb} here~~
}
  \item{ob.wt}{
%%     ~~Describe \code{ob.wt} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{b.init}{
%%     ~~Describe \code{b.init} here~~
}
  \item{intercept}{
%%     ~~Describe \code{intercept} here~~
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
  \item{area.int}{
%%     ~~Describe \code{area.int} here~~
}
  \item{interactions}{
%%     ~~Describe \code{interactions} here~~
}
  \item{max.it}{
%%     ~~Describe \code{max.it} here~~
}
  \item{standardise}{
%%     ~~Describe \code{standardise} here~~
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
function (y, X, lamb, ob.wt = rep(1, length(y)), alpha = 1, b.init = NA, 
    intercept = NA, family = "gaussian", tol = 1e-09, gamma = 0, 
    init.coef = NA, mu.min = 1e-16, mu.max = 1/mu.min, area.int = FALSE, 
    interactions, max.it = 25, standardise = TRUE) 
{
    error.flag = FALSE
    eta.min = family$linkfun(mu.min)
    eta.max = family$linkfun(mu.max)
    if (gamma == 0) {
        adapt.weights = rep(1, dim(X)[2])
    }
    if (gamma != 0) {
        adapt.weights = 1/abs(init.coef)^gamma
    }
    b.glm = rep(NA, dim(X)[2])
    b.lasso = b.glm
    if (length(lamb) == 1) {
        lambda.start = rep(lamb, dim(X)[2])
        lambda.start[1] = 0
        lambda = as.array(lambda.start)
        lambda = lambda * abs(adapt.weights[1:length(lambda)])
    }
    if (length(lamb) > 1) {
        lambda.start = rep(0, dim(X)[2])
        lambda = as.array(lambda.start)
        lambda = lambda * abs(adapt.weights[1:length(lambda)])
    }
    if (area.int == TRUE) {
        lambda = c(lambda, 0)
        X = cbind(X, interactions)
        is.in = rep(TRUE, dim(X)[2])
    }
    killed = is.infinite(lambda)
    X = X[, killed == FALSE]
    b.lasso = b.lasso[killed == FALSE]
    lambda = lambda[killed == FALSE]
    if (length(lamb) == 1 & lamb[1] == 0) {
        b.init = NA
        intercept = NA
    }
    if (any(is.na(b.init)) & is.na(intercept)) {
        mu.old = rep(mean(y), dim(X)[1])
        eta.old = etaFromMu(mu.old, family = family, mu.min, 
            mu.max, eta.min, eta.max)
        l.old = calculateLikelihood(X, family, ob.wt, mu = mu.old, 
            y, alpha, lambda, beta = rep(0, dim(X)[2]), penalty = FALSE)
        b.old = rep(1, dim(X)[2])
        diff = 10
        while (diff > tol) {
            update = irlsUpdate(y, X, ob.wt, is.in, signs, eta = eta.old, 
                mu = mu.old, alpha = alpha, lambda = lambda, 
                beta.old = b.old, penalty = FALSE, family = family, 
                mu.min, mu.max, eta.min, eta.max, tol)
            l.new = calculateLikelihood(X, family, ob.wt, mu = update$mu, 
                y, alpha, lambda, beta = update$beta, penalty = FALSE)
            diff = abs(l.new - l.old)
            mu.old = update$mu
            eta.old = update$eta
            b.old = update$beta
            l.old = l.new
        }
        b.glm = update$beta
        b.lasso = b.glm
    }
    if (any(is.na(b.init)) == FALSE) {
        b.lasso = b.init[killed == FALSE]
    }
    if (is.na(intercept) == FALSE) {
        if (family$family == "poisson") {
            b.lasso = c(log(mean(y)), rep(0, (dim(X)[2] - 1 - 
                area.int)), rep(1, area.int))
        }
        if (family$family == "binomial") {
            b.lasso = c(log(mean(y)/(1 - mean(y))), rep(0, (dim(X)[2] - 
                1 - area.int)), rep(1, area.int))
        }
        if (family$family == "gaussian") {
            b.lasso = c(mean(y), rep(0, (dim(X)[2] - 1 - area.int)), 
                rep(1, area.int))
        }
    }
    is.in = abs(b.lasso) > 1e-07
    sign.change = rep(-1, dim(X)[2])
    is.different = is.in
    diff = 10
    num.it = 0
    signs = sign(b.lasso)
    betas = c(b.lasso)
    scores = c()
    viols = c()
    likes = c()
    actions = c()
    eta.old = X \%*\% b.lasso
    eta.old[eta.old < eta.min] = eta.min
    eta.old[eta.old > eta.max] = eta.max
    mu.old = muFromEta(eta.old, family = family, mu.min, mu.max)
    like = calculateLikelihood(X, family, ob.wt, mu = mu.old, 
        y, alpha, lambda, beta = b.lasso, penalty = TRUE)
    likes = c(likes, like)
    if (length(lamb) == 1 & lamb[1] == 0) {
        diff = 0
        num.it = max.it + 1
    }
    while (diff > tol & num.it < max.it) {
        mult = 0
        update = irlsUpdate(y, X, ob.wt, is.in, signs, eta = eta.old, 
            mu = mu.old, alpha = alpha, lambda = lambda, beta.old = b.lasso, 
            penalty = TRUE, family = family, mu.min, mu.max, 
            eta.min, eta.max, tol)
        if (update$error == "Singular matrix") {
            error.flag = TRUE
            break
        }
        rawdiff = b.lasso[is.in] - update$beta
        sign.change = signs[is.in] * sign(update$beta)
        sign.change[1] = 1
        sign.change[sign.change == 0] = 1
        score.beta = rep(0, length(b.lasso))
        score.beta[is.in] = update$beta
        if (any(sign.change != 1) == TRUE) {
            delta = update$beta - b.lasso[is.in]
            prop = min(abs(b.lasso[is.in]/delta)[sign.change != 
                1])
            b.lasso[is.in] = b.lasso[is.in] + prop * delta
            is.in[abs(b.lasso) < 1e-07] = FALSE
            b.lasso[is.in == FALSE] = 0
            signs = sign(b.lasso)
            mult = 1
            score.beta = b.lasso
            update$eta = as.matrix(X[, is.in]) \%*\% b.lasso[is.in]
            update$eta[update$eta < eta.min] = eta.min
            update$eta[update$eta > eta.max] = eta.max
            update$mu = muFromEta(update$eta, family = family, 
                mu.min, mu.max)
        }
        score = t(as.vector(update$deriv) * t(update$s.wt)) \%*\% 
            (y - update$mu)
        score.lamb = alpha * lambda + (1 - alpha) * lambda * 
            as.vector(score.beta)
        score.lamb[lambda == 0] = 1e+05
        viol = as.vector(abs(score))/as.vector(score.lamb)
        bigviol = max(viol)
        if (any(sign.change != 1) != TRUE) {
            b.lasso[is.in] = update$beta
            signs = sign(b.lasso)
            if (bigviol > 1 + 1e-06 & is.in[viol == bigviol] == 
                FALSE) {
                is.in[viol == bigviol][1] = TRUE
                is.different[viol == bigviol][1] = TRUE
                signs[viol == bigviol][1] = sign(score[viol == 
                  bigviol])
                mult = 1
            }
        }
        like = calculateLikelihood(X, family, ob.wt, mu = update$mu, 
            y, alpha, lambda, beta = b.lasso, penalty = TRUE)
        mu.old = update$mu
        eta.old = update$eta
        for (act in 1:length(sign(b.lasso))) {
            if (sign(b.lasso)[act] == 0 & sign(as.matrix(betas)[act, 
                dim(as.matrix(betas))[2]]) != 0) {
                actions = rbind(actions, paste("Step ", dim(as.matrix(betas))[2] + 
                  1, ": Delete variable ", act, sep = ""))
            }
            if (sign(b.lasso)[act] != 0 & sign(as.matrix(betas)[act, 
                dim(as.matrix(betas))[2]]) == 0) {
                actions = rbind(actions, paste("Step ", dim(as.matrix(betas))[2] + 
                  1, ": Add variable ", act, sep = ""))
            }
        }
        betas = cbind(betas, b.lasso)
        scores = cbind(scores, score)
        viols = cbind(viols, viol)
        likes = cbind(likes, like)
        diff = mult + abs(likes[length(likes)] - likes[length(likes) - 
            1])
        num.it = num.it + 1
    }
    if (error.flag == TRUE) {
        return(list(b = NA, mu = NA, e.df = NA, deviance = NA, 
            likelihood = NA, GCV = NA, AIC = NA, BIC = NA, HQC = NA, 
            AICc = NA, ls = NA, pen.likelihood = NA, bs = NA, 
            s = NA, v = NA, actions = NA, flag = "Singular matrix"))
        cat(paste("Singular matrix error. No model fit.", "\n"))
        flush.console()
        stop
    }
    if (error.flag != TRUE) {
        eta = as.matrix(X[, is.in]) \%*\% as.matrix(b.lasso[is.in])
        eta[eta < eta.min] = eta.min
        eta[eta > eta.max] = eta.max
        mu = muFromEta(eta, family = family, mu.min, mu.max)
        if (area.int != FALSE) {
            is.in[dim(X)[2]] = FALSE
        }
        wt = irlsUpdate(y, X, ob.wt, is.in, signs, eta, mu, alpha, 
            lambda, beta.old = b.lasso, penalty = TRUE, family = family, 
            mu.min, mu.max, eta.min, eta.max, tol)$wt
        xwx = wt \%*\% X[, is.in]
        bpinv = rep(0, length(b.lasso))
        bpinv[is.in] = 1/abs(b.lasso[is.in])
        ginv = diag(bpinv)
        n.param = length(b.lasso[is.in]) - 1 + area.int
        if (sum(is.in) > 1) {
            eff.df = sum(diag(solve(xwx + diag(as.vector(lambda[is.in])) \%*\% 
                ginv[is.in, is.in]) \%*\% xwx)) + area.int
        }
        if (sum(is.in) == 1) {
            eff.df = sum(diag(solve(xwx + diag(as.matrix(lambda[is.in])) \%*\% 
                ginv[is.in, is.in]) \%*\% xwx)) + area.int
        }
        n.pres = sum(y > 0)
        unp.likelihood = calculateLikelihood(X, family, ob.wt, 
            mu, y, alpha, lambda, beta = b.lasso, penalty = FALSE)
        pen.likelihood = calculateLikelihood(X, family, ob.wt, 
            mu, y, alpha, lambda, beta = b.lasso, penalty = TRUE)
        if (family$family == "poisson") {
            dev = -2 * (sum(ob.wt * (y * log(mu) - mu)) + sum(y > 
                0))
        }
        if (family$family == "binomial") {
            dev = -2 * unp.likelihood
        }
        if (family$family == "gaussian") {
            dev = -2 * (unp.likelihood - calculateLikelihood(X, 
                family, ob.wt[y > 0], y[y > 0], y, alpha, lambda, 
                beta = b.lasso, penalty = FALSE))
        }
        if (eff.df < n.pres) {
            gcv = dev/(n.pres * (1 - (eff.df)/n.pres)^2)
        }
        if (eff.df >= n.pres) {
            gcv = NA
        }
        aic = -2 * unp.likelihood + 2 * (n.param + 1)
        bic = -2 * unp.likelihood + log(n.pres) * (n.param + 
            1)
        hqc = -2 * unp.likelihood + 2 * (n.param + 1) * log(log(n.pres))
        aicc = aic + 2 * (n.param + 1) * (n.param + 2)/(n.pres - 
            n.param - 2)
        b.out = rep(NA, length(killed))
        b.out[which(killed == FALSE)] = b.lasso
        b.out[which(killed == TRUE)] = 0
        return(list(b = b.out, mu = mu, e.df = eff.df, deviance = dev, 
            likelihood = unp.likelihood, GCV = gcv, AIC = aic, 
            BIC = bic, HQC = hqc, AICc = aicc, ls = likes, pen.likelihood = likes[length(likes)], 
            bs = betas, s = scores, v = viols, actions = actions))
    }
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
