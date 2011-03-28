HSROCSummary <-
function (data, burn_in = 0, i = NULL, Thin = 1, sub_rs = NULL, 
    point_estimate = c("median", "mean"), path = getwd(), chain = NULL, 
    tv = NULL, digit = 6, print_plot = FALSE) 
{
    setwd(path)
    if (missing(data)) 
        stop("You must provide a valid 'data' argument", call. = FALSE)
    N = length(data[, 1])
    if (burn_in < 0) {
        cat("The 'burn_in' argument must be greater or equal than zero. \n")
        stop("Please respecify and call HSROCSummary() again.\n")
    }
    if (is.null(i) == FALSE) {
        if (i < 0) {
            cat("The 'i' argument must be greater or equal than zero. \n")
            stop("Please respecify and call HSROCSummary() again.\n")
        }
    }
    if (Thin < 1) {
        cat("The 'Thin' argument must be greater or equal than 1. \n")
        stop("Please respecify and call HSROCSummary() again.\n")
    }
    if (is.null(sub_rs) == TRUE) {
        sub_rs = list(1, 1:N)
    }
    if (sub_rs[[1]] != (length(sub_rs) - 1)) {
        cat(paste("The value of the first element of 'sub_rs' (sub_rs[[1]] = ", 
            sub_rs[[1]], " ) does not match the number of remaining elements (length(sub_rs[[2:", 
            length(sub_rs), "]])) = ", length(2:length(sub_rs)), 
            "\n", sep = ""))
        stop("Please respecify and call HSROCSummary() again.\n")
    }
    if (is.logical(print_plot) == FALSE) {
        cat("The 'print_plot' argument must be a logical object. \n")
        stop("Please respecify and call HSROCSummary() again.\n")
    }
    file.prior = "Prior.information.txt"
    prior_dist_PI = "beta"
    if (is.null(chain) == FALSE) {
        setwd(chain[[1]])
        nb_chains = length(chain)
    }
    else {
        nb_chains = 1
    }
    point_estimate = match.arg(point_estimate)
    prior.exist = file.exists(file.prior)
    if (prior.exist == FALSE) {
        stop(paste("Make sure the \"", file.prior, "\" file created by the 'HSROC' function is in the \"", 
            path, "\" working directory. \n", sep = ""))
    }
    prior = read.table(file.prior, header = TRUE)
    model = read.table("model.txt", header = FALSE)
    Gold_se = (read.table("S2.txt", header = FALSE) == 1)
    Gold_sp = (read.table("C2.txt", header = FALSE) == 1)
    if (length(prior[, 1]) == 7) {
        Gold_Std = TRUE
        condInd = TRUE
    }
    else {
        if (length(prior[, 1]) > 7 & length(prior[, 1]) <= 7 + 
            2 * sub_rs[[1]]) {
            Gold_Std = FALSE
            condInd = TRUE
        }
        else {
            if (length(prior[, 1]) > 7 + 2 * sub_rs[[1]]) {
                Gold_Std = FALSE
                condInd = FALSE
            }
        }
    }
    if (is.null(tv) == FALSE) {
        real_life = FALSE
        if (sum(dim(tv[[1]])) != N + 5) {
            cat(paste("The true value for the within-study parameters were misspecified. Make sure the ordering described in the help file is preserved. \n"))
            stop("Please respecify and call HSROCSummary() again.\n")
        }
        if (length(tv[[2]]) != 7) {
            cat(paste("The true value for the between-study parameters were misspecified. Make sure the ordering described in the help file is preserved. \n"))
            stop("Please respecify and call HSROCSummary() again.\n")
        }
        if (Gold_Std == FALSE) {
            if (sum(dim(tv[[3]])) != sub_rs[[1]] + 2) {
                cat(paste("The true value for the test under evaluation were misspecified. Make sure the ordering described in the help file is preserved. \n"))
                stop("Please respecify and call HSROCSummary() again.\n")
            }
        }
    }
    else {
        real_life = TRUE
    }
    beta.a = prior[1, 1]
    beta.b = prior[1, 2]
    prior.THETA.lower = prior[2, 1]
    prior.THETA.upper = prior[2, 2]
    prior.LAMBDA.lower = prior[3, 1]
    prior.LAMBDA.upper = prior[3, 2]
    l.disp.alpha = prior[4, 1]
    u.disp.alpha = prior[4, 2]
    l.disp.theta = prior[5, 1]
    u.disp.theta = prior[5, 2]
    low.pi = prior[6, 1]
    up.pi = prior[6, 2]
    low.rj = prior[7, 1]
    up.rj = prior[7, 2]
    if (l.disp.theta == u.disp.theta) {
        SCO = TRUE
    }
    else {
        SCO = FALSE
    }
    rs.length = sub_rs[[1]]
    if (condInd == TRUE & Gold_Std == FALSE & model == 1) {
        if (Gold_se == TRUE & Gold_sp == FALSE) {
            low.sp = prior[8:(8 + (rs.length - 1)), 1]
            up.sp = prior[8:(8 + (rs.length - 1)), 2]
            prior_dist_S2 = NULL
            prior_dist_C2 = "beta"
        }
        else {
            if (Gold_sp == TRUE & Gold_se == FALSE) {
                low.se = prior[8:(8 + (rs.length - 1)), 1]
                up.se = prior[8:(8 + (rs.length - 1)), 2]
                prior_dist_C2 = NULL
                prior_dist_S2 = "beta"
            }
            else {
                low.se = prior[8:(8 + (rs.length - 1)), 1]
                up.se = prior[8:(8 + (rs.length - 1)), 2]
                low.sp = prior[(8 + rs.length):(8 + (2 * rs.length - 
                  1)), 1]
                up.sp = prior[(8 + rs.length):(8 + (2 * rs.length - 
                  1)), 2]
                prior_dist_S2 = "beta"
                prior_dist_C2 = "beta"
            }
        }
    }
    else {
        if (condInd == TRUE & Gold_Std == FALSE & model == 2) {
            mean.a1 = prior[8:(8 + (rs.length - 1)), 1]
            sd.a1 = prior[8:(8 + (rs.length - 1)), 2]
            mean.a0 = prior[(8 + rs.length):(8 + (2 * rs.length - 
                1)), 1]
            sd.a0 = prior[(8 + rs.length):(8 + (2 * rs.length - 
                1)), 2]
        }
        else {
            if (condInd == FALSE) {
                low.d1 = prior[8, 1]
                up.d1 = prior[8, 2]
                low.d0 = prior[9, 1]
                up.d0 = prior[9, 2]
                mean.a1 = prior[10:(10 + (rs.length - 1)), 1]
                sd.a1 = prior[10:(10 + (rs.length - 1)), 2]
                mean.a0 = prior[(10 + rs.length):(10 + (2 * rs.length - 
                  1)), 1]
                sd.a0 = prior[(10 + rs.length):(10 + (2 * rs.length - 
                  1)), 2]
                low.b1 = prior[(10 + (2 * rs.length)):(10 + (3 * 
                  rs.length - 1)), 1]
                up.b1 = prior[(10 + (2 * rs.length)):(10 + (3 * 
                  rs.length - 1)), 2]
                low.b0 = prior[(10 + (3 * rs.length)):(10 + (4 * 
                  rs.length - 1)), 1]
                up.b0 = prior[(10 + (3 * rs.length)):(10 + (4 * 
                  rs.length - 1)), 2]
            }
        }
    }
    if (prior_dist_PI == "beta") {
        alpha.PI = beta.parameter(low = low.pi, up = up.pi)[1, 
            ]
        beta.PI = beta.parameter(low = low.pi, up = up.pi)[2, 
            ]
    }
    else {
        if (prior_dist_PI == "uniform") {
            alpha.PI = low.pi
            beta.PI = up.pi
        }
    }
    if (model == 1) {
        if (Gold_se == TRUE) {
            Sens2.alpha = Sens2.beta = NULL
        }
        else {
            Sens2.alpha = beta.parameter(low = low.se, up = up.se)[1, 
                ]
            Sens2.beta = beta.parameter(low = low.se, up = up.se)[2, 
                ]
        }
    }
    if (model == 1) {
        if (Gold_sp == TRUE) {
            Spec2.alpha = Spec2.beta = NULL
        }
        else {
            Spec2.alpha = beta.parameter(low = low.sp, up = up.sp)[1, 
                ]
            Spec2.beta = beta.parameter(low = low.sp, up = up.sp)[2, 
                ]
        }
    }
    file.theta = "theta.txt"
    file.alpha = "alpha.txt"
    file.capital.THETA = "capital.THETA.txt"
    file.LAMBDA = "LAMBDA.txt"
    file.beta = "beta.txt"
    file.PI = "PI.txt"
    file.sigma.alpha = "sigma.alpha.txt"
    file.sigma.theta = "sigma.theta.txt"
    file.Sens2 = "Sens2.txt"
    file.Spec2 = "Spec2.txt"
    file.Sens1 = "Sens1.txt"
    file.Spec1 = "Spec1.txt"
    file.result = "estimate.txt"
    file.C_overall = "C_overall.txt"
    file.S_overall = "S_overall.txt"
    file.d1 = "d1.txt"
    file.d0 = "d0.txt"
    file.a1 = "a1.txt"
    file.a0 = "a0.txt"
    file.b1 = "b1.txt"
    file.b0 = "b0.txt"
    if (is.null(chain) == TRUE) {
        THETA = read.table(file.capital.THETA)
        LAMBDA = read.table(file.LAMBDA)
        beta = read.table(file.beta)
        PI = read.table(file.PI)
        S1 = read.table(file.Sens1)
        C1 = read.table(file.Spec1)
        sigma.alpha = read.table(file.sigma.alpha)
        alpha = read.table(file.alpha)
        S_overall = read.table(file.S_overall)
        C_overall = read.table(file.C_overall)
        if (is.null(i) == T) {
            iter.num = length(THETA[, 1])
        }
        else {
            iter.num = i
        }
        if ((iter.num - burn_in)/Thin < 100) 
            stop("You don't have enough iterations to estimate the MC error.  After taking into account the \"burn in\" and \"thinning interval\", you need at least 100 iterations to proceed.")
        total = iter.num
        q = burn_in
        alpha = alpha[(q + 1):total, ]
        THETA = THETA[(q + 1):total, 1]
        LAMBDA = LAMBDA[(q + 1):total, 1]
        beta = beta[(q + 1):total, 1]
        PI = PI[(q + 1):total, ]
        sigma.alpha = sigma.alpha[(q + 1):total, 1]
        S1 = S1[(q + 1):total, ]
        C1 = C1[(q + 1):total, ]
        C_overall = C_overall[(q + 1):total, 1]
        S_overall = S_overall[(q + 1):total, 1]
        taille = length((q + 1):total)
        thin.interval = Thin
        thin = seq(1, taille, by = thin.interval)
        alpha = as.matrix(alpha)
        alpha = alpha[thin, ]
        THETA = THETA[thin]
        LAMBDA = LAMBDA[thin]
        beta = beta[thin]
        PI = as.matrix(PI)
        PI = PI[thin, ]
        sigma.alpha = sigma.alpha[thin]
        S1 = as.matrix(S1)
        S1 = S1[thin, ]
        C1 = as.matrix(C1)
        C1 = C1[thin, ]
        C_overall = C_overall[thin]
        S_overall = S_overall[thin]
        if (SCO == FALSE) {
            sigma.theta = read.table(file.sigma.theta)
            theta = read.table(file.theta)
            theta = theta[(q + 1):total, ]
            sigma.theta = sigma.theta[(q + 1):total, 1]
            theta = as.matrix(theta)
            theta = theta[thin, ]
            sigma.theta = sigma.theta[thin]
        }
        if (condInd == TRUE & Gold_Std == FALSE & model == 1) {
            if (Gold_se == TRUE & Gold_sp == FALSE) {
                C2 = read.table(file.Spec2)
                C2 = C2[q:total, ]
                C2 = as.matrix(C2)
                C2 = C2[thin, ]
            }
            else {
                if (Gold_sp == TRUE & Gold_se == FALSE) {
                  S2 = read.table(file.Sens2)
                  S2 = S2[(q + 1):total, ]
                  S2 = as.matrix(S2)
                  S2 = S2[thin, ]
                }
                else {
                  S2 = read.table(file.Sens2)
                  C2 = read.table(file.Spec2)
                  S2 = S2[(q + 1):total, ]
                  C2 = C2[(q + 1):total, ]
                  S2 = as.matrix(S2)
                  C2 = as.matrix(C2)
                  S2 = S2[thin, ]
                  C2 = C2[thin, ]
                }
            }
        }
        else {
            if (condInd == TRUE & Gold_Std == FALSE & model == 
                2) {
                a1 = read.table(file.a1)
                a0 = read.table(file.a0)
                S2 = read.table(file.Sens2)
                C2 = read.table(file.Spec2)
                a1 = a1[(q + 1):total, ]
                a0 = a0[(q + 1):total, ]
                S2 = S2[(q + 1):total, ]
                C2 = C2[(q + 1):total, ]
                a1 = as.matrix(a1)
                a0 = as.matrix(a0)
                a1 = a1[thin, ]
                a0 = a0[thin, ]
                S2 = as.matrix(S2)
                C2 = as.matrix(C2)
                S2 = S2[thin, ]
                C2 = C2[thin, ]
            }
            else {
                if (condInd == FALSE) {
                  d1 = read.table(file.d1)
                  d0 = read.table(file.d0)
                  a1 = read.table(file.a1)
                  a0 = read.table(file.a0)
                  b1 = read.table(file.b1)
                  b0 = read.table(file.b0)
                  S2 = read.table(file.Sens2)
                  C2 = read.table(file.Spec2)
                  d1 = d1[(q + 1):total, ]
                  d0 = d0[(q + 1):total, ]
                  a1 = a1[(q + 1):total, ]
                  a0 = a0[(q + 1):total, ]
                  b1 = b1[(q + 1):total, ]
                  b0 = b0[(q + 1):total, ]
                  S2 = S2[(q + 1):total, ]
                  C2 = C2[(q + 1):total, ]
                  d1 = as.matrix(d1)
                  d0 = as.matrix(d0)
                  a1 = as.matrix(a1)
                  a0 = as.matrix(a0)
                  b1 = as.matrix(b1)
                  b0 = as.matrix(b0)
                  S2 = as.matrix(S2)
                  C2 = as.matrix(C2)
                  d1 = d1[thin, ]
                  d0 = d0[thin, ]
                  a1 = a1[thin, ]
                  a0 = a0[thin, ]
                  b1 = b1[thin, ]
                  b0 = b0[thin, ]
                  S2 = S2[thin, ]
                  C2 = C2[thin, ]
                }
            }
        }
    }
    else {
        if (is.null(chain) == FALSE) {
            K = length(chain)
            theta = alpha = THETA = LAMBDA = beta = PI = sigma.alpha = sigma.theta = S1 = C1 = S2 = C2 = S_overall = C_overall = a1 = a0 = b1 = b0 = d1 = d0 = numeric()
            for (k in 1:K) {
                setwd(chain[[k]])
                T = read.table(file.capital.THETA)
                L = read.table(file.LAMBDA)
                b = read.table(file.beta)
                p = read.table(file.PI)
                S.1 = read.table(file.Sens1)
                C.1 = read.table(file.Spec1)
                sig.a = read.table(file.sigma.alpha)
                a = read.table(file.alpha)
                S.ov = read.table(file.S_overall)
                C.ov = read.table(file.C_overall)
                if (is.null(i) == TRUE) {
                  iter.num = length(T[, 1])
                }
                else {
                  iter.num = i
                }
                if ((iter.num - burn_in)/Thin < 100) 
                  stop("You don't have enough iterations to estimate the MC error.  After taking into account the \"burn in\" and \"thinning interval\", you need at least 100 iterations to proceed.")
                total = iter.num
                q = burn_in
                a = a[(q + 1):total, ]
                T = T[(q + 1):total, 1]
                L = L[(q + 1):total, 1]
                b = b[(q + 1):total, 1]
                p = p[(q + 1):total, ]
                sig.a = sig.a[(q + 1):total, 1]
                S.1 = S.1[(q + 1):total, ]
                C.1 = C.1[(q + 1):total, ]
                C.ov = C.ov[(q + 1):total, 1]
                S.ov = S.ov[(q + 1):total, 1]
                taille = length((q + 1):total)
                thin.interval = Thin
                thin = seq(1, taille, by = thin.interval)
                a = as.matrix(a)
                a = a[thin, ]
                T = T[thin]
                L = L[thin]
                b = b[thin]
                p = as.matrix(p)
                p = p[thin, ]
                sig.a = sig.a[thin]
                S.1 = as.matrix(S.1)
                S.1 = S.1[thin, ]
                C.1 = as.matrix(C.1)
                C.1 = C.1[thin, ]
                C.ov = C.ov[thin]
                S.ov = S.ov[thin]
                alpha = rbind(alpha, as.matrix(a))
                THETA = rbind(THETA, as.matrix(T))
                LAMBDA = rbind(LAMBDA, as.matrix(L))
                beta = rbind(beta, as.matrix(b))
                PI = rbind(PI, as.matrix(p))
                sigma.alpha = rbind(sigma.alpha, as.matrix(sig.a))
                S1 = rbind(S1, as.matrix(S.1))
                C1 = rbind(C1, as.matrix(C.1))
                C_overall = rbind(C_overall, as.matrix(C.ov))
                S_overall = rbind(S_overall, as.matrix(S.ov))
                if (SCO == FALSE) {
                  t = read.table(file.theta)
                  sig.t = read.table(file.sigma.theta)
                  t = t[(q + 1):total, ]
                  sig.t = sig.t[(q + 1):total, 1]
                  t = as.matrix(t)
                  t = t[thin, ]
                  sig.t = sig.t[thin]
                  theta = rbind(theta, as.matrix(t))
                  sigma.theta = rbind(sigma.theta, as.matrix(sig.t))
                }
                if (condInd == TRUE & Gold_Std == FALSE & model == 
                  1) {
                  if (Gold_se == TRUE & Gold_sp == FALSE) {
                    C.2 = read.table(file.Spec2)
                    C.2 = C.2[(q + 1):total, ]
                    C.2 = as.matrix(C.2)
                    C.2 = C.2[thin, ]
                    C2 = rbind(C2, as.matrix(C.2))
                  }
                  else {
                    if (Gold_sp == TRUE & Gold_se == FALSE) {
                      S.2 = read.table(file.Sens2)
                      S.2 = S.2[(q + 1):total, ]
                      S.2 = as.matrix(S.2)
                      S.2 = S.2[thin, ]
                      S2 = rbind(S2, as.matrix(S.2))
                    }
                    else {
                      S.2 = read.table(file.Sens2)
                      C.2 = read.table(file.Spec2)
                      S.2 = S.2[(q + 1):total, ]
                      C.2 = C.2[(q + 1):total, ]
                      S.2 = as.matrix(S.2)
                      C.2 = as.matrix(C.2)
                      S.2 = S.2[thin, ]
                      C.2 = C.2[thin, ]
                      S2 = rbind(S2, as.matrix(S.2))
                      C2 = rbind(C2, as.matrix(C.2))
                    }
                  }
                }
                else {
                  if (condInd == TRUE & Gold_Std == FALSE & model == 
                    2) {
                    a.1 = read.table(file.a1)
                    a.0 = read.table(file.a0)
                    S.2 = read.table(file.Sens2)
                    C.2 = read.table(file.Spec2)
                    a.1 = a.1[(q + 1):total, ]
                    a.0 = a.0[(q + 1):total, ]
                    S.2 = S.2[(q + 1):total, ]
                    C.2 = C.2[(q + 1):total, ]
                    a.1 = as.matrix(a.1)
                    a.0 = as.matrix(a.0)
                    a.1 = a.1[thin, ]
                    a.0 = a.0[thin, ]
                    S.2 = as.matrix(S.2)
                    C.2 = as.matrix(C.2)
                    S.2 = S.2[thin, ]
                    C.2 = C.2[thin, ]
                    a1 = rbind(a1, as.matrix(a.1))
                    a0 = rbind(a0, as.matrix(a.0))
                    S2 = rbind(S2, as.matrix(S.2))
                    C2 = rbind(C2, as.matrix(C.2))
                  }
                  else {
                    if (condInd == FALSE) {
                      d.1 = read.table(file.d1)
                      d.0 = read.table(file.d0)
                      a.1 = read.table(file.a1)
                      a.0 = read.table(file.a0)
                      b.1 = read.table(file.b1)
                      b.0 = read.table(file.b0)
                      S.2 = read.table(file.Sens2)
                      C.2 = read.table(file.Spec2)
                      d.1 = d.1[(q + 1):total, ]
                      d.0 = d.0[(q + 1):total, ]
                      a.1 = a.1[(q + 1):total, ]
                      a.0 = a.0[(q + 1):total, ]
                      b.1 = b.1[(q + 1):total, ]
                      b.0 = b.0[(q + 1):total, ]
                      S.2 = S.2[(q + 1):total, ]
                      C.2 = C.2[(q + 1):total, ]
                      d.1 = as.matrix(d.1)
                      d.0 = as.matrix(d.0)
                      a.1 = as.matrix(a.1)
                      a.0 = as.matrix(a.0)
                      b.1 = as.matrix(b.1)
                      b.0 = as.matrix(b.0)
                      d.1 = d.1[thin, ]
                      d.0 = d.0[thin, ]
                      a.1 = a.1[thin, ]
                      a.0 = a.0[thin, ]
                      b.1 = b.1[thin, ]
                      b.0 = b.0[thin, ]
                      S.2 = as.matrix(S.2)
                      C.2 = as.matrix(C.2)
                      S.2 = S.2[thin, ]
                      C.2 = C.2[thin, ]
                      a1 = rbind(a1, as.matrix(a.1))
                      a0 = rbind(a0, as.matrix(a.0))
                      b1 = rbind(b1, as.matrix(b.1))
                      b0 = rbind(b0, as.matrix(b.0))
                      d1 = rbind(d1, as.matrix(d.1))
                      d0 = rbind(d0, as.matrix(d.0))
                      S2 = rbind(S2, as.matrix(S.2))
                      C2 = rbind(C2, as.matrix(C.2))
                    }
                  }
                }
            }
        }
    }
    setwd(path)
    alpha = as.mcmc(alpha)
    THETA = as.mcmc(THETA)
    LAMBDA = as.mcmc(LAMBDA)
    beta = as.mcmc(beta)
    PI = as.mcmc(PI)
    sigma.alpha = as.mcmc(sigma.alpha)
    S1 = as.mcmc(S1)
    C1 = as.mcmc(C1)
    C_overall = as.mcmc(C_overall)
    S_overall = as.mcmc(S_overall)
    if (SCO == FALSE) {
        theta = as.mcmc(theta)
        sigma.theta = as.mcmc(sigma.theta)
    }
    if (point_estimate == "mean") {
        mean_OR_med = 2
    }
    else {
        if (point_estimate == "median") {
            mean_OR_med = 3
        }
    }
    iter.size = length(beta)
    if (SCO == FALSE) {
        theta.est = apply(as.matrix(theta), 2, point_estimate)
        theta.HPD = HPDinterval(theta)
        theta.sd = sd(theta)
        sigma.theta.est = apply(as.matrix(sigma.theta), 2, point_estimate)
        sigma.theta.HPD = HPDinterval(sigma.theta)
        sigma.theta.sd = sd(sigma.theta)
    }
    alpha.est = apply(as.matrix(alpha), 2, point_estimate)
    alpha.HPD = HPDinterval(alpha)
    alpha.sd = sd(alpha)
    THETA.est = apply(as.matrix(THETA), 2, point_estimate)
    THETA.HPD = HPDinterval(THETA)
    THETA.sd = sd(THETA)
    LAMBDA.est = apply(as.matrix(LAMBDA), 2, point_estimate)
    LAMBDA.HPD = HPDinterval(LAMBDA)
    LAMBDA.sd = sd(LAMBDA)
    beta.est = apply(as.matrix(beta), 2, point_estimate)
    beta.HPD = HPDinterval(beta)
    beta.sd = sd(beta)
    PI.est = apply(as.matrix(PI), 2, point_estimate)
    PI.HPD = HPDinterval(PI)
    PI.sd = sd(PI)
    sigma.alpha.est = apply(as.matrix(sigma.alpha), 2, point_estimate)
    sigma.alpha.HPD = HPDinterval(sigma.alpha)
    sigma.alpha.sd = sd(sigma.alpha)
    S1.est = apply(as.matrix(S1), 2, point_estimate)
    S1.HPD = HPDinterval(S1)
    S1.sd = sd(S1)
    C1.est = apply(as.matrix(C1), 2, point_estimate)
    C1.HPD = HPDinterval(C1)
    C1.sd = sd(C1)
    C_overall.est = apply(as.matrix(C_overall), 2, point_estimate)
    C_overall.HPD = HPDinterval(C_overall)
    C_overall.sd = sd(C_overall)
    S_overall.est = apply(as.matrix(S_overall), 2, point_estimate)
    S_overall.HPD = HPDinterval(S_overall)
    S_overall.sd = sd(S_overall)
    if (condInd == TRUE & Gold_Std == FALSE & model == 1) {
        if (Gold_se == TRUE & Gold_sp == FALSE) {
            C2 = as.mcmc(C2)
            C2.est = apply(as.matrix(C2), 2, point_estimate)
            C2.HPD = HPDinterval(C2)
            C2.sd = sd(C2)
        }
        else {
            if (Gold_sp == TRUE & Gold_se == FALSE) {
                S2 = as.mcmc(S2)
                S2.est = apply(as.matrix(S2), 2, point_estimate)
                S2.HPD = HPDinterval(S2)
                S2.sd = sd(S2)
            }
            else {
                S2 = as.mcmc(S2)
                C2 = as.mcmc(C2)
                S2.est = apply(as.matrix(S2), 2, point_estimate)
                S2.HPD = HPDinterval(S2)
                S2.sd = sd(S2)
                C2.est = apply(as.matrix(C2), 2, point_estimate)
                C2.HPD = HPDinterval(C2)
                C2.sd = sd(C2)
            }
        }
    }
    else {
        if (condInd == TRUE & Gold_Std == FALSE & model == 2) {
            a1 = as.mcmc(a1)
            a0 = as.mcmc(a0)
            S2 = as.mcmc(S2)
            C2 = as.mcmc(C2)
            a1.est = apply(as.matrix(a1), 2, point_estimate)
            a1.HPD = HPDinterval(a1)
            a1.sd = sd(a1)
            a0.est = apply(as.matrix(a0), 2, point_estimate)
            a0.HPD = HPDinterval(a0)
            a0.sd = sd(a0)
            S2.est = apply(as.matrix(S2), 2, point_estimate)
            S2.HPD = HPDinterval(S2)
            S2.sd = sd(S2)
            C2.est = apply(as.matrix(C2), 2, point_estimate)
            C2.HPD = HPDinterval(C2)
            C2.sd = sd(C2)
        }
        else {
            if (condInd == FALSE) {
                a1 = as.mcmc(a1)
                a0 = as.mcmc(a0)
                b1 = as.mcmc(b1)
                b0 = as.mcmc(b0)
                d1 = as.mcmc(d1)
                d0 = as.mcmc(d0)
                S2 = as.mcmc(S2)
                C2 = as.mcmc(C2)
                a1.est = apply(as.matrix(a1), 2, point_estimate)
                a1.HPD = HPDinterval(a1)
                a1.sd = sd(a1)
                a0.est = apply(as.matrix(a0), 2, point_estimate)
                a0.HPD = HPDinterval(a0)
                a0.sd = sd(a0)
                b1.est = apply(as.matrix(b1), 2, point_estimate)
                b1.HPD = HPDinterval(b1)
                b1.sd = sd(b1)
                b0.est = apply(as.matrix(b0), 2, point_estimate)
                b0.HPD = HPDinterval(b0)
                b0.sd = sd(b0)
                d1.est = apply(as.matrix(d1), 2, point_estimate)
                d1.HPD = HPDinterval(d1)
                d1.sd = sd(d1)
                d0.est = apply(as.matrix(d0), 2, point_estimate)
                d0.HPD = HPDinterval(d0)
                d0.sd = sd(d0)
                S2.est = apply(as.matrix(S2), 2, point_estimate)
                S2.HPD = HPDinterval(S2)
                S2.sd = sd(S2)
                C2.est = apply(as.matrix(C2), 2, point_estimate)
                C2.HPD = HPDinterval(C2)
                C2.sd = sd(C2)
            }
        }
    }
    batch = 50
    ssize = iter.size/batch
    moy.t = moy.a = moy.T = moy.L = moy.b = moy.p = moy.sig.a = moy.sig.t = moy.s2 = moy.c2 = moy.s1 = moy.c1 = moy.s.ov = moy.c.ov = moy.a1 = moy.a0 = moy.b1 = moy.b0 = moy.d1 = moy.d0 = numeric()
    if (N == 1) {
        for (i in 1:batch) {
            moy.a = c(moy.a, mean(alpha[round((1 + ssize * (i - 
                1)), 0):round((ssize * i), 0)])/ssize)
            moy.p = c(moy.p, mean(PI[round((1 + ssize * (i - 
                1)), 0):round((ssize * i), 0)])/ssize)
            moy.s1 = c(moy.s1, mean(as.matrix(S1[round((1 + ssize * 
                (i - 1)), 0):round((ssize * i), 0)]))/ssize)
            moy.c1 = c(moy.c1, mean(as.matrix(C1[round((1 + ssize * 
                (i - 1)), 0):round((ssize * i), 0)]))/ssize)
            moy.L = c(moy.L, mean(LAMBDA[round((1 + ssize * (i - 
                1)), 0):round((ssize * i), 0)]))
            moy.T = c(moy.T, mean(THETA[round((1 + ssize * (i - 
                1)), 0):round((ssize * i), 0)]))
            moy.b = c(moy.b, mean(beta[round((1 + ssize * (i - 
                1)), 0):round((ssize * i), 0)]))
            moy.sig.a = c(moy.sig.a, mean(sigma.alpha[round((1 + 
                ssize * (i - 1)), 0):round((ssize * i), 0)]))
            moy.s.ov = c(moy.s.ov, mean(S_overall[round((1 + 
                ssize * (i - 1)), 0):round((ssize * i), 0)]))
            moy.c.ov = c(moy.c.ov, mean(C_overall[round((1 + 
                ssize * (i - 1)), 0):round((ssize * i), 0)]))
            if (SCO == FALSE) {
                moy.t = c(moy.t, mean(theta[round((1 + ssize * 
                  (i - 1)), 0):round((ssize * i), 0)])/ssize)
                moy.sig.t = c(moy.sig.t, mean(sigma.theta[round((1 + 
                  ssize * (i - 1)), 0):round((ssize * i), 0)]))
            }
        }
    }
    else {
        for (i in 1:batch) {
            moy.a = cbind(moy.a, colSums(alpha[round((1 + ssize * 
                (i - 1)), 0):round((ssize * i), 0), ])/ssize)
            moy.p = cbind(moy.p, colSums(PI[round((1 + ssize * 
                (i - 1)), 0):round((ssize * i), 0), ])/ssize)
            moy.s1 = cbind(moy.s1, colSums(as.matrix(S1[round((1 + 
                ssize * (i - 1)), 0):round((ssize * i), 0), ]))/ssize)
            moy.c1 = cbind(moy.c1, colSums(as.matrix(C1[round((1 + 
                ssize * (i - 1)), 0):round((ssize * i), 0), ]))/ssize)
            moy.L = c(moy.L, mean(LAMBDA[round((1 + ssize * (i - 
                1)), 0):round((ssize * i), 0)]))
            moy.T = c(moy.T, mean(THETA[round((1 + ssize * (i - 
                1)), 0):round((ssize * i), 0)]))
            moy.b = c(moy.b, mean(beta[round((1 + ssize * (i - 
                1)), 0):round((ssize * i), 0)]))
            moy.sig.a = c(moy.sig.a, mean(sigma.alpha[round((1 + 
                ssize * (i - 1)), 0):round((ssize * i), 0)]))
            moy.s.ov = c(moy.s.ov, mean(S_overall[round((1 + 
                ssize * (i - 1)), 0):round((ssize * i), 0)]))
            moy.c.ov = c(moy.c.ov, mean(C_overall[round((1 + 
                ssize * (i - 1)), 0):round((ssize * i), 0)]))
            if (SCO == FALSE) {
                moy.t = cbind(moy.t, colSums(theta[round((1 + 
                  ssize * (i - 1)), 0):round((ssize * i), 0), 
                  ])/ssize)
                moy.sig.t = c(moy.sig.t, mean(sigma.theta[round((1 + 
                  ssize * (i - 1)), 0):round((ssize * i), 0)]))
            }
        }
    }
    if (N == 1) {
        alpha.MCerror = sqrt(sum((moy.a - mean(alpha)/iter.size)^2)/(batch - 
            1))/sqrt(batch)
        PI.MCerror = sqrt(sum((moy.p - mean(PI)/iter.size)^2)/(batch - 
            1))/sqrt(batch)
        S1.MCerror = sqrt(sum((moy.s1 - mean(as.matrix(S1))/iter.size)^2)/(batch - 
            1))/sqrt(batch)
        C1.MCerror = sqrt(sum((moy.c1 - mean(as.matrix(C1))/iter.size)^2)/(batch - 
            1))/sqrt(batch)
        if (SCO == FALSE) {
            theta.MCerror = sqrt(sum((moy.t - mean(theta)/iter.size)^2)/(batch - 
                1))/sqrt(batch)
            sigma.theta.MCerror = sqrt(sum((moy.sig.t - mean(sigma.theta))^2)/(batch - 
                1))/sqrt(batch)
        }
    }
    else {
        alpha.MCerror = sqrt(rowSums((moy.a - colSums(alpha)/iter.size)^2)/(batch - 
            1))/sqrt(batch)
        PI.MCerror = sqrt(rowSums((moy.p - colSums(PI)/iter.size)^2)/(batch - 
            1))/sqrt(batch)
        S1.MCerror = sqrt(rowSums((moy.s1 - colSums(as.matrix(S1))/iter.size)^2)/(batch - 
            1))/sqrt(batch)
        C1.MCerror = sqrt(rowSums((moy.c1 - colSums(as.matrix(C1))/iter.size)^2)/(batch - 
            1))/sqrt(batch)
        if (SCO == FALSE) {
            theta.MCerror = sqrt(rowSums((moy.t - colSums(theta)/iter.size)^2)/(batch - 
                1))/sqrt(batch)
            sigma.theta.MCerror = sqrt(sum((moy.sig.t - mean(sigma.theta))^2)/(batch - 
                1))/sqrt(batch)
        }
    }
    THETA.MCerror = sqrt(sum((moy.T - mean(THETA))^2)/(batch - 
        1))/sqrt(batch)
    LAMBDA.MCerror = sqrt(sum((moy.L - mean(LAMBDA))^2)/(batch - 
        1))/sqrt(batch)
    beta.MCerror = sqrt(sum((moy.b - mean(beta))^2)/(batch - 
        1))/sqrt(batch)
    sigma.alpha.MCerror = sqrt(sum((moy.sig.a - mean(sigma.alpha))^2)/(batch - 
        1))/sqrt(batch)
    S_overall.MCerror = sqrt(sum((moy.s.ov - mean(S_overall))^2)/(batch - 
        1))/sqrt(batch)
    C_overall.MCerror = sqrt(sum((moy.c.ov - mean(C_overall))^2)/(batch - 
        1))/sqrt(batch)
    if (condInd == TRUE & Gold_Std == FALSE & model == 1) {
        if (Gold_se == TRUE & Gold_sp == FALSE) {
            if (sub_rs[[1]] == 1) {
                for (i in 1:batch) {
                  moy.c2 = c(moy.c2, mean(C2[round((1 + ssize * 
                    (i - 1)), 0):round((ssize * i), 0)]))
                }
                C2.MCerror = sqrt(sum((moy.c2 - mean(C2))^2)/(batch - 
                  1))/sqrt(batch)
            }
            else {
                if (sub_rs[[1]] > 1) {
                  for (i in 1:batch) {
                    moy.c2 = cbind(moy.c2, colSums(as.matrix(C2[round((1 + 
                      ssize * (i - 1)), 0):round((ssize * i), 
                      0), ]))/ssize)
                  }
                  C2.MCerror = sqrt(rowSums((moy.c2 - colSums(as.matrix(C2))/iter.size)^2)/(batch - 
                    1))/sqrt(batch)
                }
            }
        }
        else {
            if (Gold_sp == TRUE & Gold_se == FALSE) {
                if (sub_rs[[1]] == 1) {
                  for (i in 1:batch) {
                    moy.s2 = c(moy.s2, mean(S2[round((1 + ssize * 
                      (i - 1)), 0):round((ssize * i), 0)]))
                  }
                  S2.MCerror = sqrt(sum((moy.s2 - mean(S2))^2)/(batch - 
                    1))/sqrt(batch)
                }
                else {
                  if (sub_rs[[1]] > 1) {
                    for (i in 1:batch) {
                      moy.s2 = cbind(moy.s2, colSums(as.matrix(S2[round((1 + 
                        ssize * (i - 1)), 0):round((ssize * i), 
                        0), ]))/ssize)
                    }
                    S2.MCerror = sqrt(rowSums((moy.s2 - colSums(as.matrix(S2))/iter.size)^2)/(batch - 
                      1))/sqrt(batch)
                  }
                }
            }
            else {
                if (sub_rs[[1]] == 1) {
                  for (i in 1:batch) {
                    moy.s2 = c(moy.s2, mean(S2[round((1 + ssize * 
                      (i - 1)), 0):round((ssize * i), 0)]))
                    moy.c2 = c(moy.c2, mean(C2[round((1 + ssize * 
                      (i - 1)), 0):round((ssize * i), 0)]))
                  }
                  S2.MCerror = sqrt(sum((moy.s2 - mean(S2))^2)/(batch - 
                    1))/sqrt(batch)
                  C2.MCerror = sqrt(sum((moy.c2 - mean(C2))^2)/(batch - 
                    1))/sqrt(batch)
                }
                else {
                  if (sub_rs[[1]] > 1) {
                    for (i in 1:batch) {
                      moy.s2 = cbind(moy.s2, colSums(as.matrix(S2[round((1 + 
                        ssize * (i - 1)), 0):round((ssize * i), 
                        0), ]))/ssize)
                      moy.c2 = cbind(moy.c2, colSums(as.matrix(C2[round((1 + 
                        ssize * (i - 1)), 0):round((ssize * i), 
                        0), ]))/ssize)
                    }
                    S2.MCerror = sqrt(rowSums((moy.s2 - colSums(as.matrix(S2))/iter.size)^2)/(batch - 
                      1))/sqrt(batch)
                    C2.MCerror = sqrt(rowSums((moy.c2 - colSums(as.matrix(C2))/iter.size)^2)/(batch - 
                      1))/sqrt(batch)
                  }
                }
            }
        }
    }
    else {
        if (condInd == TRUE & Gold_Std == FALSE & model == 2) {
            if (sub_rs[[1]] == 1) {
                for (i in 1:batch) {
                  moy.a1 = c(moy.a1, mean(a1[round((1 + ssize * 
                    (i - 1)), 0):round((ssize * i), 0)]))
                  moy.a0 = c(moy.a0, mean(a0[round((1 + ssize * 
                    (i - 1)), 0):round((ssize * i), 0)]))
                  moy.s2 = c(moy.s2, mean(S2[round((1 + ssize * 
                    (i - 1)), 0):round((ssize * i), 0)]))
                  moy.c2 = c(moy.c2, mean(C2[round((1 + ssize * 
                    (i - 1)), 0):round((ssize * i), 0)]))
                }
                a1.MCerror = sqrt(sum((moy.a1 - mean(a1))^2)/(batch - 
                  1))/sqrt(batch)
                a0.MCerror = sqrt(sum((moy.a0 - mean(a0))^2)/(batch - 
                  1))/sqrt(batch)
                S2.MCerror = sqrt(sum((moy.s2 - mean(S2))^2)/(batch - 
                  1))/sqrt(batch)
                C2.MCerror = sqrt(sum((moy.c2 - mean(C2))^2)/(batch - 
                  1))/sqrt(batch)
            }
            else {
                if (sub_rs[[1]] > 1) {
                  for (i in 1:batch) {
                    moy.a1 = cbind(moy.a1, colSums(as.matrix(a1[round((1 + 
                      ssize * (i - 1)), 0):round((ssize * i), 
                      0), ]))/ssize)
                    moy.a0 = cbind(moy.a0, colSums(as.matrix(a0[round((1 + 
                      ssize * (i - 1)), 0):round((ssize * i), 
                      0), ]))/ssize)
                    moy.s2 = cbind(moy.s2, colSums(as.matrix(S2[round((1 + 
                      ssize * (i - 1)), 0):round((ssize * i), 
                      0), ]))/ssize)
                    moy.c2 = cbind(moy.c2, colSums(as.matrix(C2[round((1 + 
                      ssize * (i - 1)), 0):round((ssize * i), 
                      0), ]))/ssize)
                  }
                  a1.MCerror = sqrt(rowSums((moy.a1 - colSums(as.matrix(a1))/iter.size)^2)/(batch - 
                    1))/sqrt(batch)
                  a0.MCerror = sqrt(rowSums((moy.a0 - colSums(as.matrix(a0))/iter.size)^2)/(batch - 
                    1))/sqrt(batch)
                  S2.MCerror = sqrt(rowSums((moy.s2 - colSums(as.matrix(S2))/iter.size)^2)/(batch - 
                    1))/sqrt(batch)
                  C2.MCerror = sqrt(rowSums((moy.c2 - colSums(as.matrix(C2))/iter.size)^2)/(batch - 
                    1))/sqrt(batch)
                }
            }
        }
        else {
            if (condInd == FALSE) {
                if (sub_rs[[1]] == 1) {
                  for (i in 1:batch) {
                    moy.a1 = c(moy.a1, mean(a1[round((1 + ssize * 
                      (i - 1)), 0):round((ssize * i), 0)]))
                    moy.a0 = c(moy.a0, mean(a0[round((1 + ssize * 
                      (i - 1)), 0):round((ssize * i), 0)]))
                    moy.b1 = c(moy.b1, mean(b1[round((1 + ssize * 
                      (i - 1)), 0):round((ssize * i), 0)]))
                    moy.b0 = c(moy.b0, mean(b0[round((1 + ssize * 
                      (i - 1)), 0):round((ssize * i), 0)]))
                    moy.d1 = c(moy.d1, mean(d1[round((1 + ssize * 
                      (i - 1)), 0):round((ssize * i), 0)]))
                    moy.d0 = c(moy.d0, mean(d0[round((1 + ssize * 
                      (i - 1)), 0):round((ssize * i), 0)]))
                    moy.s2 = c(moy.s2, mean(S2[round((1 + ssize * 
                      (i - 1)), 0):round((ssize * i), 0)]))
                    moy.c2 = c(moy.c2, mean(C2[round((1 + ssize * 
                      (i - 1)), 0):round((ssize * i), 0)]))
                  }
                  a1.MCerror = sqrt(sum((moy.a1 - mean(a1))^2)/(batch - 
                    1))/sqrt(batch)
                  a0.MCerror = sqrt(sum((moy.a0 - mean(a0))^2)/(batch - 
                    1))/sqrt(batch)
                  b1.MCerror = sqrt(sum((moy.b1 - mean(b1))^2)/(batch - 
                    1))/sqrt(batch)
                  b0.MCerror = sqrt(sum((moy.b0 - mean(b0))^2)/(batch - 
                    1))/sqrt(batch)
                  d1.MCerror = sqrt(sum((moy.d1 - mean(d1))^2)/(batch - 
                    1))/sqrt(batch)
                  d0.MCerror = sqrt(sum((moy.d0 - mean(d0))^2)/(batch - 
                    1))/sqrt(batch)
                  S2.MCerror = sqrt(sum((moy.s2 - mean(S2))^2)/(batch - 
                    1))/sqrt(batch)
                  C2.MCerror = sqrt(sum((moy.c2 - mean(C2))^2)/(batch - 
                    1))/sqrt(batch)
                }
                else {
                  if (sub_rs[[1]] > 1) {
                    for (i in 1:batch) {
                      moy.a1 = cbind(moy.a1, colSums(as.matrix(a1[round((1 + 
                        ssize * (i - 1)), 0):round((ssize * i), 
                        0), ]))/ssize)
                      moy.a0 = cbind(moy.a0, colSums(as.matrix(a0[round((1 + 
                        ssize * (i - 1)), 0):round((ssize * i), 
                        0), ]))/ssize)
                      moy.b1 = cbind(moy.b1, colSums(as.matrix(b1[round((1 + 
                        ssize * (i - 1)), 0):round((ssize * i), 
                        0), ]))/ssize)
                      moy.b0 = cbind(moy.b0, colSums(as.matrix(b0[round((1 + 
                        ssize * (i - 1)), 0):round((ssize * i), 
                        0), ]))/ssize)
                      moy.d1 = cbind(moy.d1, colSums(as.matrix(d1[round((1 + 
                        ssize * (i - 1)), 0):round((ssize * i), 
                        0), ]))/ssize)
                      moy.d0 = cbind(moy.d0, colSums(as.matrix(d0[round((1 + 
                        ssize * (i - 1)), 0):round((ssize * i), 
                        0), ]))/ssize)
                      moy.s2 = cbind(moy.s2, colSums(as.matrix(S2[round((1 + 
                        ssize * (i - 1)), 0):round((ssize * i), 
                        0), ]))/ssize)
                      moy.c2 = cbind(moy.c2, colSums(as.matrix(C2[round((1 + 
                        ssize * (i - 1)), 0):round((ssize * i), 
                        0), ]))/ssize)
                    }
                    a1.MCerror = sqrt(rowSums((moy.a1 - colSums(as.matrix(a1))/iter.size)^2)/(batch - 
                      1))/sqrt(batch)
                    a0.MCerror = sqrt(rowSums((moy.a0 - colSums(as.matrix(a0))/iter.size)^2)/(batch - 
                      1))/sqrt(batch)
                    b1.MCerror = sqrt(rowSums((moy.b1 - colSums(as.matrix(b1))/iter.size)^2)/(batch - 
                      1))/sqrt(batch)
                    b0.MCerror = sqrt(rowSums((moy.b0 - colSums(as.matrix(b0))/iter.size)^2)/(batch - 
                      1))/sqrt(batch)
                    d1.MCerror = sqrt(rowSums((moy.d1 - colSums(as.matrix(d1))/iter.size)^2)/(batch - 
                      1))/sqrt(batch)
                    d0.MCerror = sqrt(rowSums((moy.d0 - colSums(as.matrix(d0))/iter.size)^2)/(batch - 
                      1))/sqrt(batch)
                    S2.MCerror = sqrt(rowSums((moy.s2 - colSums(as.matrix(S2))/iter.size)^2)/(batch - 
                      1))/sqrt(batch)
                    C2.MCerror = sqrt(rowSums((moy.c2 - colSums(as.matrix(C2))/iter.size)^2)/(batch - 
                      1))/sqrt(batch)
                  }
                }
            }
        }
    }
    data = list(data)
    if (real_life == FALSE) {
        pp = data[[1]][, 1]
        pn = data[[1]][, 2]
        np = data[[1]][, 3]
        nn = data[[1]][, 4]
        Sample.size = pp + pn + np + nn
        if (SCO == FALSE) {
            true.alpha = tv[[1]][, 1]
            true.theta = tv[[1]][, 2]
            true.S1 = tv[[1]][, 3]
            true.C1 = tv[[1]][, 4]
            true.PI = tv[[1]][, 5]
            true.THETA = tv[[2]][1]
            true.sigma.theta = tv[[2]][2]
            true.LAMBDA = tv[[2]][3]
            true.sigma.alpha = tv[[2]][4]
            true.beta = tv[[2]][5]
            if (Gold_Std == TRUE) {
                true.S_overall = tv[[2]][6]
                true.C_overall = tv[[2]][7]
            }
            else {
                if (condInd == TRUE & Gold_Std == FALSE & model == 
                  1) {
                  true.S2 = tv[[3]][1, ]
                  true.C2 = tv[[3]][2, ]
                  true.S_overall = tv[[2]][6]
                  true.C_overall = tv[[2]][7]
                }
                else {
                  if (condInd == TRUE & Gold_Std == FALSE & model == 
                    2) {
                    true.a1 = tv[[3]][3, ]
                    true.a0 = tv[[3]][4, ]
                    true.S2 = tv[[3]][1, ]
                    true.C2 = tv[[3]][2, ]
                    true.S_overall = tv[[2]][8]
                    true.C_overall = tv[[2]][9]
                  }
                  else {
                    if (condInd == FALSE) {
                      true.S2 = tv[[3]][1, ]
                      true.C2 = tv[[3]][2, ]
                      true.a1 = tv[[3]][3, ]
                      true.a0 = tv[[3]][4, ]
                      true.b1 = tv[[3]][5, ]
                      true.b0 = tv[[3]][6, ]
                      true.d1 = tv[[2]][6]
                      true.d0 = tv[[2]][7]
                      true.S_overall = tv[[2]][8]
                      true.C_overall = tv[[2]][9]
                    }
                  }
                }
            }
        }
        else {
            if (SCO == TRUE) {
                true.alpha = tv[[1]][, 1]
                true.S1 = tv[[1]][, 2]
                true.C1 = tv[[1]][, 3]
                true.PI = tv[[1]][, 4]
                true.THETA = tv[[2]][1]
                true.LAMBDA = tv[[2]][2]
                true.sigma.alpha = tv[[2]][3]
                true.beta = tv[[2]][4]
                if (Gold_Std == TRUE) {
                  true.S_overall = tv[[2]][5]
                  true.C_overall = tv[[2]][6]
                }
                else {
                  if (condInd == TRUE & Gold_Std == FALSE & model == 
                    1) {
                    true.S2 = tv[[3]][1, ]
                    true.C2 = tv[[3]][2, ]
                    true.S_overall = tv[[2]][5]
                    true.C_overall = tv[[2]][6]
                  }
                  else {
                    if (condInd == TRUE & Gold_Std == FALSE & 
                      model == 2) {
                      true.a1 = tv[[3]][3, ]
                      true.a0 = tv[[3]][4, ]
                      true.S2 = tv[[3]][1, ]
                      true.C2 = tv[[3]][2, ]
                      true.S_overall = tv[[2]][7]
                      true.C_overall = tv[[2]][8]
                    }
                    else {
                      if (condInd == FALSE) {
                        true.S2 = tv[[3]][1, ]
                        true.C2 = tv[[3]][2, ]
                        true.a1 = tv[[3]][3, ]
                        true.a0 = tv[[3]][4, ]
                        true.b1 = tv[[3]][5, ]
                        true.b0 = tv[[3]][6, ]
                        true.d1 = tv[[2]][6]
                        true.d0 = tv[[2]][7]
                        true.S_overall = tv[[2]][7]
                        true.C_overall = tv[[2]][8]
                      }
                    }
                  }
                }
            }
        }
        test.file = paste("Summary for N =", round((iter.num * 
            nb_chains - (burn_in) * nb_chains)/Thin, 0), ".txt")
        write(paste("______________________________________________________"), 
            file = test.file, append = TRUE)
        write(paste("Number of chains =", nb_chains), file = test.file, 
            append = TRUE)
        write(paste("Number of iteration within a chain =", iter.num, 
            "    Burn in within each chain =", burn_in), file = test.file, 
            append = TRUE)
        write(paste("Thinning interval =", Thin), file = test.file, 
            append = TRUE)
        write(paste("Total number of iteration kept =", round((iter.num * 
            nb_chains - (burn_in) * nb_chains)/Thin, 0)), file = test.file, 
            append = TRUE)
        write(paste(""), file = test.file, append = TRUE)
        write(paste("File location : ", path), file = test.file, 
            append = TRUE)
        write(paste(""), file = test.file, append = TRUE)
        write(paste("Date :", Sys.time()), file = test.file, 
            append = TRUE)
        write(paste("______________________________________________________"), 
            file = test.file, append = TRUE)
        write(paste(""), file = test.file, append = TRUE)
        if (Gold_Std == TRUE) {
            write("Perfect reference standard", file = test.file, 
                append = TRUE)
        }
        else {
            write("Imperfect reference standard", file = test.file, 
                append = TRUE)
        }
        write(paste(""), file = test.file, append = TRUE)
        write(paste("______________________________________________________"), 
            file = test.file, append = TRUE)
        write(paste("\tSAMPLE SIZE \t "), file = test.file, append = TRUE)
        write(paste("______________________________________________________"), 
            file = test.file, append = TRUE)
        write(paste(""), file = test.file, append = TRUE)
        write(paste("         Total ++ +- -+ --"), file = test.file, 
            append = TRUE)
        for (i in 1:N) {
            write(paste("Study ", i, "", Sample.size[i], "", 
                pp[i], "", pn[i], "", np[i], "", nn[i]), file = test.file, 
                append = TRUE)
        }
        write(paste(""), file = test.file, append = TRUE)
        write(paste("______________________________________________________"), 
            file = test.file, append = TRUE)
        write(paste("\tPRIOR INFORMATION \t "), file = test.file, 
            append = TRUE)
        write(paste("______________________________________________________"), 
            file = test.file, append = TRUE)
        write(paste(""), file = test.file, append = TRUE)
        write(paste(""), file = test.file, append = TRUE)
        write(paste(""), file = test.file, append = TRUE)
        write(paste("Prior of prevalence (pi) is ", prior_dist_PI, 
            "(", round(alpha.PI, digits = 4), ",", round(beta.PI, 
                digits = 4), "), <=> pi in [", low.pi, ",", up.pi, 
            "]"), file = test.file, append = TRUE)
        write(paste(""), file = test.file, append = TRUE)
        write(paste("Prior of beta is Uniform(", round(beta.a, 
            4), ",", round(beta.b, 4), ")"), file = test.file, 
            append = TRUE)
        write(paste("Prior of THETA is Uniform(", prior.THETA.lower, 
            ",", prior.THETA.upper, ")"), file = test.file, append = TRUE)
        write(paste("Prior of LAMBDA is Uniform(", prior.LAMBDA.lower, 
            ",", prior.LAMBDA.upper, ")"), file = test.file, 
            append = TRUE)
        write(paste("Prior of sigma_alpha is uniform(", l.disp.alpha, 
            ",", u.disp.alpha, ")"), file = test.file, append = TRUE)
        if (SCO == FALSE) {
            write(paste("Prior of sigma_theta is uniform(", l.disp.theta, 
                ",", u.disp.theta, ")"), file = test.file, append = TRUE)
        }
        if (condInd == TRUE & Gold_Std == FALSE & model == 1) {
            write(paste(""), file = test.file, append = TRUE)
            if (Gold_se == TRUE & Gold_sp == FALSE) {
                write(paste("Prior of S2 (Sensitivity of reference test) is "), 
                  file = test.file, append = TRUE)
                for (i in 1:rs.length) {
                  write(paste("Study(ies) ", sub_rs[[i + 1]][1], 
                    "to", sub_rs[[i + 1]][length(sub_rs[[i + 
                      1]])], "assumed to be perfect."), file = test.file, 
                    append = TRUE)
                }
                write(paste(""), file = test.file, append = TRUE)
                write(paste("Prior of C2 (Specificity of reference test) is "), 
                  file = test.file, append = TRUE)
                for (i in 1:rs.length) {
                  write(paste("Study(ies) ", sub_rs[[i + 1]][1], 
                    "to", sub_rs[[i + 1]][length(sub_rs[[i + 
                      1]])], "", prior_dist_C2, "(", round(Spec2.alpha[i], 
                      digits = 4), ",", round(Spec2.beta[i], 
                      digits = 4), "), <=> C2 in [", low.sp[i], 
                    ",", up.sp[i], "]"), file = test.file, append = TRUE)
                }
            }
            else {
                if (Gold_sp == TRUE & Gold_se == FALSE) {
                  write(paste("Prior of S2 (Sensitivity of reference test) is "), 
                    file = test.file, append = TRUE)
                  for (i in 1:rs.length) {
                    write(paste("Study(ies) ", sub_rs[[i + 1]][1], 
                      "to", sub_rs[[i + 1]][length(sub_rs[[i + 
                        1]])], "", prior_dist_S2, "(", round(Sens2.alpha[i], 
                        digits = 4), ",", round(Sens2.beta[i], 
                        digits = 4), "), <=> S2 in [", low.se[i], 
                      ",", up.se[i], "]"), file = test.file, 
                      append = TRUE)
                  }
                  write(paste(""), file = test.file, append = TRUE)
                  write(paste("Prior of C2 (Specificity of reference test) is "), 
                    file = test.file, append = TRUE)
                  for (i in 1:rs.length) {
                    write(paste("Study(ies) ", sub_rs[[i + 1]][1], 
                      "to", sub_rs[[i + 1]][length(sub_rs[[i + 
                        1]])], "assumed to be perfect."), file = test.file, 
                      append = TRUE)
                  }
                }
                else {
                  write(paste("Prior of S2 (Sensitivity of reference test) is "), 
                    file = test.file, append = TRUE)
                  for (i in 1:rs.length) {
                    write(paste("Study(ies) ", sub_rs[[i + 1]][1], 
                      "to", sub_rs[[i + 1]][length(sub_rs[[i + 
                        1]])], "", prior_dist_S2, "(", round(Sens2.alpha[i], 
                        digits = 4), ",", round(Sens2.beta[i], 
                        digits = 4), "), <=> S2 in [", low.se[i], 
                      ",", up.se[i], "]"), file = test.file, 
                      append = TRUE)
                  }
                  write(paste(""), file = test.file, append = TRUE)
                  write(paste("Prior of C2 (Specificity of reference test) is "), 
                    file = test.file, append = TRUE)
                  for (i in 1:rs.length) {
                    write(paste("Study(ies) ", sub_rs[[i + 1]][1], 
                      "to", sub_rs[[i + 1]][length(sub_rs[[i + 
                        1]])], "", prior_dist_C2, "(", round(Spec2.alpha[i], 
                        digits = 4), ",", round(Spec2.beta[i], 
                        digits = 4), "), <=> C2 in [", low.sp[i], 
                      ",", up.sp[i], "]"), file = test.file, 
                      append = TRUE)
                  }
                }
            }
        }
        else {
            if (condInd == TRUE & Gold_Std == FALSE & model == 
                2) {
                write(paste(""), file = test.file, append = TRUE)
                write(paste("Prior of a1 is "), file = test.file, 
                  append = TRUE)
                for (i in 1:rs.length) {
                  write(paste("Study(ies) ", sub_rs[[i + 1]][1], 
                    "to", sub_rs[[i + 1]][length(sub_rs[[i + 
                      1]])], "Normal (", round(mean.a1[i], digits = 4), 
                    ",", round(sd.a1[i], digits = 4), "), <=> S2 in []"), 
                    file = test.file, append = TRUE)
                }
                write(paste(""), file = test.file, append = TRUE)
                write(paste("Prior of a0 is "), file = test.file, 
                  append = TRUE)
                for (i in 1:rs.length) {
                  write(paste("Study(ies) ", sub_rs[[i + 1]][1], 
                    "to", sub_rs[[i + 1]][length(sub_rs[[i + 
                      1]])], "Normal (", round(mean.a0[i], digits = 4), 
                    ",", round(sd.a0[i], digits = 4), "), <=> C2 in []"), 
                    file = test.file, append = TRUE)
                }
            }
            else {
                if (condInd == FALSE) {
                  write(paste(""), file = test.file, append = TRUE)
                  write(paste("Prior of a1 is "), file = test.file, 
                    append = TRUE)
                  for (i in 1:rs.length) {
                    write(paste("Study(ies) ", sub_rs[[i + 1]][1], 
                      "to", sub_rs[[i + 1]][length(sub_rs[[i + 
                        1]])], "Normal (", round(mean.a1[i], 
                        digits = 4), ",", round(sd.a1[i], digits = 4), 
                      "), <=> S2 in []"), file = test.file, append = TRUE)
                  }
                  write(paste(""), file = test.file, append = TRUE)
                  write(paste("Prior of a0 is "), file = test.file, 
                    append = TRUE)
                  for (i in 1:rs.length) {
                    write(paste("Study(ies) ", sub_rs[[i + 1]][1], 
                      "to", sub_rs[[i + 1]][length(sub_rs[[i + 
                        1]])], "Normal (", round(mean.a0[i], 
                        digits = 4), ",", round(sd.a0[i], digits = 4), 
                      "), <=> C2 in []"), file = test.file, append = TRUE)
                  }
                  write(paste(""), file = test.file, append = TRUE)
                  write(paste("Prior of b1 is "), file = test.file, 
                    append = TRUE)
                  for (i in 1:rs.length) {
                    write(paste("Study(ies) ", sub_rs[[i + 1]][1], 
                      "to", sub_rs[[i + 1]][length(sub_rs[[i + 
                        1]])], "Uniform (", round(low.b1[i], 
                        digits = 4), ",", round(up.b1[i], digits = 4), 
                      "), <=> S2 in []"), file = test.file, append = TRUE)
                  }
                  write(paste(""), file = test.file, append = TRUE)
                  write(paste("Prior of b0 is "), file = test.file, 
                    append = TRUE)
                  for (i in 1:rs.length) {
                    write(paste("Study(ies) ", sub_rs[[i + 1]][1], 
                      "to", sub_rs[[i + 1]][length(sub_rs[[i + 
                        1]])], "Uniform (", round(low.b0[i], 
                        digits = 4), ",", round(up.b0[i], digits = 4), 
                      "), <=> C2 in []"), file = test.file, append = TRUE)
                  }
                  write(paste(""), file = test.file, append = TRUE)
                  write(paste("Prior of d1 is Uniform(", round(low.d1, 
                    4), ",", round(up.d1, 4), ")"), file = test.file, 
                    append = TRUE)
                  write(paste("Prior of d0 is Uniform(", round(low.d0, 
                    4), ",", round(up.d0, 4), ")"), file = test.file, 
                    append = TRUE)
                }
            }
        }
        write(paste(), file = test.file, append = TRUE)
        write(paste("______________________________________________________"), 
            file = test.file, append = TRUE)
        write(paste("______________________________________________________"), 
            file = test.file, append = TRUE)
        write(paste("\tBETWEEN-STUDY parameters (Point estimate =", 
            point_estimate, ")\t "), file = test.file, append = TRUE)
        write(paste("______________________________________________________"), 
            file = test.file, append = TRUE)
        write(paste("______________________________________________________"), 
            file = test.file, append = TRUE)
        write(paste(""), file = test.file, append = TRUE)
        write(paste("         True_value Estimate Standard_dev MC_error C.I._lower C.I._upper"), 
            file = test.file, append = TRUE)
        write(paste(""), file = test.file, append = TRUE)
        write(paste("THETA       ", round(true.THETA, digits = digit), 
            "", round(THETA.est, digits = digit), "", round(THETA.sd, 
                digits = digit), "", round(THETA.MCerror, digits = digit), 
            "", round(THETA.HPD[1], digits = digit), "", round(THETA.HPD[2], 
                digits = digit)), file = test.file, append = TRUE)
        write(paste("LAMBDA      ", round(true.LAMBDA, digits = digit), 
            "", round(LAMBDA.est, digits = digit), "", round(LAMBDA.sd, 
                digits = digit), "", round(LAMBDA.MCerror, digits = digit), 
            "", round(LAMBDA.HPD[1], digits = digit), "", round(LAMBDA.HPD[2], 
                digits = digit)), file = test.file, append = TRUE)
        write(paste("beta        ", round(true.beta, digits = digit), 
            "", round(beta.est, digits = digit), "", round(beta.sd, 
                digits = digit), "", round(beta.MCerror, digits = digit), 
            "", round(beta.HPD[1], digits = digit), "", round(beta.HPD[2], 
                digits = digit)), file = test.file, append = TRUE)
        write(paste("sigma.alpha ", round(true.sigma.alpha, digits = digit), 
            "", round(sigma.alpha.est, digits = digit), "", round(sigma.alpha.sd, 
                digits = digit), "", round(sigma.alpha.MCerror, 
                digits = digit), "", round(sigma.alpha.HPD[1], 
                digits = digit), "", round(sigma.alpha.HPD[2], 
                digits = digit)), file = test.file, append = TRUE)
        if (SCO == FALSE) {
            write(paste("sigma.theta ", round(true.sigma.theta, 
                digits = digit), "", round(sigma.theta.est, digits = digit), 
                "", round(sigma.theta.sd, digits = digit), "", 
                round(sigma.theta.MCerror, digits = digit), "", 
                round(sigma.theta.HPD[1], digits = digit), "", 
                round(sigma.theta.HPD[2], digits = digit)), file = test.file, 
                append = TRUE)
        }
        write(paste("S overall   ", round(true.S_overall, digits = digit), 
            "", round(S_overall.est, digits = digit), "", round(S_overall.sd, 
                digits = digit), "", round(S_overall.MCerror, 
                digits = digit), "", round(S_overall.HPD[1], 
                digits = digit), "", round(S_overall.HPD[2], 
                digits = digit)), file = test.file, append = TRUE)
        write(paste("C overall   ", round(true.C_overall, digits = digit), 
            "", round(C_overall.est, digits = digit), "", round(C_overall.sd, 
                digits = digit), "", round(C_overall.MCerror, 
                digits = digit), "", round(C_overall.HPD[1], 
                digits = digit), "", round(C_overall.HPD[2], 
                digits = digit)), file = test.file, append = TRUE)
        if (condInd == TRUE & Gold_Std == FALSE & model == 1) {
            write(paste(""), file = test.file, append = TRUE)
            write(paste("______________________________________________________"), 
                file = test.file, append = TRUE)
            write(paste("\tReference standard (Point estimate =", 
                point_estimate, ")\t "), file = test.file, append = TRUE)
            write(paste("______________________________________________________"), 
                file = test.file, append = TRUE)
            write(paste(""), file = test.file, append = TRUE)
            write(paste("         True_value Estimate Standard_dev MC_error C.I._lower C.I._upper"), 
                file = test.file, append = TRUE)
            write(paste(""), file = test.file, append = TRUE)
            if (Gold_se == TRUE & Gold_sp == FALSE) {
                if (rs.length != 1) {
                  for (i in 1:rs.length) {
                    write(paste("S2 of Study(ies) ", sub_rs[[i + 
                      1]][1], "to", sub_rs[[i + 1]][length(sub_rs[[i + 
                      1]])], "", round(true.S2[i], digits = digit), 
                      "(It was assumed to be perfect)"), file = test.file, 
                      append = TRUE)
                    write(paste("C2 of Study(ies) ", sub_rs[[i + 
                      1]][1], "to", sub_rs[[i + 1]][length(sub_rs[[i + 
                      1]])], "", round(true.C2[i], digits = digit), 
                      "", round(C2.est[i], digits = digit), "", 
                      round(C2.sd[i], digits = digit), "", round(C2.MCerror[i], 
                        digits = digit), "", round(C2.HPD[i, 
                        1], digits = digit), "", round(C2.HPD[i, 
                        2], digits = digit)), file = test.file, 
                      append = TRUE)
                  }
                }
                else {
                  write(paste("S2       ", round(true.S2, digits = digit), 
                    "(It was assumed to be perfect)"), file = test.file, 
                    append = TRUE)
                  write(paste("C2       ", round(true.C2, digits = digit), 
                    "", round(C2.est, digits = digit), "", round(C2.sd, 
                      digits = digit), "", round(C2.MCerror, 
                      digits = digit), "", round(C2.HPD[1], digits = digit), 
                    "", round(C2.HPD[2], digits = digit)), file = test.file, 
                    append = TRUE)
                }
            }
            else {
                if (Gold_sp == TRUE & Gold_se == FALSE) {
                  if (rs.length != 1) {
                    for (i in 1:rs.length) {
                      write(paste("S2 of Study(ies) ", sub_rs[[i + 
                        1]][1], "to", sub_rs[[i + 1]][length(sub_rs[[i + 
                        1]])], "", round(true.S2[i], digits = digit), 
                        "", round(S2.est[i], digits = digit), 
                        "", round(S2.sd[i], digits = digit), 
                        "", round(S2.MCerror[i], digits = digit), 
                        "", round(S2.HPD[i, 1], digits = digit), 
                        "", round(S2.HPD[i, 2], digits = digit)), 
                        file = test.file, append = TRUE)
                      write(paste("C2 of Study(ies) ", sub_rs[[i + 
                        1]][1], "to", sub_rs[[i + 1]][length(sub_rs[[i + 
                        1]])], "", round(true.C2[i], digits = digit), 
                        "(It was assumed to be perfect)"), file = test.file, 
                        append = TRUE)
                    }
                  }
                  else {
                    write(paste("S2       ", round(true.S2, digits = digit), 
                      "", round(S2.est, digits = digit), "", 
                      round(S2.sd, digits = digit), "", round(S2.MCerror, 
                        digits = digit), "", round(S2.HPD[1], 
                        digits = digit), "", round(S2.HPD[2], 
                        digits = digit)), file = test.file, append = TRUE)
                    write(paste("C2       ", round(true.C2, digits = digit), 
                      "(It was assumed to be perfect)"), file = test.file, 
                      append = TRUE)
                  }
                }
                else {
                  if (rs.length != 1) {
                    for (i in 1:rs.length) {
                      write(paste("S2 of Study(ies) ", sub_rs[[i + 
                        1]][1], "to", sub_rs[[i + 1]][length(sub_rs[[i + 
                        1]])], "", round(true.S2[i], digits = digit), 
                        "", round(S2.est[i], digits = digit), 
                        "", round(S2.sd[i], digits = digit), 
                        "", round(S2.MCerror[i], digits = digit), 
                        "", round(S2.HPD[i, 1], digits = digit), 
                        "", round(S2.HPD[i, 2], digits = digit)), 
                        file = test.file, append = TRUE)
                    }
                    for (i in 1:rs.length) {
                      write(paste("C2 of Study(ies) ", sub_rs[[i + 
                        1]][1], "to", sub_rs[[i + 1]][length(sub_rs[[i + 
                        1]])], "", round(true.C2[i], digits = digit), 
                        "", round(C2.est[i], digits = digit), 
                        "", round(C2.sd[i], digits = digit), 
                        "", round(C2.MCerror[i], digits = digit), 
                        "", round(C2.HPD[i, 1], digits = digit), 
                        "", round(C2.HPD[i, 2], digits = digit)), 
                        file = test.file, append = TRUE)
                    }
                  }
                  else {
                    write(paste("S2       ", round(true.S2, digits = digit), 
                      "", round(S2.est, digits = digit), "", 
                      round(S2.sd, digits = digit), "", round(S2.MCerror, 
                        digits = digit), "", round(S2.HPD[1], 
                        digits = digit), "", round(S2.HPD[2], 
                        digits = digit)), file = test.file, append = TRUE)
                    write(paste("C2       ", round(true.C2, digits = digit), 
                      "", round(C2.est, digits = digit), "", 
                      round(C2.sd, digits = digit), "", round(C2.MCerror, 
                        digits = digit), "", round(C2.HPD[1], 
                        digits = digit), "", round(C2.HPD[2], 
                        digits = digit)), file = test.file, append = TRUE)
                  }
                }
            }
        }
        else {
            if (condInd == TRUE & Gold_Std == FALSE & model == 
                2) {
                write(paste(""), file = test.file, append = TRUE)
                write(paste("______________________________________________________"), 
                  file = test.file, append = TRUE)
                write(paste("\tReference standard (Point estimate =", 
                  point_estimate, ")\t "), file = test.file, 
                  append = TRUE)
                write(paste("______________________________________________________"), 
                  file = test.file, append = TRUE)
                write(paste(""), file = test.file, append = TRUE)
                write(paste("         True_value Estimate Standard_dev MC_error C.I._lower C.I._upper"), 
                  file = test.file, append = TRUE)
                write(paste(""), file = test.file, append = TRUE)
                if (rs.length != 1) {
                  for (i in 1:rs.length) {
                    write(paste("a1 of Study(ies) ", sub_rs[[i + 
                      1]][1], "to", sub_rs[[i + 1]][length(sub_rs[[i + 
                      1]])], "", round(true.a1[i], digits = digit), 
                      "", round(a1.est[i], digits = digit), "", 
                      round(a1.sd[i], digits = digit), "", round(a1.MCerror[i], 
                        digits = digit), "", round(a1.HPD[i, 
                        1], digits = digit), "", round(a1.HPD[i, 
                        2], digits = digit)), file = test.file, 
                      append = TRUE)
                  }
                  for (i in 1:rs.length) {
                    write(paste("a0 of Study(ies) ", sub_rs[[i + 
                      1]][1], "to", sub_rs[[i + 1]][length(sub_rs[[i + 
                      1]])], "", round(true.a0[i], digits = digit), 
                      "", round(a0.est[i], digits = digit), "", 
                      round(a0.sd[i], digits = digit), "", round(a0.MCerror[i], 
                        digits = digit), "", round(a0.HPD[i, 
                        1], digits = digit), "", round(a0.HPD[i, 
                        2], digits = digit)), file = test.file, 
                      append = TRUE)
                  }
                  write(paste(""), file = test.file, append = TRUE)
                  for (i in 1:rs.length) {
                    write(paste("S2 of Study(ies) ", sub_rs[[i + 
                      1]][1], "to", sub_rs[[i + 1]][length(sub_rs[[i + 
                      1]])], "", round(true.S2[i], digits = digit), 
                      "", round(S2.est[i], digits = digit), "", 
                      round(S2.sd[i], digits = digit), "", round(S2.MCerror[i], 
                        digits = digit), "", round(S2.HPD[i, 
                        1], digits = digit), "", round(S2.HPD[i, 
                        2], digits = digit)), file = test.file, 
                      append = TRUE)
                  }
                  for (i in 1:rs.length) {
                    write(paste("C2 of Study(ies) ", sub_rs[[i + 
                      1]][1], "to", sub_rs[[i + 1]][length(sub_rs[[i + 
                      1]])], "", round(true.C2[i], digits = digit), 
                      "", round(C2.est[i], digits = digit), "", 
                      round(C2.sd[i], digits = digit), "", round(C2.MCerror[i], 
                        digits = digit), "", round(C2.HPD[i, 
                        1], digits = digit), "", round(C2.HPD[i, 
                        2], digits = digit)), file = test.file, 
                      append = TRUE)
                  }
                }
                else {
                  write(paste("a1       ", round(true.a1, digits = digit), 
                    "", round(a1.est, digits = digit), "", round(a1.sd, 
                      digits = digit), "", round(a1.MCerror, 
                      digits = digit), "", round(a1.HPD[1], digits = digit), 
                    "", round(a1.HPD[2], digits = digit)), file = test.file, 
                    append = TRUE)
                  write(paste("a0       ", round(true.a0, digits = digit), 
                    "", round(a0.est, digits = digit), "", round(a0.sd, 
                      digits = digit), "", round(a0.MCerror, 
                      digits = digit), "", round(a0.HPD[1], digits = digit), 
                    "", round(a0.HPD[2], digits = digit)), file = test.file, 
                    append = TRUE)
                  write(paste(""), file = test.file, append = TRUE)
                  write(paste("S2       ", round(true.S2, digits = digit), 
                    "", round(S2.est, digits = digit), "", round(S2.sd, 
                      digits = digit), "", round(S2.MCerror, 
                      digits = digit), "", round(S2.HPD[1], digits = digit), 
                    "", round(S2.HPD[2], digits = digit)), file = test.file, 
                    append = TRUE)
                  write(paste("C2       ", round(true.C2, digits = digit), 
                    "", round(C2.est, digits = digit), "", round(C2.sd, 
                      digits = digit), "", round(C2.MCerror, 
                      digits = digit), "", round(C2.HPD[1], digits = digit), 
                    "", round(C2.HPD[2], digits = digit)), file = test.file, 
                    append = TRUE)
                }
            }
            else {
                if (condInd == FALSE) {
                  write(paste(""), file = test.file, append = TRUE)
                  write(paste("______________________________________________________"), 
                    file = test.file, append = TRUE)
                  write(paste("\tReference standard (Point estimate =", 
                    point_estimate, ")\t "), file = test.file, 
                    append = TRUE)
                  write(paste("______________________________________________________"), 
                    file = test.file, append = TRUE)
                  write(paste(""), file = test.file, append = TRUE)
                  write(paste("         True_value Estimate Standard_dev MC_error C.I._lower C.I._upper"), 
                    file = test.file, append = TRUE)
                  write(paste(""), file = test.file, append = TRUE)
                  write(paste("d1       ", round(true.d1, digits = digit), 
                    "", round(d1.est, digits = digit), "", round(d1.sd, 
                      digits = digit), "", round(d1.MCerror, 
                      digits = digit), "", round(d1.HPD[1], digits = digit), 
                    "", round(d1.HPD[2], digits = digit)), file = test.file, 
                    append = TRUE)
                  write(paste("d0       ", round(true.d0, digits = digit), 
                    "", round(d0.est, digits = digit), "", round(d0.sd, 
                      digits = digit), "", round(d0.MCerror, 
                      digits = digit), "", round(d0.HPD[1], digits = digit), 
                    "", round(d0.HPD[2], digits = digit)), file = test.file, 
                    append = TRUE)
                  if (rs.length != 1) {
                    for (i in 1:rs.length) {
                      write(paste("a1 of Study(ies) ", sub_rs[[i + 
                        1]][1], "to", sub_rs[[i + 1]][length(sub_rs[[i + 
                        1]])], "", round(true.a1[i], digits = digit), 
                        "", round(a1.est[i], digits = digit), 
                        "", round(a1.sd[i], digits = digit), 
                        "", round(a1.MCerror[i], digits = digit), 
                        "", round(a1.HPD[i, 1], digits = digit), 
                        "", round(a1.HPD[i, 2], digits = digit)), 
                        file = test.file, append = TRUE)
                    }
                    for (i in 1:rs.length) {
                      write(paste("a0 of Study(ies) ", sub_rs[[i + 
                        1]][1], "to", sub_rs[[i + 1]][length(sub_rs[[i + 
                        1]])], "", round(true.a0[i], digits = digit), 
                        "", round(a0.est[i], digits = digit), 
                        "", round(a0.sd[i], digits = digit), 
                        "", round(a0.MCerror[i], digits = digit), 
                        "", round(a0.HPD[i, 1], digits = digit), 
                        "", round(a0.HPD[i, 2], digits = digit)), 
                        file = test.file, append = TRUE)
                    }
                    write(paste(""), file = test.file, append = TRUE)
                    for (i in 1:rs.length) {
                      write(paste("S2 of Study(ies) ", sub_rs[[i + 
                        1]][1], "to", sub_rs[[i + 1]][length(sub_rs[[i + 
                        1]])], "", round(true.S2[i], digits = digit), 
                        "", round(S2.est[i], digits = digit), 
                        "", round(S2.sd[i], digits = digit), 
                        "", round(S2.MCerror[i], digits = digit), 
                        "", round(S2.HPD[i, 1], digits = digit), 
                        "", round(S2.HPD[i, 2], digits = digit)), 
                        file = test.file, append = TRUE)
                    }
                    for (i in 1:rs.length) {
                      write(paste("C2 of Study(ies) ", sub_rs[[i + 
                        1]][1], "to", sub_rs[[i + 1]][length(sub_rs[[i + 
                        1]])], "", round(true.C2[i], digits = digit), 
                        "", round(C2.est[i], digits = digit), 
                        "", round(C2.sd[i], digits = digit), 
                        "", round(C2.MCerror[i], digits = digit), 
                        "", round(C2.HPD[i, 1], digits = digit), 
                        "", round(C2.HPD[i, 2], digits = digit)), 
                        file = test.file, append = TRUE)
                    }
                    for (i in 1:rs.length) {
                      write(paste("b1 of Study(ies) ", sub_rs[[i + 
                        1]][1], "to", sub_rs[[i + 1]][length(sub_rs[[i + 
                        1]])], "", round(true.b1[i], digits = digit), 
                        "", round(b1.est[i], digits = digit), 
                        "", round(b1.sd[i], digits = digit), 
                        "", round(b1.MCerror[i], digits = digit), 
                        "", round(b1.HPD[i, 1], digits = digit), 
                        "", round(b1.HPD[i, 2], digits = digit)), 
                        file = test.file, append = TRUE)
                    }
                    for (i in 1:rs.length) {
                      write(paste("b0 of Study(ies) ", sub_rs[[i + 
                        1]][1], "to", sub_rs[[i + 1]][length(sub_rs[[i + 
                        1]])], "", round(true.b0[i], digits = digit), 
                        "", round(b0.est[i], digits = digit), 
                        "", round(b0.sd[i], digits = digit), 
                        "", round(b0.MCerror[i], digits = digit), 
                        "", round(b0.HPD[i, 1], digits = digit), 
                        "", round(b0.HPD[i, 2], digits = digit)), 
                        file = test.file, append = TRUE)
                    }
                  }
                  else {
                    write(paste("a1       ", round(true.a1, digits = digit), 
                      "", round(a1.est, digits = digit), "", 
                      round(a1.sd, digits = digit), "", round(a1.MCerror, 
                        digits = digit), "", round(a1.HPD[1], 
                        digits = digit), "", round(a1.HPD[2], 
                        digits = digit)), file = test.file, append = TRUE)
                    write(paste("a0       ", round(true.a0, digits = digit), 
                      "", round(a0.est, digits = digit), "", 
                      round(a0.sd, digits = digit), "", round(a0.MCerror, 
                        digits = digit), "", round(a0.HPD[1], 
                        digits = digit), "", round(a0.HPD[2], 
                        digits = digit)), file = test.file, append = TRUE)
                    write(paste("b1       ", round(true.b1, digits = digit), 
                      "", round(b1.est, digits = digit), "", 
                      round(b1.sd, digits = digit), "", round(b1.MCerror, 
                        digits = digit), "", round(b1.HPD[1], 
                        digits = digit), "", round(b1.HPD[2], 
                        digits = digit)), file = test.file, append = TRUE)
                    write(paste("b0       ", round(true.b0, digits = digit), 
                      "", round(b0.est, digits = digit), "", 
                      round(b0.sd, digits = digit), "", round(b0.MCerror, 
                        digits = digit), "", round(b0.HPD[1], 
                        digits = digit), "", round(b0.HPD[2], 
                        digits = digit)), file = test.file, append = TRUE)
                    write(paste("S2       ", round(true.S2, digits = digit), 
                      "", round(S2.est, digits = digit), "", 
                      round(S2.sd, digits = digit), "", round(S2.MCerror, 
                        digits = digit), "", round(S2.HPD[1], 
                        digits = digit), "", round(S2.HPD[2], 
                        digits = digit)), file = test.file, append = TRUE)
                    write(paste("C2       ", round(true.C2, digits = digit), 
                      "", round(C2.est, digits = digit), "", 
                      round(C2.sd, digits = digit), "", round(C2.MCerror, 
                        digits = digit), "", round(C2.HPD[1], 
                        digits = digit), "", round(C2.HPD[2], 
                        digits = digit)), file = test.file, append = TRUE)
                  }
                }
            }
        }
        write(paste(""), file = test.file, append = TRUE)
        write(paste("______________________________________________________"), 
            file = test.file, append = TRUE)
        write(paste("______________________________________________________"), 
            file = test.file, append = TRUE)
        write(paste("\tWITHIN-STUDY PARAMETERS \t "), file = test.file, 
            append = TRUE)
        write(paste("______________________________________________________"), 
            file = test.file, append = TRUE)
        write(paste("______________________________________________________"), 
            file = test.file, append = TRUE)
        if (SCO == FALSE) {
            write(paste("______________________________________________________"), 
                file = test.file, append = TRUE)
            write(paste("\ttheta \t "), file = test.file, append = TRUE)
            write(paste("______________________________________________________"), 
                file = test.file, append = TRUE)
            write(paste(""), file = test.file, append = TRUE)
            write(paste("         True_value Estimate Standard_dev MC_error C.I._lower C.I._upper"), 
                file = test.file, append = TRUE)
            write(paste(""), file = test.file, append = TRUE)
            for (i in 1:N) {
                write(paste("Study ", i, "", round(true.theta[i], 
                  digits = digit), "", round(theta.est[i], digits = digit), 
                  "", round(theta.sd[i], digits = digit), "", 
                  round(theta.MCerror[i], digits = digit), "", 
                  round(theta.HPD[i, 1], digits = digit), "", 
                  round(theta.HPD[i, 2], digits = digit)), file = test.file, 
                  append = TRUE)
            }
        }
        write(paste("______________________________________________________"), 
            file = test.file, append = TRUE)
        write(paste("\talpha \t "), file = test.file, append = TRUE)
        write(paste("______________________________________________________"), 
            file = test.file, append = TRUE)
        write(paste(""), file = test.file, append = TRUE)
        write(paste("         True_value Estimate Standard_Dev MC_error C.I._lower C.I._upper"), 
            file = test.file, append = TRUE)
        write(paste(""), file = test.file, append = TRUE)
        for (i in 1:N) {
            write(paste("Study ", i, "", round(true.alpha[i], 
                digits = digit), "", round(alpha.est[i], digits = digit), 
                "", round(alpha.sd[i], digits = digit), "", round(alpha.MCerror[i], 
                  digits = digit), "", round(alpha.HPD[i, 1], 
                  digits = digit), "", round(alpha.HPD[i, 2], 
                  digits = digit)), file = test.file, append = TRUE)
        }
        write(paste("______________________________________________________"), 
            file = test.file, append = TRUE)
        write(paste("\tPrevalence  \t "), file = test.file, append = TRUE)
        write(paste("______________________________________________________"), 
            file = test.file, append = TRUE)
        write(paste(""), file = test.file, append = TRUE)
        write(paste("         True_value Estimate Standard_dev MC_error C.I._lower C.I._upper"), 
            file = test.file, append = TRUE)
        write(paste(""), file = test.file, append = TRUE)
        for (i in 1:N) {
            write(paste("Study ", i, "", round(true.PI[i], digits = digit), 
                "", round(PI.est[i], digits = digit), "", round(PI.sd[i], 
                  digits = digit), "", round(PI.MCerror[i], digits = digit), 
                "", round(PI.HPD[i, 1], digits = digit), "", 
                round(PI.HPD[i, 2], digits = digit)), file = test.file, 
                append = TRUE)
        }
        write(paste("______________________________________________________"), 
            file = test.file, append = TRUE)
        write(paste("\tSensitivity of test 1 (S1) \t "), file = test.file, 
            append = TRUE)
        write(paste("______________________________________________________"), 
            file = test.file, append = TRUE)
        write(paste(""), file = test.file, append = TRUE)
        write(paste("         True_value Estimate Standard_dev MC_error C.I._lower C.I._upper"), 
            file = test.file, append = TRUE)
        write(paste(""), file = test.file, append = TRUE)
        for (i in 1:N) {
            write(paste("Study ", i, "", round(true.S1[i], digits = digit), 
                "", round(S1.est[i], digits = digit), "", round(S1.sd[i], 
                  digits = digit), "", round(S1.MCerror[i], digits = digit), 
                "", round(S1.HPD[i, 1], digits = digit), "", 
                round(S1.HPD[i, 2], digits = digit)), file = test.file, 
                append = TRUE)
        }
        write(paste("______________________________________________________"), 
            file = test.file, append = TRUE)
        write(paste("\tSpecificity of test 1 (C1) \t "), file = test.file, 
            append = TRUE)
        write(paste("______________________________________________________"), 
            file = test.file, append = TRUE)
        write(paste(""), file = test.file, append = TRUE)
        write(paste("         True_value Estimate Standard_dev MC_error C.I._lower C.I._upper"), 
            file = test.file, append = TRUE)
        write(paste(""), file = test.file, append = TRUE)
        for (i in 1:N) {
            write(paste("Study ", i, "", round(true.C1[i], digits = digit), 
                "", round(C1.est[i], digits = digit), "", round(C1.sd[i], 
                  digits = digit), "", round(C1.MCerror[i], digits = digit), 
                "", round(C1.HPD[i, 1], digits = digit), "", 
                round(C1.HPD[i, 2], digits = digit)), file = test.file, 
                append = TRUE)
        }
        Num_study = c()
        for (i in 1:N) {
            Num_study = c(Num_study, paste("Study", i))
        }
        if (condInd == TRUE & Gold_Std == FALSE & model == 1) {
            if (Gold_se == TRUE & Gold_sp == FALSE) {
                if (rs.length != 1) {
                  refstd_parameters = array(0, dim = c(rs.length, 
                    4, 1), dimnames = list(1:rs.length, c("True Value", 
                    paste(point_estimate, "estimate"), "HPD lower", 
                    "HPD upper"), "C2"))
                  refstd_parameters[, 1, 1] <- true.C2
                  refstd_parameters[, 2, 1] <- C2.est
                  refstd_parameters[, 3, 1] <- C2.HPD[, 1]
                  refstd_parameters[, 4, 1] <- C2.HPD[, 2]
                  long = length(alpha[, 1])
                  refstd_Parameters = array(0, c(long, rs.length, 
                    1))
                  refstd_Parameters[, , 1] <- C2
                  if (print_plot == TRUE) {
                    if (is.null(chain) == FALSE) {
                      file.pdf_RS = paste("RefStd trace plots for N =", 
                        round((iter.num * nb_chains - (burn_in) * 
                          nb_chains)/Thin, 0), ".pdf")
                      pdf(file.pdf_RS, paper = "a4", height = 20)
                      param = "C2"
                      par(mfcol = c(5, 2))
                      no_chains = length(chain)
                      iter_chain = round((iter.num * nb_chains - 
                        (burn_in) * nb_chains)/Thin, 0)/no_chains
                      longueur = 1:iter_chain
                      for (i in 1:rs.length) {
                        plot(x = longueur, y = refstd_Parameters[longueur, 
                          i, 1], type = "n", col = 1, ylab = paste(param, 
                          " of reference test ", i), xlab = "iteration number", 
                          main = paste("Thinning interval = ", 
                            thin.interval, "\n Total samplesize kept = ", 
                            (iter.num * nb_chains - (burn_in) * 
                              nb_chains)/Thin))
                        for (l in 1:length(chain)) {
                          lines(x = longueur, y = refstd_Parameters[l * 
                            longueur, i, 1], col = l)
                        }
                      }
                    }
                    else {
                      file.pdf_RS = paste("RefStd trace plots for N =", 
                        round((iter.num * nb_chains - (burn_in) * 
                          nb_chains)/Thin, 0), ".pdf")
                      pdf(file.pdf_RS, paper = "a4", height = 20)
                      param = "C2"
                      par(mfcol = c(5, 2))
                      longueur = 1:long
                      for (i in 1:rs.length) {
                        plot(x = longueur, y = refstd_Parameters[, 
                          i, 1], type = "l", col = "grey", ylab = paste(param, 
                          " of study ", i), xlab = "iteration number", 
                          main = paste("Thinning interval = ", 
                            thin.interval, "\n Total samplesize kept = ", 
                            (iter.num * nb_chains - (burn_in) * 
                              nb_chains)/Thin))
                      }
                    }
                    dev.off()
                    file.pdf_RS2 = paste("RefStd density plots for N =", 
                      round((iter.num * nb_chains - (burn_in) * 
                        nb_chains)/Thin, 0), ".pdf")
                    pdf(file.pdf_RS2, paper = "a4", height = 20)
                    param = "C2"
                    par(mfcol = c(5, 2))
                    longueur = 1:long
                    for (i in 1:rs.length) {
                      plot(density(refstd_Parameters[, i, 1]), 
                        lwd = 4, type = "l", col = "grey", main = paste(param, 
                          " of study ", i, " \n Thinning interval = ", 
                          thin.interval, "\n Total samplesize kept = ", 
                          (iter.num * nb_chains - (burn_in) * 
                            nb_chains)/Thin))
                    }
                    dev.off()
                  }
                }
                else {
                  refstd_parameters = matrix(0, ncol = 4, nrow = 1)
                  rownames(refstd_parameters) = c("C2")
                  colnames(refstd_parameters) = c("True Value", 
                    paste(point_estimate, "estimate"), "HPD.low", 
                    "HPD.high")
                  refstd_parameters[1, 1] <- true.C2
                  refstd_parameters[1, 2] <- C2.est
                  refstd_parameters[1, 3] <- C2.HPD[1]
                  refstd_parameters[1, 4] <- C2.HPD[2]
                  long = length(THETA)
                  refstd_Parameters = matrix(0, nrow = long, 
                    ncol = 1)
                  colnames(refstd_Parameters) = c("C2")
                  refstd_Parameters[, 1] <- C2
                  if (print_plot == TRUE) {
                    if (is.null(chain) == FALSE) {
                      file.pdf_RS = paste("RefStd trace plots for N =", 
                        round((iter.num * nb_chains - (burn_in) * 
                          nb_chains)/Thin, 0), ".pdf")
                      pdf(file.pdf_RS, paper = "a4", height = 20)
                      param = c("C2")
                      par(mfcol = c(5, 2))
                      no_chains = length(chain)
                      iter_chain = round((iter.num * nb_chains - 
                        (burn_in) * nb_chains)/Thin, 0)/no_chains
                      longueur = 1:iter_chain
                      plot(x = longueur, y = refstd_Parameters[longueur, 
                        1], type = "n", col = 1, ylab = paste(param), 
                        xlab = "iteration number", main = paste("Thinning interval = ", 
                          thin.interval, "\n Total samplesize kept = ", 
                          (iter.num * nb_chains - (burn_in) * 
                            nb_chains)/Thin))
                      for (l in 1:length(chain)) {
                        lines(x = longueur, y = refstd_Parameters[l * 
                          longueur, 1], col = l)
                      }
                    }
                    else {
                      file.pdf_RS = paste("RefStd trace plots for N =", 
                        round((iter.num * nb_chains - (burn_in) * 
                          nb_chains)/Thin, 0), ".pdf")
                      pdf(file.pdf_RS, paper = "a4", height = 20)
                      param = c("C2")
                      par(mfcol = c(5, 2))
                      longueur = 1:long
                      plot(x = longueur, y = refstd_Parameters[, 
                        1], type = "l", col = "grey", ylab = paste(param), 
                        xlab = "iteration number", main = paste("Thinning interval = ", 
                          thin.interval, "\n Total samplesize kept = ", 
                          (iter.num * nb_chains - (burn_in) * 
                            nb_chains)/Thin))
                    }
                    dev.off()
                    file.pdf_RS2 = paste("RefStd density plots for N =", 
                      round((iter.num * nb_chains - (burn_in) * 
                        nb_chains)/Thin, 0), ".pdf")
                    pdf(file.pdf_RS2, paper = "a4", height = 20)
                    param = c("C2")
                    par(mfcol = c(5, 2))
                    longueur = 1:long
                    plot(density(refstd_Parameters[, 1]), lwd = 4, 
                      type = "l", col = "grey", main = paste(param, 
                        " of reference standard \n Thinning interval = ", 
                        thin.interval, "\n Total samplesize kept = ", 
                        (iter.num * nb_chains - (burn_in) * nb_chains)/Thin))
                    dev.off()
                  }
                }
            }
            else {
                if (Gold_sp == TRUE & Gold_se == FALSE) {
                  if (rs.length != 1) {
                    refstd_parameters = array(0, dim = c(rs.length, 
                      4, 1), dimnames = list(1:rs.length, c("True Value", 
                      paste(point_estimate, "estimate"), "HPD lower", 
                      "HPD upper"), "S2"))
                    refstd_parameters[, 1, 1] <- true.S2
                    refstd_parameters[, 2, 1] <- S2.est
                    refstd_parameters[, 3, 1] <- S2.HPD[, 1]
                    refstd_parameters[, 4, 1] <- S2.HPD[, 2]
                    long = length(alpha[, 1])
                    refstd_Parameters = array(0, c(long, rs.length, 
                      1))
                    refstd_Parameters[, , 1] <- S2
                    if (print_plot == TRUE) {
                      if (is.null(chain) == FALSE) {
                        file.pdf_RS = paste("RefStd trace plots for N =", 
                          round((iter.num * nb_chains - (burn_in) * 
                            nb_chains)/Thin, 0), ".pdf")
                        pdf(file.pdf_RS, paper = "a4", height = 20)
                        param = c("S2")
                        par(mfcol = c(5, 2))
                        no_chains = length(chain)
                        iter_chain = round((iter.num * nb_chains - 
                          (burn_in) * nb_chains)/Thin, 0)/no_chains
                        longueur = 1:iter_chain
                        for (i in 1:rs.length) {
                          plot(x = longueur, y = refstd_Parameters[longueur, 
                            i, 1], type = "n", col = 1, ylab = paste(param, 
                            " of reference test ", i), xlab = "iteration number", 
                            main = paste("Thinning interval = ", 
                              thin.interval, "\n Total samplesize kept = ", 
                              (iter.num * nb_chains - (burn_in) * 
                                nb_chains)/Thin))
                          for (l in 1:length(chain)) {
                            lines(x = longueur, y = refstd_Parameters[l * 
                              longueur, i, 1], col = l)
                          }
                        }
                      }
                      else {
                        file.pdf_RS = paste("RefStd trace plots for N =", 
                          round((iter.num * nb_chains - (burn_in) * 
                            nb_chains)/Thin, 0), ".pdf")
                        pdf(file.pdf_RS, paper = "a4", height = 20)
                        param = c("S2")
                        par(mfcol = c(5, 2))
                        longueur = 1:long
                        for (i in 1:rs.length) {
                          plot(x = longueur, y = refstd_Parameters[, 
                            i, 1], type = "l", col = "grey", 
                            ylab = paste(param, " of reference test ", 
                              i), xlab = "iteration number", 
                            main = paste("Thinning interval = ", 
                              thin.interval, "\n Total samplesize kept = ", 
                              (iter.num * nb_chains - (burn_in) * 
                                nb_chains)/Thin))
                        }
                      }
                      dev.off()
                      file.pdf_RS2 = paste("RefStd density plots for N =", 
                        round((iter.num * nb_chains - (burn_in) * 
                          nb_chains)/Thin, 0), ".pdf")
                      pdf(file.pdf_RS2, paper = "a4", height = 20)
                      param = c("S2")
                      par(mfcol = c(5, 2))
                      longueur = 1:long
                      for (i in 1:rs.length) {
                        plot(density(refstd_Parameters[, i, 1]), 
                          lwd = 4, type = "l", col = "grey", 
                          main = paste(param, " of reference standard ", 
                            i, " \n Thinning interval = ", thin.interval, 
                            "\n Total samplesize kept = ", (iter.num * 
                              nb_chains - (burn_in) * nb_chains)/Thin))
                      }
                      dev.off()
                    }
                  }
                  else {
                    refstd_parameters = matrix(0, ncol = 4, nrow = 1)
                    rownames(refstd_parameters) = c("S2")
                    colnames(refstd_parameters) = c("True Value", 
                      paste(point_estimate, "estimate"), "HPD.low", 
                      "HPD.high")
                    refstd_parameters[1, 1] <- true.S2
                    refstd_parameters[1, 2] <- S2.est
                    refstd_parameters[1, 3] <- S2.HPD[1]
                    refstd_parameters[1, 4] <- S2.HPD[2]
                    long = length(THETA)
                    refstd_Parameters = matrix(0, nrow = long, 
                      ncol = 1)
                    colnames(refstd_Parameters) = c("S2")
                    refstd_Parameters[, 1] <- S2
                    if (print_plot == TRUE) {
                      if (is.null(chain) == FALSE) {
                        file.pdf_RS = paste("RefStd trace plots for N =", 
                          round((iter.num * nb_chains - (burn_in) * 
                            nb_chains)/Thin, 0), ".pdf")
                        pdf(file.pdf_RS, paper = "a4", height = 20)
                        param = c("S2")
                        par(mfcol = c(5, 2))
                        no_chains = length(chain)
                        iter_chain = round((iter.num * nb_chains - 
                          (burn_in) * nb_chains)/Thin, 0)/no_chains
                        longueur = 1:iter_chain
                        plot(x = longueur, y = refstd_Parameters[longueur, 
                          1], type = "n", col = 1, ylab = paste(param), 
                          xlab = "iteration number", main = paste("Thinning interval = ", 
                            thin.interval, "\n Total samplesize kept = ", 
                            (iter.num * nb_chains - (burn_in) * 
                              nb_chains)/Thin))
                        for (l in 1:length(chain)) {
                          lines(x = longueur, y = refstd_Parameters[l * 
                            longueur, 1], col = l)
                        }
                      }
                      else {
                        file.pdf_RS = paste("RefStd trace plots for N =", 
                          round((iter.num * nb_chains - (burn_in) * 
                            nb_chains)/Thin, 0), ".pdf")
                        pdf(file.pdf_RS, paper = "a4", height = 20)
                        param = c("S2")
                        par(mfcol = c(5, 2))
                        longueur = 1:long
                        plot(x = longueur, y = refstd_Parameters[, 
                          1], type = "l", col = "grey", ylab = paste(param), 
                          xlab = "iteration number", main = paste("Thinning interval = ", 
                            thin.interval, "\n Total samplesize kept = ", 
                            (iter.num * nb_chains - (burn_in) * 
                              nb_chains)/Thin))
                      }
                      dev.off()
                      file.pdf_RS2 = paste("RefStd density plots for N =", 
                        round((iter.num * nb_chains - (burn_in) * 
                          nb_chains)/Thin, 0), ".pdf")
                      pdf(file.pdf_RS2, paper = "a4", height = 20)
                      param = c("S2")
                      par(mfcol = c(5, 2))
                      longueur = 1:long
                      plot(density(refstd_Parameters[, 1]), lwd = 4, 
                        type = "l", col = "grey", main = paste(param, 
                          " \n Thinning interval = ", thin.interval, 
                          "\n Total samplesize kept = ", (iter.num * 
                            nb_chains - (burn_in) * nb_chains)/Thin))
                      dev.off()
                    }
                  }
                }
                else {
                  if (rs.length != 1) {
                    refstd_parameters = array(0, dim = c(rs.length, 
                      4, 2), dimnames = list(1:rs.length, c("True Value", 
                      paste(point_estimate, "estimate"), "HPD lower", 
                      "HPD upper"), c("S2", "C2")))
                    refstd_parameters[, 1, 1] <- true.S2
                    refstd_parameters[, 2, 1] <- S2.est
                    refstd_parameters[, 3, 1] <- S2.HPD[, 1]
                    refstd_parameters[, 4, 1] <- S2.HPD[, 2]
                    refstd_parameters[, 1, 2] <- true.C2
                    refstd_parameters[, 2, 2] <- C2.est
                    refstd_parameters[, 3, 2] <- C2.HPD[, 1]
                    refstd_parameters[, 4, 2] <- C2.HPD[, 2]
                    long = length(alpha[, 1])
                    refstd_Parameters = array(0, c(long, rs.length, 
                      2))
                    refstd_Parameters[, , 1] <- S2
                    refstd_Parameters[, , 2] <- C2
                    if (print_plot == TRUE) {
                      if (is.null(chain) == FALSE) {
                        file.pdf_RS = paste("RefStd trace plots for N =", 
                          round((iter.num * nb_chains - (burn_in) * 
                            nb_chains)/Thin, 0), ".pdf")
                        pdf(file.pdf_RS, paper = "a4", height = 20)
                        param = c("S2", "C2")
                        par(mfcol = c(5, 2))
                        no_chains = length(chain)
                        iter_chain = round((iter.num * nb_chains - 
                          (burn_in) * nb_chains)/Thin, 0)/no_chains
                        longueur = 1:iter_chain
                        for (j in 1:2) {
                          for (i in 1:rs.length) {
                            plot(x = longueur, y = refstd_Parameters[longueur, 
                              i, j], type = "n", col = 1, ylab = paste(param[j], 
                              " of reference test ", i), xlab = "iteration number", 
                              main = paste("Thinning interval = ", 
                                thin.interval, "\n Total samplesize kept = ", 
                                (iter.num * nb_chains - (burn_in) * 
                                  nb_chains)/Thin))
                            for (l in 1:length(chain)) {
                              lines(x = longueur, y = refstd_Parameters[l * 
                                longueur, i, j], col = l)
                            }
                          }
                        }
                      }
                      else {
                        file.pdf_RS = paste("RefStd trace plots for N =", 
                          round((iter.num * nb_chains - (burn_in) * 
                            nb_chains)/Thin, 0), ".pdf")
                        pdf(file.pdf_RS, paper = "a4", height = 20)
                        param = c("S2", "C2")
                        par(mfcol = c(5, 2))
                        for (j in 1:2) {
                          longueur = 1:long
                          for (i in 1:rs.length) {
                            plot(x = longueur, y = refstd_Parameters[, 
                              i, j], type = "l", col = "grey", 
                              ylab = paste(param[j], " of study ", 
                                sub_rs[[i + 1]]), xlab = "iteration number", 
                              main = paste("Thinning interval = ", 
                                thin.interval, "\n Total samplesize kept = ", 
                                (iter.num * nb_chains - (burn_in) * 
                                  nb_chains)/Thin))
                          }
                        }
                      }
                      dev.off()
                      file.pdf_RS2 = paste("RefStd density plots for N =", 
                        round((iter.num * nb_chains - (burn_in) * 
                          nb_chains)/Thin, 0), ".pdf")
                      pdf(file.pdf_RS2, paper = "a4", height = 20)
                      param = c("S2", "C2")
                      par(mfcol = c(5, 2))
                      longueur = 1:long
                      for (j in 1:2) {
                        for (i in 1:rs.length) {
                          plot(density(refstd_Parameters[, i, 
                            j]), lwd = 4, type = "l", col = "grey", 
                            main = paste(param[j], " of study ", 
                              i, " \n Thinning interval = ", 
                              thin.interval, "\n Total samplesize kept = ", 
                              (iter.num * nb_chains - (burn_in) * 
                                nb_chains)/Thin))
                        }
                      }
                      dev.off()
                    }
                  }
                  else {
                    refstd_parameters = matrix(0, ncol = 4, nrow = 2)
                    rownames(refstd_parameters) = c("S2", "C2")
                    colnames(refstd_parameters) = c("True Value", 
                      paste(point_estimate, "estimate"), "HPD.low", 
                      "HPD.high")
                    refstd_parameters[1, 1] <- true.S2
                    refstd_parameters[1, 2] <- S2.est
                    refstd_parameters[1, 3] <- S2.HPD[1]
                    refstd_parameters[1, 4] <- S2.HPD[2]
                    refstd_parameters[2, 1] <- true.C2
                    refstd_parameters[2, 2] <- C2.est
                    refstd_parameters[2, 3] <- C2.HPD[1]
                    refstd_parameters[2, 4] <- C2.HPD[2]
                    long = length(THETA)
                    refstd_Parameters = matrix(0, nrow = long, 
                      ncol = 2)
                    colnames(refstd_Parameters) = c("S2", "C2")
                    refstd_Parameters[, 1] <- S2
                    refstd_Parameters[, 2] <- C2
                    if (print_plot == TRUE) {
                      if (is.null(chain) == FALSE) {
                        file.pdf_RS = paste("RefStd trace plots for N =", 
                          round((iter.num * nb_chains - (burn_in) * 
                            nb_chains)/Thin, 0), ".pdf")
                        pdf(file.pdf_RS, paper = "a4", height = 20)
                        no_chains = length(chain)
                        iter_chain = round((iter.num * nb_chains - 
                          (burn_in) * nb_chains)/Thin, 0)/no_chains
                        longueur = 1:iter_chain
                        param = c("S2", "C2")
                        par(mfcol = c(5, 2))
                        for (j in 1:2) {
                          plot(x = longueur, y = refstd_Parameters[longueur, 
                            j], type = "n", col = 1, ylab = paste(param[j]), 
                            xlab = "iteration number", main = paste("Thinning interval = ", 
                              thin.interval, "\n Total samplesize kept = ", 
                              (iter.num * nb_chains - (burn_in) * 
                                nb_chains)/Thin))
                          for (l in 1:length(chain)) {
                            lines(x = longueur, y = refstd_Parameters[l * 
                              longueur, j], col = l)
                          }
                        }
                      }
                      else {
                        file.pdf_RS = paste("RefStd trace plots for N =", 
                          round((iter.num * nb_chains - (burn_in) * 
                            nb_chains)/Thin, 0), ".pdf")
                        pdf(file.pdf_RS, paper = "a4", height = 20)
                        Param = c("S2", "C2")
                        par(mfcol = c(5, 2))
                        for (j in 1:2) {
                          longueur = 1:long
                          plot(x = longueur, y = refstd_Parameters[, 
                            j], type = "l", col = "grey", ylab = paste(Param[j]), 
                            xlab = "iteration number", main = paste("Thinning interval = ", 
                              thin.interval, "\n Total samplesize kept = ", 
                              (iter.num * nb_chains - (burn_in) * 
                                nb_chains)/Thin))
                        }
                      }
                      dev.off()
                      file.pdf_RS2 = paste("RefStd density plots for N =", 
                        round((iter.num * nb_chains - (burn_in) * 
                          nb_chains)/Thin, 0), ".pdf")
                      pdf(file.pdf_RS2, paper = "a4", height = 20)
                      param = c("S2", "C2")
                      par(mfcol = c(5, 2))
                      for (j in 1:2) {
                        longueur = 1:long
                        plot(density(refstd_Parameters[, j]), 
                          lwd = 4, type = "l", col = "grey", 
                          main = paste(param[j], " of study ", 
                            j, " \n Thinning interval = ", thin.interval, 
                            "\n Total samplesize kept = ", (iter.num * 
                              nb_chains - (burn_in) * nb_chains)/Thin))
                      }
                      dev.off()
                    }
                  }
                }
            }
        }
        else {
            if (condInd == TRUE & Gold_Std == FALSE & model == 
                2) {
                if (rs.length != 1) {
                  refstd_parameters = array(0, dim = c(rs.length, 
                    4, 4), dimnames = list(1:rs.length, c("True Value", 
                    paste(point_estimate, "estimate"), "HPD lower", 
                    "HPD upper"), c("S2", "C2", "a1", "a0")))
                  refstd_parameters[, 1, 1] <- true.S2
                  refstd_parameters[, 2, 1] <- S2.est
                  refstd_parameters[, 3, 1] <- S2.HPD[, 1]
                  refstd_parameters[, 4, 1] <- S2.HPD[, 2]
                  refstd_parameters[, 1, 2] <- true.C2
                  refstd_parameters[, 2, 2] <- C2.est
                  refstd_parameters[, 3, 2] <- C2.HPD[, 1]
                  refstd_parameters[, 4, 2] <- C2.HPD[, 2]
                  refstd_parameters[, 1, 3] <- true.a1
                  refstd_parameters[, 2, 3] <- a1.est
                  refstd_parameters[, 3, 3] <- a1.HPD[, 1]
                  refstd_parameters[, 4, 3] <- a1.HPD[, 2]
                  refstd_parameters[, 1, 4] <- true.a0
                  refstd_parameters[, 2, 4] <- a0.est
                  refstd_parameters[, 3, 4] <- a0.HPD[, 1]
                  refstd_parameters[, 4, 4] <- a0.HPD[, 2]
                  long = length(alpha[, 1])
                  refstd_Parameters = array(0, c(long, rs.length, 
                    4))
                  refstd_Parameters[, , 1] <- S2
                  refstd_Parameters[, , 2] <- C2
                  refstd_Parameters[, , 3] <- a1
                  refstd_Parameters[, , 4] <- a0
                  if (print_plot == TRUE) {
                    file.pdf_RS = paste("RefStd trace plots for N =", 
                      round((iter.num * nb_chains - (burn_in) * 
                        nb_chains)/Thin, 0), ".pdf")
                    pdf(file.pdf_RS, paper = "a4", height = 20)
                    param = c("S2", "C2", "a1", "a0")
                    par(mfcol = c(5, 2))
                    for (j in 1:4) {
                      longueur = 1:long
                      for (i in 1:rs.length) {
                        plot(x = longueur, y = refstd_Parameters[, 
                          i, j], type = "l", col = "grey", ylab = paste(param[j], 
                          " of study ", i), xlab = "iteration number", 
                          main = paste("Thinning interval = ", 
                            thin.interval, "\n Total samplesize kept = ", 
                            (iter.num * nb_chains - (burn_in) * 
                              nb_chains)/Thin))
                        abline(a = refstd_parameters[i, 2, j], 
                          b = 0, col = "black", lwd = 3)
                        abline(a = refstd_parameters[i, 1, j], 
                          b = 0, col = "red", lwd = 3)
                        abline(a = refstd_parameters[i, 3, j], 
                          b = 0, col = "green", lwd = 3)
                        abline(a = refstd_parameters[i, 4, j], 
                          b = 0, col = "green", lwd = 3)
                      }
                    }
                    dev.off()
                    file.pdf_RS2 = paste("RefStd density plots for N =", 
                      round((iter.num * nb_chains - (burn_in) * 
                        nb_chains)/Thin, 0), ".pdf")
                    pdf(file.pdf_RS2, paper = "a4", height = 20)
                    param = c("S2", "C2", "a1", "a0")
                    par(mfcol = c(5, 2))
                    longueur = 1:long
                    for (j in 1:4) {
                      for (i in 1:rs.length) {
                        plot(density(refstd_Parameters[, i, j]), 
                          lwd = 4, type = "l", col = "grey", 
                          main = paste(param[j], " of study ", 
                            i, " \n Thinning interval = ", thin.interval, 
                            "\n Total samplesize kept = ", (iter.num * 
                              nb_chains - (burn_in) * nb_chains)/Thin))
                      }
                    }
                    dev.off()
                  }
                }
                else {
                  refstd_parameters = matrix(0, ncol = 4, nrow = 4)
                  rownames(refstd_parameters) = c("S2", "C2", 
                    "a1", "a0")
                  colnames(refstd_parameters) = c("True Value", 
                    paste(point_estimate, "estimate"), "HPD.low", 
                    "HPD.high")
                  refstd_parameters[1, 1] <- true.S2
                  refstd_parameters[1, 2] <- S2.est
                  refstd_parameters[1, 3] <- S2.HPD[1]
                  refstd_parameters[1, 4] <- S2.HPD[2]
                  refstd_parameters[2, 1] <- true.C2
                  refstd_parameters[2, 2] <- C2.est
                  refstd_parameters[2, 3] <- C2.HPD[1]
                  refstd_parameters[2, 4] <- C2.HPD[2]
                  refstd_parameters[3, 1] <- true.a1
                  refstd_parameters[3, 2] <- a1.est
                  refstd_parameters[3, 3] <- a1.HPD[1]
                  refstd_parameters[3, 4] <- a1.HPD[2]
                  refstd_parameters[4, 1] <- true.a0
                  refstd_parameters[4, 2] <- a0.est
                  refstd_parameters[4, 3] <- a0.HPD[1]
                  refstd_parameters[4, 4] <- a0.HPD[2]
                  long = length(THETA)
                  refstd_Parameters = matrix(0, nrow = long, 
                    ncol = 4)
                  colnames(refstd_Parameters) = c("S2", "C2", 
                    "a1", "a0")
                  refstd_Parameters[, 1] <- S2
                  refstd_Parameters[, 2] <- C2
                  refstd_Parameters[, 3] <- a1
                  refstd_Parameters[, 4] <- a0
                  if (print_plot == TRUE) {
                    file.pdf_RS = paste("RefStd trace plots for N =", 
                      round((iter.num * nb_chains - (burn_in) * 
                        nb_chains)/Thin, 0), ".pdf")
                    pdf(file.pdf_RS, paper = "a4", height = 20)
                    Param = c("S2", "C2", "a1", "a0")
                    par(mfcol = c(5, 2))
                    for (j in 1:4) {
                      longueur = 1:long
                      plot(x = longueur, y = refstd_Parameters[, 
                        j], type = "l", col = "grey", ylab = paste(Param[j]), 
                        xlab = "iteration number", main = paste("Thinning interval = ", 
                          thin.interval, "\n Total samplesize kept = ", 
                          (iter.num * nb_chains - (burn_in) * 
                            nb_chains)/Thin))
                      abline(a = refstd_parameters[j, 2], b = 0, 
                        col = "black", lwd = 3)
                      abline(a = refstd_parameters[j, 1], b = 0, 
                        col = "red", lwd = 3)
                      abline(a = refstd_parameters[j, 3], b = 0, 
                        col = "green", lwd = 3)
                      abline(a = refstd_parameters[j, 4], b = 0, 
                        col = "green", lwd = 3)
                    }
                    dev.off()
                    file.pdf_RS2 = paste("RefStd density plots for N =", 
                      round((iter.num * nb_chains - (burn_in) * 
                        nb_chains)/Thin, 0), ".pdf")
                    pdf(file.pdf_RS2, paper = "a4", height = 20)
                    param = c("S2", "C2", "a1", "a0")
                    par(mfcol = c(5, 2))
                    for (j in 1:4) {
                      longueur = 1:long
                      plot(density(refstd_Parameters[, j]), lwd = 4, 
                        type = "l", col = "grey", main = paste(param[j], 
                          " of study ", j, " \n Thinning interval = ", 
                          thin.interval, "\n Total samplesize kept = ", 
                          (iter.num * nb_chains - (burn_in) * 
                            nb_chains)/Thin))
                    }
                    dev.off()
                  }
                }
            }
        }
        if (SCO == FALSE) {
            parameters = array(0, dim = c(N, 4, 5), dimnames = list(Num_study, 
                c("True value", paste(point_estimate, "estimate"), 
                  "HPD lower", "HPD upper"), c("theta", "alpha", 
                  "pi", "S1", "C1")))
            parameters[, 1, 1] <- true.theta
            parameters[, 2, 1] <- theta.est
            parameters[, 3, 1] <- theta.HPD[, 1]
            parameters[, 4, 1] <- theta.HPD[, 2]
            parameters[, 1, 2] <- true.alpha
            parameters[, 2, 2] <- alpha.est
            parameters[, 3, 2] <- alpha.HPD[, 1]
            parameters[, 4, 2] <- alpha.HPD[, 2]
            parameters[, 1, 3] <- true.PI
            parameters[, 2, 3] <- PI.est
            parameters[, 3, 3] <- PI.HPD[, 1]
            parameters[, 4, 3] <- PI.HPD[, 2]
            parameters[, 1, 4] <- true.S1
            parameters[, 2, 4] <- S1.est
            parameters[, 3, 4] <- S1.HPD[, 1]
            parameters[, 4, 4] <- S1.HPD[, 2]
            parameters[, 1, 5] <- true.C1
            parameters[, 2, 5] <- C1.est
            parameters[, 3, 5] <- C1.HPD[, 1]
            parameters[, 4, 5] <- C1.HPD[, 2]
            long = length(alpha[, 1])
            Parameters = array(0, c(long, N, 5))
            Parameters[, , 1] <- theta
            Parameters[, , 2] <- alpha
            Parameters[, , 3] <- PI
            Parameters[, , 4] <- S1
            Parameters[, , 5] <- C1
            parameter = matrix(0, ncol = 4, nrow = 7)
            rownames(parameter) = c("THETA", "LAMBDA", "beta", 
                "sigma.alpha", "sigma.theta", "S Overall", "C Overall")
            colnames(parameter) = c("True.value", paste(point_estimate, 
                "estimate"), "HPD.low", "HPD.high")
            parameter[1, 1] <- true.THETA
            parameter[1, 2] <- THETA.est
            parameter[1, 3] <- THETA.HPD[1]
            parameter[1, 4] <- THETA.HPD[2]
            parameter[2, 1] <- true.LAMBDA
            parameter[2, 2] <- LAMBDA.est
            parameter[2, 3] <- LAMBDA.HPD[1]
            parameter[2, 4] <- LAMBDA.HPD[2]
            parameter[3, 1] <- true.beta
            parameter[3, 2] <- beta.est
            parameter[3, 3] <- beta.HPD[1]
            parameter[3, 4] <- beta.HPD[2]
            parameter[4, 1] <- true.sigma.alpha
            parameter[4, 2] <- sigma.alpha.est
            parameter[4, 3] <- sigma.alpha.HPD[1]
            parameter[4, 4] <- sigma.alpha.HPD[2]
            parameter[5, 1] <- true.sigma.theta
            parameter[5, 2] <- sigma.theta.est
            parameter[5, 3] <- sigma.theta.HPD[1]
            parameter[5, 4] <- sigma.theta.HPD[2]
            parameter[6, 1] <- true.S_overall
            parameter[6, 2] <- S_overall.est
            parameter[6, 3] <- S_overall.HPD[1]
            parameter[6, 4] <- S_overall.HPD[2]
            parameter[7, 1] <- true.C_overall
            parameter[7, 2] <- C_overall.est
            parameter[7, 3] <- C_overall.HPD[1]
            parameter[7, 4] <- C_overall.HPD[2]
            long = length(THETA)
            Parameter = matrix(0, nrow = long, ncol = 7)
            colnames(Parameter) = c("THETA", "LAMBDA", "beta", 
                "sigma.alpha", "sigma.theta", "S overall", "C overall")
            Parameter[, 1] <- THETA
            Parameter[, 2] <- LAMBDA
            Parameter[, 3] <- beta
            Parameter[, 4] <- sigma.alpha
            Parameter[, 5] <- sigma.theta
            Parameter[, 6] <- S_overall
            Parameter[, 7] <- C_overall
        }
        else {
            if (SCO == TRUE) {
                parameters = array(0, dim = c(N, 4, 4), dimnames = list(Num_study, 
                  c("True value", paste(point_estimate, "estimate"), 
                    "HPD lower", "HPD upper"), c("alpha", "pi", 
                    "S1", "C1")))
                parameters[, 1, 1] <- true.alpha
                parameters[, 2, 1] <- alpha.est
                parameters[, 3, 1] <- alpha.HPD[, 1]
                parameters[, 4, 1] <- alpha.HPD[, 2]
                parameters[, 1, 2] <- true.PI
                parameters[, 2, 2] <- PI.est
                parameters[, 3, 2] <- PI.HPD[, 1]
                parameters[, 4, 2] <- PI.HPD[, 2]
                parameters[, 1, 3] <- true.S1
                parameters[, 2, 3] <- S1.est
                parameters[, 3, 3] <- S1.HPD[, 1]
                parameters[, 4, 3] <- S1.HPD[, 2]
                parameters[, 1, 4] <- true.C1
                parameters[, 2, 4] <- C1.est
                parameters[, 3, 4] <- C1.HPD[, 1]
                parameters[, 4, 4] <- C1.HPD[, 2]
                long = length(alpha[, 1])
                Parameters = array(0, c(long, N, 4))
                Parameters[, , 1] <- alpha
                Parameters[, , 2] <- PI
                Parameters[, , 3] <- S1
                Parameters[, , 4] <- C1
                parameter = matrix(0, ncol = 4, nrow = 6)
                rownames(parameter) = c("THETA", "LAMBDA", "beta", 
                  "sigma.alpha", "S Overall", "C Overall")
                colnames(parameter) = c("True.value", paste(point_estimate, 
                  "estimate"), "HPD.low", "HPD.high")
                parameter[1, 1] <- true.THETA
                parameter[1, 2] <- THETA.est
                parameter[1, 3] <- THETA.HPD[1]
                parameter[1, 4] <- THETA.HPD[2]
                parameter[2, 1] <- true.LAMBDA
                parameter[2, 2] <- LAMBDA.est
                parameter[2, 3] <- LAMBDA.HPD[1]
                parameter[2, 4] <- LAMBDA.HPD[2]
                parameter[3, 1] <- true.beta
                parameter[3, 2] <- beta.est
                parameter[3, 3] <- beta.HPD[1]
                parameter[3, 4] <- beta.HPD[2]
                parameter[4, 1] <- true.sigma.alpha
                parameter[4, 2] <- sigma.alpha.est
                parameter[4, 3] <- sigma.alpha.HPD[1]
                parameter[4, 4] <- sigma.alpha.HPD[2]
                parameter[5, 1] <- true.S_overall
                parameter[5, 2] <- S_overall.est
                parameter[5, 3] <- S_overall.HPD[1]
                parameter[5, 4] <- S_overall.HPD[2]
                parameter[6, 1] <- true.C_overall
                parameter[6, 2] <- C_overall.est
                parameter[6, 3] <- C_overall.HPD[1]
                parameter[6, 4] <- C_overall.HPD[2]
                long = length(THETA)
                Parameter = matrix(0, nrow = long, ncol = 6)
                colnames(Parameter) = c("THETA", "LAMBDA", "beta", 
                  "sigma.alpha", "S overall", "C overall")
                Parameter[, 1] <- THETA
                Parameter[, 2] <- LAMBDA
                Parameter[, 3] <- beta
                Parameter[, 4] <- sigma.alpha
                Parameter[, 5] <- S_overall
                Parameter[, 6] <- C_overall
            }
        }
        if (print_plot == TRUE) {
            if (SCO == FALSE) {
                if (is.null(chain) == FALSE) {
                  file.pdf5 = paste("Trace plots for N =", round((iter.num * 
                    nb_chains - (burn_in) * nb_chains)/Thin, 
                    0), ".pdf")
                  pdf(file.pdf5, paper = "a4", height = 20)
                  param = c("theta", "alpha", "PI", "S1", "C1")
                  Param = c("Capital Theta", "Capital Lambda", 
                    "beta", "~sigma[alpha]", "~sigma[theta]", 
                    "S Overall", "C Overall")
                  no_chains = length(chain)
                  iter_chain = round((iter.num * nb_chains - 
                    (burn_in) * nb_chains)/Thin, 0)/no_chains
                  min_param = c(min(Parameters[, , 1]), min(Parameters[, 
                    , 2]), min(Parameters[, , 3]), min(Parameters[, 
                    , 4]), min(Parameters[, , 5]))
                  max_param = c(max(Parameters[, , 1]), max(Parameters[, 
                    , 2]), max(Parameters[, , 3]), max(Parameters[, 
                    , 4]), max(Parameters[, , 5]))
                  dlag = (max_param - min_param)/100
                  range_param = numeric()
                  for (j in 1:5) {
                    range_param = cbind(range_param, seq(min_param[j] + 
                      dlag[j]/2, max_param[j] - dlag[j]/2, by = dlag[j]))
                  }
                  par(mfcol = c(5, 2))
                  longueur = 1:iter_chain
                  for (j in 1:5) {
                    for (i in 1:N) {
                      plot(x = longueur, y = Parameters[longueur, 
                        i, j], type = "n", col = 1, ylab = paste(param[j], 
                        " of study ", i), xlab = "iteration number", 
                        main = paste("Thinning interval = ", 
                          thin.interval, "\n Total samplesize kept = ", 
                          (iter.num * nb_chains - (burn_in) * 
                            nb_chains)/Thin), ylim = range(range_param[, 
                          j]))
                      for (l in 1:length(chain)) {
                        lines(x = longueur, y = Parameters[l * 
                          longueur, i, j], col = l)
                      }
                    }
                  }
                  min_Param = c(min(Parameter[, 1]), min(Parameter[, 
                    2]), min(Parameter[, 3]), min(Parameter[, 
                    4]), min(Parameter[, 5]), min(Parameter[, 
                    6]), min(Parameter[, 7]))
                  max_Param = c(max(Parameter[, 1]), max(Parameter[, 
                    2]), max(Parameter[, 3]), max(Parameter[, 
                    4]), max(Parameter[, 5]), max(Parameter[, 
                    6]), max(Parameter[, 7]))
                  dlag = (max_Param - min_Param)/100
                  range_Param = numeric()
                  for (j in 1:7) {
                    range_Param = cbind(range_Param, seq(min_Param[j] + 
                      dlag[j]/2, max_Param[j] - dlag[j]/2, by = dlag[j]))
                  }
                  for (j in 1:7) {
                    plot(x = longueur, y = Parameter[longueur, 
                      j], type = "n", col = 1, ylab = paste(Param[j]), 
                      xlab = "iteration number", main = paste("Thinning interval = ", 
                        thin.interval, "\n Total samplesize kept = ", 
                        (iter.num * nb_chains - (burn_in) * nb_chains)/Thin), 
                      ylim = range(range_Param[, j]))
                    for (l in 1:length(chain)) {
                      lines(x = longueur, y = Parameter[l * longueur, 
                        j], col = l)
                    }
                  }
                  dev.off()
                }
                else {
                  file.pdf2 = paste("Trace plots for N =", round((iter.num * 
                    nb_chains - (burn_in) * nb_chains)/Thin, 
                    0), ".pdf")
                  pdf(file.pdf2, paper = "a4", height = 20)
                  param = c("theta", "alpha", "PI", "S1", "C1")
                  Param = c("Capital Theta", "Capital Lambda", 
                    "beta", "~sigma[alpha]", "~sigma[theta]", 
                    "S Overall", "C Overall")
                  min_param = c(min(Parameters[, , 1]), min(Parameters[, 
                    , 2]), min(Parameters[, , 3]), min(Parameters[, 
                    , 4]), min(Parameters[, , 5]))
                  max_param = c(max(Parameters[, , 1]), max(Parameters[, 
                    , 2]), max(Parameters[, , 3]), max(Parameters[, 
                    , 4]), max(Parameters[, , 5]))
                  dlag = (max_param - min_param)/100
                  range_param = numeric()
                  for (j in 1:5) {
                    range_param = cbind(range_param, seq(min_param[j] + 
                      dlag[j]/2, max_param[j] - dlag[j]/2, by = dlag[j]))
                  }
                  par(mfcol = c(5, 2))
                  longueur = 1:long
                  for (j in 1:5) {
                    for (i in 1:N) {
                      plot(x = longueur, y = Parameters[, i, 
                        j], type = "l", col = "grey", ylab = paste(param[j], 
                        " of study ", i), xlab = "iteration number", 
                        main = paste("Thinning interval = ", 
                          thin.interval, "\n Total samplesize kept = ", 
                          (iter.num * nb_chains - (burn_in) * 
                            nb_chains)/Thin), ylim = range(range_param[, 
                          j]))
                      abline(a = parameters[i, 3, j], b = 0, 
                        col = "green", lwd = 3)
                      abline(a = parameters[i, 4, j], b = 0, 
                        col = "green", lwd = 3)
                    }
                  }
                  min_Param = c(min(Parameter[, 1]), min(Parameter[, 
                    2]), min(Parameter[, 3]), min(Parameter[, 
                    4]), min(Parameter[, 5]), min(Parameter[, 
                    6]), min(Parameter[, 7]))
                  max_Param = c(max(Parameter[, 1]), max(Parameter[, 
                    2]), max(Parameter[, 3]), max(Parameter[, 
                    4]), max(Parameter[, 5]), max(Parameter[, 
                    6]), max(Parameter[, 7]))
                  dlag = (max_Param - min_Param)/100
                  range_Param = numeric()
                  for (j in 1:7) {
                    range_Param = cbind(range_Param, seq(min_Param[j] + 
                      dlag[j]/2, max_Param[j] - dlag[j]/2, by = dlag[j]))
                  }
                  for (j in 1:7) {
                    plot(x = longueur, y = Parameter[, j], type = "l", 
                      col = "grey", ylab = paste(Param[j]), xlab = "iteration number", 
                      main = paste("Thinning interval = ", thin.interval, 
                        "\n Total samplesize kept = ", (iter.num * 
                          nb_chains - (burn_in) * nb_chains)/Thin), 
                      ylim = range(range_Param[, j]))
                    abline(a = parameter[j, 3], b = 0, col = "green", 
                      lwd = 3)
                    abline(a = parameter[j, 4], b = 0, col = "green", 
                      lwd = 3)
                  }
                  dev.off()
                }
                file.pdf3 = paste("Density plots for N =", round((iter.num * 
                  nb_chains - (burn_in) * nb_chains)/Thin, 0), 
                  ".pdf")
                pdf(file.pdf3, paper = "a4", height = 20)
                param = c("theta", "alpha", "PI", "S1", "C1")
                Param = c("Capital Theta", "Capital Lambda", 
                  "beta", "~sigma[alpha]", "~sigma[theta]", "S Overall", 
                  "C Overall")
                par(mfcol = c(5, 2))
                longueur = 1:long
                for (j in 1:5) {
                  for (i in 1:N) {
                    plot(density(Parameters[, i, j]), lwd = 4, 
                      type = "l", col = "grey", main = paste(param[j], 
                        " of study ", i, " \n Thinning interval = ", 
                        thin.interval, "\n Total samplesize kept = ", 
                        (iter.num * nb_chains - (burn_in) * nb_chains)/Thin))
                  }
                }
                for (j in 1:7) {
                  plot(density(Parameter[, j]), lwd = 4, type = "l", 
                    col = "grey", main = paste(Param[j], " \n Thinning interval = ", 
                      thin.interval, "\n Total samplesize kept = ", 
                      (iter.num * nb_chains - (burn_in) * nb_chains)/Thin))
                }
                dev.off()
                pdf("Summary ROC curve.pdf")
                default.x = range(1, 0)
                default.y = range(0, 1)
                plot(x = default.x, y = default.y, type = "n", 
                  xlim = rev(range(default.x)), xlab = "", ylab = "")
                title(xlab = "Specificity", ylab = "Sensitivity", 
                  cex.lab = 1.5, main = "Summary ROC curve")
                Sensi1 = apply(as.matrix(Parameters[, , 4]), 
                  2, median)
                Speci1 = apply(as.matrix(Parameters[, , 5]), 
                  2, median)
                Scale_factor = 10
                symbols(Speci1, Sensi1, circles = rowSums(as.matrix(data[[1]])), 
                  inches = 0.1 * Scale_factor/7, add = TRUE)
                Ov_Se = 1 - pnorm((median(Parameter[, 1]) - median(Parameter[, 
                  2])/2)/exp(median(Parameter[, 3])/2))
                Ov_Sp = pnorm((median(Parameter[, 1]) + median(Parameter[, 
                  2])/2)/exp(-median(Parameter[, 3])/2))
                points(Ov_Sp, Ov_Se, pch = 19, cex = 2)
                thet = qnorm((1 - as.matrix(Parameters[, , 4])) + 
                  1e-14) * exp(Parameter[, 3]/2) + Parameter[, 
                  2]/2
                min_TH = quantile(thet, 0.05)
                max_TH = quantile(thet, 0.95)
                dTH = 5e-05
                TH_range = seq(min_TH + dTH/2, max_TH - dTH/2, 
                  dTH)
                S_sroc = 1 - pnorm((TH_range - median(Parameter[, 
                  2])/2)/exp(median(Parameter[, 3])/2))
                C_sroc = pnorm((TH_range + median(Parameter[, 
                  2])/2)/exp(-median(Parameter[, 3])/2))
                lines(C_sroc, S_sroc, lwd = 3, col = "black", 
                  lty = 1)
                dev.off()
            }
            else {
                if (SCO == TRUE) {
                  if (is.null(chain) == FALSE) {
                    file.pdf5 = paste("Trace plots for N =", 
                      round((iter.num * nb_chains - (burn_in) * 
                        nb_chains)/Thin, 0), ".pdf")
                    pdf(file.pdf5, paper = "a4", height = 20)
                    param = c("alpha", "PI", "S1", "C1")
                    Param = c("Capital Theta", "Capital Lambda", 
                      "beta", "~sigma[alpha]", "S Overall", "C Overall")
                    no_chains = length(chain)
                    iter_chain = round((iter.num * nb_chains - 
                      (burn_in) * nb_chains)/Thin, 0)/no_chains
                    min_param = c(min(Parameters[, , 1]), min(Parameters[, 
                      , 2]), min(Parameters[, , 3]), min(Parameters[, 
                      , 4]))
                    max_param = c(max(Parameters[, , 1]), max(Parameters[, 
                      , 2]), max(Parameters[, , 3]), max(Parameters[, 
                      , 4]))
                    dlag = (max_param - min_param)/100
                    range_param = numeric()
                    for (j in 1:4) {
                      range_param = cbind(range_param, seq(min_param[j] + 
                        dlag[j]/2, max_param[j] - dlag[j]/2, 
                        by = dlag[j]))
                    }
                    par(mfcol = c(5, 2))
                    longueur = 1:iter_chain
                    for (j in 1:4) {
                      for (i in 1:N) {
                        plot(x = longueur, y = Parameters[longueur, 
                          i, j], type = "n", col = 1, ylab = paste(param[j], 
                          " of study ", i), xlab = "iteration number", 
                          main = paste("Thinning interval = ", 
                            thin.interval, "\n Total samplesize kept = ", 
                            (iter.num * nb_chains - (burn_in) * 
                              nb_chains)/Thin), ylim = range(range_param[, 
                            j]))
                        for (l in 1:length(chain)) {
                          lines(x = longueur, y = Parameters[l * 
                            longueur, i, j], col = l)
                        }
                      }
                    }
                    min_Param = c(min(Parameter[, 1]), min(Parameter[, 
                      2]), min(Parameter[, 3]), min(Parameter[, 
                      4]), min(Parameter[, 5]), min(Parameter[, 
                      6]))
                    max_Param = c(max(Parameter[, 1]), max(Parameter[, 
                      2]), max(Parameter[, 3]), max(Parameter[, 
                      4]), max(Parameter[, 5]), max(Parameter[, 
                      6]))
                    dlag = (max_Param - min_Param)/100
                    range_Param = numeric()
                    for (j in 1:6) {
                      range_Param = cbind(range_Param, seq(min_Param[j] + 
                        dlag[j]/2, max_Param[j] - dlag[j]/2, 
                        by = dlag[j]))
                    }
                    for (j in 1:6) {
                      plot(x = longueur, y = Parameter[longueur, 
                        j], type = "n", col = 1, ylab = paste(Param[j]), 
                        xlab = "iteration number", main = paste("Thinning interval = ", 
                          thin.interval, "\n Total samplesize kept = ", 
                          (iter.num * nb_chains - (burn_in) * 
                            nb_chains)/Thin), ylim = range(range_Param[, 
                          j]))
                      for (l in 1:length(chain)) {
                        lines(x = longueur, y = Parameter[l * 
                          longueur, j], col = l)
                      }
                    }
                    dev.off()
                  }
                  else {
                    file.pdf2 = paste("Trace plots for N =", 
                      round((iter.num * nb_chains - (burn_in) * 
                        nb_chains)/Thin, 0), ".pdf")
                    pdf(file.pdf2, paper = "a4", height = 20)
                    param = c("alpha", "PI", "S1", "C1")
                    Param = c("Capital Theta", "Capital Lambda", 
                      "beta", "~sigma[alpha]", "S Overall", "C Overall")
                    min_param = c(min(Parameters[, , 1]), min(Parameters[, 
                      , 2]), min(Parameters[, , 3]), min(Parameters[, 
                      , 4]))
                    max_param = c(max(Parameters[, , 1]), max(Parameters[, 
                      , 2]), max(Parameters[, , 3]), max(Parameters[, 
                      , 4]))
                    dlag = (max_param - min_param)/100
                    range_param = numeric()
                    for (j in 1:4) {
                      range_param = cbind(range_param, seq(min_param[j] + 
                        dlag[j]/2, max_param[j] - dlag[j]/2, 
                        by = dlag[j]))
                    }
                    par(mfcol = c(5, 2))
                    longueur = 1:long
                    for (j in 1:4) {
                      for (i in 1:N) {
                        plot(x = longueur, y = Parameters[, i, 
                          j], type = "l", col = "grey", ylab = paste(param[j], 
                          " of study ", i), xlab = "iteration number", 
                          main = paste("Thinning interval = ", 
                            thin.interval, "\n Total samplesize kept = ", 
                            (iter.num * nb_chains - (burn_in) * 
                              nb_chains)/Thin), ylim = range(range_param[, 
                            j]))
                        abline(a = parameters[i, 3, j], b = 0, 
                          col = "green", lwd = 3)
                        abline(a = parameters[i, 4, j], b = 0, 
                          col = "green", lwd = 3)
                      }
                    }
                    min_Param = c(min(Parameter[, 1]), min(Parameter[, 
                      2]), min(Parameter[, 3]), min(Parameter[, 
                      4]), min(Parameter[, 5]), min(Parameter[, 
                      6]))
                    max_Param = c(max(Parameter[, 1]), max(Parameter[, 
                      2]), max(Parameter[, 3]), max(Parameter[, 
                      4]), max(Parameter[, 5]), max(Parameter[, 
                      6]))
                    dlag = (max_Param - min_Param)/100
                    range_Param = numeric()
                    for (j in 1:6) {
                      range_Param = cbind(range_Param, seq(min_Param[j] + 
                        dlag[j]/2, max_Param[j] - dlag[j]/2, 
                        by = dlag[j]))
                    }
                    for (j in 1:6) {
                      plot(x = longueur, y = Parameter[, j], 
                        type = "l", col = "grey", ylab = paste(Param[j]), 
                        xlab = "iteration number", main = paste("Thinning interval = ", 
                          thin.interval, "\n Total samplesize kept = ", 
                          (iter.num * nb_chains - (burn_in) * 
                            nb_chains)/Thin), ylim = range(range_Param[, 
                          j]))
                      abline(a = parameter[j, 3], b = 0, col = "green", 
                        lwd = 3)
                      abline(a = parameter[j, 4], b = 0, col = "green", 
                        lwd = 3)
                    }
                    dev.off()
                  }
                  file.pdf3 = paste("Density plots for N =", 
                    round((iter.num * nb_chains - (burn_in) * 
                      nb_chains)/Thin, 0), ".pdf")
                  pdf(file.pdf3, paper = "a4", height = 20)
                  param = c("alpha", "PI", "S1", "C1")
                  Param = c("Capital Theta", "Capital Lambda", 
                    "beta", "~sigma[alpha]", "S Overall", "C Overall")
                  par(mfcol = c(5, 2))
                  longueur = 1:long
                  for (j in 1:4) {
                    for (i in 1:N) {
                      plot(density(Parameters[, i, j]), lwd = 4, 
                        type = "l", col = "grey", main = paste(param[j], 
                          " of study ", i, " \n Thinning interval = ", 
                          thin.interval, "\n Total samplesize kept = ", 
                          (iter.num * nb_chains - (burn_in) * 
                            nb_chains)/Thin))
                    }
                  }
                  for (j in 1:6) {
                    plot(density(Parameter[, j]), lwd = 4, type = "l", 
                      col = "grey", main = paste(Param[j], " \n Thinning interval = ", 
                        thin.interval, "\n Total samplesize kept = ", 
                        (iter.num * nb_chains - (burn_in) * nb_chains)/Thin))
                  }
                  dev.off()
                  pdf("Summary ROC curve.pdf")
                  default.x = range(1, 0)
                  default.y = range(0, 1)
                  plot(x = default.x, y = default.y, type = "n", 
                    xlim = rev(range(default.x)), xlab = "", 
                    ylab = "")
                  title(xlab = "Specificity", ylab = "Sensitivity", 
                    cex.lab = 1.5, main = "Summary ROC curve")
                  Sensi1 = apply(as.matrix(Parameters[, , 3]), 
                    2, median)
                  Speci1 = apply(as.matrix(Parameters[, , 4]), 
                    2, median)
                  Scale_factor = 10
                  symbols(Speci1, Sensi1, circles = rowSums(as.matrix(data[[1]])), 
                    inches = 0.1 * Scale_factor/7, add = TRUE)
                  Ov_Se = 1 - pnorm((median(Parameter[, 1]) - 
                    median(Parameter[, 2])/2)/exp(median(Parameter[, 
                    3])/2))
                  Ov_Sp = pnorm((median(Parameter[, 1]) + median(Parameter[, 
                    2])/2)/exp(-median(Parameter[, 3])/2))
                  points(Ov_Sp, Ov_Se, pch = 19, cex = 2)
                  thet = qnorm((1 - as.matrix(Parameters[, , 
                    3])) + 1e-14) * exp(Parameter[, 3]/2) + Parameter[, 
                    2]/2
                  min_TH = quantile(thet, 0.05)
                  max_TH = quantile(thet, 0.95)
                  dTH = 5e-05
                  TH_range = seq(min_TH + dTH/2, max_TH - dTH/2, 
                    dTH)
                  S_sroc = 1 - pnorm((TH_range - median(Parameter[, 
                    2])/2)/exp(median(Parameter[, 3])/2))
                  C_sroc = pnorm((TH_range + median(Parameter[, 
                    2])/2)/exp(-median(Parameter[, 3])/2))
                  lines(C_sroc, S_sroc, lwd = 3, col = "black", 
                    lty = 1)
                  dev.off()
                }
            }
        }
    }
    else {
        if (real_life == TRUE) {
            d = as.matrix(data[[1]])
            Sample.size = d[, 1] + d[, 2] + d[, 3] + d[, 4]
            pp = d[, 1]
            pn = d[, 2]
            np = d[, 3]
            nn = d[, 4]
            test.file = paste("Summary for N =", round((iter.num * 
                nb_chains - (burn_in) * nb_chains)/Thin, 0), 
                ".txt")
            write(paste("Number of chains =", nb_chains), file = test.file, 
                append = TRUE)
            write(paste("Number of iteration within a chain =", 
                iter.num, "    Burn in within each chain =", 
                burn_in), file = test.file, append = TRUE)
            write(paste("Thinning interval =", Thin), file = test.file, 
                append = TRUE)
            write(paste("Total number of iteration kept =", round((iter.num * 
                nb_chains - (burn_in) * nb_chains)/Thin, 0)), 
                file = test.file, append = TRUE)
            write(paste(""), file = test.file, append = TRUE)
            write(paste("File location : ", path), file = test.file, 
                append = TRUE)
            write(paste(""), file = test.file, append = TRUE)
            write(paste("Date :", Sys.time()), file = test.file, 
                append = TRUE)
            write(paste("______________________________________________________"), 
                file = test.file, append = TRUE)
            write(paste(""), file = test.file, append = TRUE)
            if (Gold_Std == TRUE) {
                write("Perfect reference standard", file = test.file, 
                  append = TRUE)
            }
            else {
                write("Imperfect reference standard", file = test.file, 
                  append = TRUE)
            }
            write(paste(""), file = test.file, append = TRUE)
            write(paste("______________________________________________________"), 
                file = test.file, append = TRUE)
            write(paste("\tSAMPLE SIZE \t "), file = test.file, 
                append = TRUE)
            write(paste("______________________________________________________"), 
                file = test.file, append = TRUE)
            write(paste(""), file = test.file, append = TRUE)
            write(paste("         Total ++ +- -+ --"), file = test.file, 
                append = TRUE)
            for (i in 1:N) {
                write(paste("Study ", i, "", Sample.size[i], 
                  "", pp[i], "", pn[i], "", np[i], "", nn[i]), 
                  file = test.file, append = TRUE)
            }
            write(paste(""), file = test.file, append = TRUE)
            write(paste("______________________________________________________"), 
                file = test.file, append = TRUE)
            write(paste("\tPRIOR INFORMATION \t "), file = test.file, 
                append = TRUE)
            write(paste("______________________________________________________"), 
                file = test.file, append = TRUE)
            write(paste(""), file = test.file, append = TRUE)
            write(paste(""), file = test.file, append = TRUE)
            write(paste(""), file = test.file, append = TRUE)
            write(paste("Prior of prevalence (pi) is ", prior_dist_PI, 
                "(", round(alpha.PI, digits = 4), ",", round(beta.PI, 
                  digits = 4), "), <=> pi in [", low.pi, ",", 
                up.pi, "]"), file = test.file, append = TRUE)
            write(paste(""), file = test.file, append = TRUE)
            write(paste("Prior of beta is Uniform(", round(beta.a, 
                4), ",", round(beta.b, 4), ")"), file = test.file, 
                append = TRUE)
            write(paste("Prior of THETA is Uniform(", prior.THETA.lower, 
                ",", prior.THETA.upper, ")"), file = test.file, 
                append = TRUE)
            write(paste("Prior of LAMBDA is Uniform(", prior.LAMBDA.lower, 
                ",", prior.LAMBDA.upper, ")"), file = test.file, 
                append = TRUE)
            write(paste("Prior of sigma_alpha is uniform(", l.disp.alpha, 
                ",", u.disp.alpha, ")"), file = test.file, append = TRUE)
            if (SCO == FALSE) {
                write(paste("Prior of sigma_theta is uniform(", 
                  l.disp.theta, ",", u.disp.theta, ")"), file = test.file, 
                  append = TRUE)
            }
            if (condInd == TRUE & Gold_Std == FALSE & model == 
                1) {
                write(paste(""), file = test.file, append = TRUE)
                if (Gold_se == TRUE & Gold_sp == FALSE) {
                  write(paste("Prior of S2 (Sensitivity of reference test) is "), 
                    file = test.file, append = TRUE)
                  for (i in 1:rs.length) {
                    write(paste("Study(ies) ", sub_rs[[i + 1]][1], 
                      "to", sub_rs[[i + 1]][length(sub_rs[[i + 
                        1]])], "assumed to be perfect."), file = test.file, 
                      append = TRUE)
                  }
                  write(paste(""), file = test.file, append = TRUE)
                  write(paste("Prior of C2 (Specificity of reference test) is "), 
                    file = test.file, append = TRUE)
                  for (i in 1:rs.length) {
                    write(paste("Study(ies) ", sub_rs[[i + 1]][1], 
                      "to", sub_rs[[i + 1]][length(sub_rs[[i + 
                        1]])], "", prior_dist_C2, "(", round(Spec2.alpha[i], 
                        digits = 4), ",", round(Spec2.beta[i], 
                        digits = 4), "), <=> C2 in [", low.sp[i], 
                      ",", up.sp[i], "]"), file = test.file, 
                      append = TRUE)
                  }
                }
                else {
                  if (Gold_sp == TRUE & Gold_se == FALSE) {
                    write(paste("Prior of S2 (Sensitivity of reference test) is "), 
                      file = test.file, append = TRUE)
                    for (i in 1:rs.length) {
                      write(paste("Study(ies) ", sub_rs[[i + 
                        1]][1], "to", sub_rs[[i + 1]][length(sub_rs[[i + 
                        1]])], "", prior_dist_S2, "(", round(Sens2.alpha[i], 
                        digits = 4), ",", round(Sens2.beta[i], 
                        digits = 4), "), <=> S2 in [", low.se[i], 
                        ",", up.se[i], "]"), file = test.file, 
                        append = TRUE)
                    }
                    write(paste(""), file = test.file, append = TRUE)
                    write(paste("Prior of C2 (Specificity of reference test) is "), 
                      file = test.file, append = TRUE)
                    for (i in 1:rs.length) {
                      write(paste("Study(ies) ", sub_rs[[i + 
                        1]][1], "to", sub_rs[[i + 1]][length(sub_rs[[i + 
                        1]])], "assumed to be perfect."), file = test.file, 
                        append = TRUE)
                    }
                  }
                  else {
                    write(paste("Prior of S2 (Sensitivity of reference test) is "), 
                      file = test.file, append = TRUE)
                    for (i in 1:rs.length) {
                      write(paste("Study(ies) ", sub_rs[[i + 
                        1]][1], "to", sub_rs[[i + 1]][length(sub_rs[[i + 
                        1]])], "", prior_dist_S2, "(", round(Sens2.alpha[i], 
                        digits = 4), ",", round(Sens2.beta[i], 
                        digits = 4), "), <=> S2 in [", low.se[i], 
                        ",", up.se[i], "]"), file = test.file, 
                        append = TRUE)
                    }
                    write(paste(""), file = test.file, append = TRUE)
                    write(paste("Prior of C2 (Specificity of reference test) is "), 
                      file = test.file, append = TRUE)
                    for (i in 1:rs.length) {
                      write(paste("Study(ies) ", sub_rs[[i + 
                        1]][1], "to", sub_rs[[i + 1]][length(sub_rs[[i + 
                        1]])], "", prior_dist_C2, "(", round(Spec2.alpha[i], 
                        digits = 4), ",", round(Spec2.beta[i], 
                        digits = 4), "), <=> C2 in [", low.sp[i], 
                        ",", up.sp[i], "]"), file = test.file, 
                        append = TRUE)
                    }
                  }
                }
            }
            else {
                if (condInd == TRUE & Gold_Std == FALSE & model == 
                  2) {
                  write(paste(""), file = test.file, append = TRUE)
                  write(paste("Prior of a1 is "), file = test.file, 
                    append = TRUE)
                  for (i in 1:rs.length) {
                    write(paste("Study(ies) ", sub_rs[[i + 1]][1], 
                      "to", sub_rs[[i + 1]][length(sub_rs[[i + 
                        1]])], "Normal (", round(mean.a1[i], 
                        digits = 4), ",", round(sd.a1[i], digits = 4), 
                      "), <=> S2 in []"), file = test.file, append = TRUE)
                  }
                  write(paste(""), file = test.file, append = TRUE)
                  write(paste("Prior of a0 is "), file = test.file, 
                    append = TRUE)
                  for (i in 1:rs.length) {
                    write(paste("Study(ies) ", sub_rs[[i + 1]][1], 
                      "to", sub_rs[[i + 1]][length(sub_rs[[i + 
                        1]])], "Normal (", round(mean.a0[i], 
                        digits = 4), ",", round(sd.a0[i], digits = 4), 
                      "), <=> C2 in []"), file = test.file, append = TRUE)
                  }
                }
                else {
                  if (condInd == FALSE) {
                    write(paste(""), file = test.file, append = TRUE)
                    write(paste("Prior of a1 is "), file = test.file, 
                      append = TRUE)
                    for (i in 1:rs.length) {
                      write(paste("Study(ies) ", sub_rs[[i + 
                        1]][1], "to", sub_rs[[i + 1]][length(sub_rs[[i + 
                        1]])], "Normal (", round(mean.a1[i], 
                        digits = 4), ",", round(sd.a1[i], digits = 4), 
                        "), <=> S2 in []"), file = test.file, 
                        append = TRUE)
                    }
                    write(paste(""), file = test.file, append = TRUE)
                    write(paste("Prior of a0 is "), file = test.file, 
                      append = TRUE)
                    for (i in 1:rs.length) {
                      write(paste("Study(ies) ", sub_rs[[i + 
                        1]][1], "to", sub_rs[[i + 1]][length(sub_rs[[i + 
                        1]])], "Normal (", round(mean.a0[i], 
                        digits = 4), ",", round(sd.a0[i], digits = 4), 
                        "), <=> C2 in []"), file = test.file, 
                        append = TRUE)
                    }
                    write(paste(""), file = test.file, append = TRUE)
                    write(paste("Prior of b1 is "), file = test.file, 
                      append = TRUE)
                    for (i in 1:rs.length) {
                      write(paste("Study(ies) ", sub_rs[[i + 
                        1]][1], "to", sub_rs[[i + 1]][length(sub_rs[[i + 
                        1]])], "Uniform (", round(low.b1[i], 
                        digits = 4), ",", round(up.b1[i], digits = 4), 
                        "), <=> S2 in []"), file = test.file, 
                        append = TRUE)
                    }
                    write(paste(""), file = test.file, append = TRUE)
                    write(paste("Prior of b0 is "), file = test.file, 
                      append = TRUE)
                    for (i in 1:rs.length) {
                      write(paste("Study(ies) ", sub_rs[[i + 
                        1]][1], "to", sub_rs[[i + 1]][length(sub_rs[[i + 
                        1]])], "Uniform (", round(low.b0[i], 
                        digits = 4), ",", round(up.b0[i], digits = 4), 
                        "), <=> C2 in []"), file = test.file, 
                        append = TRUE)
                    }
                    write(paste(""), file = test.file, append = TRUE)
                    write(paste("Prior of d1 is Uniform(", round(low.d1, 
                      4), ",", round(up.d1, 4), ")"), file = test.file, 
                      append = TRUE)
                    write(paste("Prior of d0 is Uniform(", round(low.d0, 
                      4), ",", round(up.d0, 4), ")"), file = test.file, 
                      append = TRUE)
                  }
                }
            }
            write(paste(), file = test.file, append = TRUE)
            write(paste("______________________________________________________"), 
                file = test.file, append = TRUE)
            write(paste("______________________________________________________"), 
                file = test.file, append = TRUE)
            write(paste("\tBETWEEN_STUDY parameters (Point estimate =", 
                point_estimate, ")\t "), file = test.file, append = TRUE)
            write(paste("______________________________________________________"), 
                file = test.file, append = TRUE)
            write(paste("______________________________________________________"), 
                file = test.file, append = TRUE)
            write(paste(""), file = test.file, append = TRUE)
            write(paste("         Estimate Standard_Dev MC_error C.I._lower C.I._upper"), 
                file = test.file, append = TRUE)
            write(paste(""), file = test.file, append = TRUE)
            write(paste("THETA       ", round(THETA.est, digits = digit), 
                "", round(THETA.sd, digits = digit), "", round(THETA.MCerror, 
                  digits = digit), "", round(THETA.HPD[1], digits = digit), 
                "", round(THETA.HPD[2], digits = digit)), file = test.file, 
                append = TRUE)
            write(paste("LAMBDA      ", round(LAMBDA.est, digits = digit), 
                "", round(LAMBDA.sd, digits = digit), "", round(LAMBDA.MCerror, 
                  digits = digit), "", round(LAMBDA.HPD[1], digits = digit), 
                "", round(LAMBDA.HPD[2], digits = digit)), file = test.file, 
                append = TRUE)
            write(paste("beta        ", round(beta.est, digits = digit), 
                "", round(beta.sd, digits = digit), "", round(beta.MCerror, 
                  digits = digit), "", round(beta.HPD[1], digits = digit), 
                "", round(beta.HPD[2], digits = digit)), file = test.file, 
                append = TRUE)
            write(paste("sigma.alpha ", round(sigma.alpha.est, 
                digits = digit), "", round(sigma.alpha.sd, digits = digit), 
                "", round(sigma.alpha.MCerror, digits = digit), 
                "", round(sigma.alpha.HPD[1], digits = digit), 
                "", round(sigma.alpha.HPD[2], digits = digit)), 
                file = test.file, append = TRUE)
            if (SCO == FALSE) {
                write(paste("sigma.theta ", round(sigma.theta.est, 
                  digits = digit), "", round(sigma.theta.sd, 
                  digits = digit), "", round(sigma.theta.MCerror, 
                  digits = digit), "", round(sigma.theta.HPD[1], 
                  digits = digit), "", round(sigma.theta.HPD[2], 
                  digits = digit)), file = test.file, append = TRUE)
            }
            write(paste("S overall          ", round(S_overall.est, 
                digits = digit), "", round(S_overall.sd, digits = digit), 
                "", round(S_overall.MCerror, digits = digit), 
                "", round(S_overall.HPD[1], digits = digit), 
                "", round(S_overall.HPD[2], digits = digit)), 
                file = test.file, append = TRUE)
            write(paste("C overall          ", round(C_overall.est, 
                digits = digit), "", round(C_overall.sd, digits = digit), 
                "", round(C_overall.MCerror, digits = digit), 
                "", round(C_overall.HPD[1], digits = digit), 
                "", round(C_overall.HPD[2], digits = digit)), 
                file = test.file, append = TRUE)
            write(paste(""), file = test.file, append = TRUE)
            if (condInd == TRUE & Gold_Std == FALSE & model == 
                1) {
                write(paste(""), file = test.file, append = TRUE)
                write(paste("______________________________________________________"), 
                  file = test.file, append = TRUE)
                write(paste("\tReference standard (Point estimate =", 
                  point_estimate, ")\t "), file = test.file, 
                  append = TRUE)
                write(paste("______________________________________________________"), 
                  file = test.file, append = TRUE)
                write(paste(""), file = test.file, append = TRUE)
                write(paste("         Estimate Standard_Dev MC_error C.I._lower C.I._upper"), 
                  file = test.file, append = TRUE)
                write(paste(""), file = test.file, append = TRUE)
                if (Gold_se == TRUE & Gold_sp == FALSE) {
                  if (rs.length != 1) {
                    for (i in 1:rs.length) {
                      write(paste("S2 of Study(ies) ", sub_rs[[i + 
                        1]][1], "to", sub_rs[[i + 1]][length(sub_rs[[i + 
                        1]])], "(It was assumed to be perfect)"), 
                        file = test.file, append = TRUE)
                      write(paste("C2 of Study(ies) ", sub_rs[[i + 
                        1]][1], "to", sub_rs[[i + 1]][length(sub_rs[[i + 
                        1]])], "", round(C2.est[i], digits = digit), 
                        "", round(C2.sd[i], digits = digit), 
                        "", round(C2.MCerror[i], digits = digit), 
                        "", round(C2.HPD[i, 1], digits = digit), 
                        "", round(C2.HPD[i, 2], digits = digit)), 
                        file = test.file, append = TRUE)
                    }
                  }
                  else {
                    write(paste("S2           (It was assumed to be perfect)"), 
                      file = test.file, append = TRUE)
                    write(paste("C2       ", round(C2.est, digits = digit), 
                      "", round(C2.sd, digits = digit), "", round(C2.MCerror, 
                        digits = digit), "", round(C2.HPD[1], 
                        digits = digit), "", round(C2.HPD[2], 
                        digits = digit)), file = test.file, append = TRUE)
                  }
                }
                else {
                  if (Gold_sp == TRUE & Gold_se == FALSE) {
                    if (rs.length != 1) {
                      for (i in 1:rs.length) {
                        write(paste("S2 of Study(ies) ", sub_rs[[i + 
                          1]][1], "to", sub_rs[[i + 1]][length(sub_rs[[i + 
                          1]])], "", round(S2.est[i], digits = digit), 
                          "", round(S2.sd[i], digits = digit), 
                          "", round(S2.MCerror[i], digits = digit), 
                          "", round(S2.HPD[i, 1], digits = digit), 
                          "", round(S2.HPD[i, 2], digits = digit)), 
                          file = test.file, append = TRUE)
                        write(paste("C2 of Study(ies) ", sub_rs[[i + 
                          1]][1], "to", sub_rs[[i + 1]][length(sub_rs[[i + 
                          1]])], "(It was assumed to be perfect)"), 
                          file = test.file, append = TRUE)
                      }
                    }
                    else {
                      write(paste("S2       ", round(S2.est, 
                        digits = digit), "", round(S2.sd, digits = digit), 
                        "", round(S2.MCerror, digits = digit), 
                        "", round(S2.HPD[1], digits = digit), 
                        "", round(S2.HPD[2], digits = digit)), 
                        file = test.file, append = TRUE)
                      write(paste("C2           (It was assumed to be perfect)"), 
                        file = test.file, append = TRUE)
                    }
                  }
                  else {
                    if (rs.length != 1) {
                      for (i in 1:rs.length) {
                        write(paste("S2 of Study(ies) ", sub_rs[[i + 
                          1]][1], "to", sub_rs[[i + 1]][length(sub_rs[[i + 
                          1]])], "", round(S2.est[i], digits = digit), 
                          "", round(S2.sd[i], digits = digit), 
                          "", round(S2.MCerror[i], digits = digit), 
                          "", round(S2.HPD[i, 1], digits = digit), 
                          "", round(S2.HPD[i, 2], digits = digit)), 
                          file = test.file, append = TRUE)
                      }
                      for (i in 1:rs.length) {
                        write(paste("C2 of Study(ies) ", sub_rs[[i + 
                          1]][1], "to", sub_rs[[i + 1]][length(sub_rs[[i + 
                          1]])], "", round(C2.est[i], digits = digit), 
                          "", round(C2.sd[i], digits = digit), 
                          "", round(C2.MCerror[i], digits = digit), 
                          "", round(C2.HPD[i, 1], digits = digit), 
                          "", round(C2.HPD[i, 2], digits = digit)), 
                          file = test.file, append = TRUE)
                      }
                    }
                    else {
                      write(paste("S2       ", round(S2.est, 
                        digits = digit), "", round(S2.sd, digits = digit), 
                        "", round(S2.MCerror, digits = digit), 
                        "", round(S2.HPD[1], digits = digit), 
                        "", round(S2.HPD[2], digits = digit)), 
                        file = test.file, append = TRUE)
                      write(paste("C2       ", round(C2.est, 
                        digits = digit), "", round(C2.sd, digits = digit), 
                        "", round(C2.MCerror, digits = digit), 
                        "", round(C2.HPD[1], digits = digit), 
                        "", round(C2.HPD[2], digits = digit)), 
                        file = test.file, append = TRUE)
                    }
                  }
                }
            }
            else {
                if (condInd == TRUE & Gold_Std == FALSE & model == 
                  2) {
                  write(paste(""), file = test.file, append = TRUE)
                  write(paste("______________________________________________________"), 
                    file = test.file, append = TRUE)
                  write(paste("\tReference standard (Point estimate =", 
                    point_estimate, ")\t "), file = test.file, 
                    append = TRUE)
                  write(paste("______________________________________________________"), 
                    file = test.file, append = TRUE)
                  write(paste(""), file = test.file, append = TRUE)
                  write(paste("         Estimate Standard_Dev MC_error C.I._lower C.I._upper"), 
                    file = test.file, append = TRUE)
                  write(paste(""), file = test.file, append = TRUE)
                  if (rs.length != 1) {
                    for (i in 1:rs.length) {
                      write(paste("a1 of Study(ies) ", sub_rs[[i + 
                        1]][1], "to", sub_rs[[i + 1]][length(sub_rs[[i + 
                        1]])], "", round(a1.est[i], digits = digit), 
                        "", round(a1.sd[i], digits = digit), 
                        "", round(a1.MCerror[i], digits = digit), 
                        "", round(a1.HPD[i, 1], digits = digit), 
                        "", round(a1.HPD[i, 2], digits = digit)), 
                        file = test.file, append = TRUE)
                    }
                    for (i in 1:rs.length) {
                      write(paste("a0 of Study(ies) ", sub_rs[[i + 
                        1]][1], "to", sub_rs[[i + 1]][length(sub_rs[[i + 
                        1]])], "", round(a0.est[i], digits = digit), 
                        "", round(a0.sd[i], digits = digit), 
                        "", round(a0.MCerror[i], digits = digit), 
                        "", round(a0.HPD[i, 1], digits = digit), 
                        "", round(a0.HPD[i, 2], digits = digit)), 
                        file = test.file, append = TRUE)
                    }
                    write(paste(""), file = test.file, append = TRUE)
                    for (i in 1:rs.length) {
                      write(paste("S2 of Study(ies) ", sub_rs[[i + 
                        1]][1], "to", sub_rs[[i + 1]][length(sub_rs[[i + 
                        1]])], "", round(S2.est[i], digits = digit), 
                        "", round(S2.sd[i], digits = digit), 
                        "", round(S2.MCerror[i], digits = digit), 
                        "", round(S2.HPD[i, 1], digits = digit), 
                        "", round(S2.HPD[i, 2], digits = digit)), 
                        file = test.file, append = TRUE)
                    }
                    for (i in 1:rs.length) {
                      write(paste("C2 of Study(ies) ", sub_rs[[i + 
                        1]][1], "to", sub_rs[[i + 1]][length(sub_rs[[i + 
                        1]])], "", round(C2.est[i], digits = digit), 
                        "", round(C2.sd[i], digits = digit), 
                        "", round(C2.MCerror[i], digits = digit), 
                        "", round(C2.HPD[i, 1], digits = digit), 
                        "", round(C2.HPD[i, 2], digits = digit)), 
                        file = test.file, append = TRUE)
                    }
                  }
                  else {
                    write(paste("a1       ", round(a1.est, digits = digit), 
                      "", round(a1.sd, digits = digit), "", round(a1.MCerror, 
                        digits = digit), "", round(a1.HPD[1], 
                        digits = digit), "", round(a1.HPD[2], 
                        digits = digit)), file = test.file, append = TRUE)
                    write(paste("a0       ", round(a0.est, digits = digit), 
                      "", round(a0.sd, digits = digit), "", round(a0.MCerror, 
                        digits = digit), "", round(a0.HPD[1], 
                        digits = digit), "", round(a0.HPD[2], 
                        digits = digit)), file = test.file, append = TRUE)
                    write(paste(""), file = test.file, append = TRUE)
                    write(paste("S2       ", round(S2.est, digits = digit), 
                      "", round(S2.sd, digits = digit), "", round(S2.MCerror, 
                        digits = digit), "", round(S2.HPD[1], 
                        digits = digit), "", round(S2.HPD[2], 
                        digits = digit)), file = test.file, append = TRUE)
                    write(paste("C2       ", round(C2.est, digits = digit), 
                      "", round(C2.sd, digits = digit), "", round(C2.MCerror, 
                        digits = digit), "", round(C2.HPD[1], 
                        digits = digit), "", round(C2.HPD[2], 
                        digits = digit)), file = test.file, append = TRUE)
                  }
                }
                else {
                  if (condInd == FALSE) {
                    write(paste(""), file = test.file, append = TRUE)
                    write(paste("______________________________________________________"), 
                      file = test.file, append = TRUE)
                    write(paste("\tReference standard (Point estimate =", 
                      point_estimate, ")\t "), file = test.file, 
                      append = TRUE)
                    write(paste("______________________________________________________"), 
                      file = test.file, append = TRUE)
                    write(paste(""), file = test.file, append = TRUE)
                    write(paste("         Estimate Standard_Dev MC_error C.I._lower C.I._upper"), 
                      file = test.file, append = TRUE)
                    write(paste(""), file = test.file, append = TRUE)
                    write(paste("d1       ", round(d1.est, digits = digit), 
                      "", round(d1.sd, digits = digit), "", round(d1.MCerror, 
                        digits = digit), "", round(d1.HPD[1], 
                        digits = digit), "", round(d1.HPD[2], 
                        digits = digit)), file = test.file, append = TRUE)
                    write(paste("d0       ", round(d0.est, digits = digit), 
                      "", round(d0.sd, digits = digit), "", round(d0.MCerror, 
                        digits = digit), "", round(d0.HPD[1], 
                        digits = digit), "", round(d0.HPD[2], 
                        digits = digit)), file = test.file, append = TRUE)
                    if (rs.length != 1) {
                      for (i in 1:rs.length) {
                        write(paste("a1 of Study(ies) ", sub_rs[[i + 
                          1]][1], "to", sub_rs[[i + 1]][length(sub_rs[[i + 
                          1]])], "", round(a1.est[i], digits = digit), 
                          "", round(a1.sd[i], digits = digit), 
                          "", round(a1.MCerror[i], digits = digit), 
                          "", round(a1.HPD[i, 1], digits = digit), 
                          "", round(a1.HPD[i, 2], digits = digit)), 
                          file = test.file, append = TRUE)
                      }
                      for (i in 1:rs.length) {
                        write(paste("a0 of Study(ies) ", sub_rs[[i + 
                          1]][1], "to", sub_rs[[i + 1]][length(sub_rs[[i + 
                          1]])], "", round(a0.est[i], digits = digit), 
                          "", round(a0.sd[i], digits = digit), 
                          "", round(a0.MCerror[i], digits = digit), 
                          "", round(a0.HPD[i, 1], digits = digit), 
                          "", round(a0.HPD[i, 2], digits = digit)), 
                          file = test.file, append = TRUE)
                      }
                      write(paste(""), file = test.file, append = TRUE)
                      for (i in 1:rs.length) {
                        write(paste("S2 of Study(ies) ", sub_rs[[i + 
                          1]][1], "to", sub_rs[[i + 1]][length(sub_rs[[i + 
                          1]])], "", round(S2.est[i], digits = digit), 
                          "", round(S2.sd[i], digits = digit), 
                          "", round(S2.MCerror[i], digits = digit), 
                          "", round(S2.HPD[i, 1], digits = digit), 
                          "", round(S2.HPD[i, 2], digits = digit)), 
                          file = test.file, append = TRUE)
                      }
                      for (i in 1:rs.length) {
                        write(paste("C2 of Study(ies) ", sub_rs[[i + 
                          1]][1], "to", sub_rs[[i + 1]][length(sub_rs[[i + 
                          1]])], "", round(C2.est[i], digits = digit), 
                          "", round(C2.sd[i], digits = digit), 
                          "", round(C2.MCerror[i], digits = digit), 
                          "", round(C2.HPD[i, 1], digits = digit), 
                          "", round(C2.HPD[i, 2], digits = digit)), 
                          file = test.file, append = TRUE)
                      }
                      for (i in 1:rs.length) {
                        write(paste("b1 of Study(ies) ", sub_rs[[i + 
                          1]][1], "to", sub_rs[[i + 1]][length(sub_rs[[i + 
                          1]])], "", round(b1.est[i], digits = digit), 
                          "", round(b1.sd[i], digits = digit), 
                          "", round(b1.MCerror[i], digits = digit), 
                          "", round(b1.HPD[i, 1], digits = digit), 
                          "", round(b1.HPD[i, 2], digits = digit)), 
                          file = test.file, append = TRUE)
                      }
                      for (i in 1:rs.length) {
                        write(paste("b0 of Study(ies) ", sub_rs[[i + 
                          1]][1], "to", sub_rs[[i + 1]][length(sub_rs[[i + 
                          1]])], "", round(b0.est[i], digits = digit), 
                          "", round(b0.sd[i], digits = digit), 
                          "", round(b0.MCerror[i], digits = digit), 
                          "", round(b0.HPD[i, 1], digits = digit), 
                          "", round(b0.HPD[i, 2], digits = digit)), 
                          file = test.file, append = TRUE)
                      }
                    }
                    else {
                      write(paste("a1       ", round(a1.est, 
                        digits = digit), "", round(a1.sd, digits = digit), 
                        "", round(a1.MCerror, digits = digit), 
                        "", round(a1.HPD[1], digits = digit), 
                        "", round(a1.HPD[2], digits = digit)), 
                        file = test.file, append = TRUE)
                      write(paste("a0       ", round(a0.est, 
                        digits = digit), "", round(a0.sd, digits = digit), 
                        "", round(a0.MCerror, digits = digit), 
                        "", round(a0.HPD[1], digits = digit), 
                        "", round(a0.HPD[2], digits = digit)), 
                        file = test.file, append = TRUE)
                      write(paste("b1       ", round(b1.est, 
                        digits = digit), "", round(b1.sd, digits = digit), 
                        "", round(b1.MCerror, digits = digit), 
                        "", round(b1.HPD[1], digits = digit), 
                        "", round(b1.HPD[2], digits = digit)), 
                        file = test.file, append = TRUE)
                      write(paste("b0       ", round(b0.est, 
                        digits = digit), "", round(b0.sd, digits = digit), 
                        "", round(b0.MCerror, digits = digit), 
                        "", round(b0.HPD[1], digits = digit), 
                        "", round(b0.HPD[2], digits = digit)), 
                        file = test.file, append = TRUE)
                      write(paste("S2       ", round(S2.est, 
                        digits = digit), "", round(S2.sd, digits = digit), 
                        "", round(S2.MCerror, digits = digit), 
                        "", round(S2.HPD[1], digits = digit), 
                        "", round(S2.HPD[2], digits = digit)), 
                        file = test.file, append = TRUE)
                      write(paste("C2       ", round(C2.est, 
                        digits = digit), "", round(C2.sd, digits = digit), 
                        "", round(C2.MCerror, digits = digit), 
                        "", round(C2.HPD[1], digits = digit), 
                        "", round(C2.HPD[2], digits = digit)), 
                        file = test.file, append = TRUE)
                    }
                  }
                }
            }
            write(paste(""), file = test.file, append = TRUE)
            write(paste("______________________________________________________"), 
                file = test.file, append = TRUE)
            write(paste("______________________________________________________"), 
                file = test.file, append = TRUE)
            write(paste("\tWITHIN-STUDY PARAMETERS \t "), file = test.file, 
                append = TRUE)
            write(paste("______________________________________________________"), 
                file = test.file, append = TRUE)
            write(paste("______________________________________________________"), 
                file = test.file, append = TRUE)
            if (SCO == FALSE) {
                write(paste("______________________________________________________"), 
                  file = test.file, append = TRUE)
                write(paste("\ttheta \t "), file = test.file, 
                  append = TRUE)
                write(paste("______________________________________________________"), 
                  file = test.file, append = TRUE)
                write(paste(""), file = test.file, append = TRUE)
                write(paste("         Estimate Standard_Dev MC_error C.I._lower C.I._upper"), 
                  file = test.file, append = TRUE)
                write(paste(""), file = test.file, append = TRUE)
                for (i in 1:N) {
                  write(paste("Study ", i, "", round(theta.est[i], 
                    digits = digit), "", round(theta.sd[i], digits = digit), 
                    "", round(theta.MCerror[i], digits = digit), 
                    "", round(theta.HPD[i, 1], digits = digit), 
                    "", round(theta.HPD[i, 2], digits = digit)), 
                    file = test.file, append = TRUE)
                }
            }
            write(paste("______________________________________________________"), 
                file = test.file, append = TRUE)
            write(paste("\talpha \t "), file = test.file, append = TRUE)
            write(paste("______________________________________________________"), 
                file = test.file, append = TRUE)
            write(paste(""), file = test.file, append = TRUE)
            write(paste("         Estimate Standard_Dev MC_error C.I._lower C.I._upper"), 
                file = test.file, append = TRUE)
            write(paste(""), file = test.file, append = TRUE)
            for (i in 1:N) {
                write(paste("Study ", i, "", round(alpha.est[i], 
                  digits = digit), "", round(alpha.sd[i], digits = digit), 
                  "", round(alpha.MCerror[i], digits = digit), 
                  "", round(alpha.HPD[i, 1], digits = digit), 
                  "", round(alpha.HPD[i, 2], digits = digit)), 
                  file = test.file, append = TRUE)
            }
            write(paste("______________________________________________________"), 
                file = test.file, append = TRUE)
            write(paste("\tPrevalence \t "), file = test.file, 
                append = TRUE)
            write(paste("______________________________________________________"), 
                file = test.file, append = TRUE)
            write(paste(""), file = test.file, append = TRUE)
            write(paste("         Estimate Standard_Dev MC_error C.I._lower C.I._upper"), 
                file = test.file, append = TRUE)
            write(paste(""), file = test.file, append = TRUE)
            for (i in 1:N) {
                write(paste("Study ", i, "", round(PI.est[i], 
                  digits = digit), "", round(PI.sd[i], digits = digit), 
                  "", round(PI.MCerror[i], digits = digit), "", 
                  round(PI.HPD[i, 1], digits = digit), "", round(PI.HPD[i, 
                    2], digits = digit)), file = test.file, append = TRUE)
            }
            write(paste("______________________________________________________"), 
                file = test.file, append = TRUE)
            write(paste("\tSensitivity of test 1 (S1) \t "), 
                file = test.file, append = TRUE)
            write(paste("______________________________________________________"), 
                file = test.file, append = TRUE)
            write(paste(""), file = test.file, append = TRUE)
            write(paste("         Estimate Standard_Dev MC_error C.I._lower C.I._upper"), 
                file = test.file, append = TRUE)
            write(paste(""), file = test.file, append = TRUE)
            for (i in 1:N) {
                write(paste("Study ", i, "", round(S1.est[i], 
                  digits = digit), "", round(S1.sd[i], digits = digit), 
                  "", round(S1.MCerror[i], digits = digit), "", 
                  round(S1.HPD[i, 1], digits = digit), "", round(S1.HPD[i, 
                    2], digits = digit)), file = test.file, append = TRUE)
            }
            write(paste("______________________________________________________"), 
                file = test.file, append = TRUE)
            write(paste("\tSpecificity of test 1 (C1) \t "), 
                file = test.file, append = TRUE)
            write(paste("______________________________________________________"), 
                file = test.file, append = TRUE)
            write(paste(""), file = test.file, append = TRUE)
            write(paste("         Estimate Standard_Dev MC_error C.I._lower C.I._upper"), 
                file = test.file, append = TRUE)
            write(paste(""), file = test.file, append = TRUE)
            for (i in 1:N) {
                write(paste("Study ", i, "", round(C1.est[i], 
                  digits = digit), "", round(C1.sd[i], digits = digit), 
                  "", round(C1.MCerror[i], digits = digit), "", 
                  round(C1.HPD[i, 1], digits = digit), "", round(C1.HPD[i, 
                    2], digits = digit)), file = test.file, append = TRUE)
            }
            Num_study = c()
            for (i in 1:N) {
                Num_study = c(Num_study, paste("Study", i))
            }
            if (condInd == TRUE & Gold_Std == FALSE & model == 
                1) {
                if (Gold_se == TRUE & Gold_sp == FALSE) {
                  if (rs.length != 1) {
                    refstd_parameters = array(0, dim = c(rs.length, 
                      3, 1), dimnames = list(1:rs.length, c(paste(point_estimate, 
                      "estimate"), "HPD lower", "HPD upper"), 
                      "C2"))
                    refstd_parameters[, 1, 1] <- C2.est
                    refstd_parameters[, 2, 1] <- C2.HPD[, 1]
                    refstd_parameters[, 3, 1] <- C2.HPD[, 2]
                    long = length(alpha[, 1])
                    refstd_Parameters = array(0, c(long, rs.length, 
                      1))
                    refstd_Parameters[, , 1] <- C2
                    if (print_plot == TRUE) {
                      if (is.null(chain) == FALSE) {
                        file.pdf_RS = paste("RefStd trace plots for N =", 
                          round((iter.num * nb_chains - (burn_in) * 
                            nb_chains)/Thin, 0), ".pdf")
                        pdf(file.pdf_RS, paper = "a4", height = 20)
                        param = "C2"
                        par(mfcol = c(5, 2))
                        no_chains = length(chain)
                        iter_chain = round((iter.num * nb_chains - 
                          (burn_in) * nb_chains)/Thin, 0)/no_chains
                        longueur = 1:iter_chain
                        for (i in 1:rs.length) {
                          plot(x = longueur, y = refstd_Parameters[longueur, 
                            i, 1], type = "n", col = 1, ylab = paste(param, 
                            " of reference test ", i), xlab = "iteration number", 
                            main = paste("Thinning interval = ", 
                              thin.interval, "\n Total samplesize kept = ", 
                              (iter.num * nb_chains - (burn_in) * 
                                nb_chains)/Thin))
                          for (l in 1:length(chain)) {
                            lines(x = longueur, y = refstd_Parameters[l * 
                              longueur, i, 1], col = l)
                          }
                        }
                      }
                      else {
                        file.pdf_RS = paste("RefStd trace plots for N =", 
                          round((iter.num * nb_chains - (burn_in) * 
                            nb_chains)/Thin, 0), ".pdf")
                        pdf(file.pdf_RS, paper = "a4", height = 20)
                        param = "C2"
                        par(mfcol = c(5, 2))
                        longueur = 1:long
                        for (i in 1:rs.length) {
                          plot(x = longueur, y = refstd_Parameters[, 
                            i, 1], type = "l", col = "grey", 
                            ylab = paste(param, " of study ", 
                              i), xlab = "iteration number", 
                            main = paste("Thinning interval = ", 
                              thin.interval, "\n Total samplesize kept = ", 
                              (iter.num * nb_chains - (burn_in) * 
                                nb_chains)/Thin))
                        }
                      }
                      dev.off()
                      file.pdf_RS2 = paste("RefStd density plots for N =", 
                        round((iter.num * nb_chains - (burn_in) * 
                          nb_chains)/Thin, 0), ".pdf")
                      pdf(file.pdf_RS2, paper = "a4", height = 20)
                      param = "C2"
                      par(mfcol = c(5, 2))
                      longueur = 1:long
                      for (i in 1:rs.length) {
                        plot(density(refstd_Parameters[, i, 1]), 
                          lwd = 4, type = "l", col = "grey", 
                          main = paste(param, " of study ", i, 
                            " \n Thinning interval = ", thin.interval, 
                            "\n Total samplesize kept = ", (iter.num * 
                              nb_chains - (burn_in) * nb_chains)/Thin))
                      }
                      dev.off()
                    }
                  }
                  else {
                    refstd_parameters = matrix(0, ncol = 3, nrow = 1)
                    rownames(refstd_parameters) = c("C2")
                    colnames(refstd_parameters) = c(paste(point_estimate, 
                      "estimate"), "HPD.low", "HPD.high")
                    refstd_parameters[1, 1] <- C2.est
                    refstd_parameters[1, 2] <- C2.HPD[1]
                    refstd_parameters[1, 3] <- C2.HPD[2]
                    long = length(THETA)
                    refstd_Parameters = matrix(0, nrow = long, 
                      ncol = 1)
                    colnames(refstd_Parameters) = c("C2")
                    refstd_Parameters[, 1] <- C2
                    if (print_plot == TRUE) {
                      if (is.null(chain) == FALSE) {
                        file.pdf_RS = paste("RefStd trace plots for N =", 
                          round((iter.num * nb_chains - (burn_in) * 
                            nb_chains)/Thin, 0), ".pdf")
                        pdf(file.pdf_RS, paper = "a4", height = 20)
                        param = c("C2")
                        par(mfcol = c(5, 2))
                        no_chains = length(chain)
                        iter_chain = round((iter.num * nb_chains - 
                          (burn_in) * nb_chains)/Thin, 0)/no_chains
                        longueur = 1:iter_chain
                        plot(x = longueur, y = refstd_Parameters[longueur, 
                          1], type = "n", col = 1, ylab = paste(param), 
                          xlab = "iteration number", main = paste("Thinning interval = ", 
                            thin.interval, "\n Total samplesize kept = ", 
                            (iter.num * nb_chains - (burn_in) * 
                              nb_chains)/Thin))
                        for (l in 1:length(chain)) {
                          lines(x = longueur, y = refstd_Parameters[l * 
                            longueur, 1], col = l)
                        }
                      }
                      else {
                        file.pdf_RS = paste("RefStd trace plots for N =", 
                          round((iter.num * nb_chains - (burn_in) * 
                            nb_chains)/Thin, 0), ".pdf")
                        pdf(file.pdf_RS, paper = "a4", height = 20)
                        param = c("C2")
                        par(mfcol = c(5, 2))
                        longueur = 1:long
                        plot(x = longueur, y = refstd_Parameters[, 
                          1], type = "l", col = "grey", ylab = paste(param), 
                          xlab = "iteration number", main = paste("Thinning interval = ", 
                            thin.interval, "\n Total samplesize kept = ", 
                            (iter.num * nb_chains - (burn_in) * 
                              nb_chains)/Thin))
                      }
                      dev.off()
                      file.pdf_RS2 = paste("RefStd density plots for N =", 
                        round((iter.num * nb_chains - (burn_in) * 
                          nb_chains)/Thin, 0), ".pdf")
                      pdf(file.pdf_RS2, paper = "a4", height = 20)
                      param = c("C2")
                      par(mfcol = c(5, 2))
                      longueur = 1:long
                      plot(density(refstd_Parameters[, 1]), lwd = 4, 
                        type = "l", col = "grey", main = paste(param, 
                          " of reference standard \n Thinning interval = ", 
                          thin.interval, "\n Total samplesize kept = ", 
                          (iter.num * nb_chains - (burn_in) * 
                            nb_chains)/Thin))
                      dev.off()
                    }
                  }
                }
                else {
                  if (Gold_sp == TRUE & Gold_se == FALSE) {
                    if (rs.length != 1) {
                      refstd_parameters = array(0, dim = c(rs.length, 
                        3, 1), dimnames = list(1:rs.length, c(paste(point_estimate, 
                        "estimate"), "HPD lower", "HPD upper"), 
                        "S2"))
                      refstd_parameters[, 1, 1] <- S2.est
                      refstd_parameters[, 2, 1] <- S2.HPD[, 1]
                      refstd_parameters[, 3, 1] <- S2.HPD[, 2]
                      long = length(alpha[, 1])
                      refstd_Parameters = array(0, c(long, rs.length, 
                        1))
                      refstd_Parameters[, , 1] <- S2
                      if (print_plot == TRUE) {
                        if (is.null(chain) == FALSE) {
                          file.pdf_RS = paste("RefStd trace plots for N =", 
                            round((iter.num * nb_chains - (burn_in) * 
                              nb_chains)/Thin, 0), ".pdf")
                          pdf(file.pdf_RS, paper = "a4", height = 20)
                          param = c("S2")
                          par(mfcol = c(5, 2))
                          no_chains = length(chain)
                          iter_chain = round((iter.num * nb_chains - 
                            (burn_in) * nb_chains)/Thin, 0)/no_chains
                          longueur = 1:iter_chain
                          for (i in 1:rs.length) {
                            plot(x = longueur, y = refstd_Parameters[longueur, 
                              i, 1], type = "n", col = 1, ylab = paste(param, 
                              " of reference test ", i), xlab = "iteration number", 
                              main = paste("Thinning interval = ", 
                                thin.interval, "\n Total samplesize kept = ", 
                                (iter.num * nb_chains - (burn_in) * 
                                  nb_chains)/Thin))
                            for (l in 1:length(chain)) {
                              lines(x = longueur, y = refstd_Parameters[l * 
                                longueur, i, 1], col = l)
                            }
                          }
                        }
                        else {
                          file.pdf_RS = paste("RefStd trace plots for N =", 
                            round((iter.num * nb_chains - (burn_in) * 
                              nb_chains)/Thin, 0), ".pdf")
                          pdf(file.pdf_RS, paper = "a4", height = 20)
                          param = c("S2")
                          par(mfcol = c(5, 2))
                          longueur = 1:long
                          for (i in 1:rs.length) {
                            plot(x = longueur, y = refstd_Parameters[, 
                              i, 1], type = "l", col = "grey", 
                              ylab = paste(param, " of reference test ", 
                                i), xlab = "iteration number", 
                              main = paste("Thinning interval = ", 
                                thin.interval, "\n Total samplesize kept = ", 
                                (iter.num * nb_chains - (burn_in) * 
                                  nb_chains)/Thin))
                          }
                        }
                        dev.off()
                        file.pdf_RS2 = paste("RefStd density plots for N =", 
                          round((iter.num * nb_chains - (burn_in) * 
                            nb_chains)/Thin, 0), ".pdf")
                        pdf(file.pdf_RS2, paper = "a4", height = 20)
                        param = c("S2")
                        par(mfcol = c(5, 2))
                        longueur = 1:long
                        for (i in 1:rs.length) {
                          plot(density(refstd_Parameters[, i, 
                            1]), lwd = 4, type = "l", col = "grey", 
                            main = paste(param, " of reference standard ", 
                              i, " \n Thinning interval = ", 
                              thin.interval, "\n Total samplesize kept = ", 
                              (iter.num * nb_chains - (burn_in) * 
                                nb_chains)/Thin))
                        }
                        dev.off()
                      }
                    }
                    else {
                      refstd_parameters = matrix(0, ncol = 3, 
                        nrow = 1)
                      rownames(refstd_parameters) = c("S2")
                      colnames(refstd_parameters) = c(paste(point_estimate, 
                        "estimate"), "HPD.low", "HPD.high")
                      refstd_parameters[1, 1] <- S2.est
                      refstd_parameters[1, 2] <- S2.HPD[1]
                      refstd_parameters[1, 3] <- S2.HPD[2]
                      long = length(THETA)
                      refstd_Parameters = matrix(0, nrow = long, 
                        ncol = 1)
                      colnames(refstd_Parameters) = c("S2")
                      refstd_Parameters[, 1] <- S2
                      if (print_plot == TRUE) {
                        if (is.null(chain) == FALSE) {
                          file.pdf_RS = paste("RefStd trace plots for N =", 
                            round((iter.num * nb_chains - (burn_in) * 
                              nb_chains)/Thin, 0), ".pdf")
                          pdf(file.pdf_RS, paper = "a4", height = 20)
                          param = c("S2")
                          par(mfcol = c(5, 2))
                          no_chains = length(chain)
                          iter_chain = round((iter.num * nb_chains - 
                            (burn_in) * nb_chains)/Thin, 0)/no_chains
                          longueur = 1:iter_chain
                          plot(x = longueur, y = refstd_Parameters[longueur, 
                            1], type = "n", col = 1, ylab = paste(param), 
                            xlab = "iteration number", main = paste("Thinning interval = ", 
                              thin.interval, "\n Total samplesize kept = ", 
                              (iter.num * nb_chains - (burn_in) * 
                                nb_chains)/Thin))
                          for (l in 1:length(chain)) {
                            lines(x = longueur, y = refstd_Parameters[l * 
                              longueur, 1], col = l)
                          }
                        }
                        else {
                          file.pdf_RS = paste("RefStd trace plots for N =", 
                            round((iter.num * nb_chains - (burn_in) * 
                              nb_chains)/Thin, 0), ".pdf")
                          pdf(file.pdf_RS, paper = "a4", height = 20)
                          param = c("S2")
                          par(mfcol = c(5, 2))
                          longueur = 1:long
                          plot(x = longueur, y = refstd_Parameters[, 
                            1], type = "l", col = "grey", ylab = paste(param), 
                            xlab = "iteration number", main = paste("Thinning interval = ", 
                              thin.interval, "\n Total samplesize kept = ", 
                              (iter.num * nb_chains - (burn_in) * 
                                nb_chains)/Thin))
                        }
                        dev.off()
                        file.pdf_RS2 = paste("RefStd density plots for N =", 
                          round((iter.num * nb_chains - (burn_in) * 
                            nb_chains)/Thin, 0), ".pdf")
                        pdf(file.pdf_RS2, paper = "a4", height = 20)
                        param = c("S2")
                        par(mfcol = c(5, 2))
                        longueur = 1:long
                        plot(density(refstd_Parameters[, 1]), 
                          lwd = 4, type = "l", col = "grey", 
                          main = paste(param, " \n Thinning interval = ", 
                            thin.interval, "\n Total samplesize kept = ", 
                            (iter.num * nb_chains - (burn_in) * 
                              nb_chains)/Thin))
                        dev.off()
                      }
                    }
                  }
                  else {
                    if (rs.length != 1) {
                      refstd_parameters = array(0, dim = c(rs.length, 
                        3, 2), dimnames = list(1:rs.length, c(paste(point_estimate, 
                        "estimate"), "HPD lower", "HPD upper"), 
                        c("S2", "C2")))
                      refstd_parameters[, 1, 1] <- S2.est
                      refstd_parameters[, 2, 1] <- S2.HPD[, 1]
                      refstd_parameters[, 3, 1] <- S2.HPD[, 2]
                      refstd_parameters[, 1, 2] <- C2.est
                      refstd_parameters[, 2, 2] <- C2.HPD[, 1]
                      refstd_parameters[, 3, 2] <- C2.HPD[, 2]
                      long = length(alpha[, 1])
                      refstd_Parameters = array(0, c(long, rs.length, 
                        2))
                      refstd_Parameters[, , 1] <- S2
                      refstd_Parameters[, , 2] <- C2
                      if (print_plot == TRUE) {
                        if (is.null(chain) == FALSE) {
                          file.pdf_RS = paste("RefStd trace plots for N =", 
                            round((iter.num * nb_chains - (burn_in) * 
                              nb_chains)/Thin, 0), ".pdf")
                          pdf(file.pdf_RS, paper = "a4", height = 20)
                          param = c("S2", "C2")
                          par(mfcol = c(5, 2))
                          no_chains = length(chain)
                          iter_chain = round((iter.num * nb_chains - 
                            (burn_in) * nb_chains)/Thin, 0)/no_chains
                          longueur = 1:iter_chain
                          for (j in 1:2) {
                            for (i in 1:rs.length) {
                              plot(x = longueur, y = refstd_Parameters[longueur, 
                                i, j], type = "n", col = 1, ylab = paste(param[j], 
                                " of reference test ", i), xlab = "iteration number", 
                                main = paste("Thinning interval = ", 
                                  thin.interval, "\n Total samplesize kept = ", 
                                  (iter.num * nb_chains - (burn_in) * 
                                    nb_chains)/Thin))
                              for (l in 1:length(chain)) {
                                lines(x = longueur, y = refstd_Parameters[l * 
                                  longueur, i, j], col = l)
                              }
                            }
                          }
                        }
                        else {
                          file.pdf_RS = paste("RefStd trace plots for N =", 
                            round((iter.num * nb_chains - (burn_in) * 
                              nb_chains)/Thin, 0), ".pdf")
                          pdf(file.pdf_RS, paper = "a4", height = 20)
                          param = c("S2", "C2")
                          par(mfcol = c(5, 2))
                          longueur = 1:long
                          for (j in 1:2) {
                            for (i in 1:rs.length) {
                              plot(x = longueur, y = refstd_Parameters[, 
                                i, j], type = "l", col = "grey", 
                                ylab = paste(param[j], " of reference test ", 
                                  i), xlab = "iteration number", 
                                main = paste("Thinning interval = ", 
                                  thin.interval, "\n Total samplesize kept = ", 
                                  (iter.num * nb_chains - (burn_in) * 
                                    nb_chains)/Thin))
                            }
                          }
                        }
                        dev.off()
                        file.pdf_RS2 = paste("RefStd density plots for N =", 
                          round((iter.num * nb_chains - (burn_in) * 
                            nb_chains)/Thin, 0), ".pdf")
                        pdf(file.pdf_RS2, paper = "a4", height = 20)
                        param = c("S2", "C2")
                        par(mfcol = c(5, 2))
                        longueur = 1:long
                        for (j in 1:2) {
                          for (i in 1:rs.length) {
                            plot(density(refstd_Parameters[, 
                              i, j]), lwd = 4, type = "l", col = "grey", 
                              main = paste(param[j], " of study ", 
                                i, " \n Thinning interval = ", 
                                thin.interval, "\n Total samplesize kept = ", 
                                (iter.num * nb_chains - (burn_in) * 
                                  nb_chains)/Thin))
                          }
                        }
                        dev.off()
                      }
                    }
                    else {
                      refstd_parameters = matrix(0, ncol = 3, 
                        nrow = 2)
                      rownames(refstd_parameters) = c("S2", "C2")
                      colnames(refstd_parameters) = c(paste(point_estimate, 
                        "estimate"), "HPD.low", "HPD.high")
                      refstd_parameters[1, 1] <- S2.est
                      refstd_parameters[1, 2] <- S2.HPD[1]
                      refstd_parameters[1, 3] <- S2.HPD[2]
                      refstd_parameters[2, 1] <- C2.est
                      refstd_parameters[2, 2] <- C2.HPD[1]
                      refstd_parameters[2, 3] <- C2.HPD[2]
                      long = length(THETA)
                      refstd_Parameters = matrix(0, nrow = long, 
                        ncol = 2)
                      colnames(refstd_Parameters) = c("S2", "C2")
                      refstd_Parameters[, 1] <- S2
                      refstd_Parameters[, 2] <- C2
                      if (print_plot == TRUE) {
                        if (is.null(chain) == FALSE) {
                          file.pdf_RS = paste("RefStd trace plots for N =", 
                            round((iter.num * nb_chains - (burn_in) * 
                              nb_chains)/Thin, 0), ".pdf")
                          pdf(file.pdf_RS, paper = "a4", height = 20)
                          param = c("S2", "C2")
                          par(mfcol = c(5, 2))
                          no_chains = length(chain)
                          iter_chain = round((iter.num * nb_chains - 
                            (burn_in) * nb_chains)/Thin, 0)/no_chains
                          longueur = 1:iter_chain
                          for (j in 1:2) {
                            plot(x = longueur, y = refstd_Parameters[longueur, 
                              j], type = "n", col = 1, ylab = paste(param[j]), 
                              xlab = "iteration number", main = paste("Thinning interval = ", 
                                thin.interval, "\n Total samplesize kept = ", 
                                (iter.num * nb_chains - (burn_in) * 
                                  nb_chains)/Thin))
                            for (l in 1:length(chain)) {
                              lines(x = longueur, y = refstd_Parameters[l * 
                                longueur, j], col = l)
                            }
                          }
                        }
                        else {
                          file.pdf_RS = paste("RefStd trace plots for N =", 
                            round((iter.num * nb_chains - (burn_in) * 
                              nb_chains)/Thin, 0), ".pdf")
                          pdf(file.pdf_RS, paper = "a4", height = 20)
                          param = c("S2", "C2")
                          par(mfcol = c(5, 2))
                          longueur = 1:long
                          for (j in 1:2) {
                            plot(x = longueur, y = refstd_Parameters[, 
                              j], type = "l", col = "grey", ylab = paste(param[j]), 
                              xlab = "iteration number", main = paste("Thinning interval = ", 
                                thin.interval, "\n Total samplesize kept = ", 
                                (iter.num * nb_chains - (burn_in) * 
                                  nb_chains)/Thin))
                          }
                        }
                        dev.off()
                        file.pdf_RS2 = paste("RefStd density plots for N =", 
                          round((iter.num * nb_chains - (burn_in) * 
                            nb_chains)/Thin, 0), ".pdf")
                        pdf(file.pdf_RS2, paper = "a4", height = 20)
                        param = c("S2", "C2")
                        par(mfcol = c(5, 2))
                        for (j in 1:2) {
                          longueur = 1:long
                          plot(density(refstd_Parameters[, j]), 
                            lwd = 4, type = "l", col = "grey", 
                            main = paste(param[j], " \n Thinning interval = ", 
                              thin.interval, "\n Total samplesize kept = ", 
                              (iter.num * nb_chains - (burn_in) * 
                                nb_chains)/Thin))
                        }
                        dev.off()
                      }
                    }
                  }
                }
            }
            else {
                if (condInd == TRUE & Gold_Std == FALSE & model == 
                  2) {
                  if (rs.length != 1) {
                    refstd_parameters = array(0, dim = c(rs.length, 
                      3, 4), dimnames = list(1:rs.length, c(paste(point_estimate, 
                      "estimate"), "HPD lower", "HPD upper"), 
                      c("S2", "C2", "a1", "a0")))
                    refstd_parameters[, 1, 1] <- S2.est
                    refstd_parameters[, 2, 1] <- S2.HPD[, 1]
                    refstd_parameters[, 3, 1] <- S2.HPD[, 2]
                    refstd_parameters[, 1, 2] <- C2.est
                    refstd_parameters[, 2, 2] <- C2.HPD[, 1]
                    refstd_parameters[, 3, 2] <- C2.HPD[, 2]
                    refstd_parameters[, 1, 3] <- a1.est
                    refstd_parameters[, 2, 3] <- a1.HPD[, 1]
                    refstd_parameters[, 3, 3] <- a1.HPD[, 2]
                    refstd_parameters[, 1, 4] <- a0.est
                    refstd_parameters[, 2, 4] <- a0.HPD[, 1]
                    refstd_parameters[, 3, 4] <- a0.HPD[, 2]
                    long = length(alpha[, 1])
                    refstd_Parameters = array(0, dim = c(long, 
                      rs.length, 4))
                    refstd_Parameters[, , 1] <- S2
                    refstd_Parameters[, , 2] <- C2
                    refstd_Parameters[, , 3] <- a1
                    refstd_Parameters[, , 4] <- a0
                    if (print_plot == TRUE) {
                      file.pdf_RS = paste("RefStd trace plots for N =", 
                        round((iter.num * nb_chains - (burn_in) * 
                          nb_chains)/Thin, 0), ".pdf")
                      pdf(file.pdf_RS, paper = "a4", height = 20)
                      param = c("S2", "C2", "a1", "a0")
                      par(mfcol = c(5, 2))
                      for (j in 1:4) {
                        longueur = 1:long
                        for (i in 1:rs.length) {
                          plot(x = longueur, y = refstd_Parameters[, 
                            i, j], type = "l", col = "grey", 
                            ylab = paste(param[j], " of study ", 
                              i), xlab = "iteration number", 
                            main = paste("Thinning interval = ", 
                              thin.interval, "\n Total samplesize kept = ", 
                              (iter.num * nb_chains - (burn_in) * 
                                nb_chains)/Thin))
                          abline(a = refstd_parameters[i, 2, 
                            j], b = 0, col = "green", lwd = 3)
                          abline(a = refstd_parameters[i, 3, 
                            j], b = 0, col = "green", lwd = 3)
                        }
                      }
                      dev.off()
                      file.pdf_RS2 = paste("RefStd density plots for N =", 
                        round((iter.num * nb_chains - (burn_in) * 
                          nb_chains)/Thin, 0), ".pdf")
                      pdf(file.pdf_RS2, paper = "a4", height = 20)
                      param = c("S2", "C2", "a1", "a0")
                      par(mfcol = c(5, 2))
                      longueur = 1:long
                      for (j in 1:4) {
                        for (i in 1:rs.length) {
                          plot(density(refstd_Parameters[, i, 
                            j]), lwd = 4, type = "l", col = "grey", 
                            main = paste(param[j], " of study ", 
                              i, " \n Thinning interval = ", 
                              thin.interval, "\n Total samplesize kept = ", 
                              (iter.num * nb_chains - (burn_in) * 
                                nb_chains)/Thin))
                        }
                      }
                      dev.off()
                    }
                  }
                  else {
                    refstd_parameters = matrix(0, ncol = 3, nrow = 4)
                    rownames(refstd_parameters) = c("S2", "C2", 
                      "a1", "a0")
                    colnames(refstd_parameters) = c(paste(point_estimate, 
                      "estimate"), "HPD.low", "HPD.high")
                    refstd_parameters[1, 1] <- S2.est
                    refstd_parameters[1, 2] <- S2.HPD[1]
                    refstd_parameters[1, 3] <- S2.HPD[2]
                    refstd_parameters[2, 1] <- C2.est
                    refstd_parameters[2, 2] <- C2.HPD[1]
                    refstd_parameters[2, 3] <- C2.HPD[2]
                    refstd_parameters[3, 1] <- a1.est
                    refstd_parameters[3, 2] <- a1.HPD[1]
                    refstd_parameters[3, 3] <- a1.HPD[2]
                    refstd_parameters[4, 1] <- a0.est
                    refstd_parameters[4, 2] <- a0.HPD[1]
                    refstd_parameters[4, 3] <- a0.HPD[2]
                    long = length(THETA)
                    refstd_Parameters = matrix(0, nrow = long, 
                      ncol = 4)
                    colnames(refstd_Parameters) = c("S2", "C2", 
                      "a1", "a0")
                    refstd_Parameters[, 1] <- S2
                    refstd_Parameters[, 2] <- C2
                    refstd_Parameters[, 3] <- a1
                    refstd_Parameters[, 4] <- a0
                    if (print_plot == TRUE) {
                      file.pdf_RS = paste("RefStd trace plots for N =", 
                        round((iter.num * nb_chains - (burn_in) * 
                          nb_chains)/Thin, 0), ".pdf")
                      pdf(file.pdf_RS, paper = "a4", height = 20)
                      Param = c("S2", "C2", "a1", "a0")
                      par(mfcol = c(5, 2))
                      for (j in 1:4) {
                        longueur = 1:long
                        plot(x = longueur, y = refstd_Parameters[, 
                          j], type = "l", col = "grey", ylab = paste(Param[j]), 
                          xlab = "iteration number", main = paste("Thinning interval = ", 
                            thin.interval, "\n Total samplesize kept = ", 
                            (iter.num * nb_chains - (burn_in) * 
                              nb_chains)/Thin))
                        abline(a = refstd_parameters[j, 2], b = 0, 
                          col = "green", lwd = 3)
                        abline(a = refstd_parameters[j, 3], b = 0, 
                          col = "green", lwd = 3)
                      }
                      dev.off()
                      file.pdf_RS2 = paste("Density plots for N =", 
                        round((iter.num * nb_chains - (burn_in) * 
                          nb_chains)/Thin, 0), ".pdf")
                      pdf(file.pdf_RS2, paper = "a4", height = 20)
                      param = c("S2", "C2", "a1", "a0")
                      par(mfcol = c(5, 2))
                      for (j in 1:4) {
                        longueur = 1:long
                        plot(density(refstd_Parameters[, j]), 
                          lwd = 4, type = "l", col = "grey", 
                          main = paste(param[j], " of study ", 
                            j, " \n Thinning interval = ", thin.interval, 
                            "\n Total samplesize kept = ", (iter.num * 
                              nb_chains - (burn_in) * nb_chains)/Thin))
                      }
                      dev.off()
                    }
                  }
                }
            }
            if (SCO == FALSE) {
                parameters = array(0, dim = c(N, 3, 5), dimnames = list(Num_study, 
                  c(paste(point_estimate, "estimate"), "HPD lower", 
                    "HPD upper"), c("theta", "alpha", "pi", "S1", 
                    "C1")))
                parameters[, 1, 1] <- theta.est
                parameters[, 2, 1] <- theta.HPD[, 1]
                parameters[, 3, 1] <- theta.HPD[, 2]
                parameters[, 1, 2] <- alpha.est
                parameters[, 2, 2] <- alpha.HPD[, 1]
                parameters[, 3, 2] <- alpha.HPD[, 2]
                parameters[, 1, 3] <- PI.est
                parameters[, 2, 3] <- PI.HPD[, 1]
                parameters[, 3, 3] <- PI.HPD[, 2]
                parameters[, 1, 4] <- S1.est
                parameters[, 2, 4] <- S1.HPD[, 1]
                parameters[, 3, 4] <- S1.HPD[, 2]
                parameters[, 1, 5] <- C1.est
                parameters[, 2, 5] <- C1.HPD[, 1]
                parameters[, 3, 5] <- C1.HPD[, 2]
                long = length(alpha[, 1])
                Parameters = array(0, c(long, N, 5))
                Parameters[, , 1] <- theta
                Parameters[, , 2] <- alpha
                Parameters[, , 3] <- PI
                Parameters[, , 4] <- S1
                Parameters[, , 5] <- C1
                parameter = matrix(0, ncol = 3, nrow = 7)
                rownames(parameter) = c("THETA", "LAMBDA", "beta", 
                  "sigma.alpha", "sigma.theta", "S Overall", 
                  "C Overall")
                colnames(parameter) = c(paste(point_estimate, 
                  "estimate"), "HPD.low", "HPD.high")
                parameter[1, 1] <- THETA.est
                parameter[1, 2] <- THETA.HPD[1]
                parameter[1, 3] <- THETA.HPD[2]
                parameter[2, 1] <- LAMBDA.est
                parameter[2, 2] <- LAMBDA.HPD[1]
                parameter[2, 3] <- LAMBDA.HPD[2]
                parameter[3, 1] <- beta.est
                parameter[3, 2] <- beta.HPD[1]
                parameter[3, 3] <- beta.HPD[2]
                parameter[4, 1] <- sigma.alpha.est
                parameter[4, 2] <- sigma.alpha.HPD[1]
                parameter[4, 3] <- sigma.alpha.HPD[2]
                parameter[5, 1] <- sigma.theta.est
                parameter[5, 2] <- sigma.theta.HPD[1]
                parameter[5, 3] <- sigma.theta.HPD[2]
                parameter[6, 1] <- S_overall.est
                parameter[6, 2] <- S_overall.HPD[1]
                parameter[6, 3] <- S_overall.HPD[2]
                parameter[7, 1] <- C_overall.est
                parameter[7, 2] <- C_overall.HPD[1]
                parameter[7, 3] <- C_overall.HPD[2]
                long = length(THETA)
                Parameter = matrix(0, nrow = long, ncol = 7)
                colnames(Parameter) = c("THETA", "LAMBDA", "beta", 
                  "sigma.alpha", "sigma.theta", "S Overall", 
                  "C Overall")
                Parameter[, 1] <- THETA
                Parameter[, 2] <- LAMBDA
                Parameter[, 3] <- beta
                Parameter[, 4] <- sigma.alpha
                Parameter[, 5] <- sigma.theta
                Parameter[, 6] <- S_overall
                Parameter[, 7] <- C_overall
            }
            else {
                if (SCO == TRUE) {
                  parameters = array(0, dim = c(N, 3, 4), dimnames = list(Num_study, 
                    c(paste(point_estimate, "estimate"), "HPD lower", 
                      "HPD upper"), c("alpha", "pi", "S1", "C1")))
                  parameters[, 1, 1] <- alpha.est
                  parameters[, 2, 1] <- alpha.HPD[, 1]
                  parameters[, 3, 1] <- alpha.HPD[, 2]
                  parameters[, 1, 2] <- PI.est
                  parameters[, 2, 2] <- PI.HPD[, 1]
                  parameters[, 3, 2] <- PI.HPD[, 2]
                  parameters[, 1, 3] <- S1.est
                  parameters[, 2, 3] <- S1.HPD[, 1]
                  parameters[, 3, 3] <- S1.HPD[, 2]
                  parameters[, 1, 4] <- C1.est
                  parameters[, 2, 4] <- C1.HPD[, 1]
                  parameters[, 3, 4] <- C1.HPD[, 2]
                  long = length(alpha[, 1])
                  Parameters = array(0, c(long, N, 4))
                  Parameters[, , 1] <- alpha
                  Parameters[, , 2] <- PI
                  Parameters[, , 3] <- S1
                  Parameters[, , 4] <- C1
                  parameter = matrix(0, ncol = 3, nrow = 6)
                  rownames(parameter) = c("THETA", "LAMBDA", 
                    "beta", "sigma.alpha", "S Overall", "C Overall")
                  colnames(parameter) = c(paste(point_estimate, 
                    "estimate"), "HPD.low", "HPD.high")
                  parameter[1, 1] <- THETA.est
                  parameter[1, 2] <- THETA.HPD[1]
                  parameter[1, 3] <- THETA.HPD[2]
                  parameter[2, 1] <- LAMBDA.est
                  parameter[2, 2] <- LAMBDA.HPD[1]
                  parameter[2, 3] <- LAMBDA.HPD[2]
                  parameter[3, 1] <- beta.est
                  parameter[3, 2] <- beta.HPD[1]
                  parameter[3, 3] <- beta.HPD[2]
                  parameter[4, 1] <- sigma.alpha.est
                  parameter[4, 2] <- sigma.alpha.HPD[1]
                  parameter[4, 3] <- sigma.alpha.HPD[2]
                  parameter[5, 1] <- S_overall.est
                  parameter[5, 2] <- S_overall.HPD[1]
                  parameter[5, 3] <- S_overall.HPD[2]
                  parameter[6, 1] <- C_overall.est
                  parameter[6, 2] <- C_overall.HPD[1]
                  parameter[6, 3] <- C_overall.HPD[2]
                  long = length(THETA)
                  Parameter = matrix(0, nrow = long, ncol = 6)
                  colnames(Parameter) = c("THETA", "LAMBDA", 
                    "beta", "sigma.alpha", "S Overall", "C Overall")
                  Parameter[, 1] <- THETA
                  Parameter[, 2] <- LAMBDA
                  Parameter[, 3] <- beta
                  Parameter[, 4] <- sigma.alpha
                  Parameter[, 5] <- S_overall
                  Parameter[, 6] <- C_overall
                }
            }
            if (print_plot == TRUE) {
                if (SCO == FALSE) {
                  if (is.null(chain) == FALSE) {
                    file.pdf5 = paste("Trace plots for N =", 
                      round((iter.num * nb_chains - (burn_in) * 
                        nb_chains)/Thin, 0), ".pdf")
                    pdf(file.pdf5, paper = "a4", height = 20)
                    param = c("theta", "alpha", "PI", "S1", "C1")
                    Param = c("Capital Theta", "Capital Lambda", 
                      "beta", "~sigma[alpha]", "~sigma[theta]", 
                      "S Overall", "C Overall")
                    no_chains = length(chain)
                    iter_chain = round((iter.num * nb_chains - 
                      (burn_in) * nb_chains)/Thin, 0)/no_chains
                    min_param = c(min(Parameters[, , 1]), min(Parameters[, 
                      , 2]), min(Parameters[, , 3]), min(Parameters[, 
                      , 4]), min(Parameters[, , 5]))
                    max_param = c(max(Parameters[, , 1]), max(Parameters[, 
                      , 2]), max(Parameters[, , 3]), max(Parameters[, 
                      , 4]), max(Parameters[, , 5]))
                    dlag = (max_param - min_param)/100
                    range_param = numeric()
                    for (j in 1:5) {
                      range_param = cbind(range_param, seq(min_param[j] + 
                        dlag[j]/2, max_param[j] - dlag[j]/2, 
                        by = dlag[j]))
                    }
                    par(mfcol = c(5, 2))
                    longueur = 1:iter_chain
                    for (j in 1:5) {
                      for (i in 1:N) {
                        plot(x = longueur, y = Parameters[longueur, 
                          i, j], type = "n", col = 1, ylab = paste(param[j], 
                          " of study ", i), xlab = "iteration number", 
                          main = paste("Thinning interval = ", 
                            thin.interval, "\n Total samplesize kept = ", 
                            (iter.num * nb_chains - (burn_in) * 
                              nb_chains)/Thin), ylim = range(range_param[, 
                            j]))
                        for (l in 1:length(chain)) {
                          lines(x = longueur, y = Parameters[l * 
                            longueur, i, j], col = l)
                        }
                      }
                    }
                    min_Param = c(min(Parameter[, 1]), min(Parameter[, 
                      2]), min(Parameter[, 3]), min(Parameter[, 
                      4]), min(Parameter[, 5]), min(Parameter[, 
                      6]), min(Parameter[, 7]))
                    max_Param = c(max(Parameter[, 1]), max(Parameter[, 
                      2]), max(Parameter[, 3]), max(Parameter[, 
                      4]), max(Parameter[, 5]), max(Parameter[, 
                      6]), max(Parameter[, 7]))
                    dlag = (max_Param - min_Param)/100
                    range_Param = numeric()
                    for (j in 1:7) {
                      range_Param = cbind(range_Param, seq(min_Param[j] + 
                        dlag[j]/2, max_Param[j] - dlag[j]/2, 
                        by = dlag[j]))
                    }
                    for (j in 1:7) {
                      plot(x = longueur, y = Parameter[longueur, 
                        j], type = "n", col = 1, ylab = paste(Param[j]), 
                        xlab = "iteration number", main = paste("Thinning interval = ", 
                          thin.interval, "\n Total samplesize kept = ", 
                          (iter.num * nb_chains - (burn_in) * 
                            nb_chains)/Thin), ylim = range(range_Param[, 
                          j]))
                      for (l in 1:length(chain)) {
                        lines(x = longueur, y = Parameter[l * 
                          longueur, j], col = l)
                      }
                    }
                    dev.off()
                  }
                  else {
                    file.pdf2 = paste("Trace plots for N =", 
                      round((iter.num * nb_chains - (burn_in) * 
                        nb_chains)/Thin, 0), ".pdf")
                    pdf(file.pdf2, paper = "a4", height = 20)
                    param = c("theta", "alpha", "PI", "S1", "C1")
                    Param = c("Capital Theta", "Capital Lambda", 
                      "beta", "~sigma[alpha]", "~sigma[theta]", 
                      "S Overall", "C Overall")
                    min_param = c(min(Parameters[, , 1]), min(Parameters[, 
                      , 2]), min(Parameters[, , 3]), min(Parameters[, 
                      , 4]), min(Parameters[, , 5]))
                    max_param = c(max(Parameters[, , 1]), max(Parameters[, 
                      , 2]), max(Parameters[, , 3]), max(Parameters[, 
                      , 4]), max(Parameters[, , 5]))
                    dlag = (max_param - min_param)/100
                    range_param = numeric()
                    for (j in 1:5) {
                      range_param = cbind(range_param, seq(min_param[j] + 
                        dlag[j]/2, max_param[j] - dlag[j]/2, 
                        by = dlag[j]))
                    }
                    par(mfcol = c(5, 2))
                    longueur = 1:long
                    for (j in 1:5) {
                      for (i in 1:N) {
                        plot(x = longueur, y = Parameters[, i, 
                          j], type = "l", col = "grey", ylab = paste(param[j], 
                          " of study ", i), xlab = "iteration number", 
                          main = paste("Thinning interval = ", 
                            thin.interval, "\n Total samplesize kept = ", 
                            (iter.num * nb_chains - (burn_in) * 
                              nb_chains)/Thin), ylim = range(range_param[, 
                            j]))
                        abline(a = parameters[i, 2, j], b = 0, 
                          col = "green", lwd = 3)
                        abline(a = parameters[i, 3, j], b = 0, 
                          col = "green", lwd = 3)
                      }
                    }
                    min_Param = c(min(Parameter[, 1]), min(Parameter[, 
                      2]), min(Parameter[, 3]), min(Parameter[, 
                      4]), min(Parameter[, 5]), min(Parameter[, 
                      6]), min(Parameter[, 7]))
                    max_Param = c(max(Parameter[, 1]), max(Parameter[, 
                      2]), max(Parameter[, 3]), max(Parameter[, 
                      4]), max(Parameter[, 5]), max(Parameter[, 
                      6]), max(Parameter[, 7]))
                    dlag = (max_Param - min_Param)/100
                    range_Param = numeric()
                    for (j in 1:7) {
                      range_Param = cbind(range_Param, seq(min_Param[j] + 
                        dlag[j]/2, max_Param[j] - dlag[j]/2, 
                        by = dlag[j]))
                    }
                    for (j in 1:7) {
                      plot(x = longueur, y = Parameter[, j], 
                        type = "l", col = "grey", ylab = paste(Param[j]), 
                        xlab = "iteration number", main = paste("Thinning interval = ", 
                          thin.interval, "\n Total samplesize kept = ", 
                          (iter.num * nb_chains - (burn_in) * 
                            nb_chains)/Thin), ylim = range(range_Param[, 
                          j]))
                      abline(a = parameter[j, 2], b = 0, col = "green", 
                        lwd = 3)
                      abline(a = parameter[j, 3], b = 0, col = "green", 
                        lwd = 3)
                    }
                    dev.off()
                  }
                  file.pdf3 = paste("Density plots for N =", 
                    round((iter.num * nb_chains - (burn_in) * 
                      nb_chains)/Thin, 0), ".pdf")
                  pdf(file.pdf3, paper = "a4", height = 20)
                  param = c("theta", "alpha", "PI", "S1", "C1")
                  Param = c("Capital Theta", "Capital Lambda", 
                    "beta", "~sigma[alpha]", "~sigma[theta]", 
                    "S Overall", "C Overall")
                  par(mfcol = c(5, 2))
                  longueur = 1:long
                  for (j in 1:5) {
                    for (i in 1:N) {
                      plot(density(Parameters[, i, j]), lwd = 4, 
                        type = "l", col = "grey", main = paste(param[j], 
                          " of study ", i, " \n Thinning interval = ", 
                          thin.interval, "\n Total samplesize kept = ", 
                          (iter.num * nb_chains - (burn_in) * 
                            nb_chains)/Thin))
                    }
                  }
                  for (j in 1:7) {
                    plot(density(Parameter[, j]), lwd = 4, type = "l", 
                      col = "grey", main = paste(Param[j], " \n Thinning interval = ", 
                        thin.interval, "\n Total samplesize kept = ", 
                        (iter.num * nb_chains - (burn_in) * nb_chains)/Thin))
                  }
                  dev.off()
                  pdf("Summary ROC curve.pdf")
                  default.x = range(1, 0)
                  default.y = range(0, 1)
                  plot(x = default.x, y = default.y, type = "n", 
                    xlim = rev(range(default.x)), xlab = "", 
                    ylab = "")
                  title(xlab = "Specificity", ylab = "Sensitivity", 
                    cex.lab = 1.5, main = "Summary ROC curve")
                  Sensi1 = apply(as.matrix(Parameters[, , 4]), 
                    2, median)
                  Speci1 = apply(as.matrix(Parameters[, , 5]), 
                    2, median)
                  Scale_factor = 10
                  symbols(Speci1, Sensi1, circles = rowSums(as.matrix(data[[1]])), 
                    inches = 0.1 * Scale_factor/7, add = TRUE)
                  Ov_Se = 1 - pnorm((median(Parameter[, 1]) - 
                    median(Parameter[, 2])/2)/exp(median(Parameter[, 
                    3])/2))
                  Ov_Sp = pnorm((median(Parameter[, 1]) + median(Parameter[, 
                    2])/2)/exp(-median(Parameter[, 3])/2))
                  points(Ov_Sp, Ov_Se, pch = 19, cex = 2)
                  thet = qnorm((1 - as.matrix(Parameters[, , 
                    4])) + 1e-14) * exp(Parameter[, 3]/2) + Parameter[, 
                    2]/2
                  min_TH = quantile(thet, 0.05)
                  max_TH = quantile(thet, 0.95)
                  dTH = 5e-05
                  TH_range = seq(min_TH + dTH/2, max_TH - dTH/2, 
                    dTH)
                  S_sroc = 1 - pnorm((TH_range - median(Parameter[, 
                    2])/2)/exp(median(Parameter[, 3])/2))
                  C_sroc = pnorm((TH_range + median(Parameter[, 
                    2])/2)/exp(-median(Parameter[, 3])/2))
                  lines(C_sroc, S_sroc, lwd = 3, col = "black", 
                    lty = 1)
                  dev.off()
                }
                else {
                  if (SCO == TRUE) {
                    if (is.null(chain) == FALSE) {
                      file.pdf5 = paste("Trace plots for N =", 
                        round((iter.num * nb_chains - (burn_in) * 
                          nb_chains)/Thin, 0), ".pdf")
                      pdf(file.pdf5, paper = "a4", height = 20)
                      param = c("alpha", "PI", "S1", "C1")
                      Param = c("Capital Theta", "Capital Lambda", 
                        "beta", "~sigma[alpha]", "S Overall", 
                        "C Overall")
                      no_chains = length(chain)
                      iter_chain = round((iter.num * nb_chains - 
                        (burn_in) * nb_chains)/Thin, 0)/no_chains
                      min_param = c(min(Parameters[, , 1]), min(Parameters[, 
                        , 2]), min(Parameters[, , 3]), min(Parameters[, 
                        , 4]))
                      max_param = c(max(Parameters[, , 1]), max(Parameters[, 
                        , 2]), max(Parameters[, , 3]), max(Parameters[, 
                        , 4]))
                      dlag = (max_param - min_param)/100
                      range_param = numeric()
                      for (j in 1:4) {
                        range_param = cbind(range_param, seq(min_param[j] + 
                          dlag[j]/2, max_param[j] - dlag[j]/2, 
                          by = dlag[j]))
                      }
                      par(mfcol = c(5, 2))
                      longueur = 1:iter_chain
                      for (j in 1:4) {
                        for (i in 1:N) {
                          plot(x = longueur, y = Parameters[longueur, 
                            i, j], type = "n", col = 1, ylab = paste(param[j], 
                            " of study ", i), xlab = "iteration number", 
                            main = paste("Thinning interval = ", 
                              thin.interval, "\n Total samplesize kept = ", 
                              (iter.num * nb_chains - (burn_in) * 
                                nb_chains)/Thin), ylim = range(range_param[, 
                              j]))
                          for (l in 1:length(chain)) {
                            lines(x = longueur, y = Parameters[l * 
                              longueur, i, j], col = l)
                          }
                        }
                      }
                      min_Param = c(min(Parameter[, 1]), min(Parameter[, 
                        2]), min(Parameter[, 3]), min(Parameter[, 
                        4]), min(Parameter[, 5]), min(Parameter[, 
                        6]))
                      max_Param = c(max(Parameter[, 1]), max(Parameter[, 
                        2]), max(Parameter[, 3]), max(Parameter[, 
                        4]), max(Parameter[, 5]), max(Parameter[, 
                        6]))
                      dlag = (max_Param - min_Param)/100
                      range_Param = numeric()
                      for (j in 1:6) {
                        range_Param = cbind(range_Param, seq(min_Param[j] + 
                          dlag[j]/2, max_Param[j] - dlag[j]/2, 
                          by = dlag[j]))
                      }
                      for (j in 1:6) {
                        plot(x = longueur, y = Parameter[longueur, 
                          j], type = "n", col = 1, ylab = paste(Param[j]), 
                          xlab = "iteration number", main = paste("Thinning interval = ", 
                            thin.interval, "\n Total samplesize kept = ", 
                            (iter.num * nb_chains - (burn_in) * 
                              nb_chains)/Thin), ylim = range(range_Param[, 
                            j]))
                        for (l in 1:length(chain)) {
                          lines(x = longueur, y = Parameter[l * 
                            longueur, j], col = l)
                        }
                      }
                      dev.off()
                    }
                    else {
                      file.pdf2 = paste("Trace plots for N =", 
                        round((iter.num * nb_chains - (burn_in) * 
                          nb_chains)/Thin, 0), ".pdf")
                      pdf(file.pdf2, paper = "a4", height = 20)
                      param = c("alpha", "PI", "S1", "C1")
                      Param = c("Capital Theta", "Capital Lambda", 
                        "beta", "~sigma[alpha]", "S Overall", 
                        "C Overall")
                      min_param = c(min(Parameters[, , 1]), min(Parameters[, 
                        , 2]), min(Parameters[, , 3]), min(Parameters[, 
                        , 4]))
                      max_param = c(max(Parameters[, , 1]), max(Parameters[, 
                        , 2]), max(Parameters[, , 3]), max(Parameters[, 
                        , 4]))
                      dlag = (max_param - min_param)/100
                      range_param = numeric()
                      for (j in 1:4) {
                        range_param = cbind(range_param, seq(min_param[j] + 
                          dlag[j]/2, max_param[j] - dlag[j]/2, 
                          by = dlag[j]))
                      }
                      par(mfcol = c(5, 2))
                      longueur = 1:long
                      for (j in 1:4) {
                        for (i in 1:N) {
                          plot(x = longueur, y = Parameters[, 
                            i, j], type = "l", col = "grey", 
                            ylab = paste(param[j], " of study ", 
                              i), xlab = "iteration number", 
                            main = paste("Thinning interval = ", 
                              thin.interval, "\n Total samplesize kept = ", 
                              (iter.num * nb_chains - (burn_in) * 
                                nb_chains)/Thin), ylim = range(range_param[, 
                              j]))
                        }
                      }
                      min_Param = c(min(Parameter[, 1]), min(Parameter[, 
                        2]), min(Parameter[, 3]), min(Parameter[, 
                        4]), min(Parameter[, 5]), min(Parameter[, 
                        6]))
                      max_Param = c(max(Parameter[, 1]), max(Parameter[, 
                        2]), max(Parameter[, 3]), max(Parameter[, 
                        4]), max(Parameter[, 5]), max(Parameter[, 
                        6]))
                      dlag = (max_Param - min_Param)/100
                      range_Param = numeric()
                      for (j in 1:6) {
                        range_Param = cbind(range_Param, seq(min_Param[j] + 
                          dlag[j]/2, max_Param[j] - dlag[j]/2, 
                          by = dlag[j]))
                      }
                      for (j in 1:6) {
                        plot(x = longueur, y = Parameter[, j], 
                          type = "l", col = "grey", ylab = paste(Param[j]), 
                          xlab = "iteration number", main = paste("Thinning interval = ", 
                            thin.interval, "\n Total samplesize kept = ", 
                            (iter.num * nb_chains - (burn_in) * 
                              nb_chains)/Thin), ylim = range(range_Param[, 
                            j]))
                      }
                      dev.off()
                    }
                    file.pdf3 = paste("Density plots for N =", 
                      round((iter.num * nb_chains - (burn_in) * 
                        nb_chains)/Thin, 0), ".pdf")
                    pdf(file.pdf3, paper = "a4", height = 20)
                    param = c("alpha", "PI", "S1", "C1")
                    Param = c("Capital Theta", "Capital Lambda", 
                      "beta", "~sigma[alpha]", "S Overall", "C Overall")
                    par(mfcol = c(5, 2))
                    longueur = 1:long
                    for (j in 1:4) {
                      for (i in 1:N) {
                        plot(density(Parameters[, i, j]), lwd = 4, 
                          type = "l", col = "grey", main = paste(param[j], 
                            " of study ", i, " \n Thinning interval = ", 
                            thin.interval, "\n Total samplesize kept = ", 
                            (iter.num * nb_chains - (burn_in) * 
                              nb_chains)/Thin))
                      }
                    }
                    for (j in 1:6) {
                      plot(density(Parameter[, j]), lwd = 4, 
                        type = "l", col = "grey", main = paste(Param[j], 
                          " \n Thinning interval = ", thin.interval, 
                          "\n Total samplesize kept = ", (iter.num * 
                            nb_chains - (burn_in) * nb_chains)/Thin))
                    }
                    dev.off()
                    pdf("Summary ROC curve.pdf")
                    default.x = range(1, 0)
                    default.y = range(0, 1)
                    plot(x = default.x, y = default.y, type = "n", 
                      xlim = rev(range(default.x)), xlab = "", 
                      ylab = "")
                    title(xlab = "Specificity", ylab = "Sensitivity", 
                      cex.lab = 1.5, main = "Summary ROC curve")
                    Sensi1 = apply(as.matrix(Parameters[, , 3]), 
                      2, median)
                    Speci1 = apply(as.matrix(Parameters[, , 4]), 
                      2, median)
                    Scale_factor = 10
                    symbols(Speci1, Sensi1, circles = rowSums(as.matrix(data[[1]])), 
                      inches = 0.1 * Scale_factor/7, add = TRUE)
                    Ov_Se = 1 - pnorm((median(Parameter[, 1]) - 
                      median(Parameter[, 2])/2)/exp(median(Parameter[, 
                      3])/2))
                    Ov_Sp = pnorm((median(Parameter[, 1]) + median(Parameter[, 
                      2])/2)/exp(-median(Parameter[, 3])/2))
                    points(Ov_Sp, Ov_Se, pch = 19, cex = 2)
                    thet = qnorm((1 - as.matrix(Parameters[, 
                      , 3])) + 1e-14) * exp(Parameter[, 3]/2) + 
                      Parameter[, 2]/2
                    min_TH = quantile(thet, 0.05)
                    max_TH = quantile(thet, 0.95)
                    dTH = 5e-05
                    TH_range = seq(min_TH + dTH/2, max_TH - dTH/2, 
                      dTH)
                    S_sroc = 1 - pnorm((TH_range - median(Parameter[, 
                      2])/2)/exp(median(Parameter[, 3])/2))
                    C_sroc = pnorm((TH_range + median(Parameter[, 
                      2])/2)/exp(-median(Parameter[, 3])/2))
                    lines(C_sroc, S_sroc, lwd = 3, col = "black", 
                      lty = 1)
                    dev.off()
                  }
                }
            }
        }
    }
    if (Gold_Std == FALSE) {
        Results = list(parameter, parameters, refstd_parameters, 
            paste("See '", getwd(), "' for complete results", 
                sep = ""))
        names(Results) = c("Between-study parameters", "Within-study parameters", 
            "Reference standard", "")
    }
    else {
        Results = list(parameter, parameters, paste("See '", 
            getwd(), "' for complete results", sep = ""))
        names(Results) = c("Between-study parameters", "Within-study parameters", 
            "")
    }
    return(Results)
    cat(paste("See \"", getwd(), "\" for complete results", sep = ""))
}
