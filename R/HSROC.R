HSROC <-
function (data, iter.num, init = NULL, sub_rs = NULL, first.run = TRUE, 
    path = getwd(), refresh = 100, prior.SEref = NULL, prior.SPref = NULL, 
    prior_PI = c(0, 1), prior_LAMBDA = c(-3, 3), prior_THETA = c(-1.5, 
        1.5), prior_sd_alpha = list(0, 2, "sd"), prior_sd_theta = list(0, 
        2, "sd"), prior_beta = c(-0.75, 0.75)) 
{
    if (missing(data)) 
        stop("You must provide a valid 'data' argument", call. = FALSE)
    N = length(data[, 1])
    if (missing(iter.num) | iter.num <= 0) {
        cat("The number of iteration is either missing or less than 1. \n", 
            call. = FALSE)
        stop("Please respecify and call HSROC() again.\n")
    }
    if (is.null(sub_rs) == TRUE) {
        sub_rs = list(1, 1:N)
    }
    if (sub_rs[[1]] != (length(sub_rs) - 1)) {
        cat(paste("The value of the first element of 'sub_rs' (sub_rs[[1]] = ", 
            sub_rs[[1]], " ) does not match the number of remaining elements (length(sub_rs[[2:", 
            length(sub_rs), "]])) = ", length(2:length(sub_rs)), 
            "\n", sep = ""))
        stop("Please respecify and call HSROC() again.\n")
    }
    if (is.logical(first.run) == FALSE) {
        cat("The 'first.run' argument must be a logical object. \n")
        stop("Please respecify and call HSROC() again.\n")
    }
    if (is.null(prior.SEref) == FALSE | is.null(prior.SPref) == 
        FALSE) {
        if ((length(prior.SEref)/2 + length(prior.SEref)/2)/2 != 
            sub_rs[[1]]) {
            cat("The number of reference standards in 'prior.SEref' and(or) 'prior.SPref' is not matching the one defined in the 'sub_rs[[1]]' argument. \n")
            stop("Please respecify and call HSROC() again.\n")
        }
    }
    if (is.null(prior.SEref) == TRUE & is.null(prior.SPref) == 
        TRUE) {
        Gold_Std = TRUE
    }
    else {
        Gold_Std = FALSE
    }
    if (is.null(prior.SEref) == TRUE) {
        write(1, file = "S2.txt", ncolumns = 1)
    }
    else {
        write(2, file = "S2.txt", ncolumns = 1)
    }
    if (is.null(prior.SPref) == TRUE) {
        write(1, file = "C2.txt", ncolumns = 1)
    }
    else {
        write(2, file = "C2.txt", ncolumns = 1)
    }
    if (is.null(init) == FALSE) {
        random = FALSE
        if (sum(dim(init[[1]])) != N + 5) {
            cat(paste("Initial values for the within-study parameters were misspecified. Make sure the ordering described in the help file is preserved. \n"))
            stop("Please respecify and call HSROC() again.\n")
        }
        if (length(init[[2]]) != 5) {
            cat(paste("Initial values for the between-study parameters were misspecified. Make sure the ordering described in the help file is preserved. \n"))
            stop("Please respecify and call HSROC() again.\n")
        }
        if (Gold_Std == FALSE) {
            if (sum(dim(init[[3]])) != sub_rs[[1]] + 2) {
                cat(paste("Initial values for the test under evaluation were misspecified. Make sure the ordering described in the help file is preserved. \n"))
                stop("Please respecify and call HSROC() again.\n")
            }
        }
    }
    else {
        random = TRUE
    }
    low.pi = prior_PI[1]
    up.pi = prior_PI[2]
    if (all(low.pi < up.pi) == FALSE) {
        cat("The 'prior_PI' argument is a vector with 2 components specifying a range.  Thus, the first component of the vector must be less than the second component. Type '? HSROC' for more help. \n")
        stop("Please respecify and call HSROC() again.\n")
    }
    prior.LAMBDA.lower = prior_LAMBDA[1]
    prior.LAMBDA.upper = prior_LAMBDA[2]
    if (all(prior.LAMBDA.lower < prior.LAMBDA.upper) == FALSE) {
        cat("The 'prior_LAMBDA' argument is a vector with 2 components specifying a range.  Thus, the first component of the vector must be less than the second component. Type '? HSROC' for more help. \n")
        stop("Please respecify and call HSROC() again.\n")
    }
    prior.THETA.lower = prior_THETA[1]
    prior.THETA.upper = prior_THETA[2]
    if (all(prior.THETA.lower < prior.THETA.upper) == FALSE) {
        cat("The 'prior_THETA' argument is a vector with 2 components specifying a range.  Thus, the first component of the vector must be less than the second component. Type '? HSROC' for more help. \n")
        stop("Please respecify and call HSROC() again.\n")
    }
    l.disp.alpha = prior_sd_alpha[[1]]
    u.disp.alpha = prior_sd_alpha[[2]]
    if (all(l.disp.alpha < u.disp.alpha) == FALSE & prior_sd_alpha[[3]] != 
        "p") {
        cat("The 'prior_sd_alpha' argument is a list with the first 2 components specifying a range.  Thus, the first component of the list must be less than the second component. Type '? HSROC' for more help. \n")
        stop("Please respecify and call HSROC() again.\n")
    }
    l.disp.theta = prior_sd_theta[[1]]
    u.disp.theta = prior_sd_theta[[2]]
    if (l.disp.theta != 0 & u.disp.theta != 0) {
        if (all(l.disp.theta < u.disp.theta) == FALSE & prior_sd_theta[[3]] != 
            "p") {
            cat("The 'prior_sd_theta' argument is a list with the first 2 components specifying a range.  Thus, the first component of the list must be less than the second component. Type '? HSROC' for more help. \n")
            stop("Please respecify and call HSROC() again.\n")
        }
    }
    if (prior_sd_theta[[1]] == 0 & prior_sd_theta[[2]] == 0) {
        SCO = TRUE
    }
    else {
        SCO = FALSE
    }
    if (is.null(prior_beta)) {
        beta.a = -log((prior.LAMBDA.upper/3) + 1)
        beta.b = log((prior.LAMBDA.upper/3) + 1)
    }
    else {
        beta.a = prior_beta[1]
        beta.b = prior_beta[2]
        if (all(beta.a < beta.b) == FALSE) {
            cat("The 'prior_beta' argument is a vector with 2 components specifying a range.  Thus, the first component of the vector must be less than the second component. Type '? HSROC' for more help. \n")
            stop("Please respecify and call HSROC() again.\n")
        }
    }
    write(1, file = "model.txt", ncolumns = 1)
    data = list(data)
    file.pi = "PI.txt"
    file.C2 = "Spec2.txt"
    file.S2 = "Sens2.txt"
    file.alpha = "alpha.txt"
    file.theta = "theta.txt"
    file.sig.theta = "sigma.theta.txt"
    file.sig.alpha = "sigma.alpha.txt"
    file.THETA = "capital.THETA.txt"
    file.LAMBDA = "LAMBDA.txt"
    file.beta = "beta.txt"
    file.C1 = "Spec1.txt"
    file.S1 = "Sens1.txt"
    file.C_overall = "C_overall.txt"
    file.S_overall = "S_overall.txt"
    file.choix = "choix.txt"
    file.ll = "log.likelihood.txt"
    file.Yj = "Y_j.txt"
    file.TV = "Start_values.txt"
    file.TV2 = "Start_values2.txt"
    file.TV3 = "Start_REFSTD.txt"
    file.Restart = "REstarting values.txt"
    file.Restart2 = "REstarting values 2.txt"
    file.Restart_REFSTD = "REstarting REFSTD.txt"
    file.Restart_index = "REstart_values_index.txt"
    file.A.alpha = "A.alpha.txt"
    file.B.alpha = "B.alpha.txt"
    file.mean.rij.one = "mean.rij.one.txt"
    file.mean.rij.zero = "mean.rij.zero.txt"
    setwd(path)
    condInd = TRUE
    prior_dist_PI = "beta"
    range_rij = c(-prior.LAMBDA.upper/2 - 3 * exp(-beta.a/2), 
        prior.LAMBDA.upper/2 + 3 * exp(beta.b/2))
    low.rj = range_rij[1]
    up.rj = range_rij[2]
    write(c(low.rj, up.rj), file = "range of latent variable.txt", 
        ncolumns = 2)
    alpha.PI = beta.parameter(low = low.pi, up = up.pi)[1, ]
    beta.PI = beta.parameter(low = low.pi, up = up.pi)[2, ]
    L.disp.alpha = L.disp.theta = numeric()
    if (l.disp.alpha == 0) {
        L.disp.alpha = 1e-10
    }
    else {
        if (l.disp.alpha > 0) {
            L.disp.alpha = l.disp.alpha
        }
    }
    if (prior_sd_alpha[[3]] == "sd") {
        low.disp.alpha = u.disp.alpha^(-2)
        up.disp.alpha = L.disp.alpha^(-2)
        write(1, file = "Prior on sigma_alpha.txt", ncolumns = 1)
    }
    else {
        if (prior_sd_alpha[[3]] == "v") {
            low.disp.alpha = u.disp.alpha^(-1)
            up.disp.alpha = L.disp.alpha^(-1)
            write(2, file = "Prior on sigma_alpha.txt", ncolumns = 1)
        }
        else {
            if (prior_sd_alpha[[3]] == "p") {
                low.disp.alpha = L.disp.alpha
                up.disp.alpha = u.disp.alpha
                write(3, file = "Prior on sigma_alpha.txt", ncolumns = 1)
            }
        }
    }
    if (SCO == FALSE) {
        if (l.disp.theta == 0) {
            L.disp.theta = 1e-10
        }
        else {
            if (l.disp.theta > 0) {
                L.disp.theta = l.disp.theta
            }
        }
        if (prior_sd_theta[[3]] == "sd") {
            low.disp.theta = u.disp.theta^(-2)
            up.disp.theta = L.disp.theta^(-2)
            write(1, file = "Prior on sigma_theta.txt", ncolumns = 1)
        }
        else {
            if (prior_sd_theta[[3]] == "v") {
                low.disp.theta = u.disp.theta^(-1)
                up.disp.theta = L.disp.theta^(-1)
                write(2, file = "Prior on sigma_theta.txt", ncolumns = 1)
            }
            else {
                if (prior_sd_theta[[3]] == "p") {
                  low.disp.theta = L.disp.theta
                  up.disp.theta = u.disp.theta
                  write(3, file = "Prior on sigma_theta.txt", 
                    ncolumns = 1)
                }
            }
        }
    }
    else {
        if (SCO == TRUE) {
            low.disp.theta = up.disp.theta = 1
        }
    }
    long.se = length(prior.SEref)
    low.se = prior.SEref[1:sub_rs[[1]]]
    up.se = prior.SEref[(sub_rs[[1]] + 1):long.se]
    long.sp = length(prior.SPref)
    low.sp = prior.SPref[1:sub_rs[[1]]]
    up.sp = prior.SPref[(sub_rs[[1]] + 1):long.sp]
    if (Gold_Std == FALSE) {
        if (is.null(prior.SEref) == TRUE) {
            Sens2.alpha = Sens2.beta = NULL
            Gold_se = TRUE
        }
        else {
            if (is.null(prior.SEref) == FALSE) {
                Sens2.alpha = beta.parameter(low = low.se, up = up.se)[1, 
                  ]
                Sens2.beta = beta.parameter(low = low.se, up = up.se)[2, 
                  ]
                Gold_se = FALSE
            }
        }
    }
    else {
        if (Gold_Std == TRUE) {
            Sens2.alpha = Sens2.beta = Gold_se = NULL
        }
    }
    if (Gold_Std == FALSE) {
        if (is.null(prior.SPref) == TRUE) {
            Spec2.alpha = Spec2.beta = NULL
            Gold_sp = TRUE
        }
        else {
            if (is.null(prior.SPref) == FALSE) {
                Spec2.alpha = beta.parameter(low = low.sp, up = up.sp)[1, 
                  ]
                Spec2.beta = beta.parameter(low = low.sp, up = up.sp)[2, 
                  ]
                Gold_sp = FALSE
            }
        }
    }
    else {
        if (Gold_Std == TRUE) {
            Spec2.alpha = Spec2.beta = Gold_sp = NULL
        }
    }
    if (first.run == TRUE) {
        RESTART_i = NA
        RESTART = NA
        RESTART_REFSTD = NA
        file.create(file.Restart)
        file.create(file.Restart2)
        file.create(file.Restart_REFSTD)
    }
    else {
        DATA.restart = read.table(file.Restart)
        DATA.restart2 = read.table(file.Restart2)
        DATA.restart_refstd = read.table(file.Restart_REFSTD)
        RESTART_i = t(DATA.restart)
        RESTART = DATA.restart2
        RESTART_REFSTD = DATA.restart_refstd
    }
    PRIOR.Parameters = c(beta.a, beta.b, prior.THETA.lower, prior.THETA.upper, 
        prior.LAMBDA.lower, prior.LAMBDA.upper, low.disp.alpha, 
        up.disp.alpha, low.disp.theta, up.disp.theta, alpha.PI, 
        beta.PI, Sens2.alpha, Sens2.beta, Spec2.alpha, Spec2.beta)
    test.results = data[[1]]
    Start.values = Which_data(RANDOM = random, data = data, init = init, 
        GS = Gold_Std)[[1]]
    Start.values2 = Which_data(RANDOM = random, data = data, 
        init = init, GS = Gold_Std)[[2]]
    Start.REFSTD = Which_data(RANDOM = random, data = data, init = init, 
        GS = Gold_Std)[[3]]
    INITS = Initialization(first.run = first.run, random = random, 
        param = PRIOR.Parameters, cond.Ind = condInd, rs = sub_rs, 
        GS_se = Gold_se, GS_sp = Gold_sp, Data1 = Start.values, 
        Data2 = RESTART_i, Data3 = RESTART, Data4 = Start.values2, 
        Data5 = Start.REFSTD, Data6 = RESTART_REFSTD, path = path, 
        studies = N, sco = SCO, psa = prior_sd_alpha[[3]], pst = prior_sd_theta[[3]])
    if (INITS[[1]][3] == 0 | INITS[[1]][1] == 0) {
        cat(paste("Unsuitable initial values were provided. "))
        stop("Please respecify and call HSROC() again.\n  If you're using 'init=NULL' you need just to run the 'HSROC' function again.\n")
    }
    init.sigma.alpha = INITS[[1]][3]
    prec.alpha = INITS[[1]][4]
    init.THETA = INITS[[1]][5]
    init.LAMBDA = INITS[[1]][6]
    init.beta = INITS[[1]][7]
    init.alpha = as.vector(INITS[[2]][, 1])
    init.S1 = as.vector(INITS[[2]][, 3])
    init.C1 = as.vector(INITS[[2]][, 4])
    init.PI = as.vector(INITS[[2]][, 5])
    if (SCO == TRUE) {
        init.theta = rep(init.THETA, N)
    }
    else {
        if (SCO == FALSE) {
            init.sigma.theta = INITS[[1]][1]
            prec.theta = INITS[[1]][2]
            init.theta = as.vector(INITS[[2]][, 2])
        }
    }
    if (Gold_Std == FALSE) {
        if (Gold_se == TRUE) {
            init.C2 = as.vector(INITS[[3]][2, ])
        }
        else {
            if (Gold_sp == TRUE) {
                init.S2 = as.vector(INITS[[3]][1, ])
            }
            else {
                init.S2 = as.vector(INITS[[3]][1, ])
                init.C2 = as.vector(INITS[[3]][2, ])
            }
        }
    }
    D = DATA.organizer(d = test.results, m = N)
    n = D[[1]]
    All.Studies = D[[2]]
    t1 = numeric()
    t2 = numeric()
    T = t(mapply(Test, All.Studies))
    t1 = T[, 1]
    t2 = T[, 2]
    studygroup = rep((1:N), n)
    n_rs = numeric()
    n_REFSTD = REFSTD_3(rs = sub_rs, n.sample = D[[1]], studies = N)
    studygroup_REFSTD = REFSTD_4(rs = sub_rs, n.sample = D[[1]], 
        n_rs = n_REFSTD)
    Total = sum(n)
    PRIOR.BETWEEN = rbind(c(beta.a, beta.b), c(prior.THETA.lower, 
        prior.THETA.upper), c(prior.LAMBDA.lower, prior.LAMBDA.upper), 
        c(l.disp.alpha, u.disp.alpha), c(l.disp.theta, u.disp.theta), 
        c(low.pi, up.pi), c(low.rj, up.rj))
    colnames(PRIOR.BETWEEN) = c("Lower bound", "Upper bound")
    rownames(PRIOR.BETWEEN) = c("beta", "THETA", "LAMBDA", "sigma_alpha", 
        "sigma_theta", "prevalence", "Range_rij")
    write.table(PRIOR.BETWEEN, file = "Prior.information.txt")
    if (Gold_Std == FALSE) {
        if (Gold_se == TRUE) {
            C2.p = c()
            for (i in 1:length(n_REFSTD)) {
                C2.p = c(C2.p, paste("C2", i, sep = ""))
            }
            PRIOR.C2 = cbind(low.sp, up.sp)
            rownames(PRIOR.C2) = C2.p
            colnames(PRIOR.C2) = c("lower", "upper")
            write.table(PRIOR.C2, file = "Prior.information.txt", 
                append = TRUE, col.names = FALSE)
        }
        else {
            if (Gold_sp == TRUE) {
                S2.p = c()
                for (i in 1:length(n_REFSTD)) {
                  S2.p = c(S2.p, paste("S2", i, sep = ""))
                }
                PRIOR.S2 = cbind(low.se, up.se)
                rownames(PRIOR.S2) = S2.p
                colnames(PRIOR.S2) = c("lower", "upper")
                write.table(PRIOR.S2, file = "Prior.information.txt", 
                  append = TRUE, col.names = FALSE)
            }
            else {
                S2.p = C2.p = c()
                for (i in 1:length(n_REFSTD)) {
                  S2.p = c(S2.p, paste("S2", i, sep = ""))
                  C2.p = c(C2.p, paste("C2", i, sep = ""))
                }
                PRIOR.S2 = cbind(low.se, up.se)
                rownames(PRIOR.S2) = S2.p
                colnames(PRIOR.S2) = c("lower", "upper")
                PRIOR.C2 = cbind(low.sp, up.sp)
                rownames(PRIOR.C2) = C2.p
                colnames(PRIOR.C2) = c("lower", "upper")
                write.table(PRIOR.S2, file = "Prior.information.txt", 
                  append = TRUE, col.names = FALSE)
                write.table(PRIOR.C2, file = "Prior.information.txt", 
                  append = TRUE, col.names = FALSE)
            }
        }
    }
    vec.PI = as.numeric(init.PI)
    vec.S1 = as.numeric(init.S1)
    vec.C1 = as.numeric(init.C1)
    vec.alpha = as.numeric(init.alpha)
    vec.sigma.alpha = as.numeric(init.sigma.alpha)
    vec.THETA = as.numeric(init.THETA)
    vec.LAMBDA = as.numeric(init.LAMBDA)
    vec.beta = as.numeric(init.beta)
    vec.MH = as.numeric(exp(vec.beta))
    if (SCO == TRUE) {
        vec.theta = rep(vec.THETA, N)
    }
    else {
        if (SCO == FALSE) {
            vec.sigma.theta = as.numeric(init.sigma.theta)
            vec.theta = as.numeric(init.theta)
        }
    }
    if (Gold_Std == TRUE) {
        vec.S2 = vec.C2 = 1
    }
    else {
        if (Gold_se == TRUE) {
            vec.C2 = as.numeric(init.C2)
            vec.S2 = 1
        }
        else {
            if (Gold_sp == TRUE) {
                vec.C2 = 1
                vec.S2 = as.numeric(init.S2)
            }
            else {
                vec.C2 = as.numeric(init.C2)
                vec.S2 = as.numeric(init.S2)
            }
        }
    }
    count = 1
    cat("Starting the Gibbs sampler for ", iter.num, " iterations ... \n")
    for (l in 1:iter.num) {
        if (l == refresh * count) {
            cat(paste(refresh * count, " iterations completed out of ", 
                iter.num, "... \n"))
            count = count + 1
        }
        PI.rep = rep(vec.PI, times = n)
        S1.rep = rep(vec.S1, times = n)
        C1.rep = rep(vec.C1, times = n)
        S2.rep = rep(vec.S2, REFSTD_1(rs = sub_rs, n = n_REFSTD, 
            t = Total))
        C2.rep = rep(vec.C2, REFSTD_1(rs = sub_rs, n = n_REFSTD, 
            t = Total))
        prob.Yj = Prob(t1, t2, pi = PI.rep, S1 = S1.rep, S2 = S2.rep, 
            C1 = C1.rep, C2 = C2.rep)
        Y.ij = rbinom(n = sum(n), size = 1, prob = prob.Yj)
        y = cbind(studygroup, t1, t2, Y.ij)
        Y1 = by(y[, 2:4], y[, 1], FUN = XY.function, b = 1)
        Y2 = by(y[, 2:4], y[, 1], FUN = XY.function, b = 2)
        Y3 = by(y[, 2:4], y[, 1], FUN = XY.function, b = 3)
        Y4 = by(y[, 2:4], y[, 1], FUN = XY.function, b = 4)
        PI.alpha = by(y[, 4], y[, 1], pi.alpha, b = alpha.PI)
        PI.beta = by(y[, 4], y[, 1], pi.beta, b = beta.PI)
        PI = rbeta(n = N, shape1 = PI.alpha, shape2 = PI.beta)
        vec.PI = PI
        if (Gold_Std == FALSE) {
            YY1 = REFSTD_5(rs = sub_rs, yy1 = Y1, yy2 = Y2, yy3 = Y3, 
                yy4 = Y4, stu.gr = studygroup_REFSTD, t1 = t1, 
                t2 = t2, yij = Y.ij)[[1]]
            YY2 = REFSTD_5(rs = sub_rs, yy1 = Y1, yy2 = Y2, yy3 = Y3, 
                yy4 = Y4, stu.gr = studygroup_REFSTD, t1 = t1, 
                t2 = t2, yij = Y.ij)[[2]]
            YY3 = REFSTD_5(rs = sub_rs, yy1 = Y1, yy2 = Y2, yy3 = Y3, 
                yy4 = Y4, stu.gr = studygroup_REFSTD, t1 = t1, 
                t2 = t2, yij = Y.ij)[[3]]
            YY4 = REFSTD_5(rs = sub_rs, yy1 = Y1, yy2 = Y2, yy3 = Y3, 
                yy4 = Y4, stu.gr = studygroup_REFSTD, t1 = t1, 
                t2 = t2, yij = Y.ij)[[4]]
            X.ij = 1 - Y.ij
            x = cbind(studygroup_REFSTD, t1, t2, X.ij)
            X1 = by(x[, 2:4], x[, 1], FUN = XY.function, b = 1)
            X2 = by(x[, 2:4], x[, 1], FUN = XY.function, b = 2)
            X3 = by(x[, 2:4], x[, 1], FUN = XY.function, b = 3)
            X4 = by(x[, 2:4], x[, 1], FUN = XY.function, b = 4)
            n.refstd = REFSTD_2(n_rs = sub_rs, likelihood = list(YY1, 
                YY2, YY3, YY4, X1, X2, X3, X4), prior = list(Sens2.alpha, 
                Sens2.beta, Spec2.alpha, Spec2.beta))[[1]]
            if (Gold_se == TRUE) {
                a.Sp2 = REFSTD_2(n_rs = sub_rs, likelihood = list(YY1, 
                  YY2, YY3, YY4, X1, X2, X3, X4), prior = list(Sens2.alpha, 
                  Sens2.beta, Spec2.alpha, Spec2.beta))[[4]]
                b.Sp2 = REFSTD_2(n_rs = sub_rs, likelihood = list(YY1, 
                  YY2, YY3, YY4, X1, X2, X3, X4), prior = list(Sens2.alpha, 
                  Sens2.beta, Spec2.alpha, Spec2.beta))[[5]]
                Spec2 = REFSTD_6_SP(refstd = Gold_Std, N.refstd = n.refstd, 
                  A.Sp2 = a.Sp2, B.Sp2 = b.Sp2)$SP
                vec.C2 = Spec2
                vec.S2 = 1
            }
            else {
                if (Gold_sp == TRUE) {
                  a.Se2 = REFSTD_2(n_rs = sub_rs, likelihood = list(YY1, 
                    YY2, YY3, YY4, X1, X2, X3, X4), prior = list(Sens2.alpha, 
                    Sens2.beta, Spec2.alpha, Spec2.beta))[[2]]
                  b.Se2 = REFSTD_2(n_rs = sub_rs, likelihood = list(YY1, 
                    YY2, YY3, YY4, X1, X2, X3, X4), prior = list(Sens2.alpha, 
                    Sens2.beta, Spec2.alpha, Spec2.beta))[[3]]
                  Sens2 = REFSTD_6_SE(refstd = Gold_Std, N.refstd = n.refstd, 
                    A.Se2 = a.Se2, B.Se2 = b.Se2)$SE
                  vec.C2 = 1
                  vec.S2 = Sens2
                }
                else {
                  a.Se2 = REFSTD_2(n_rs = sub_rs, likelihood = list(YY1, 
                    YY2, YY3, YY4, X1, X2, X3, X4), prior = list(Sens2.alpha, 
                    Sens2.beta, Spec2.alpha, Spec2.beta))[[2]]
                  b.Se2 = REFSTD_2(n_rs = sub_rs, likelihood = list(YY1, 
                    YY2, YY3, YY4, X1, X2, X3, X4), prior = list(Sens2.alpha, 
                    Sens2.beta, Spec2.alpha, Spec2.beta))[[3]]
                  a.Sp2 = REFSTD_2(n_rs = sub_rs, likelihood = list(YY1, 
                    YY2, YY3, YY4, X1, X2, X3, X4), prior = list(Sens2.alpha, 
                    Sens2.beta, Spec2.alpha, Spec2.beta))[[4]]
                  b.Sp2 = REFSTD_2(n_rs = sub_rs, likelihood = list(YY1, 
                    YY2, YY3, YY4, X1, X2, X3, X4), prior = list(Sens2.alpha, 
                    Sens2.beta, Spec2.alpha, Spec2.beta))[[5]]
                  Sens2 = REFSTD_6_SE(refstd = Gold_Std, N.refstd = n.refstd, 
                    A.Se2 = a.Se2, B.Se2 = b.Se2)$SE
                  Spec2 = REFSTD_6_SP(refstd = Gold_Std, N.refstd = n.refstd, 
                    A.Sp2 = a.Sp2, B.Sp2 = b.Sp2)$SP
                  vec.C2 = Spec2
                  vec.S2 = Sens2
                }
            }
        }
        else {
            if (Gold_Std == TRUE) {
                vec.S2 = 1
                vec.C2 = 1
            }
        }
        alpha.rep = rep(vec.alpha, times = n)
        theta.rep = rep(vec.theta, times = n)
        mu = alpha.rep * (Y.ij/2) - alpha.rep * ((1 - Y.ij)/2)
        s = exp(vec.beta * (Y.ij - 0.5))
        rr.j = apply(as.matrix(t1), 2, truncnorm, n = Total, 
            m = mu, sd = s, limit = cbind(theta.rep, low.rj, 
                up.rj))
        rrr.j = matrix(rr.j, ncol = 3, byrow = FALSE)
        r.ij = mapply(f.rij, rrr.j[, 2], rrr.j[, 3], rrr.j[, 
            1])
        r.ij.one = sum(r.ij * Y.ij)
        r.ij.zero = sum(r.ij * (1 - Y.ij))
        y = cbind(y, r.ij)
        if (SCO == FALSE) {
            by.theta = cbind(studygroup, t1, r.ij)
            lower1 = by(by.theta[, 3], INDICES = list(by.theta[, 
                1], by.theta[, 2]), FUN = max)[1:N]
            upper1 = by(by.theta[, 3], INDICES = list(by.theta[, 
                1], by.theta[, 2]), FUN = min)[(N + 1):(2 * N)]
            lower = apply(as.matrix(lower1), 1, limit.theta, 
                low.rj)
            upper = apply(as.matrix(upper1), 1, limit.theta, 
                up.rj)
            theta = mapply(truncnorm2, lower, upper, MoreArgs = list(vec.THETA, 
                vec.sigma.theta, 1))
            vec.theta = theta
        }
        else {
            if (SCO == TRUE) {
                vec.theta = rep(vec.THETA, N)
            }
        }
        by.alpha = cbind(studygroup, Y.ij, (1 - Y.ij), r.ij)
        A = 0.25 * by(by.alpha[, 2:3], by.alpha[, 1], A.fonction, 
            b = 2 * vec.beta) + 1/vec.sigma.alpha^2
        B = by(by.alpha[, c(2, 4)], by.alpha[, 1], B.fonction.Ind, 
            b = 2 * vec.beta) + (2 * vec.LAMBDA)/vec.sigma.alpha^2
        alpha = rnorm(n = N, mean = as.matrix(B/(2 * A)), sd = as.matrix((A^(-0.5))))
        vec.alpha = alpha
        LAMBDA = truncnorm2(n = 1, m = mean(vec.alpha), sd = (vec.sigma.alpha/sqrt(N)), 
            l = prior.LAMBDA.lower, u = prior.LAMBDA.upper)
        vec.LAMBDA = LAMBDA
        f.beta = M_H2_IND(r = r.ij, y = Y.ij, a = alpha.rep, 
            N = Total, low.beta = exp(beta.a), up.beta = exp(beta.b), 
            x = vec.MH)
        vec.beta = log(f.beta[1])
        vec.MH = f.beta[1]
        if (SCO == FALSE) {
            THETA = truncnorm2(n = 1, m = mean(vec.theta), sd = (vec.sigma.theta/sqrt(N)), 
                l = prior.THETA.lower, u = prior.THETA.upper)
            vec.THETA = THETA
        }
        else {
            if (SCO == TRUE) {
                by.THETA = cbind(t1, r.ij)
                lower1 = by(by.THETA[, 2], INDICES = list(by.THETA[, 
                  1]), FUN = max)[1]
                upper1 = by(by.THETA[, 2], INDICES = list(by.THETA[, 
                  1]), FUN = min)[2]
                lower = apply(as.matrix(lower1), 1, limit.theta, 
                  low.rj)
                upper = apply(as.matrix(upper1), 1, limit.theta, 
                  up.rj)
                THETA_bounds = truncunif(bornes = c(lower, upper), 
                  prior.l = prior.THETA.lower, prior.u = prior.THETA.upper)
                THETA = runif(1, THETA_bounds[1], THETA_bounds[2])
                vec.THETA = THETA
            }
        }
        if (prior_sd_alpha[[3]] == "sd") {
            prec.alpha.shape = (N/2) - (1/2)
            prec.alpha.scale = (0.5 * sum((vec.alpha - vec.LAMBDA)^2))^(-1)
            prec.alpha = truncgamma(n = 1, shape = prec.alpha.shape, 
                scale = prec.alpha.scale, l = low.disp.alpha, 
                u = up.disp.alpha)
        }
        else {
            if (prior_sd_alpha[[3]] == "v") {
                prec.alpha.shape = (N/2) - 1
                prec.alpha.scale = (0.5 * sum((vec.alpha - vec.LAMBDA)^2))^(-1)
                prec.alpha = truncgamma(n = 1, shape = prec.alpha.shape, 
                  scale = prec.alpha.scale, l = low.disp.alpha, 
                  u = up.disp.alpha)
            }
            else {
                if (prior_sd_alpha[[3]] == "p") {
                  prec.alpha.shape = (N/2) + low.disp.alpha
                  prec.alpha.scale = (up.disp.alpha + 0.5 * sum((vec.alpha - 
                    vec.LAMBDA)^2))^(-1)
                  prec.alpha = rgamma(n = 1, shape = prec.alpha.shape, 
                    scale = prec.alpha.scale)
                }
            }
        }
        sigma.alpha = sqrt(1/prec.alpha)
        vec.sigma.alpha = sigma.alpha
        if (SCO == FALSE) {
            if (prior_sd_theta[[3]] == "sd") {
                prec.theta.shape = (N/2) - (1/2)
                prec.theta.scale = (0.5 * sum((vec.theta - vec.THETA)^2))^(-1)
                prec.theta = truncgamma(n = 1, shape = prec.theta.shape, 
                  scale = prec.theta.scale, l = low.disp.theta, 
                  u = up.disp.theta)
            }
            else {
                if (prior_sd_theta[[3]] == "v") {
                  prec.theta.shape = (N/2) - (1/2)
                  prec.theta.scale = (0.5 * sum((vec.theta - 
                    vec.THETA)^2))^(-1)
                  prec.theta = truncgamma(n = 1, shape = prec.theta.shape, 
                    scale = prec.theta.scale, l = low.disp.theta, 
                    u = up.disp.theta)
                }
                else {
                  if (prior_sd_alpha[[3]] == "p") {
                    prec.theta.shape = (N/2) + low.disp.theta
                    prec.theta.scale = (up.disp.theta + 0.5 * 
                      sum((vec.theta - vec.THETA)^2))^(-1)
                    prec.theta = rgamma(n = 1, shape = prec.theta.shape, 
                      scale = prec.theta.scale)
                  }
                }
            }
            sigma.theta = sqrt(1/prec.theta)
            vec.sigma.theta = sigma.theta
        }
        Sens1 = pnorm(exp(-vec.beta/2) * (vec.theta - vec.alpha/2), 
            mean = 0, sd = 1, lower.tail = FALSE)
        Spec1 = pnorm(exp(vec.beta/2) * (vec.theta + vec.alpha/2), 
            mean = 0, sd = 1, lower.tail = TRUE)
        S_overall = pnorm((vec.THETA - vec.LAMBDA/2) * exp(-vec.beta/2), 
            mean = 0, sd = 1, lower.tail = FALSE)
        C_overall = pnorm((vec.THETA + vec.LAMBDA/2) * exp(vec.beta/2), 
            mean = 0, sd = 1, lower.tail = TRUE)
        vec.S1 = Sens1
        vec.C1 = Spec1
        beta_new = runif(1, beta.a, beta.b)
        theta_new = rnorm(1, mean = vec.THETA, sd = vec.sigma.theta)
        alpha_new = rnorm(1, mean = vec.LAMBDA, sd = vec.sigma.alpha)
        Sens1_new = pnorm(exp(-beta_new/2) * (theta_new - alpha_new/2), 
            mean = 0, sd = 1, lower.tail = FALSE)
        Spec1_new = pnorm(exp(beta_new/2) * (theta_new + alpha_new/2), 
            mean = 0, sd = 1, lower.tail = TRUE)
        if (SCO == FALSE) {
            write(cbind(vec.alpha, vec.theta, vec.S1, vec.C1, 
                vec.PI), file = file.Restart, ncolumns = N)
            write(c(vec.LAMBDA, vec.sigma.alpha, vec.THETA, vec.sigma.theta, 
                vec.beta), file = file.Restart2)
            write(t(rbind(vec.S2, vec.C2)), file = file.Restart_REFSTD, 
                ncolumns = length(n_REFSTD))
        }
        else {
            if (SCO == TRUE) {
                write(cbind(vec.alpha, vec.S1, vec.C1, vec.PI), 
                  file = file.Restart, ncolumns = N)
                write(c(vec.LAMBDA, vec.sigma.alpha, vec.THETA, 
                  vec.beta), file = file.Restart2)
                write(t(rbind(vec.S2, vec.C2)), file = file.Restart_REFSTD, 
                  ncolumns = length(n_REFSTD))
            }
        }
        write(Sens1_new, file = "Sens1_new.txt", append = TRUE, 
            ncolumns = N)
        write(Spec1_new, file = "Spec1_new.txt", append = TRUE, 
            ncolumns = N)
        write(vec.PI, file = file.pi, append = TRUE, ncolumns = N)
        write(vec.alpha, file = file.alpha, append = TRUE, ncolumns = N)
        write(vec.sigma.alpha, file = file.sig.alpha, append = TRUE, 
            ncolumns = N)
        write(vec.THETA, file = file.THETA, append = TRUE, ncolumns = N)
        write(vec.LAMBDA, file = file.LAMBDA, append = TRUE, 
            ncolumns = N)
        write(vec.beta, file = file.beta, append = TRUE, ncolumns = N)
        write(vec.S1, file = file.S1, append = TRUE, ncolumns = N)
        write(vec.C1, file = file.C1, append = TRUE, ncolumns = N)
        write(f.beta, file = file.choix, append = TRUE, ncolumns = 11)
        write(C_overall, file = file.C_overall, append = TRUE, 
            ncolumns = N)
        write(S_overall, file = file.S_overall, append = TRUE, 
            ncolumns = N)
        if (SCO == FALSE) {
            write(vec.theta, file = file.theta, append = TRUE, 
                ncolumns = N)
            write(vec.sigma.theta, file = file.sig.theta, append = TRUE, 
                ncolumns = N)
        }
        if (Gold_Std == FALSE) {
            if (Gold_se == TRUE) {
                write(vec.C2, file = file.C2, append = TRUE, 
                  ncolumns = length(n_REFSTD))
            }
            else {
                if (Gold_sp == TRUE) {
                  write(vec.S2, file = file.S2, append = TRUE, 
                    ncolumns = length(n_REFSTD))
                }
                else {
                  write(vec.C2, file = file.C2, append = TRUE, 
                    ncolumns = length(n_REFSTD))
                  write(vec.S2, file = file.S2, append = TRUE, 
                    ncolumns = length(n_REFSTD))
                }
            }
        }
    }
    file.remove(file.Restart2)
    file.remove(file.Restart)
    file.remove(file.Restart_REFSTD)
    if (SCO == FALSE) {
        write(cbind(vec.alpha, vec.theta, vec.S1, vec.C1, vec.PI), 
            file = file.Restart, append = TRUE, ncolumns = N)
        write(c(vec.LAMBDA, vec.sigma.alpha, vec.THETA, vec.sigma.theta, 
            vec.beta), file = file.Restart2, append = TRUE)
        write(t(rbind(vec.S2, vec.C2)), file = file.Restart_REFSTD, 
            append = TRUE, ncolumns = length(n_REFSTD))
        write(paste("______________________________________________________"), 
            file = file.Restart_index, append = TRUE)
        write(paste("\t   REstarting values.txt "), file = file.Restart_index, 
            append = TRUE)
        write(paste("Row 1 : alpha parameters for all M = ", 
            N, " study(ies)\t    "), file = file.Restart_index, 
            append = TRUE)
        write(paste("Row 2 : theta parameters for all M = ", 
            N, " study(ies)\t    "), file = file.Restart_index, 
            append = TRUE)
        write(paste("Row 3 : sensitivity of test under evaluation (S1) for all M = ", 
            N, " study(ies)\t    "), file = file.Restart_index, 
            append = TRUE)
        write(paste("Row 4 : specificity of test under evaluation (C1) for all M = ", 
            N, " study(ies)\t    "), file = file.Restart_index, 
            append = TRUE)
        write(paste("Row 5 : prevalence for all M = ", N, " study(ies)\t    "), 
            file = file.Restart_index, append = TRUE)
        write(paste("______________________________________________________"), 
            file = file.Restart_index, append = TRUE)
        write(paste("\t   REstarting values 2.txt "), file = file.Restart_index, 
            append = TRUE)
        write(paste("Column 1 : LAMBDA parameter\t    "), file = file.Restart_index, 
            append = TRUE)
        write(paste("Column 2 : sigma alpha parameter\t    "), 
            file = file.Restart_index, append = TRUE)
        write(paste("Column 3 : THETA parameter\t    "), file = file.Restart_index, 
            append = TRUE)
        write(paste("Column 4 : sigma theta parameter\t    "), 
            file = file.Restart_index, append = TRUE)
        write(paste("Column 5 : beta parameter\t    "), file = file.Restart_index, 
            append = TRUE)
        write(paste("______________________________________________________"), 
            file = file.Restart_index, append = TRUE)
    }
    else {
        if (SCO == TRUE) {
            write(cbind(vec.alpha, vec.S1, vec.C1, vec.PI), file = file.Restart, 
                append = TRUE, ncolumns = N)
            write(c(vec.LAMBDA, vec.sigma.alpha, vec.THETA, vec.beta), 
                file = file.Restart2, append = TRUE)
            write(t(rbind(vec.S2, vec.C2)), file = file.Restart_REFSTD, 
                append = TRUE, ncolumns = length(n_REFSTD))
            write(paste("______________________________________________________"), 
                file = file.Restart_index, append = TRUE)
            write(paste("\t   REstarting values.txt "), file = file.Restart_index, 
                append = TRUE)
            write(paste("Row 1 : alpha parameters for all M = ", 
                N, " study(ies)\t    "), file = file.Restart_index, 
                append = TRUE)
            write(paste("Row 2 : sensitivity of test under evaluation (S1) for all M = ", 
                N, " study(ies)\t    "), file = file.Restart_index, 
                append = TRUE)
            write(paste("Row 3 : specificity of test under evaluation (C1) for all M = ", 
                N, " study(ies)\t    "), file = file.Restart_index, 
                append = TRUE)
            write(paste("Row 4 : prevalence for all M = ", N, 
                " study(ies)\t    "), file = file.Restart_index, 
                append = TRUE)
            write(paste("______________________________________________________"), 
                file = file.Restart_index, append = TRUE)
            write(paste("\t   REstarting values 2.txt "), file = file.Restart_index, 
                append = TRUE)
            write(paste("Column 1 : LAMBDA parameter\t    "), 
                file = file.Restart_index, append = TRUE)
            write(paste("Column 2 : sigma alpha parameter\t    "), 
                file = file.Restart_index, append = TRUE)
            write(paste("Column 3 : THETA parameter\t    "), 
                file = file.Restart_index, append = TRUE)
            write(paste("Column 4 : beta parameter\t    "), file = file.Restart_index, 
                append = TRUE)
            write(paste("______________________________________________________"), 
                file = file.Restart_index, append = TRUE)
        }
    }
    if (Gold_Std == FALSE) {
        if (Gold_se == TRUE) {
            write(paste("\t   REstarting REFSTD.txt "), file = file.Restart_index, 
                append = TRUE)
            write(paste("Row 1 : specificity of reference test (C2) \t    "), 
                file = file.Restart_index, append = TRUE)
        }
        else {
            if (Gold_sp == TRUE) {
                write(paste("\t   REstarting REFSTD.txt "), file = file.Restart_index, 
                  append = TRUE)
                write(paste("Row 1 : sensitivity of reference test (S2) \t    "), 
                  file = file.Restart_index, append = TRUE)
            }
            else {
                write(paste("\t   REstarting REFSTD.txt "), file = file.Restart_index, 
                  append = TRUE)
                write(paste("Row 1 : sensitivity of reference test (S2) \t    "), 
                  file = file.Restart_index, append = TRUE)
                write(paste("Row 2 : specificity of reference test (C2) \t    "), 
                  file = file.Restart_index, append = TRUE)
            }
        }
    }
    cat(paste("The files created during the Gibbs sampler process are in  \"", 
        getwd(), "\" ", sep = ""))
}
