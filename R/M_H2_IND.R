M_H2_IND <-
function (r, y, a, N, low.beta, up.beta, x) 
{
    acc_rate = 0
    alpha.G = sum(0.5 * (1 - y))
    beta.G = 0.5 * sum((1 - y) * (r + 0.5 * a)^2)
    alpha.IG = sum(0.5 * y)
    beta.IG = 0.5 * sum(y * (r - 0.5 * a)^2)
    can <- rinvgamma(1, alpha.IG, beta.IG)
    ratio = dgamma(can, alpha.G, rate = beta.G)/dgamma(x, alpha.G, 
        rate = beta.G)
    if (is.na(ratio) == TRUE) {
        files.remove()
    }
    if (is.na(ratio) == TRUE) {
        cat(paste("Unsuitable initial values were provided. "))
        stop("Please respecify and call HSROC() again.\n")
    }
    aprob = min(1, ratio)
    u = runif(1)
    if (u < aprob) {
        x = can
        acc_rate = 1
    }
    if (x < low.beta) {
        xx = low.beta
    }
    else {
        if (x > up.beta) {
            xx = up.beta
        }
        else {
            xx = x
        }
    }
    return(c(xx, acc_rate, aprob))
}
