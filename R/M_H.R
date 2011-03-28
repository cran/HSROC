M_H <-
function (r, y, a, N, d1, d0, e.i.j, low.beta, up.beta, x) 
{
    Y.pos = sum(y)
    Y.neg = N - Y.pos
    prob.pos = mean(y)
    Y = rbinom(n = 1, size = 1, prob.pos)
    alpha.IG = sum(0.5 * y)
    beta.IG = (0.5 * sum(y * (r - 0.5 * (a + d1 * e.i.j))^2))
    alpha.G = sum(0.5 * (1 - y))
    beta.G = 0.5 * sum((1 - y) * (r + 0.5 * (a + d0 * e.i.j))^2)
    if (Y == 1) {
        can <- rgamma(1, alpha.G, scale = 1/beta.G)
        ratio = dinvgamma(can, alpha.IG, beta.IG)/dinvgamma(x, 
            alpha.IG, beta.IG)
        aprob = min(1, ratio)
        u = runif(1)
        if (u < aprob) {
            x = can
        }
    }
    else {
        can <- rinvgamma(1, alpha.IG, beta.IG)
        ratio = dgamma(can, alpha.G, rate = beta.G)/dgamma(x, 
            alpha.G, rate = beta.G)
        aprob = min(1, ratio)
        u = runif(1)
        if (u < aprob) {
            x = can
        }
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
    results = c(x, Y)
    return(results)
}
