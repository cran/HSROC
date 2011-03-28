M_H2 <-
function (r, y, a, N, d1, d0, e.i.j, low.beta, up.beta, x) 
{
    alpha.IG = sum(0.5 * y)
    beta.IG = (0.5 * sum(y * (r - 0.5 * (a + d1 * e.i.j))^2))
    alpha.G = sum(0.5 * (1 - y))
    beta.G = 0.5 * sum((1 - y) * (r + 0.5 * (a + d0 * e.i.j))^2)
    can <- rgamma(1, alpha.G, scale = 1/beta.G)
    ratio = dinvgamma(can, alpha.IG, beta.IG)/dinvgamma(x, alpha.IG, 
        beta.IG)
    aprob = min(1, ratio)
    u = runif(1)
    if (u < aprob) {
        x = can
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
    return(xx)
}
