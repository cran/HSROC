f.rij <-
function (l, u, r) 
{
    if (r < l) {
        b = l
    }
    else {
        if (r > u) {
            b = u
        }
        else {
            b = r
        }
    }
    return(b)
}
