.onAttach <-
function (...) 
{
    date <- date()
    x <- regexpr("[0-9]{4}", date)
    this.year <- substr(date, x[1], x[1] + attr(x, "match.length") - 
        1)
    cat("##\n## Hierarchical Summary Receiver Operating Characteristic package (HSROC)\n")
    cat("## Copyright (C) 2010-", this.year, " Ian Schiller and Nandini Dendukuri \n", 
        sep = "")
    cat("##\n## Development of HSROC package was supported by grants from the \n")
    cat("## Canadian Institutes of Health Research (MOP #89857) \n")
    cat("## and the Fonds de la Recherche en Santé Québec. \n##\n")
    require(lattice, quietly = TRUE)
    require(coda, quietly = TRUE)
    require(MASS, quietly = TRUE)
    require(MCMCpack, quietly = TRUE)
}
