.onAttach <-
function (...)
{
    date <- date()
    x <- regexpr("[0-9]{4}", date)
    this.year <- substr(date, x[1], x[1] + attr(x, "match.length") -
        1)
    packageStartupMessage("##\n## Hierarchical Summary Receiver Operating Characteristic package (HSROC)")
    packageStartupMessage("## Copyright (C) 2010-", this.year,
        " Ian Schiller and Nandini Dendukuri ", sep = "")
    packageStartupMessage("##\n## Development of HSROC package was supported by grants from the ")
    packageStartupMessage("## Canadian Institutes of Health Research (MOP #89857) ")
    packageStartupMessage("## and the Fonds de la Recherche en Sant\u00E9 Qu\u00E9bec. ")
}
