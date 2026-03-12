#' FDID example dataset
#'
#' A long-format panel dataset for demonstrating the fdid package.
#'
#' @format A data frame with 11973 rows and 17 columns:
#' \describe{
#'   \item{provid}{Province ID}
#'   \item{countyid}{County ID}
#'   \item{zupu}{Genealogy book count}
#'   \item{pczupu}{Genealogy book density (per capita); 45\% of counties have zero}
#'   \item{lnpczupu}{Log-transformed genealogy density: \code{log(pczupu + 1)}; used as a continuous treatment in Xu, Zhao, and Ding (2026)}
#'   \item{anyzupu}{Indicator: any genealogy book present}
#'   \item{avggrain}{Average grain output}
#'   \item{nograin}{Indicator: no grain data}
#'   \item{urban}{Urban population share}
#'   \item{dis_bj}{Distance to Beijing}
#'   \item{dis_pc}{Distance to provincial capital}
#'   \item{rice}{Rice cultivation indicator}
#'   \item{minority}{Minority population share}
#'   \item{edu}{Education level}
#'   \item{lnpop}{Log population}
#'   \item{year}{Year (1954--1966)}
#'   \item{mortality}{Mortality rate}
#' }
"mortality"