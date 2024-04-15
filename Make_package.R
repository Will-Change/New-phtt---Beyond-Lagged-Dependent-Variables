## ################################
##
## Make Package Codes
##
## ################################

## Remove pkg
#remove.packages("phtt")

## Create/update documentation and (re-)write NAMESPACE
devtools::document()

## CRAN-check pkg
devtools::check()       # check the package

## Install
devtools::install_local("phtt", force = TRUE)

devtools::install_github("lidom/R-package-phtt", force = TRUE)
##
## #################################

## usethis
# Run once to configure package to use pkgdown
#usethis::use_pkgdown()
# Run to build the website
pkgdown::build_site()
# 
## usethis::use_pkgdown_github_pages()



data(Cigar)
## Panel-Dimensions:
print(dim(Cigar))
N <- 46
T <- 30
## Dependent variable:
## Cigarette-Sales per Capita
l.Consumption    <- log(matrix(Cigar$sales, T,N))
## Independent variables:
## Consumer Price Index
cpi        <- matrix(Cigar$cpi, T,N)
## Real Price per Pack of Cigarettes 
l.Price  <- log(matrix(Cigar$price, T,N)/cpi)
## Real Disposable Income per Capita  
l.Income    <- log(matrix(Cigar$ndi,   T,N)/cpi)

## Estimation:
KSS.fit <- KSS(l.Consumption ~ l.Price + l.Income, CV = TRUE)
KSS.summary <- summary(KSS.fit)
print(KSS.summary)
