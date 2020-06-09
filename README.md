# comoOdeCpp

RCpp version of CoMo Consortium's COVID-19 transmission model


by Ricardo Aguas, Sompob Saralamba, Bo Gao


## Install

- GUI option:

 - Download the latest release from [https://github.com/ocelhay/comoOdeCpp/releases](https://github.com/ocelhay/comoOdeCpp/releases), e.g. `comoOdeCpp_13.13.1.tar.gz`

 - RStudio >> Tools >> Install Packages >> Install from: Package Archive File (.tar.gz) >> Package archive: Browse to select `comoOdeCpp_13.13.1.tar.gz`


- R commandline option:

 - Download the latest release from [https://github.com/ocelhay/comoOdeCpp/releases](https://github.com/ocelhay/comoOdeCpp/releases), e.g. `comoOdeCpp_13.13.1.tar.gz`

 - `install.packages("path_to_download/comoOdeCpp_13.13.1.tar.gz", repos = NULL, type="source")`


- Using `devtools` package

```r

install.packages("devtools")
library("devtools")
install_github("ocelhay/comoOdeCpp", subdir="comoOdeCpp")

```


## Usage:

Add `library("comoOdeCpp")` at the beginning of your script, or `require("comoOdeCpp")` if you are writing a package


Examples use with ode function:

```r

covidOdeCpp_reset()
out0 <- ode(y = Y, times = times, func = covidOdeCpp, parms = parameters,
            input=vectors0, A=A,
            contact_home=contact_home, contact_school=contact_school,
            contact_work=contact_work, contact_other=contact_other,
            popbirth_col2=popbirth[,2], popstruc_col2=popstruc[,2],
            ageing=ageing,
            ifr_col2=ifr[,2], ihr_col2=ihr[,2], mort_col=mort)


covidOdeCpp_reset()
out0 <- ode(y = Y, times = times, method = "euler", hini = 0.05, func = covidOdeCpp, parms = parameters,
            input=vectors0, A=A,
            contact_home=contact_home, contact_school=contact_school,
            contact_work=contact_work, contact_other=contact_other,
            popbirth_col2=popbirth[,2], popstruc_col2=popstruc[,2],
            ageing=ageing,
            ifr_col2=ifr[,2], ihr_col2=ihr[,2], mort_col=mort)
)


```
