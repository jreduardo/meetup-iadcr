## ---- include = FALSE----------------------------------------------------

##-------------------------------------------
## Definições knitr
library(knitr)

knitr::knit_hooks$set(
    ## Para fontsize de chunk R
    mysize = function(before, options, envir) {
        if (before) return(options$size)
    }
)

opts_chunk$set(
    mysize = TRUE,
    tidy = FALSE,
    cache = TRUE,
    echo = FALSE,
    size = "\\normalsize",
    out.width = "1\\textwidth",
    fig.path = "figures/",
    cache.path = "cache/",
    fig.align = "center",
    fig.height = 6,
    fig.width = 9,
    dev.args = list(family = "Palatino")
    )

##-------------------------------------------
## Packages
library(lattice)
library(latticeExtra)


## ---- fig.cap = "Probabilidades para modelos Poisson"--------------------

## Definindo parametros e calculando as probabilidas
lambdas <- c(3, 8, 15)
x <- 0:30; xx <- rep(x, 3)
px <- NULL
for(i in 1:3) px <- c(px, dpois(x, lambdas[i]))

## Criando categorias para split da lattice
caso <- rep(c("1", "2", "3"), each = length(x))

## Definindo nome para os splits da lattice
fl <- as.expression(lapply(lambdas,
    function(x){ bquote(lambda==.(x)) }))
    
xyplot(px ~ xx | caso, type = c("h", "g"),
       lwd = 3, xlab = "y", ylab = expression(P(Y==y)),
       layout = c(NA, 1), col = 1,
       between = list(x = 0.2, y = 0.3),
       strip = strip.custom(bg = "gray90",
                          factor.levels = fl))


## ---- eval = FALSE, echo = TRUE------------------------------------------
## 
## model <- glm(y ~ preditor, family = poisson)
## 

## ---- fig.cap = "Probabilidades para a não verificação de equidispersão"----

library(compoisson)

## Parametros da distribuição
lambdas <- c(1.1, 8, 915); nus <- c(0.25, 1, 2.5)
medias <- mapply(com.mean, lambda = lambdas, nu = nus)
variancias <- mapply(com.var, lambda = lambdas, nu = nus)

## Calculando as probabilidades
y <- 0:30; yy <- rep(x, 3)
py.com <- py.pois <- NULL
for(i in 1:3) py.com <- c(py.com, dcom(x, lambdas[i], nus[i]))
for(i in 1:3) py.pois <- c(py.pois, dpois(x, medias[i]))

## Criando categorias para split da lattice
caso <- rep(c("1", "2", "3"), each = length(x))
fl <- expression(Superdispersão%->%mu[Y]>sigma[Y]^2,
                 Equidispersão%->%mu[Y]==sigma[Y]^2,
                 Subdispersão%->%mu[Y]<sigma[Y]^2
                 )

xyplot(py.com ~ c(yy - 0.14) | caso, type = c("h", "g"),
       lwd = 2.5, xlab = "y", ylab = expression(P(Y==y)),
       col = "forestgreen", ylim = c(-0.005, 0.24),
       xlim = extendrange(y), layout = c(NA, 1),
       between = list(x = 0.2,y = 0.3),
       strip = strip.custom(bg = "gray90",
                            factor.levels = fl)) + 
  as.layer(xyplot(py.pois ~ c(yy + 0.14) | caso, 
                  lwd = 2.5, col = 1,
                  type = "h"))


## ---- eval = FALSE, echo = TRUE------------------------------------------
## 
## ## Support Functions and Datasets for Venables and Ripley's MASS
## library(MASS)
## model <- glm.nb(y ~ preditor)
## 

## ---- eval = FALSE, echo = TRUE------------------------------------------
## 
## model <- glm(y ~ preditor, family = quasipoisson)
## 

## ---- fig.cap = "Contagens que apresentam excesso de zeros"--------------

library(CompGLM)

set.seed(100)
n <- 500
rb <- rbinom(n, 1, 0.98)
for (i in 1:n) {
    rb[i] <- ifelse(rb[i] == 0, 0, rcomp(1, 10, 2))
}

y <- rb
py_real <- prop.table(table(y))

m0 <- glm(y ~ 1, family = poisson)
py_pois<- dpois(sort(unique(y)), exp(m0$coef))

m1 <- glm.comp(y ~ 1)
coefs <- sapply(coef(m1), exp)
## py_quasi <- dpois(sort(unique(y)), exp(m0$coef))
py_cmp <- dcom(sort(unique(y)), coefs[1], coefs[2])

yu <- sort(unique(y))
xyplot(py_real ~ yu, type = "h",
       lwd = 3, grid = T, col = 1,
       xlab = "y",
       ylab = expression(P(Y==y)),
       key = list(space = "right",
                  lines = list(
                      lty=1, col = c(1, 4, "forestgreen"), lwd = 3),
                  text = list(c("Real", "Poisson", "COM-Poisson"))),
       panel = function(x, y, ...) {
           panel.xyplot(x, y, ...)
           panel.lines(x = x - 0.1, y = py_pois, type = "h",
                    col = 4, lwd = 3)
           panel.lines(x = x + 0.1, y = py_cmp, type = "h",
                       col = "forestgreen", lwd = 3)
       })


## ---- eval = FALSE, echo = TRUE------------------------------------------
## 
## ## Political Science Computational Laboratory, Stanford University
## library(pscl)
## 
## hurdle(resp ~ py_preditor | pi_preditor, dist = "poisson")
## 

## ---- eval = FALSE, echo = TRUE------------------------------------------
## 
## ## Political Science Computational Laboratory, Stanford University
## library(pscl)
## 
## zeroinfl(resp ~ py_preditor | pi_preditor, dist = "poisson")
## 

