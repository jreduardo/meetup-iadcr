## ---- include = FALSE----------------------------------------------------

##-------------------------------------------
## Definições knitr
library(knitr)

opts_chunk$set(
    cache = FALSE,
    fig.path = "figures/",
    cache.path = "cache/",
    fig.align = "center",
    dev.args = list(family = "Palatino")
    )

##-------------------------------------------
## Definições lattice
library(lattice)
library(latticeExtra)

ps <- list(
    box.rectangle = list(fill = c("gray70")),
    box.umbrella = list(lty = 1),
    dot.symbol = list(pch = 19),
    dot.line = list(col = "gray50", lty = 3),
    plot.symbol = list(pch = 19),
    strip.background = list(col = c("gray80", "gray50"))
    )
trellis.par.set(ps)


## ------------------------------------------------------------------------

## Carrega o conjunto de dados
load("ninfas.rda")
str(ninfas)


## ------------------------------------------------------------------------

## Visualizando graficamente
xyplot(ntot ~ dias | cult, data = ninfas,
       jitter.x = TRUE, grid = TRUE)

## Verificando a relação média-variância
mv <- aggregate(ntot ~ dias + cult, data = ninfas,
                FUN = function(x)
                    c(mean = mean(x), var = var(x)))

xlim <- ylim <- extendrange(c(mv$ntot), f = 0.03)
xyplot(ntot[, "var"] ~ ntot[, "mean"],
       data = mv, grid = TRUE,
       type = c("p", "smooth"),
       xlab = "Média amostral",
       ylab = "Variância amostral") +
    layer(panel.abline(a = 0, b = 1, lty = 2))



## ------------------------------------------------------------------------

## Preditores Considerados
f0 <- ntot ~ bloc
f1 <- ntot ~ bloc + cult
f2 <- ntot ~ bloc + cult + dias
f3 <- ntot ~ bloc + cult * dias

## Ajustando os modelos Poisson
m0P <- glm(f0, data = ninfas, family = poisson)
m1P <- glm(f1, data = ninfas, family = poisson)
m2P <- glm(f2, data = ninfas, family = poisson)
m3P <- glm(f3, data = ninfas, family = poisson)

## Ajustando os modelos Quasi-Poisson
m0Q <- glm(f0, data = ninfas, family = quasipoisson)
m1Q <- glm(f1, data = ninfas, family = quasipoisson)
m2Q <- glm(f2, data = ninfas, family = quasipoisson)
m3Q <- glm(f3, data = ninfas, family = quasipoisson)

## Ajustando os modelos Binomial Negativo
library(MASS)
m0B <- glm.nb(f0, data = ninfas)
m1B <- glm.nb(f1, data = ninfas)
m2B <- glm.nb(f2, data = ninfas)
m3B <- glm.nb(f3, data = ninfas)


## ------------------------------------------------------------------------

##-------------------------------------------
## Testes de razão de verossimilhanças
anova(m0P, m1P, m2P, m3P, test = "Chisq")
anova(m0Q, m1Q, m2Q, m3Q, test = "Chisq")
anova(m0B, m1B, m2B, m3B, test = "Chisq")


## ------------------------------------------------------------------------

## Lista dos modelos ajustados
modelList <- list("Poisson" = m2P, "Quasi" = m2Q, "BinNeg" = m2B)

## Comparando via logVerossimilhança
sapply(modelList, logLik)

## Estimativas dos parâmetros
## lapply(modelList, function(x)
##     cbind("Est" = coef(x), "SE" = sqrt(diag(vcov(x)))))
car::compareCoefs(m2P, m2Q, m2B)


## ------------------------------------------------------------------------

## Resumo do modelo Quasi-Poisson
summary(m2Q)

## Resumo do modelo Binomial Negativo
summary(m2B)


## ------------------------------------------------------------------------

## Gráficos padrão, cuidado!
par(mfrow = c(2, 3))
plot(m2Q, which = 1:6)

par(mfrow = c(2, 3))
plot(m2B, which = 1:6)

## Análise de diagnóstico
boot::glm.diag.plots(m2Q)
boot::glm.diag.plots(m2B)

## Gráficos qauntil-quantil com envelopes simulados
hnp::hnp(m2Q)
hnp::hnp(m2B)


## ------------------------------------------------------------------------

##-------------------------------------------
## Predição

## Obtendo o preditor linear que considera o efeito médio de blocos
library(doBy)
X <- LSmatrix(lm(ntot ~ bloc + cult + dias, ninfas),
              effect = c("cult", "dias"))

## Estimando as contagens médias com intervalo de confiança
library(multcomp)

## pela Poisson
aux <- exp(confint(glht(m2P, linfct = X),
               calpha = univariate_calpha())$confint)
colnames(aux) <- c("fit", "lwr", "upr")
pred <- aux

## pela Quasi-Poisson
aux <- exp(confint(glht(m2Q, linfct = X),
                   calpha = univariate_calpha())$confint)
## colnames(aux) <- c("Qfit", "Qlwr", "Qupr")
pred <- rbind(pred, aux)

## pela Binomial Negativa
aux <- family(m2B)$linkinv(confint(glht(m2B, linfct = X),
                   calpha = univariate_calpha())$confint)
## colnames(aux) <- c("Qfit", "Qlwr", "Qupr")
pred <- rbind(pred, aux)

aux <- expand.grid(cult = levels(ninfas$cult),
                   dias = unique(ninfas$dias))
pred <- cbind(rbind(aux, aux, aux), pred)
pred$modelo <- as.factor(rep(c("1.Poisson", "2.Quasi", "3.BinNeg"),
                           each = nrow(aux)))
pred <- pred[, c(6, 1:5)]
pred <- pred[order(pred$cult, pred$dias, pred$modelo), ]

## Para fazer o painel com segmentos
source(paste0("https://gitlab.c3sl.ufpr.br/leg/legTools/raw/",
              "issue%2315/R/panel.segplot.by.R"))

segplot(dias ~ lwr + upr | cult,
        centers = fit, data = pred,
        horizontal=FALSE, draw=FALSE,
        lwd = 2, grid = TRUE,
        panel = panel.segplot.by, groups = pred$model, f = 0.15,
        pch = 1:nlevels(pred$model)+3,
        as.table = TRUE,
        key = list(type = "o", divide = 1,
                   lines = list(pch = 1:nlevels(pred$model)+3,
                              lty = 1),
                   text = list(c("Poisson", "Quasi-Poisson",
                                 "Binomial Negativa"))))


