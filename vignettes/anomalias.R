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
load("anomalias.rda")
str(anomalias)


## ------------------------------------------------------------------------

## Visualizando graficamente
xyplot(nanom/ncel ~ l2taxa, groups = dose,
       data = anomalias,
       grid = TRUE, type = c("p", "spline"),
       xlab = expression(log[2]~taxa~de~liberação),
       ylab = "Proporção de anomalias",
       auto.key = list(space = "right", title = "dose",
                       cex.title = 1))


## ------------------------------------------------------------------------

## Preditores Considerados
f1 <- nanom ~ offset(log(ncel)) + dose
f2 <- nanom ~ offset(log(ncel)) + l2taxa
f3 <- nanom ~ offset(log(ncel)) + dose + l2taxa
f4 <- nanom ~ offset(log(ncel)) + dose * l2taxa

## Ajustando os modelos Poisson
m1 <- glm(f1, data = anomalias, family = poisson)
m2 <- glm(f2, data = anomalias, family = poisson)
m3 <- glm(f3, data = anomalias, family = poisson)
m4 <- glm(f4, data = anomalias, family = poisson)


## ------------------------------------------------------------------------

##-------------------------------------------
## Testes de razão de verossimilhanças

## Testando a inclusão de dose dado l2taxa
anova(m2, m3, test = "Chisq")

## Testando a inclusão de l2taxa dado dose
anova(m1, m3, test = "Chisq")

## Testando a inclusão dos efeitos para interação
anova(m3, m4, test = "Chisq")


## ------------------------------------------------------------------------

## Resumo do modelo
summary(m4)


## ------------------------------------------------------------------------

## Teste de adequação do ajuste 
with(m4, pchisq(q = deviance, df = df.residual, lower.tail = FALSE))

## Gráficos padrão, cuidado!
par(mfrow = c(2, 3))
plot(m4, which = 1:6)

## Análise de diagnóstico
boot::glm.diag.plots(m4)

## Verificar
hnp::hnp(m4)


## ------------------------------------------------------------------------

## Predição
pred <- with(anomalias,
             expand.grid(
                 ncel = 1,
                 dose = c("1", "2.5", "5"),
                 l2taxa = seq(min(l2taxa), max(l2taxa), l = 30)
             ))

qn <- qnorm(0.975)

## pela Poisson
aux <- predict(m4, newdata = pred, se.fit = TRUE)
aux <- with(aux, fit + se.fit %*%
                 cbind(lwr = -1, fit = 0, upr = 1) * qn)
pred <- cbind(pred, exp(aux))

xyplot(nanom/ncel ~ l2taxa, groups = dose,
       data = anomalias,
       grid = TRUE,
       xlab = expression(log[2]~taxa~de~liberação),
       ylab = "Proporção de anomalias",
       auto.key = list(title = "dose", columns = 3,
                       cex.title = 1, lines = TRUE)) +
    as.layer(
        xyplot(fit ~ l2taxa,
               groups = dose, data = pred,
               type = "l")
    ) +
    as.layer(
        xyplot(lwr ~ l2taxa,
               groups = dose, data = pred,
               type = "l", lty = 2)
    ) +
    as.layer(
        xyplot(upr ~ l2taxa,
               groups = dose, data = pred,
               type = "l", lty = 2)
    )
    

## ------------------------------------------------------------------------

## Visualizando as distribuições Poisson
pred <- expand.grid(ncel = 1, dose = c("1", "2.5", "5"),
                    l2taxa = c(-3.5, -1, 0, 2), x = 0:6)
lambda <- predict(m4, newdata = pred, type = "response")
pred <- cbind(pred, lambda)
pred$px <- with(pred, dpois(x, lambda))

## Definindo nome para os splits da lattice
fl <- as.expression(lapply(unique(pred$l2taxa),
    function(x){ bquote(l2taxa==.(x)) }))

xyplot(px ~ x | factor(l2taxa), groups = dose,
       data = subset(pred, dose == 2.5),
       auto.key = list(title = "dose", lines = TRUE, points = FALSE,
                       columns = 3, cex.title = 1),
       type = "h", grid = TRUE,
       as.table = TRUE,
       lwd = 2.5, xlab = "y: contagem de anomalias cromossômicas",
       ylab = expression(P(Y==y)),
       ylim = c(0, max(pred$px))*1.1,
       strip = strip.custom(bg = "gray90",
                            factor.levels = fl)) +
    as.layer(
        xyplot(px ~ (x - 0.14) | l2taxa, groups = dose, lwd = 2.5,
               data = subset(pred, dose == 1), type = "h")
        ) +
    as.layer(
        xyplot(px ~ (x + 0.14) | l2taxa, groups = dose, lwd = 2.5,
               data = subset(pred, dose == 5), type = "h")
        )


