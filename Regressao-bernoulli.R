
rm(list = ls())

## Instalação e Leitura de Pacotes______________________________________________

install.packages('ggplot2')
install.packages('stringi')
install.packages('kableExtra')
install.packages('knitr')
install.packages('dplyr')
install.packages('gridExtra')
install.packages('grid')
install.packages('reshape2')
install.packages('doBy')
install.packages('gdata')
install.packages("RColorBrewer")
install.packages("pROC")


require('ggplot2')
require('stringi')
require('kableExtra')
require('knitr')
require('dplyr')
require('gridExtra')
require('grid')
require('reshape2')
require('doBy')
require('gdata')
library("RColorBrewer")
library("pROC")


# Leitura dos dados_____________________________________________________________

# Dados de interesse: Fatores que influenciam na admissão de estudantes em Universidades Americanas

# Variáveis explicativas: 
# - GRE (Graduate Record Examination)
# - GPA (Grade Point Average)
# - Rank (Prestígio da escola do ensino médio)

# Variável resposta (admit):
# - 0 não admissão
# - 1 admissão

dados <- read.csv("https://stats.idre.ucla.edu/stat/data/binary.csv")
#transformação de var num - > factor
dados$rank = factor(dados$rank)
head(dados)

## Análise Descritiva___________________________________________________________

par(mfrow = c(1, 3))

# Gráfico Boxplot Admissão X GRE (Graduate Record Examination)
g1 = ggplot(dados, aes(x = factor(admit), y = gre)) + 
  geom_boxplot(fill = "lightblue")+
  labs(y = "GRE", x = "Admissão")+
  theme_bw()+
  theme(legend.position = "top", text = element_text(size = 8))

# Gráfico Boxplot Admissão X GPA (Grade Point Average)
g2 = ggplot(dados, aes(x=factor(admit), y=gpa)) + 
  geom_boxplot(fill = "yellow")+
  labs(y = "GPA", x = "Admissão")+
  theme_bw()+
  theme(legend.position = "top", text = element_text(size = 8))

# Gráfico Admissão X Prestígio da escola
g3 = ggplot(dados, aes(rank, ..count..)) + geom_bar(aes(fill = factor(admit)), position = "dodge")+
  labs(y = "Frequência", x = "Rank (Prestígio)", fill = "Admissão")+
  theme(legend.position = "top" , text = element_text(size = 8), legend.key.size = unit(0.4, "cm"))

grid.arrange(arrangeGrob(g1, g2, ncol = 2, nrow = 1), 
             arrangeGrob(g3, nrow = 1, ncol = 1), widths = c(1, 1))

## Ajuste do Modelo de Regressão Logística______________________________________

modelo = glm(formula = admit ~ ., family = "binomial", data = dados)
summary(modelo)

## Intervalo de confiança --------------------------------------------------------

confint(modelo)

# Com razão de chance
exp(confint(modelo))

# Gráfico -----------------------------------------------------------------------

par(mfrow = c(1, 1))

dados2 <- with(dados, data.frame(gre = rep(seq(from = 200, to = 800, length.out = 100), 4), 
                                 gpa = 3, rank = factor(rep( 1 : 4, each = 100))))
dados3 <- cbind(dados2, PredictedProb = predict(modelo, newdata = dados2, type = "response"))

# Gráfico da probabilidade de Sucesso com GPA fixado = 3
# Obs: variando o GRE visualizamos as escolas 1, 2, 3, 4 conforme admissião

ggplot(dados3, aes(x = gre, y = PredictedProb)) + 
  geom_line(aes(colour = rank), size = 1)+
  labs(x = "GRE", y = "Probilidade de Sucesso", color = "Rank")+
  theme_bw() +
  theme(text = element_text(size = 10))


# Análise do ajuste de modelo - Deviance________________________________________

deviance = with(modelo, cbind(Deviance = deviance, 
                              "Graus de Liberdade" = df.residual,
                              "P-valor" = pchisq(deviance, df.residual, 
                                                 lower.tail=FALSE)))
deviance
# Obs: 
# H0: O modelo está bem ajustado
# H1: O modelo não está bem ajustado
# com p-valor < 0.05 rejeitamos H0, o modelo não está bem ajustado
# problema: aproximação assitotica, algumas distribuições precisam de amostras muito grandes,
# para que a aproximação seja razoável, utilizaremos outros recursos.

# Tentativa: Análise da diferença de Deviance___________________________________

# modelo nulo: Apenas B0 -> modelo sem nenhuma variavel explicativa
# modelo ajustado

diferenca = modelo$null.deviance - modelo$deviance
diferenca.gl =  modelo$df.null - modelo$df.residual
pvalor = 1 - pchisq(diferenca, diferenca.gl)
pvalor
# Obs:
# com p-valor significativo < 0.05 rejeita a hipotese de igualdade do modelos

Deviance = round(c(modelo$null.deviance, modelo$deviance, diferenca),2)
GL = c(modelo$df.null, modelo$df.residual, diferenca.gl)
Modelo = c("Modelo Nulo", "Modelo Ajustado", "Diferença")
frame = data.frame(Modelo, Deviance, GL)
frame

# Análise dos Resíduos__________________________________________________________

# Função que calcula a Matriz H
retorna.H <- function(modelo){
  X = model.matrix(modelo)
  W = diag(modelo$weights)
  M = solve(t(X)%*%W%*%X)
  H = sqrt(W)%*%X%*%M%*%t(X)%*%sqrt(W)
  # Extrai os elementos da diagonal da Matriz H
  h = diag(H)
  return(h)}

h = retorna.H(modelo)
residuos.modelo = resid(modelo, type = "deviance")/sqrt((1-h))
n <- nrow(dados)
m = 100
matriz.residuos <- matrix(NA, n, m)

for(i in 1:m){
  nresp = rbinom(n, 1, fitted(modelo))
  #nresp <- rpois(n, fitted(modelo))
  fit <- glm(nresp ~ gre + gpa + rank, family = "binomial", data = dados)
  h = retorna.H(fit)
  matriz.residuos[,i] = sort(resid(fit, type = "deviance")/sqrt(1-h))
}

# Cálculo dos percentis
intervalos = apply(matriz.residuos, 1, function(x) 
  quantile(x, c(0.025, 0.975)))

# Cálculo da média
med <- apply(matriz.residuos, 1, mean)

# Amplitude de variação dos resíduos
faixa = range(residuos.modelo, intervalos)

# Gráfico de envelope
qqnorm(residuos.modelo, xlab = "Percentil da Normal Padrão", ylab = "Resíduos", 
       ylim = faixa, pch = 16, cex.lab = 0.7, cex.axis = 0.7, main = " ")

# Adicionando o envelope simulado ao gráfico
par(new = T)
qqnorm(intervalos[1,], axes = F, xlab = " ", ylab = " ", 
       type = "l", ylim = faixa, lty = 1, cex = 0.6, main = "Gráfico Envelope")
par(new = T)
qqnorm(intervalos[2,], axes = F, xlab = " ", ylab = " ", type = "l",
       ylim = faixa, lty = 1, cex = 0.6, main = " ")
par(new = T)
qqnorm(med, axes = F, xlab = "", ylab = "", type = "l",
       ylim = faixa, lty = 2, cex = 0.6, main = " ")


# Comparando os ajustes de modelo pelas funções de ligação______________________

modelo2 = glm(formula = admit ~ ., 
              family = binomial(link = "probit"), data = dados)

modelo3 = glm(formula = admit ~ ., 
              family = binomial(link = "cloglog"), data = dados)

deviance.m = c(modelo$deviance, modelo2$deviance, 
               modelo3$deviance)
aic.m = c(modelo$aic, modelo2$aic, modelo3$aic)
nome = c("Modelo Logit", "Modelo Probit", "Modelo cLog-log")
frame = data.frame(Modelo = nome, Deviance = deviance.m, AIC = aic.m)
frame

## Gráfico de comparação de funções de ligações com rank fixado = 1

novos.dados = data.frame(gre = seq(min(dados$gre), max(dados$gre), 
                                   by = 0.1), gpa = 3, rank = factor(1))

novos.dados$modelo1 = predict(modelo, newdata = novos.dados, type = "response")
novos.dados$modelo2 = predict(modelo2, newdata = novos.dados, type = "response")
novos.dados$modelo3 = predict(modelo3, newdata = novos.dados, type = "response")

dados.long = reshape(novos.dados, direction = "long", varying = list(names(novos.dados)[4:6]),
                     v.names = "Y",
                     idvar = c("gre"),
                     timevar = "Modelo",
                     times = c("Logistico", "Probit", "cLogLog"))

row.names(dados.long) = NULL

# Gráfico GRE(Graduate Record Examination) X Probabilidade de Sucesso
ggplot(dados.long, aes(x = gre, y = Y)) +
  geom_line(aes(colour = Modelo),size = 1)+
  labs(x = "GRE", y="Probilidade de Sucesso", color = "Ligação")+
  theme_bw() +
  theme(text = element_text(size = 10))

## Função de Cálculo da Sensibilidade e Especificidade__________________________

sensib_especif = function(modelo, limiar){
  valor_predito <- ifelse(predict(modelo, type = "response") > limiar, 1, 0)
  valor_predito = factor(valor_predito, levels = c(0, 1))
  valor_observado <- modelo$y
  conf_matrix <- table(valor_observado, valor_predito)
  sensibilidade = conf_matrix[2, 2]/rowSums(conf_matrix)[2]
  especificidade = conf_matrix[1, 1]/rowSums(conf_matrix)[1]
  frame = data.frame(sensibilidade, especificidade)
  return(frame)}

sensib_especif(modelo, limiar = 0.5)


# Obs: A redução do limiar de corte gera aumento da sensibilidade e redução da especificidade.
# O aumento do limiar de corte gera redução da sensibilidade e aumento da especificidade.
# Solução: uma maneira comum de escolha do limiar é justamente onde sensibilidade e especificidade
# se cruzam no gráfico, ou seja, quando elas são iguais.

limiar.vetor = seq(0.0, 1, by = 0.01)
n = length(limiar.vetor)
frame.se = data.frame(sensibidade = rep(NA, n), especificidade = rep(NA, n))
for(i in 1:n)
  frame.se[i,] = sensib_especif(modelo, limiar.vetor[i])

frame.long = data.frame(limiar = c(limiar.vetor, limiar.vetor), y = c(frame.se[,1], frame.se[,2]), 
                        identificador = c(rep("sensibilidade", n), rep("especificidade", n)))

## Gráfico de Sensibilidade / Especificidade
ggplot(frame.long, aes(x = limiar, y = y)) + 
  geom_line(aes(colour = identificador), size = 1)+
  labs(x = "Limiar", y = "Sensibilidade e Especificidade", color = " ")+
  theme_bw() +
  theme(text = element_text(size = 10))


## Curva ROC 

# DEF: A curva ROC é a representação gráfica da sensibilidade x 1 - especificidade.
# Sensibilidade = Taxa de verdadeiros positivos
# Especificidade = Taxa de falsos positivos
# Obs: Assim conseguimos a medida AUC(area under curve), muito comum utilizar a curva ROC
# e a estatística AUC para comparar diferentes modelos de classificação.

prob = predict(modelo, type = c("response"))
dados$prob = prob
#comparando o verdadeiro valor x probabilidades ajustadas
curva.roc <- roc(admit ~ prob, data = dados)
curva.roc$auc

## Gráfico da curva ROC

roc.data = data.frame(x = 1 - curva.roc$specificities, y =  curva.roc$sensitivities)

ggplot(roc.data, aes(x = x, y = y)) + 
  geom_line(col = "blue")+
  labs(x = " 1 - Especificidade", y = "Sensibilidade", color = " ") +
  theme_bw() +
  theme(text = element_text(size = 10)) +
  geom_abline(intercept = 0, slope = 1, col = "gray", lty = 2)

