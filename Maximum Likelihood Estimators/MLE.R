#install.packages("bbmle")
library(bbmle)

# 21. Двойное цензурирование. Исследовать зависимость ширины доверительного интервала от размера выборки

start.seed <- 42
set.seed(start.seed)

# Проверим на какой-то выборке адекватность оценок

# Выборка X из гамма-распределения
n <- 500
X <- rgamma(n, shape=2, scale=2)
# Промоделируем U и V с каким-то распределением 
U <- rexp(n, rate=0.5)
V <- U + runif(n, 0, 10)

# Z -- цензурированное значение X
Z <- pmax(U, pmin(X,V))
delta1 <- as.numeric(X<U)
delta2 <- as.numeric(X<V)

# Слагаемые логарифма функции правдоподобия, зависящие от неизвестных параметров.
# Другие слагаемые не зависят от параметров, но зависят от ф.р. U и V, которые неизвестны.
# (вывод см. в файле log-likelihood-derivation.pdf)
minusLL <- function(shape, scale) {
    -sum(delta1*pgamma(Z, shape=shape, scale=scale, log.p=T)+
             (1-delta2)*pgamma(Z, shape=shape, scale=scale, log.p=T, lower.tail=F)+
             (1-delta1)*delta2*(dgamma(Z, shape=shape, scale=scale, log=T)))
}

model <- mle2(minusLL, start=list(shape=3,scale=0.5),
              method = "L-BFGS-B", 
              lower=c(shape=0,scale=1e-10))
confint(model, level=0.95)
confint(model, level=0.99)

# График (корня из) профиля правдоподобия.
# То что он похож на |x| (в окрестности оценки) говорит о том, что модель адекватна (профиль хорошо приближается параболой)
# т.е. (асимптотическим) доверительным интервалам для параметров можно доверять
plot(profile(model))

# Исследуем зависимость ширины доверительных интервалов от размера выборки
# Для этого будем находить их для выборок размера от 100 до 10,000 с шагом 100 
load("confint_width")
if (recalculate) {
    shape.wconfint <- numeric(0)
    scale.wconfint <- numeric(0)
    for (n in seq(1, 1e2, 1)*100) {
        # Выборка X из гамма-распределения
        X <- rgamma(n, shape=2, scale=2)
        # Промоделируем U и V с каким-то распределением 
        U <- rexp(n, rate=0.5)
        V <- U + runif(n, 0, 10)
        
        Z <- pmax(U, pmin(X,V))
        delta1 <- as.numeric(X<U)
        delta2 <- as.numeric(X<V)
        
        # Слагаемые логарифма функции правдоподобия, зависящие от неизвестных параметров.
        # Другие слагаемые не зависят от параметров, но зависят от ф.р. U и V, которые неизвестны.
        minusLL <- function(shape, scale) {
          -sum(delta1*pgamma(Z, shape=shape, scale=scale, log.p=T)+
                   (1-delta2)*pgamma(Z, shape=shape, scale=scale, log.p=T, lower.tail=F)+(1-delta1)*delta2*((shape-1)*log(Z) - Z/scale - shape*log(scale) - lgamma(shape)))
        }
        
        model <- mle2(minusLL, start=list(shape=3,scale=0.5),
                      method = "L-BFGS-B", 
                      lower=c(shape=0,scale=1e-10))
        ci <- suppressMessages(confint(model, level=0.95))
        
        shape.wconfint <- c(shape.wconfint, ci["shape", 2] - ci["shape", 1])
        scale.wconfint <- c(scale.wconfint, ci["scale", 2] - ci["scale", 1])
    }
    recalculate <- F
    current.seed <- .Random.seed
    save(shape.wconfint, scale.wconfint, recalculate, start.seed, current.seed, file="confint_width")
}

# Построим график зависимости ширины от размера выборки
plot(x=seq(1, 1e2, 1)*100, shape.wconfint, col="blue", xlab="n", ylab="Conf. int. width")
points(seq(1, 1e2, 1)*100, scale.wconfint, col="black", pch=8)
legend("topright", legend = c("shape", "scale"), col = c("blue", "black"), lty = 1)

# После возведения ширины в степень -2, зависимость стала похожей на линейную
# Ширина доверительного интервала пропорциональна (в каком-то смысле?) n^(-1/2)
# (теоретически получается из ас. нормальности в многомерном случае)
plot(x=seq(1, 1e2, 1)*100, (shape.wconfint)^(-2), col="blue", xlab="n", ylab="Conf. int. width")
points(seq(1, 1e2, 1)*100, (scale.wconfint)^(-2), col="black", pch=8)
legend("topright", legend = c("shape", "scale"), col = c("blue", "black"), lty = 1)
