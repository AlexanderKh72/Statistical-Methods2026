### Вариант 1.
### Реализовать EM-алгоритм для оценивания параметров нормальной смеси. 
### Визуализировать полученное решение.
library(mvtnorm)

# Данные об извержениях гейзера Old Faithful: продолжительность извержения и время между извержениями
data <- read.table("faithful.txt", header = TRUE) |> as.matrix()
n <- nrow(data)

# Построим scatter-plot для данных.
# То, что точки визуально разбиваются на два облака, может означать, что данные неоднородны,
# т.е. распределение представляет собой смесь двух нормальных (2-мерных) распределений
# Плотность распределения: (1 - alpha)*phi(x; mu1, sigma1) + alpha*phi(x; mu2, sigma2)
plot(data)

# Мы не наблюдаем метку, содержащую информацию о том, к какой компоненте относится точка.
# Начальное приближение --- случайная точка, взятая равномерно из (0,1).
# (начальное приближение параметра бернуллиевского распределения метки, если быть точным)
# (можно было бы выбирать начальное приближение параметров mu_i, sigma_i, alpha, и делать E-шаг, чтобы найти Z)
set.seed(4399)
z <- runif(n)

alpha_prev <- +Inf
mu1_prev <- rep(+Inf, 2)
mu2_prev <- rep(+Inf, 2)
sigma1_prev <- rep(+Inf, 4)
sigma2_prev <- rep(+Inf, 4)

# Устанавливаем максимальное число итераций.
# Цикл завершается, когда изменение параметров достаточно мало
max.iter <- 1000
for (i in 1:max.iter) {
    # Оценка параметров (M-шаг)
    alpha <- mean(z)
    mu1 <- colSums((1-z) * data) / (n - sum(z))
    mu2 <- colSums(z * data) / (sum(z))
    sigma1 <- (t(data) - mu1) %*% ((1-z) * t(t(data) - mu1)) / sum((1-z))
    sigma2 <- (t(data) - mu2) %*% (z * t(t(data) - mu2)) / sum(z)
    # Пересчёт z (E-шаг)
    z <- alpha * dmvnorm(data, mu2, sigma2) / ((1-alpha) * dmvnorm(data, mu1, sigma1) + alpha * dmvnorm(data, mu2, sigma2))
    
    # Расчет изменений параметров и проверка критерия остановки
    delta_alpha <- abs(alpha - alpha_prev)
    delta_mu1 <- sqrt(sum((mu1 - mu1_prev)^2))
    delta_mu2 <- sqrt(sum((mu2 - mu2_prev)^2))
    delta_sigma1 <- norm(sigma1 - sigma1_prev, "F")
    delta_sigma2 <- norm(sigma2 - sigma2_prev, "F")
    
    if (all(c(delta_alpha, delta_mu1, delta_mu2, delta_sigma1, delta_sigma2) < 1e-6)) break
    
    alpha_prev <- alpha
    mu1_prev <- mu1
    mu2_prev <- mu2
    sigma1_prev <- sigma1
    sigma2_prev <- sigma2
    
    i0 <- i
}   

# Число потребовавшихся итераций
# (на разных сидах получалось около 40)
i0

# Визуализация: на фоне данных отрисованы линии уровня плотностей многомерных распределений с оценёнными параметрами.
x <- seq(min(data[,1]), max(data[,1]), by=.1)
y <- seq(min(data[,2]), max(data[,2]), by=.1)
net1 <- outer(x, y, function(x,y) dmvnorm(cbind(x,y), mu1, sigma1))
net2 <- outer(x, y, function(x,y) dmvnorm(cbind(x,y), mu2, sigma2))
plot(data)
contour(x, y, net1, add=T)
contour(x, y, net2, add=T)
# Эллипсоиды не пересекаются, т.е. два распределения действительно отличаются (мы и не сомневались, впрочем).

# Дополнительно проверим "качество классификации", а именно, насколько полученные значения z близки к 0 или 1.
z[z > 0.01 & z < 0.99]
# Всего  два наблюдения не попали в 0.01 окрестность 0 или 1, да и они достаточно близки к 0 или 1.
# Значит можно достаточно уверенно говорить о том, к какому типу извержения относится каждое



