#### Моделирование распределения Бирнбаума-Саундерса адаптивным методом отбора (вариант 15) ####

# Плотность
dBirnSaund <- function(x) {
  gamma <- 1
  (sqrt(x) + sqrt(1/x)) / (2*gamma*x) * dnorm((sqrt(x) - sqrt(1/x)) / gamma)
}

# ВАЖНО! Если считать log(dBirnSaund(x)), то будут -Inf, хотя всего-то -1000
log.dBirnSaund<- function(x) {
  gamma <- 1
  log(sqrt(x)+1/sqrt(x)) - log(2*gamma*x) + 1/log(2*pi) - (x^2 - 2*x + 1) / (2*x*gamma^2)
}

# Функция распределения
pBirnSaund <- function(x) {
  gamma <- 1
  pnorm((sqrt(x) - sqrt(1/x)) / gamma)
}

# Логарифм плотности НЕ выпуклый вверх: у него одна точка перегиба
# Поэтому параметр \gamma фиксирован и захардкоден правый хвост мажоранты


# Точки z --- точки пересечения прямых, точки изломов мажоранты
# (формула выведена на листочке и переписана в код)
z <- function(x,h) {
  j <- 1:(length(x)-3)
  ((h[j+2] - h[j]) + ((h[j+1] - h[j])/(x[j+1] - x[j]) * x[j] - (h[j+3] - h[j+2])/(x[j+3]-x[j+2])*x[j+2])) / ((h[j+1] - h[j])/(x[j+1] - x[j]) - (h[j+3] - h[j+2])/(x[j+3] - x[j+2]))
}

# Находит кусочно-линейную мажоранту по наборам x и h
# Возвращает список из 3х элементов:
# $corners --- вектор точек излома мажоранты. Если x состоит из k точек, то точек излома 2*(k-2) - 1 
# $slope --- вектор длины 2*(k-2), состоящий из угловых коэффициентов мажоранты на участках линейности
# $intercept --- то же самое со свободными коэффициентами
u <- function(x, h) {
  k <- length(x)
  z.k <- z(x, h)
  
  corners <- c(rbind(x[2:(k-2)], z.k), x[k-1])
  
  # индексы точек, через которые проходят прямые на каждом участке
  # (заметим, что...)
  j <- cumsum(rep(c(+2, -1), k - 2))
  
  slope <- (h[j+1] - h[j]) / (x[j+1] - x[j])

  # Момент хардкодинга из-за того, что логарифм плотности не выпуклый вверх
  slope[length(slope)] <- -1/2  # h(x) ~ -1/(2*gamma^2) * x, при x -> Inf
  
  intercept <- h[j] - x[j] * slope
  
  list(corners=corners, slope=slope, intercept=intercept)
}

# интеграл мажоранты (вектор интегралов по отрезкам, где мажоранта линейна)
integrate.u <- function(u.k) {
  corners <- c(0, u.k$corners, +Inf)
  j <- 1:(length(corners)-1)
  (exp(u.k$slope[j]*corners[j+1] + u.k$intercept[j]) - exp(u.k$slope[j]*corners[j] + u.k$intercept[j])) / u.k$slope[j]
}

# Моделирование n случайных величин из распределения Бирнбаума-Саундерса
# Алгоритм по статье (Gilks & Wild, 1990)
# (с другим способом построения мажоранты)
rBirnSaund <- function(n) {
  k <- 4  # текущее количество точек, по которым строится мажоранта
  x <- c(0.1, 0.25, 0.7, 0.75)  # точки, в которых был посчитан логарифм плотности
  h <- log.dBirnSaund(x)  # значение логарифма плотности в точках из x
  
  u.k <- u(x,h)
  integrals <- integrate.u(u.k)
  corners <- c(0, u.k$corners)
  
  res.sample <- numeric(0)
  while (length(res.sample) < n) {
    #seed <- .Random.seed
    accepted <- F
    # Sampling step
    ## Sample x* from u_k
    j <- sample(1:length(corners), size=1, prob=integrals)  # выбираем компоненту
    x.star <- (log(integrals[j] * u.k$slope[j] * runif(1) + 
                     exp(u.k$slope[j]*corners[j] + u.k$intercept[j])) - u.k$intercept[j]) / u.k$slope[j] # моделируем методом обратной функции
    u.k.star <- u.k$slope[j]*x.star + u.k$intercept[j]  # значение мажоранты в смоделирванной точке 
    ## Sample w from Uniform(0,1)
    w <- runif(1)
    ## Squeezing test
    i <- max(which(x < x.star)) |> suppressWarnings()  # между какими иксами попала x* (можно было через j посчитать)
    if ((i < k) && (i != -Inf)) {  # отрезки, где определена миноранта
      l.k <- h[i]+((h[i+1] - h[i]) / (x[i+1] - x[i])) * (x.star - x[i])
      if (w <= exp(l.k - u.k.star)) {
        res.sample <- c(res.sample, x.star)
        accepted <- T
      }
    }
    ## Rejection test
    if(!accepted) {
      h.star <- log.dBirnSaund(x.star)
      if (w <= exp(h.star - u.k.star)) {
        res.sample <- c(res.sample, x.star)
        accepted <- T
      }
      # Updating step
      if (k < 1000 && i < k) {
        if (i == -Inf) {
          x <- c(x.star, x)
          h <- c(h.star, h)
        }
        # else if (i == k) {
        #   x <- c(x, x.star)  # Мы не добавляем точки, больше 0.75 
        #   h <- c(h, h.star)  # (точка перегиба чуть правее 0.75)
        # }
        else if (i < k) {
          x <- c(x[1:i], x.star, x[(i+1):k])
          h <- c(h[1:i], h.star, h[(i+1):k])
        }
        k <- k + 1
        u.k <- u(x,h)
        integrals <- integrate.u(u.k)
        corners <- c(0, u.k$corners)
      }
    }
  }
  res.sample
}

# Для проверки корректности моделирования, построим выборку из p-value и проверим её распределение на равномерность
# 1000 раз промоделируем выборку размера 1000, на каждом шаге воспользуемся критерием Колмогорова-Смирнова, для получения p-value
load("pvalue_sample")
#recalculate <- T
if (recalculate){
  set.seed(start.seed <- 42)
  pval.sample <- numeric(0)
  for (counter in (1:1000)) {
    BS.sample <- rBirnSaund(1000)
    pval.sample <- c(pval.sample, ks.test(BS.sample, y=pBirnSaund)$p.value)
  }
  recalculate <- F
  curr.seed <- .Random.seed
  save(pval.sample, start.seed, curr.seed, recalculate, file="pvalue_sample")
}
# Опять применим критерий Колмогорова-Смирнова для выборки из p-value
# (если гипотеза верна, то p-value должно иметь распределение U(0,1))
ks.test(pval.sample, y="punif")
hist(pval.sample)


