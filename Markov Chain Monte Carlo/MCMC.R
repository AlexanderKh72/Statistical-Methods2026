# Markov chain Monte Carlo
# Вариант 4. Моделирование B(0.7, 0.3)

# Предлагающее распределение. 
# Функции записаны в общем виде, хотя в нашем случае распределение не зависит от текущего состоряния.
dproposal <- function(x.new, x.curr) (dbeta(x.new, 1, 1))
rproposal <- function(n, x.curr) (rbeta(n, 1, 1))

# Плотность моделируемого распределения. Вообще, B считать не обязательно, оно далее будет сокращаться
B <- beta(0.7, 0.3)
dsimul <- function(x) {
    x^(0.7-1) * (1-x)^(0.3 - 1)  / B
}

# Функция вероятности принятия нового состояния. По формуле
acceptance <- function(x.new, x.curr) {
    min(1, (dproposal(x.curr, x.new) * dsimul(x.new)) / (dproposal(x.new, x.curr) * dsimul(x.curr)))
}

# В силу специфики метода и "неприятности" распределения (плотность бесконечная на краях),
    # нам нужен очень большой объём выборки, чтобы оставить из неё только часть.
load("MCMC_B11")
if (recalculate){
    set.seed(start.seed <- 4399)
    # начальное сосотояние
    sample.simul <- runif(1)
    x.curr <- sample.simul
    # будем считать, сколько раз мы приняли состояние
    accepted <- numeric(0)
    # большая выборка
    N <- 100000
    for (i in 1:(N-1)) {
        # предлагаемое состояние
        x.new <- rproposal(1, x.curr)
        # уровень принятия
        alpha <- acceptance(x.new, x.curr)
        if (runif(1) > alpha) {
            # либо не принимаем, и новое состояние совпадает с текущим
            accepted <- c(accepted, 0)
            x.new  <- x.curr
        }
        else {
            # либо принимаем и обновляем текущее состояние
            accepted <- c(accepted, 1)
            x.curr  <- x.new
        }
        sample.simul <- c(sample.simul, x.new)
    }
    cur.seed <- .Random.seed
    recalculate <- F
    save(sample.simul, accepted, start.seed, cur.seed, recalculate, N, file="MCMC_B11")
}

# Для начала, надо оценить, когда цепь вышла на стационарное состояние.
# Будем оценивать это по значению среднего. При стационарности, оно должно достаточно быстро выйти на значение 0.7
plot((cumsum(sample.simul)/seq_along(sample.simul)), ylim=c(0,1))

# Попробуем обрубить первые 10,000 значений (burn-in)
burnin <- 10000
ind <- burnin:N
# Выход на прямую стал резче, значит цепь к этому моменту вышла близко к стационарному распределению
plot((cumsum(sample.simul[ind])/seq_along(sample.simul[ind])), ylim=c(0,1))

# В качестве дополнительной проверки, посмотрим среднее значение на последовательных непересекающихся участках
block_size <- length(ind) / 10
block_means <- sapply(split(sample.simul[ind], ceiling(seq_along(sample.simul[ind])/block_size)), mean)
# Значение, в целом, близко к нужному
block_means

# Теперь надо справится с автокорреляцией при помощи прореживания
# Из-за неприятности распределения нужно быть готовым прореживать достаточно часто
ind <- burnin:N
# Построим АКФ для цепи после выхода на стационарность. 
# Как видно по графику, она достаточно ярко выражена даже при больших лагах.
acf(sample.simul[ind], lag.max = 500)

# Попробуем подобрать лаг достаточно небольшой, но при этом, чтобы АКФ резко падала при увеличении лага
# За кадром остаются эксперименты, значение, на котором я остановился -- 50
lag.s <- 50
ind <- seq(burnin, N, by=lag.s)
acf(sample.simul[ind], lag.max = 50)

# Окончательно, мы отбросили первые 10% выборки и проредили с шагом 50. Оставшийся объём выборки -- 1801
# Проверим её соответствие распределению B(0.7, 0.3)
ind <- seq(burnin, N, by=lag.s)
sample.final <- sample.simul[ind]
# Гистограмма на фоне плотности очень хорошо совпадают
hist(sample.final, breaks = 30, prob = TRUE, 
     main = "Выборка vs Beta(0.7,0.3)", col = "lightblue")
curve(dbeta(x, 0.7, 0.3), add = TRUE, col = "red", lwd = 2)

# Построим график qq. Если мы нигде не ошиблись, точки должны лежать на прямой y=x.
qqplot(qbeta(ppoints(length(sample.final)), 0.7, 0.3), 
       sample.final, 
       main = "QQ-plot против Beta(0.7,0.3)",
       xlab = "Теоретические квантили", ylab = "Эмпирические квантили")
# Результат радует мой глаз
abline(0, 1, col = "red")

# В задании требуется оценить вероятность принятия, при разных дисперсиях
mean(accepted[burnin:N - 1])

# Попробуем другие распределения с большей дисперсией.
# Выбор параметров не очень логичный с практической точки зрения,
    # но при попытке выбрать те параметры, при которых моделирование B(alpha, beta) было бы простой задачей
    # (например, B(1, n)), даже просто среднее значение не успевало выйти на нужный уровень (оставалось около 0.6, например).
    # Видимо, это связно с формой распределения B(0.7, 0.3), что выйти на него занимает очень много шагов.
# Распределение выбиралось так, чтобы хотя бы среднее proposal distribution совпадало с предельным
dproposal <- function(x.new, x.curr) (dbeta(x.new, 3.5, 1.5))
rproposal <- function(n, x.curr) (rbeta(n, 3.5, 1.5))

load("MCMC_B3515")
# recalculate <- T
if (recalculate){
    set.seed(start.seed <- cur.seed)
    sample.simul <- runif(1)
    x.curr <- sample.simul
    accepted <- numeric(0)
    N <- 100000
    for (i in 1:(N-1)) {
        x.new <- rproposal(1, x.curr)
        alpha <- acceptance(x.new, x.curr)
        if (runif(1) > alpha) {
            accepted <- c(accepted, 0)
            x.new  <- x.curr
        }
        else {
            accepted <- c(accepted, 1)
            x.curr  <- x.new
        }
        sample.simul <- c(sample.simul, x.new)
    }
    cur.seed <- .Random.seed
    recalculate <- F
    save(sample.simul, accepted, start.seed, cur.seed, recalculate, N, file="MCMC_B3515")
}

#plot(sample.simul, ylim=c(0,1))
plot((cumsum(sample.simul)/seq_along(sample.simul)), ylim=c(0,1))
burnin <- 20000
ind <- burnin:N
plot((cumsum(sample.simul[ind])/seq_along(sample.simul[ind])), ylim=c(0,1))

block_size <- length(ind) / 10
block_means <- sapply(split(sample.simul[ind], ceiling(seq_along(sample.simul[ind])/block_size)), mean)
print(block_means)

ind <- burnin:N
acf(sample.simul[ind], lag.max = 500)

# Лаг сильно больше
lag.s <- 200
ind <- seq(burnin, N, by=lag.s)
acf(sample.simul[ind], lag.max = 50)

ind <- seq(burnin, N, by=lag.s)
sample.final <- sample.simul[ind]
hist(sample.final, breaks = 30, prob = TRUE, 
     main = "Выборка vs Beta(0.7,0.3)", col = "lightblue")
curve(dbeta(x, 0.7, 0.3), add = TRUE, col = "red", lwd = 2)

qqplot(qbeta(ppoints(length(sample.final)), 0.7, 0.3), 
       sample.final, 
       main = "QQ-plot против Beta(0.7,0.3)",
       xlab = "Теоретические квантили", ylab = "Эмпирические квантили")
abline(0, 1, col = "red")

# Вероятность принятия уменьшилась. 
# Оно и понятно - proposal ditribution имеет маленькую дисперсиcию (0.035) и собрано у среднего значения 0.7
# А искомое - имеет большую дисперсию (0.105) и концентрируется на краях
mean(accepted[burnin:N - 1])


dproposal <- function(x.new, x.curr) (dbeta(x.new, 7, 3))
rproposal <- function(n, x.curr) (rbeta(n, 7, 3))

load("MCMC_B73")
# recalculate <- T
if (recalculate){
    set.seed(start.seed <- cur.seed)
    sample.simul <- runif(1)
    x.curr <- sample.simul
    accepted <- numeric(0)
    N <- 100000
    for (i in 1:(N-1)) {
        x.new <- rproposal(1, x.curr)
        alpha <- acceptance(x.new, x.curr)
        if (runif(1) > alpha) {
            accepted <- c(accepted, 0)
            x.new  <- x.curr
        }
        else {
            accepted <- c(accepted, 1)
            x.curr  <- x.new
        }
        sample.simul <- c(sample.simul, x.new)
    }
    cur.seed <- .Random.seed
    recalculate <- F
    save(sample.simul, accepted, start.seed, cur.seed, recalculate, N, file="MCMC_B73")
}

plot(sample.simul, ylim=c(0,1))
plot((cumsum(sample.simul)/seq_along(sample.simul)), ylim=c(0,1))
burnin <- 15000
ind <- burnin:N
plot((cumsum(sample.simul[ind])/seq_along(sample.simul[ind])), ylim=c(0,1))

block_size <- length(ind) / 10
block_means <- sapply(split(sample.simul[ind], ceiling(seq_along(sample.simul[ind])/block_size)), mean)
print(block_means)

ind <- burnin:N
acf(sample.simul[ind], lag.max = 500)

lag.s <- 1000
ind <- seq(burnin, N, by=lag.s)
acf(sample.simul[ind], lag.max = 50)

ind <- seq(burnin, N, by=lag.s)
sample.final <- sample.simul[ind]
hist(sample.final, breaks = 10, prob = TRUE, 
     main = "Выборка vs Beta(0.7,0.3)", col = "lightblue")
curve(dbeta(x, 0.7, 0.3), add = TRUE, col = "red", lwd = 2)

qqplot(qbeta(ppoints(length(sample.final)), 0.7, 0.3), 
       sample.final, 
       main = "QQ-plot против Beta(0.7,0.3)",
       xlab = "Теоретические квантили", ylab = "Эмпирические квантили")
abline(0, 1, col = "red")

mean(accepted[burnin:N - 1])


