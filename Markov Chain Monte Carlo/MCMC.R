dproposal <- function(x.new, x.curr) (dbeta(x.new, 1, 1))
rproposal <- function(n, x.curr) (rbeta(n, 1, 1))

B <- beta(0.7, 0.3)
dsimul <- function(x) {
    x^(0.7-1) * (1-x)^(0.3 - 1)  / B
}

acceptance <- function(x.new, x.curr) {
    min(1, (dproposal(x.curr, x.new) * dsimul(x.new)) / (dproposal(x.new, x.curr) * dsimul(x.curr)))
}

load("MCMC_B11")
if (recalculate){
    set.seed(start.seed <- 4399)
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
    save(sample.simul, accepted, start.seed, cur.seed, recalculate, file="MCMC_B11")
}

burnin <- 10000
ind <- burnin:N
plot((cumsum(sample.simul[ind])/seq_along(sample.simul[ind])), ylim=c(0,1))

block_size <- length(ind) / 10
block_means <- sapply(split(sample.simul[ind], ceiling(seq_along(sample.simul[ind])/block_size)), mean)
print(block_means)

ind <- burnin:N
acf(sample.simul[ind], lag.max = 500)

lag.s <- 50
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

mean(accepted[burnin:N - 1])


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
    save(sample.simul, accepted, start.seed, cur.seed, recalculate, file="MCMC_B3515")
}

plot(sample.simul, ylim=c(0,1))
plot((cumsum(sample.simul)/seq_along(sample.simul)), ylim=c(0,1))
burnin <- 20000
ind <- burnin:N
plot((cumsum(sample.simul[ind])/seq_along(sample.simul[ind])), ylim=c(0,1))

block_size <- length(ind) / 10
block_means <- sapply(split(sample.simul[ind], ceiling(seq_along(sample.simul[ind])/block_size)), mean)
print(block_means)

ind <- burnin:N
acf(sample.simul[ind], lag.max = 500)

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
    save(sample.simul, accepted, start.seed, cur.seed, recalculate, file="MCMC_B73")
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
hist(sample.final, breaks = 30, prob = TRUE, 
     main = "Выборка vs Beta(0.7,0.3)", col = "lightblue")
curve(dbeta(x, 0.7, 0.3), add = TRUE, col = "red", lwd = 2)

qqplot(qbeta(ppoints(length(sample.final)), 0.7, 0.3), 
       sample.final, 
       main = "QQ-plot против Beta(0.7,0.3)",
       xlab = "Теоретические квантили", ylab = "Эмпирические квантили")
abline(0, 1, col = "red")

mean(accepted[burnin:N - 1])


