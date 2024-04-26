#BERRAHO Aya & MARGHT Abde-rrahman

rm(list=ls())
library(ggplot2)

#Suppositions de l'énoncé
r <- 0
S0 <- 1
t0 <- 0
T <- 1
K <- 1
sigma <-0.25

#Constantes a fixer
n <- 1000
N <- 1000

#Question 1#

#Grille de discrétisation  
t <- seq(0,T,T/n)
delta <- T/n

#Simulation du mouvement brownien 
mouvement_brownien <- function(n){
  W <- numeric(n+1)
  W[1] <- 0
  for (i in 1:n) {
    Z <- rnorm(1)
    W[i+1] <- W[i] + sqrt(t[i+1] - t[i]) * Z
  }
  return(W)
}

W <- mouvement_brownien(n)

# Tracé du graphique
df <- data.frame(t = t, W = W)
ggplot(df, aes(x = t, y = W)) +
  geom_line() +
  labs(x = "t", y = "W(t)", title = "Graphique du processus W")

prix_actif <- function(n){
  S <- numeric(n)
  W <- mouvement_brownien(n)
  S <- S0*exp((r-sigma^2/2)*t+sigma*W)
  return(S)
}

S <- prix_actif(n)

# Tracé du graphique
df <- data.frame(t = t, S = S)
ggplot(df, aes(x = t, y = S)) +
  geom_line() +
  labs(x = "t", y = "S(t)", title = "Graphique du prix S")

#Calcul du prix E[G]
phi <- function(X){
  res <- X-K
  res[res<0] <- 0
  return (res)
}

A_n<- function(n,N){ #Simule un vecteur de taille N idd A-n^T
  A <- numeric(N)
  for (i in 1:N){
    S <- prix_actif(n)
    A[i] <- mean(S)
  }
  return (A)
}

MC <- function(n,N){
  A <- A_n(n,N)
  h <- phi(A)
  return(c(mean(h),var(h)/N))
}
prix_MC <- MC(n,N)

#2

S_V2 <- function(n,W){
  S <- numeric(n+1)
  for(i in 1:n+1){
    S[i] <- S0*exp( (r-(sigma^2)/2)*t[i] + sigma * W[i])
  }
  return(S)
}
A_n_V2 <- function(n,W){ #Simule un vecteur de taille 1 de loi A-n^T
  S <- S_V2(n,W)
  return (mean(S))
}

MC_ANT <- function(n,N){
  A <- numeric(N)
  B <- numeric(N)
  for(i in 1:N){
    W <- mouvement_brownien(n)
    A[i] <- A_n_V2(n,W) - K
    B[i] <- A_n_V2(n,(-1)*W) - K
  }
  A[A<0] <- 0
  B[B<0] <- 0
  h <- 0.5*(A+B)
  return(c(mean(h),var(h)/(4*N)))
}
  
prix_MC_ANTI <- MC_ANT(n,N)

#3
l <- (sqrt(6*(n+1))*sigma)/(4*n*sqrt(2*n+1))
m <- sigma/n * sqrt((n+1)*(2*n+1)/6)
F2 <- pnorm(l)
F1 <- pnorm(l,m,1)
mu <- exp((m**2)/2 - (sigma**2)*(n+1)/(4*n))*(1-F1)-(1-F2)


b_optimal <- function(m){
  b <- numeric(m)
  A <- numeric(m)
  Y <- numeric(m)
  for (i in 1:m){
    S <- prix_actif(n)
    A[i] <- phi(mean(S))
    Y[i] <- max(0,(prod(S)^(1/n))-K)
  }
  for(i in 1:m){
    num <- sum( (A[1:i]-mean(A[1:i]))*(Y[1:i]-mean(Y[1:i])))
    denom <- sum((Y[1:i]-mean(Y[1:i]))**2)
    b[i] <- -num/denom
  }
  return(b)
}

m <- 200
b <- b_optimal(m)
plot(b) #on prend donc n0 = 50

n0 <- 100
b_n0 <- b[n0]

MC_C <- function(n,N){
  A <- numeric(N)
  Y <- numeric(N)
  for (i in 1:N){
    S <- prix_actif(n)
    A[i] <- mean(S)
    Y[i] <- max(0,(prod(S)^(1/n))-K)
  }
  h <- phi(A)+b_n0*(Y-mu)
  return(c(mean(h), var(h)/N))
}
prix_MC_C <- MC_C(n,N)


#4
#a

Z<- rnorm(n)

g <- function(Z){
  res <- max(0,g_tilde(Z)-K)
  return(res)
}

g_tilde <- function(Z){
  G <- numeric(n)
  for(i in 1:n){
    G[i] <- exp( (r-(sigma**2)/2)*t[i] + sigma*sum(Z[1:i])/sqrt(n))
  }
  res <- S0* mean(G)
  return(res)
}




#4
#c

vecteur_z_S <- function(y,n){
  z <- numeric(n)
  S <- numeric(n)
  z[1] <- sigma*sqrt(delta)*(y+K)/y
  for(j in 2:n){
    S[j] <- S0*exp(-(sigma**2)*t[j]/2 + sigma*sqrt(delta)*sum(z[1:j]))
    z[j] <- z[j-1] - sigma*sqrt(delta)*S[j-1]/(n*y)
  }
  return(c(z,S))
}

zS <- vecteur_z_S(2,n)
z <- zS[1:n]
S <- zS[n:length(zS)]
S

abs <- seq(0.1550196035,0.1550196037,length.out=100)
res <- numeric(length(abs))
for(i in 1:length(abs)){
  res[i] <- mean(vecteur_z_S(abs[i],N)[n:length(zS)])-K-abs[i]
}


df <- data.frame(x = abs, y = res)
ggplot(df, aes(x = abs, y = res)) +
  geom_line() +
  labs(x = "y", y = "res", title = "Graphique:")

# Avec la methode de bifurcation, je trouve y_hat = 0.1550196035
y_hat <- 0.1550196035
mu_hat <- vecteur_z_S(y_hat,n)[1:n]


MC_EP <- function(n,N){
  h<-numeric(N)
  for(i in 1:N){
    Z <- rnorm(n)
    h[i] <- g(Z+mu_hat)*exp(- mu_hat%*%Z - 0.5 *mu_hat%*%mu_hat)
  }
  return(c(mean(h),var(h)/N))
}
prix_MC_EP <- MC_EP(n,N)

