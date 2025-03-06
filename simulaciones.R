# SIMULACIÓN
library(VGAM)


valores = seq(0.0001,1,0.0001)

m = function(x,a,b){
  return( 6 * dkumar(x,a,b) * pkumar(x,a,b)  * (1- pkumar(x,a,b) ) )
}

M = function(x,a,b){
  return( 3 * pkumar(x,a,b)^2 - 2 * pkumar(x,a,b)^3 )
}


# definir forma de la densidad

f3 = function(x,q,alpha,a,b){
  
  t1 = - alpha * M(x,a,b)^(alpha - 1) * m(x,a,b) * (1-inverse_logit(q))
  t2 = log(inverse_logit(q)) * ( 1 - (1-inverse_logit(q))*M(x,a,b)^(alpha) )
  return(t1/t2)
}

# definir la acumulada

acum = function(x, q,alpha,a,b){
  t1 = log(1 - (1-q)*M(x,a,b)^(alpha) )
  t2 = log(q)
  return(t1/t2)
}


acum3 = function(x, q,alpha,a,b){
  t1 = log(1 - (1- inverse_logit(q)  )*M(x,a,b)^(alpha) )
  t2 = log(inverse_logit(q))
  return(t1/t2)
}

# funcion cuantil


quantile_l = function(u,q,alpha,a,b){
  c = ( (1-q^u)/(1-q) )^(1/alpha)
  z = 0.5 - sin( asin(1-2*c)/3  )
  
  valor = ( 1 - (1-z)^(1/b)  )^(1/a)
  return(valor)
}

# logits

logit_transform <- function(lambda) {
  if (lambda <= 0 || lambda >= 1) {
    stop("Lambda must be between 0 and 1.")
  }
  return(log(lambda / (1 - lambda)))  
}


inverse_logit <- function(logit_lambda) {
  return(1 / (1 + exp(-logit_lambda)))  
}






###----- definir funciones -----###



log_like = function(parametros, data) {
  q <- parametros[1]
  alpha <- parametros[2]
  a <- parametros[3]
  b <- parametros[4]
  
  -sum(log(f3(data, q, alpha, a, b)))
}





## prueba

f = function(x,q,alpha,a,b){
  t1 = - alpha * M(x,a,b)^(alpha - 1) * m(x,a,b) * (1-q)
  t2 = log(q) * ( 1 - (1-q)*M(x,a,b)^(alpha) )
  return(t1/t2)
}


log_like2 <- function(parametros, data) {
  q <- parametros[1]
  alpha <- parametros[2]
  a <- parametros[3]
  b <- parametros[4]
  
  -sum(log(f1(data, q, alpha, a, b)))
}


f1 = function(x,q,alpha,a,b){
  
  t1 = - exp(alpha) * M(x,a,b)^(exp(alpha) - 1) * m(x,a,b) * (1-inverse_logit(q))
  t2 = log(inverse_logit(q)) * ( 1 - (1-inverse_logit(q))*M(x,a,b)^( exp(alpha)) )
  return(t1/t2)
}




muestra1 = quantile_l(runif(50),0.5,1,1,1)
hist(muestra1, freq = F)
points(valores,f(valores, 0.5,1,1,1), type = "l")


initial_params <- c(inverse_logit(0.2), 0.85, 0.85, 0.85) 
result <- optim(par = initial_params, fn = log_like2, data = muestra1)
c(inverse_logit(result$par[1]), exp(result$par[2]), result$par[3],result$par[4] )
inverse_logit(result$par[1])
exp(result$par[2])

hist(muestra1, freq = F)
points(valores,f3(valores, result$par[1] , result$par[2], result$par[3],result$par[4] ), type = "l", col = "blue")
points(valores,f1(valores, inverse_logit(result$par[1]), exp(result$par[2]), result$par[3],result$par[4] ), type = "l", col = "blue")

c(inverse_logit(result$par[1]), exp(result$par[2]), result$par[3],result$par[4] )

log_like(c(inverse_logit(result$par[1]), exp(result$par[2]), result$par[3],result$par[4] ),muestra1)

hist(data$FIRMCOST, freq = F, breaks = 20)
result <- optim(par = c(5, 1.54, 0.5, 0.5) , fn = log_like, data = data$FIRMCOST)
result$par
points(valores,f3(valores, inverse_logit(result$par[1]),result$par[2], result$par[3],result$par[4] ), type = "l", col = "blue")
log_like(c(inverse_logit(result$par[1]),result$par[2], result$par[3],result$par[4] ), data$FIRMCOST) 
ks.test(
  data$FIRMCOST, 
  function(x) acum3(data$FIRMCOST, inverse_logit(result$par[1]),result$par[2], result$par[3],result$par[4] )
)



# Escenario 1



# 0.5,1,1,1



# n = 50 ------------------------------------------------------------------





set.seed(50)

N = 10000

matrix_sim11 = matrix(NA_real_,N,4 )

for(i in 1:N){
  
  muestra = quantile_l( runif(50),  0.5,1,1,1) 
  initial_params <- c(inverse_logit(-0.2),0.5,0.5,0.5) 
  estimados <- optim(par = initial_params, fn = log_like, data = muestra)
  
  matrix_sim11[i,1] = inverse_logit(estimados$par[1])
  matrix_sim11[i,2] = estimados$par[2]
  matrix_sim11[i,3] = estimados$par[3]
  matrix_sim11[i,4] = estimados$par[4]
  
  
}



colMeans(matrix_sim11)

colMeans(matrix_sim11) - c(0.5,1,1,1)

sqrt(c(mean( (matrix_sim11[,1]  - 0.5 )^2 ),
  mean( (matrix_sim11[,2]  - 1 )^2 ),
  mean( (matrix_sim11[,3]  - 1 )^2 ),
  mean( (matrix_sim11[,4]  - 1 )^2 )))




# n = 100 ------------------------------------------------------------------





set.seed(100)

N = 10000

matrix_sim21 = matrix(NA_real_,N,4 )

for(i in 1:N){
  
  muestra = quantile_l( runif(100), 0.5,1,1,1) 
  initial_params <- c(inverse_logit(-0.2),0.5,0.5,0.5) 
  estimados <- optim(par = initial_params, fn = log_like, data = muestra)
  
  matrix_sim21[i,1] = inverse_logit(estimados$par[1])
  matrix_sim21[i,2] = estimados$par[2]
  matrix_sim21[i,3] = estimados$par[3]
  matrix_sim21[i,4] = estimados$par[4]
  
  
}


colMeans(matrix_sim21)

colMeans(matrix_sim21) - c(0.5,1,1,1)

sqrt(c(mean( (matrix_sim21[,1]  - 0.5 )^2 ),
  mean( (matrix_sim21[,2]  - 1 )^2 ),
  mean( (matrix_sim21[,3]  - 1 )^2 ),
  mean( (matrix_sim21[,4]  - 1 )^2 )))


# n = 200 ------------------------------------------------------------------





set.seed(200)

N = 10000

matrix_sim31 = matrix(NA_real_,N,4 )

for(i in 1:N){
  
  muestra = quantile_l( runif(200), 0.5,1,1,1) 
  initial_params <- c(inverse_logit(-0.2),0.5,0.5,0.5) 
  estimados <- optim(par = initial_params, fn = log_like, data = muestra)
  
  matrix_sim31[i,1] = inverse_logit(estimados$par[1])
  matrix_sim31[i,2] = estimados$par[2]
  matrix_sim31[i,3] = estimados$par[3]
  matrix_sim31[i,4] = estimados$par[4]
  
  
}


colMeans(matrix_sim31)

colMeans(matrix_sim31) - c(0.5,1,1,1)

sqrt(c(mean( (matrix_sim31[,1]  - 0.5 )^2 ),
  mean( (matrix_sim31[,2]  - 1 )^2 ),
  mean( (matrix_sim31[,3]  - 1 )^2 ),
  mean( (matrix_sim31[,4]  - 1 )^2 )))












# n = 500 ------------------------------------------------------------------





set.seed(500)

N = 10000

matrix_sim41 = matrix(NA_real_,N,4 )

for(i in 1:N){
  
  muestra = quantile_l( runif(500), 0.5,1,1,1) 
  initial_params <- c(inverse_logit(-0.2),0.5,0.5,0.5)
  estimados <- optim(par = initial_params, fn = log_like, data = muestra)
  
  matrix_sim41[i,1] = inverse_logit(estimados$par[1])
  matrix_sim41[i,2] = estimados$par[2]
  matrix_sim41[i,3] = estimados$par[3]
  matrix_sim41[i,4] = estimados$par[4]
  
  
}



colMeans(matrix_sim41)

colMeans(matrix_sim41) - c(0.5,1,1,1)

sqrt(c(mean( (matrix_sim41[,1]  - 0.5 )^2 ),
  mean( (matrix_sim41[,2]  - 1 )^2 ),
  mean( (matrix_sim41[,3]  - 1 )^2 ),
  mean( (matrix_sim41[,4]  - 1 )^2 )))







# Escenario 2



# 0.1, 2,  1.5, 2



# n = 50  - 2 -------------------------------------------------------------


set.seed(502)

N = 10000

matrix_sim12 = matrix(NA_real_,N,4 )

for(i in 1:N){
  
  muestra = quantile_l( runif(50), 0.1, 2,  1.5, 2) 
  initial_params <- c(inverse_logit(-2.8), 1, 0.5, 1.5) 
  estimados <- optim(par = initial_params, fn = log_like, data = muestra)
  
  matrix_sim12[i,1] = inverse_logit(estimados$par[1])
  matrix_sim12[i,2] = estimados$par[2]
  matrix_sim12[i,3] = estimados$par[3]
  matrix_sim12[i,4] = estimados$par[4]
  
  
}




colMeans(matrix_sim12)

colMeans(matrix_sim12) - c(0.1, 2,  1.5, 2)

sqrt(c(mean( (matrix_sim12[,1]  - 0.1 )^2 ),
  mean( (matrix_sim12[,2]  - 2 )^2 ),
  mean( (matrix_sim12[,3]  - 1.5 )^2 ),
  mean( (matrix_sim12[,4]  - 2 )^2 )))



# n = 100 -2 ------------------------------------------------------------------





set.seed(1002)

N = 10000

matrix_sim22 = matrix(NA_real_,N,4 )

for(i in 1:N){
  
  muestra = quantile_l( runif(100),   0.1, 2,  1.5, 2) 
  initial_params <- c(inverse_logit(-2.8), 1, 0.5, 1.5) 
  estimados <- optim(par = initial_params, fn = log_like, data = muestra)
  
  matrix_sim22[i,1] = inverse_logit(estimados$par[1])
  matrix_sim22[i,2] = estimados$par[2]
  matrix_sim22[i,3] = estimados$par[3]
  matrix_sim22[i,4] = estimados$par[4]
  
  
}





colMeans(matrix_sim22)

colMeans(matrix_sim22) - c(0.1, 2,  1.5, 2)

sqrt(c(mean( (matrix_sim22[,1]  - 0.1 )^2 ),
  mean( (matrix_sim22[,2]  - 2 )^2 ),
  mean( (matrix_sim22[,3]  - 1.5 )^2 ),
  mean( (matrix_sim22[,4]  - 2 )^2 )))



# n = 200 - 2 ------------------------------------------------------------------





set.seed(2002)

N = 10000

matrix_sim32 = matrix(NA_real_,N,4 )

for(i in 1:N){
  
  muestra = quantile_l( runif(200),  0.1, 2,  1.5, 2) 
  initial_params <- c(inverse_logit(-2.8), 1, 0.5, 1.5) 
  estimados <- optim(par = initial_params, fn = log_like, data = muestra)
  
  matrix_sim32[i,1] = inverse_logit(estimados$par[1])
  matrix_sim32[i,2] = estimados$par[2]
  matrix_sim32[i,3] = estimados$par[3]
  matrix_sim32[i,4] = estimados$par[4]
  
  
}




colMeans(matrix_sim32)

colMeans(matrix_sim32) - c(0.1, 2,  1.5, 2)

sqrt(c(mean( (matrix_sim32[,1]  - 0.1 )^2 ),
  mean( (matrix_sim32[,2]  - 2 )^2 ),
  mean( (matrix_sim32[,3]  - 1.5 )^2 ),
  mean( (matrix_sim32[,4]  - 2 )^2 )))





# n = 500 - 2 ------------------------------------------------------------------





set.seed(5002)

N = 10000

matrix_sim42 = matrix(NA_real_,N,4 )

for(i in 1:N){
  
  muestra = quantile_l( runif(500),  0.1, 2,  1.5, 2) 
  initial_params <- c(inverse_logit(-2.8), 1, 0.5, 1.5) 
  estimados <- optim(par = initial_params, fn = log_like, data = muestra)
  
  matrix_sim42[i,1] = inverse_logit(estimados$par[1])
  matrix_sim42[i,2] = estimados$par[2]
  matrix_sim42[i,3] = estimados$par[3]
  matrix_sim42[i,4] = estimados$par[4]
  
  
}


colMeans(matrix_sim42)

colMeans(matrix_sim42) - c(0.1, 2,  1.5, 2)

sqrt(c(mean( (matrix_sim42[,1]  - 0.1 )^2 ),
       mean( (matrix_sim42[,2]  - 2 )^2 ),
       mean( (matrix_sim42[,3]  - 1.5 )^2 ),
       mean( (matrix_sim42[,4]  - 2 )^2 )))





# Escenario 3



# 0.8, 1,  3, 3.5



# n = 50  - 3 -------------------------------------------------------------


set.seed(503)

N = 10000

matrix_sim13 = matrix(NA_real_,N,4 )

for(i in 1:N){
  
  muestra = quantile_l( runif(50), 0.8, 1,  3, 3.5) 
  initial_params <- c(inverse_logit(1.2), 0.5, 2.5, 1.5) 
  estimados <- optim(par = initial_params, fn = log_like, data = muestra)
  
  matrix_sim13[i,1] = inverse_logit(estimados$par[1])
  matrix_sim13[i,2] = estimados$par[2]
  matrix_sim13[i,3] = estimados$par[3]
  matrix_sim13[i,4] = estimados$par[4]
  
  
}



colMeans(matrix_sim13)

colMeans(matrix_sim13) - c(0.8, 1,  3, 3.5)

sqrt(c(mean( (matrix_sim13[,1]  - 0.8 )^2 ),
       mean( (matrix_sim13[,2]  - 1 )^2 ),
       mean( (matrix_sim13[,3]  - 3 )^2 ),
       mean( (matrix_sim42[,4]  - 3.5 )^2 )))






# n = 100 -3 ------------------------------------------------------------------





set.seed(1003)

N = 10000

matrix_sim23 = matrix(NA_real_,N,4 )

for(i in 1:N){
  
  muestra = quantile_l( runif(100), 0.8, 1,  3, 3.5) 
  initial_params <- c(inverse_logit(1.2), 0.5, 2.5, 1.5) 
  estimados <- optim(par = initial_params, fn = log_like, data = muestra)
  
  matrix_sim23[i,1] = inverse_logit(estimados$par[1])
  matrix_sim23[i,2] = estimados$par[2]
  matrix_sim23[i,3] = estimados$par[3]
  matrix_sim23[i,4] = estimados$par[4]
  
  
}

colMeans(matrix_sim23)

colMeans(matrix_sim23) - c(0.8, 1,  3, 3.5)

sqrt(c(mean( (matrix_sim23[,1]  - 0.8 )^2 ),
       mean( (matrix_sim23[,2]  - 1 )^2 ),
       mean( (matrix_sim23[,3]  - 3 )^2 ),
       mean( (matrix_sim23[,4]  - 3.5 )^2 )))



# n = 200 - 2 ------------------------------------------------------------------





set.seed(2002)

N = 10000

matrix_sim33= matrix(NA_real_,N,4 )

for(i in 1:N){
  
  muestra = quantile_l( runif(200), 0.8, 1,  3, 3.5) 
  initial_params <- c(inverse_logit(1.2), 1, 0.5, 1.5) 
  estimados <- optim(par = initial_params, fn = log_like, data = muestra)
  
  matrix_sim33[i,1] = inverse_logit(estimados$par[1])
  matrix_sim33[i,2] = estimados$par[2]
  matrix_sim33[i,3] = estimados$par[3]
  matrix_sim33[i,4] = estimados$par[4]
  
  
}


colMeans(matrix_sim33)

colMeans(matrix_sim33) - c(0.8, 1,  3, 3.5)

sqrt(c(mean( (matrix_sim33[,1]  - 0.8 )^2 ),
       mean( (matrix_sim33[,2]  - 1 )^2 ),
       mean( (matrix_sim33[,3]  - 3 )^2 ),
       mean( (matrix_sim33[,4]  - 3.5 )^2 )))






# n = 500 - 2 ------------------------------------------------------------------





set.seed(5002)

N = 10000

matrix_sim43 = matrix(NA_real_,N,4 )

for(i in 1:N){
  
  muestra = quantile_l( runif(500), 0.8, 1,  3, 3.5) 
  initial_params <- c(inverse_logit(1.2), 1, 0.5, 1.5) 
  estimados <- optim(par = initial_params, fn = log_like, data = muestra)
  
  matrix_sim43[i,1] = inverse_logit(estimados$par[1])
  matrix_sim43[i,2] = estimados$par[2]
  matrix_sim43[i,3] = estimados$par[3]
  matrix_sim43[i,4] = estimados$par[4]
  
  
}

colMeans(matrix_sim43)

colMeans(matrix_sim43) - c(0.8, 1,  3, 3.5)

sqrt(c(mean( (matrix_sim43[,1]  - 0.8 )^2 ),
       mean( (matrix_sim43[,2]  - 1 )^2 ),
       mean( (matrix_sim43[,3]  - 3 )^2 ),
       mean( (matrix_sim43[,4]  - 3.5 )^2 )))



