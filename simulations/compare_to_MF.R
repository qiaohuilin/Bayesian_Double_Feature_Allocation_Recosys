rm(list=ls())
library(MCMCpack)
library(invgamma)
library(truncnorm)
library(extraDistr)
library(EnvStats)


Rcpp::sourceCpp("MCMC_cpp.cpp")

set.seed(123)
load("trueR5_lam=3_b=0.2_m=100_n=150.RData")
#load("trueR5_lam=3_b=0.5_m=50.RData")

R_obs
R_true
train_bi=1*(R_obs!=0)
test_bi=(1*(R_true!=0))*(1*(R_obs==0))
R_test=R_true*test_bi


#compare to matrix factorization with a grid number of features
#in this section, lambda is regularization, alpha is step size.
unroll_Vecs <- function (params, Y, R, num_users, num_movies, num_features) {
  # Unrolls vector into X and Theta
  # Also calculates difference between preduction and actual 
  
  endIdx <- num_movies * num_features
  
  X     <- matrix(params[1:endIdx], nrow = num_movies, ncol = num_features)
  Theta <- matrix(params[(endIdx + 1): (endIdx + (num_users * num_features))], 
                  nrow = num_users, ncol = num_features)
  
  Y_dash     <-   (((X %*% t(Theta)) - Y) * R) # Prediction error
  
  return(list(X = X, Theta = Theta, Y_dash = Y_dash))
}

J_cost <-  function(params, Y, R, num_users, num_movies, num_features, lambda, alpha) {
  # Calculates the cost
  
  unrolled <- unroll_Vecs(params, Y, R, num_users, num_movies, num_features)
  X <- unrolled$X
  Theta <- unrolled$Theta
  Y_dash <- unrolled$Y_dash
  
  J <-  .5 * sum(   Y_dash ^2)  + lambda/2 * sum(Theta^2) + lambda/2 * sum(X^2)
  
  return (J)
}

grr <- function(params, Y, R, num_users, num_movies, num_features, lambda, alpha) {
  # Calculates the gradient step
  # Here lambda is the regularization parameter
  # Alpha is the step size
  
  unrolled <- unroll_Vecs(params, Y, R, num_users, num_movies, num_features)
  X <- unrolled$X
  Theta <- unrolled$Theta
  Y_dash <- unrolled$Y_dash
  
  X_grad     <- ((   Y_dash  %*% Theta) + lambda * X     )
  Theta_grad <- (( t(Y_dash) %*% X)     + lambda * Theta )
  
  grad = c(X_grad, Theta_grad)
  return(grad)
}

num_features_list=c(10,15,20,25,30,35)
cost <- NULL
for(i in 1:4){
  num_features=num_features_list[i]
  res <- optim(par = c(runif(m * num_features), runif(n * num_features)), # Random starting parameters
               fn = J_cost, gr = grr, 
               Y=t(R_obs), R=t(train_bi), 
               num_users=m, num_movies=n,num_features=num_features, 
               lambda=1, alpha = 0.2, 
               method = "L-BFGS-B", control=list(maxit=3000, trace=1))
  cost[i]=res$value
}
cost[2:4]-cost[1:3]


num_features=15
res <- optim(par = c(runif(m * num_features), runif(n * num_features)), # Random starting parameters
             fn = J_cost, gr = grr, 
             Y=t(R_obs), R=t(train_bi), 
             num_users=m, num_movies=n,num_features=num_features, 
             lambda=1, alpha = 0.2, 
             method = "L-BFGS-B", control=list(maxit=3000, trace=1))
U=matrix(res$par[1:m * num_features],m,num_features)
V=matrix(res$par[1:n * num_features],n,num_features)

dim(U)
dim(V)
R_pred=U %*% t(V)
R_pred[R_pred>=5]=5
R_pred[R_pred<=1]=1
R_pred[R_pred<5 & R_pred>1]=round(R_pred[R_pred<5 & R_pred>1])
R_pred_test=R_pred*test_bi
R_pred_test[R_pred_test!=0] 
R_test[R_test!=0]
sum(R_pred_test[R_pred_test!=0] == R_test[R_test!=0])




if(!require(recosystem)) 
  install.packages("recosystem", repos = "http://cran.us.r-project.org")
library(recosystem)
set.seed(123, sample.kind = "Rounding") # This is a randomized algorithm

user_train=vector()
item_train=vector()
rating_train=vector()
for(i in 1:m){
  for(j in 1:n){
    if(R_obs[i,j]!=0){
      user_train=c(user_train,i)
      item_train=c(item_train,j)
      rating_train=c(rating_train,R_obs[i,j])
    }
  }
}
train_set=data.frame(cbind(user_train,item_train,rating_train))

user_test=vector()
item_test=vector()
rating_test=vector()
for(i in 1:m){
  for(j in 1:n){
    if(R_test[i,j]!=0){
      user_test=c(user_test,i)
      item_test=c(item_test,j)
      rating_test=c(rating_test,R_test[i,j])
    }
  }
}
test_set=data.frame(cbind(user_test,item_test,rating_test))


# Convert the train and test sets into recosystem input format
train_data <-  with(train_set, data_memory(user_index = user_train, 
                                           item_index = item_train, 
                                           rating     = rating_train
                                           ))
test_data  <-  with(test_set,  data_memory(user_index = user_test, 
                                           item_index = item_test, 
                                           rating     = rating_test))

# Create the model object
r <-  recosystem::Reco()

# Select the best tuning parameters
opts <- r$tune(train_data, opts = list(dim = c(10, 20, 30), 
                                       lrate = c(0.1, 0.2),
                                       costp_l2 = c(0.01, 0.1), 
                                       costq_l2 = c(0.01, 0.1),
                                       nthread  = 4, niter = 10))

# Train the algorithm  
r$train(train_data, opts = c(opts$min, nthread = 4, niter = 20))
y_hat_reco <-  r$predict(test_data, out_memory())
head(y_hat_reco, 10)

y_hat_reco_orig=y_hat_reco

###
y_hat_reco2=y_hat_reco
y_hat_reco2[y_hat_reco2>4]=5
y_hat_reco2[y_hat_reco2<=1]=1
y_hat_reco2[y_hat_reco2>1 & y_hat_reco2<=2]=2
y_hat_reco2[y_hat_reco2>2 & y_hat_reco2<=3]=3
y_hat_reco2[y_hat_reco2>3 & y_hat_reco2<=4]=4
head(y_hat_reco2, 10)
sum(y_hat_reco2==rating_test)/length(rating_test)


####
y_hat_reco[y_hat_reco>5]=5
y_hat_reco[y_hat_reco<1]=1
y_hat_reco[y_hat_reco<=5 & y_hat_reco>=1]=round(y_hat_reco[y_hat_reco<=5 & y_hat_reco>=1])

sum(y_hat_reco==rating_test)/length(rating_test)
######3



# Define Mean Absolute Error (MAE)
MAE <- function(true_ratings, predicted_ratings){
  mean(abs(true_ratings - predicted_ratings))
}

# Define Mean Squared Error (MSE)
MSE <- function(true_ratings, predicted_ratings){
  mean((true_ratings - predicted_ratings)^2)
}

# Define Root Mean Squared Error (RMSE)
RMSE <- function(true_ratings, predicted_ratings){
  sqrt(mean((true_ratings - predicted_ratings)^2))
}
RMSE(test_set$rating, y_hat_reco)
MSE(test_set$rating, y_hat_reco)
MAE(test_set$rating, y_hat_reco)


