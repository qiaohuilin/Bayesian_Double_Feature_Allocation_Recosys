set.seed(123)

load("sim5_lam=3_b=0.5_m=100_n=150_withmu_rcpp.RData")

user_unseen_pred_acc=vector()
user_unseen_top_movie=vector()

for(i in 1:100){
# pred accuracy hit
unseen_true=R_true[i,R_obs[i,]==0]
unseen_pred=R_pred_vote[i,R_obs[i,]==0]
user_unseen_pred_acc[i]=sum(unseen_true==unseen_pred)/length(unseen_true)

# top 10 pred movies hit 
unseen_mu_true=mu_true[i,R_obs[i,]==0]
unseen_mu_pred=mu_pred_mat[i,R_obs[i,]==0]

ix=sort(unseen_mu_true,decreasing = T,index.return=T)$ix
ix_t=ix[1:10]
ix=sort(unseen_mu_pred,decreasing = T,index.return=T)$ix
ix_p=ix[1:10]
user_unseen_top_movie[i]= length(intersect(ix_t,ix_p))/10
}


# compare to MF
train_bi=1*(R_obs!=0)
test_bi=(1*(R_true!=0))*(1*(R_obs==0))
R_test=R_true*test_bi

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


####
y_hat_reco[y_hat_reco>5]=5
y_hat_reco[y_hat_reco<1]=1
y_hat_reco[y_hat_reco<=5 & y_hat_reco>=1]=round(y_hat_reco[y_hat_reco<=5 & y_hat_reco>=1])

sum(y_hat_reco==rating_test)/length(rating_test)


#### for the format
y_hat_reco_orig_mat=matrix(0,100,150)
y_hat_reco_mat=matrix(0,100,150)

tt=1

for(i in 1:100){
  for(j in 1:150){
    if(R_obs[i,j]==0){
      y_hat_reco_mat[i,j]=y_hat_reco[tt]
      y_hat_reco_orig_mat[i,j]=y_hat_reco_orig[tt]
      tt=tt+1
    }
  }
}


user_unseen_pred_acc_mf=vector()
user_unseen_top_movie_mf=vector()

for(i in 1:100){
  # pred accuracy hit
  unseen_true=R_true[i,R_obs[i,]==0]
  unseen_pred=y_hat_reco_mat[i,R_obs[i,]==0]
  user_unseen_pred_acc_mf[i]=sum(unseen_true==unseen_pred)/length(unseen_true)
  
  # top 10 pred movies hit 
  unseen_mu_true=mu_true[i,R_obs[i,]==0]
  unseen_mu_pred=y_hat_reco_orig_mat[i,R_obs[i,]==0]
  
  ix=sort(unseen_mu_true,decreasing = T,index.return=T)$ix
  ix_t=ix[1:10]
  ix=sort(unseen_mu_pred,decreasing = T,index.return=T)$ix
  ix_p=ix[1:10]
  user_unseen_top_movie_mf[i]= length(intersect(ix_t,ix_p))/10
}

boxplot(user_unseen_top_movie,user_unseen_top_movie_mf)

library(ggplot2)
library(dplyr)

top_movie=data.frame(cbind(c(user_unseen_top_movie,user_unseen_top_movie_mf),c(rep("BDFA",100),rep("MF",100))))
colnames(top_movie)=c("Accuracy","Methods")
top_movie$Accuracy %>% as.numeric() -> top_movie$Accuracy

ggplot(data=top_movie)+
  geom_boxplot(mapping = aes(x=Methods,y=Accuracy))+
  ylab("Top 10 Unseen Movie Accurary")+
  theme(axis.title.y = element_text(size = rel(1)))

user_rate=data.frame(cbind(c(user_unseen_pred_acc,user_unseen_pred_acc_mf),c(rep("BDFA",100),rep("MF",100))))
colnames(user_rate)=c("Accuracy","Methods")
user_rate$Accuracy %>% as.numeric() -> user_rate$Accuracy

ggplot(data=user_rate)+
  geom_boxplot(mapping = aes(x=Methods,y=Accuracy))+
  ylab("Unseen Rating(1-5) Accuracy")+
  theme(axis.title.y = element_text(size = rel(1)))


k_store2=vector()
for(i in 1:length(k_store)){
  k_store2[i]=k_store[[i]]
}
table(k_store2)

#K_map=18
K_map=43
idx=which(k_store2==K_map)

A_store_sort=list()
for(i in idx){
  Ai=A_store[[i]]
  Aisub=Ai[,2:19]
  cix=sort(colSums(Aisub),decreasing = T,index.return=T)$ix
  A_sort=cbind(A[,1],Aisub[,cix])
  A_store_sort[[i]]=A_sort
}

h_dist=vector()

check_perm_dist <- function(Ai,Aj){
  Aisub=Ai[,2:19]
  Ajsub=Aj[,2:19]
  j_visited=vector()
  Ajsub_perm=matrix(0,m,18)
  hd=vector()
  for(i in 1:18){
    d=vector()
    for(j in 1:dim(Ajsub)[2]){
      if(j %in% j_visited){
        d=c(d,1000)
        }
      else{
      d=c(d,sum(Aisub[,i]!=Ajsub[,j]))
      }
    }
    #print(which.min(d))
    Ajsub_perm[,i]=Ajsub[,which.min(d)]
    j_visited=c(j_visited,which.min(d))
    #hd=c(hd,sum(Ajsub_perm[,i]!=Aisub[,i]))
    hd=c(hd,min(d))
  }
  return(sum(hd)/18)
}

#for(i in 1:18){
#  print(sum(Ajsub_perm[,i]!=Aisub[,i]))
#}

h_dist=vector()
for(tt in 1:length(idx)){
  i=idx[tt]
  print(i)
  Ai=A_store_sort[[i]]
  di=vector()
  for(j in idx){
    if(i != j){
      print(j)
      Aj=A_store_sort[[j]]
      d=check_perm_dist(Ai,Aj)
      di=c(di,d)
    }
  }
  h_dist[tt]=sum(di)/length(idx)
}

which.min(h_dist)
idx[which.min(h_dist)]

idx_star=idx[which.min(h_dist)]

A_star=A_store[[idx_star]]
B_star=B_store[[idx_star]]
theta_star=theta_store[[idx_star]]

######## for simulation
Aisub=A_true[,2:19]
cix=sort(colSums(Aisub),decreasing = T,index.return=T)$ix
Aisub=Aisub[,cix]

Ajsub=A_star[,2:19]
Bjsub=B_star[,2:19]
j_visited=vector()
Astar_perm=matrix(0,m,18)
Bstar_perm=matrix(0,n,18)
hd=vector()
for(i in 1:18){
  d=vector()
  for(j in 1:dim(Ajsub)[2]){
    if(j %in% j_visited){
      d=c(d,1000)
    }
    else{
      d=c(d,sum(Aisub[,i]!=Ajsub[,j]))
    }
  }
  #print(which.min(d))
  Astar_perm[,i]=Ajsub[,which.min(d)]
  Bstar_perm[,i]=Bjsub[,which.min(d)]
  j_visited=c(j_visited,which.min(d))
  #hd=c(hd,sum(Ajsub_perm[,i]!=Aisub[,i]))
  hd=c(hd,min(d))
}

Ajsub=A_star[,2:19]
Bjsub=B_star[,2:19]
cix=sort(colSums(Ajsub),decreasing = T,index.return=T)$ix
Astar_perm=Ajsub[,cix]
Bstar_perm=Bjsub[,cix]


testplotx<- matrix(sample(0:1, 25, replace = TRUE), nrow = 5)
image(t(testplotx), col = c("white", "blue"), main = "base image plot")


image(t(Aisub), col = c("white", "blue"), main = "base image plot")
image(t(Astar_perm), col = c("white", "blue"))

library(ggplot2)
plotDat <- reshape::melt(t(Aisub))

#and plot

ggplot(plotDat, aes(X1, X2, fill = c("white", "blue")[value + 1 ])) +
  geom_tile() +
  scale_fill_identity() +
  theme_minimal() +
  ylab("Users") +
  xlab("Features") +
  theme(axis.title.y = element_text(size = rel(1)))


plotDat <- reshape::melt(t(Astar_perm))

#and plot
ggplot(plotDat, aes(X1, X2, fill = c("white", "blue")[ value + 1 ])) +
  geom_tile() +
  scale_fill_identity() +
  theme_minimal() +
  ylab("Users") +
  xlab("Features") +
  theme(axis.title.y = element_text(size = rel(1)))

plotDat <- reshape::melt(t(B_true[,2:19]))

#and plot
ggplot(plotDat, aes(X1, X2, fill = c("white", "blue")[ value + 1 ])) +
  geom_tile() +
  scale_fill_identity() +
  theme_minimal() +
  ylab("Movies") +
  xlab("Features") +
  theme(axis.title.y = element_text(size = rel(1)))

plotDat <- reshape::melt(t(Bstar_perm))

#and plot
ggplot(plotDat, aes(X1, X2, fill = c("white", "blue")[ value + 1 ])) +
  geom_tile() +
  scale_fill_identity() +
  theme_minimal() +
  ylab("Movies") +
  xlab("Features") +
  theme(axis.title.y = element_text(size = rel(1)))#+


theta_true_sub=theta_true[2:19]
barplot(theta_true_sub,col='blue',ylim=c(-5,3))

plotDat_theta_true=data.frame(cbind(1:18,theta_true_sub))
colnames(plotDat_theta_true)=c("features",'theta')
pdf("sim-postsum-theta-true.pdf",width=6,height=2)
ggplot(plotDat_theta_true)+
  geom_bar(mapping=aes(x=features,y=theta),stat = 'identity',fill="blue")+
  ylab("Theta") +
  xlab("Features") +
  theme(axis.title.y = element_text(size = rel(1)))
dev.off()


theta_star_sub=theta_star[2:19]
barplot(theta_star_sub,col='blue',ylim=c(-5,3))

plotDat_theta_star=data.frame(cbind(1:18,theta_star_sub))
colnames(plotDat_theta_star)=c("features",'theta')
pdf("sim-postsum-theta-est.pdf",width=6,height=2)
ggplot(plotDat_theta_star)+
  geom_bar(mapping=aes(x=features,y=theta),stat = 'identity',fill="blue")+
  ylab("Theta") +
  xlab("Features") +
  theme(axis.title.y = element_text(size = rel(1)))
dev.off()


######################################

####for real data#####################
Astar_sub=A_star[,2:44]
Bstar_sub=B_star[,2:44]


library(ggplot2)
plotDat <- reshape::melt(t(Astar_sub))

#and plot
ggplot(plotDat, aes(X1, X2, fill = c("white", "blue")[ value + 1 ])) +
  geom_tile() +
  scale_fill_identity() +
  theme_minimal() +
  ylab("Users") +
  xlab("Features") +
  theme(axis.title.y = element_text(size = rel(1)))#+
  #coord_cartesian(ylim = c(1, 400),xlim=c(1,40)) +
  #scale_x_continuous(limits=c(1, 400), labels = c(1, seq(10, 40, by=10)))+
  #scale_y_continuous(limits=c(1, 400), labels = c(1, seq(100, 400, by=100)))


plotDat <- reshape::melt(t(Bstar_sub))

#and plot
ggplot(plotDat, aes(X1, X2, fill = c("white", "blue")[ value + 1 ])) +
  geom_tile() +
  scale_fill_identity() +
  theme_minimal() +
  ylab("Movies") +
  xlab("Features") +
  theme(axis.title.y = element_text(size = rel(1)))#+
#coord_cartesian(ylim = c(1, 400),xlim=c(1,40)) +
#scale_x_continuous(limits=c(1, 400), labels = c(1, seq(10, 40, by=10)))+
#scale_y_continuous(limits=c(1, 400), labels = c(1, seq(100, 400, by=100)))

theta_star_sub=theta_star[2:44]
barplot(theta_star_sub,col='blue',ylim=c(-5,3))

plotDat_theta_star=data.frame(cbind(1:43,theta_star_sub))
colnames(plotDat_theta_star)=c("features",'theta')
pdf("postsum-theta.pdf",width=6,height=2)
ggplot(plotDat_theta_star)+
  geom_bar(mapping=aes(x=features,y=theta),stat = 'identity',fill="blue")+
  ylab("Theta") +
  xlab("Features") +
  theme(axis.title.y = element_text(size = rel(1)))
dev.off()