# Main function to conduct mediation analysis in presence of time-varying mediators, a survival outcome and competing risks in survival probability scale

med_sp=function(L=NULL, M, m, Y, treat='logcocr', control.value=a, treat.value=a_star, data, time_points, t=20){
  
  N=dim(data)[1]
  NL=length(L)
  NM=length(M)
  
  boot <- foreach(i=1:1000, .combine='rbind') %dopar% {
    ind <- sample(1:N, replace=TRUE)
    PredictL_a <- PredictL_astar <- PredictL_astar_a <- matrix(NA,nrow=N,ncol=NL)
    PredictM_a <- PredictM_astar <- PredictM_a_astar <- matrix(NA,nrow=N,ncol=NM)
    
    # Predict M1
    pred.data.astar.m1 <- pred.data.a.m1 <- model.frame(M[[1]])[ind,]
    pred.data.astar.m1[, treat] <- treat.value
    pred.data.a.m1[, treat] <- control.value
    PredictM_astar[,1] <- predict(M[[1]], pred.data.astar.m1)
    PredictM_a[,1] <- predict(M[[1]], pred.data.a.m1)
  
    # Predict L1  
    if(NL > 1){
      pred.data.astar.l1 <- pred.data.a.l1 <- pred.data.astar.a.l1 <- model.frame(L[[1]])[ind,]
      pred.data.astar.l1[, treat] <- pred.data.astar.a.l1[, treat] <- treat.value
      pred.data.a.l1[, treat] <- control.value
      pred.data.astar.l1[, m[1]] <- PredictM_astar[,1]
      pred.data.a.l1[, m[1]] <- pred.data.astar.a.l1[, m[1]] <- PredictM_a[,1]
      
      PredictL_astar_a[,1] <- rbinom(N, size=1, prob=predict(L[[1]], pred.data.astar.a.l1, type='response'))
      PredictL_a[,1] <- rbinom(N, size=1, prob=predict(L[[1]], pred.data.a.l1, type='response'))
      PredictL_astar[,1] <- rbinom(N, size=1, prob=predict(L[[1]], pred.data.astar.l1, type='response'))
    }
    
    # Predict Li and Mi (only if more than one time point)
    if (NM > 1){
        for (i in 2:NM){
          pred.data.a.m <- pred.data.astar.m <- pred.data.a.astar.m <- as.data.frame(matrix(nrow=N, ncol=(dim(model.frame(M[[i]]))[2]-1)))
          colnames(pred.data.a.m) <- colnames(pred.data.astar.m) <- colnames(pred.data.a.astar.m) <- attr(terms(M[[i]]),"term.labels")
          names <- colnames(pred.data.a.m)[which(colnames(pred.data.a.m) %in% attr(terms(M[[1]]),"term.labels"))]
  
          pred.data.a.m[, names] <- pred.data.a.astar.m[, names] <- pred.data.astar.m[, names] <- model.frame(M[[1]])[ind,names]
          pred.data.a.m[, treat] <- pred.data.a.astar.m[, treat] <- control.value
          pred.data.astar.m[, treat] <- treat.value
          
          pred.data.a.m[, m[i-1]] <- PredictM_a[,i-1]
          pred.data.astar.m[, m[i-1]] <- PredictM_astar[,i-1]
          pred.data.a.astar.m[, m[i-1]] <- PredictM_a_astar[,i-1]
          if(i==2){
            pred.data.a.astar.m[, m[1]] <- PredictM_a[,1]
          }
          if(NL > 1){
            PredictM_a[which(PredictL_a[,i-1]==0),i] <- predict(M[[i]], pred.data.a.m[which(PredictL_a[,i-1]==0),])
            PredictM_astar[which(PredictL_astar[,i-1]==0),i] <- predict(M[[i]], pred.data.astar.m[which(PredictL_astar[,i-1]==0),])
            PredictM_a_astar[which(PredictL_astar_a[,i-1]==0),i] <- predict(M[[i]], pred.data.a.astar.m[which(PredictL_astar_a[,i-1]==0),])
          } else{
            PredictM_a[,i] <- predict(M[[i]], pred.data.a.m)
            PredictM_astar[,i] <- predict(M[[i]], pred.data.astar.m)
            PredictM_a_astar[,i] <- predict(M[[i]], pred.data.a.astar.m)
          }
          
          if(NL > 1 & i<=NL){
            pred.data.a.l <- pred.data.astar.l <- pred.data.astar.a.l <- as.data.frame(matrix(nrow=N, ncol=(dim(model.frame(L[[i]]))[2]-1)))
            colnames(pred.data.a.l) <- colnames(pred.data.astar.l) <- colnames(pred.data.astar.a.l) <- attr(terms(L[[i]]),"term.labels")
            names <- colnames(pred.data.a.l)[which(colnames(pred.data.a.l) %in% attr(terms(L[[1]]),"term.labels"))]
  
            pred.data.a.l[, names] <- pred.data.astar.a.l[, names] <- pred.data.astar.l[, names] <- model.frame(L[[1]])[ind,names]
            pred.data.a.l[, treat] <- control.value
            pred.data.astar.l[, treat] <- pred.data.astar.a.l[, treat] <- treat.value
            
            pred.data.a.l[, m[i]] <- PredictM_a[,i]
            pred.data.astar.l[, m[i]] <- PredictM_astar[,i]
            pred.data.astar.a.l[, m[i]] <- PredictM_a_astar[,i]
           
            PredictL_a[which(PredictL_a[,i-1]==0),i] <- rbinom(length(which(PredictL_a[,i-1]==0)), size=1, prob=predict(L[[i]], pred.data.a.l[which(PredictL_a[,i-1]==0),], type='response'))
            PredictL_astar[which(PredictL_astar[,i-1]==0),i] <- rbinom(length(which(PredictL_astar[,i-1]==0)), size=1, prob=predict(L[[i]], pred.data.astar.l[which(PredictL_astar[,i-1]==0),], type='response'))
            PredictL_astar_a[which(PredictL_astar_a[,i-1]==0),i] <- rbinom(length(which(PredictL_astar_a[,i-1]==0)), size=1, prob=predict(L[[i]], pred.data.astar.a.l[which(PredictL_astar_a[,i-1]==0),], type='response'))
          }}}
    
    # Predict Y
    # Data augmentation method for person-time database
    # PredictY_DEIEM: a*, D1_a, M1_aD1a, D2_aD1aM1aD1a, M2_aD1aM1aD2a
    # PredictY_TEDE_2: a, D1_a, M1_aD1a, D2_aD1a M1aD1a, M2_aD1aM1aD2a
    # PredictY_IEMIED: a*, D1_a, M1_a*D1a, D2_aD1a M1a*D1a, M2_a*D1aM1a*D2a
    # PredictY_IEDTE_1: a*, D1_a*, M1_a*D1a*, D2_a*D1a*M1a*D1a*, M2_a*D1a*M1a*D2a*
  
    pred.data.a.y <- pred.data.astar.y <- pred.data.astar.a.y <- pred.data.astar.a.astar.a.y <- data[ind,c('idno',getvarnames(Y$call)$xvar[-2],m,colnames(data)[grep('time.since.first.exam', colnames(data))])]
    pred.data.a.y[, treat] <- control.value
    pred.data.astar.y[, treat] <- pred.data.astar.a.y[, treat] <- pred.data.astar.a.astar.a.y[, treat] <- treat.value
    pred.data.a.y[, m] <- pred.data.astar.a.y[, m] <- PredictM_a
    pred.data.astar.y[, m] <- PredictM_astar
    pred.data.astar.a.astar.a.y[, m] <- PredictM_a_astar
    pred.data.astar.a.astar.a.y[, m[1]] <- PredictM_a[,1]

    
    ########################
    
    # Data augmentation method for the counterfactuals
    vector_time_points <- c()
    for (i in 1:length(time_points)){
      vector_time_points <- c(vector_time_points, m[i], time_points[i])
    }
    
    # pred.data.a.y
    pred.data.a.y$id_boot <- seq(1:dim(pred.data.a.y)[1])
    df_tv <- reshape(pred.data.a.y, direction = "long", varying = vector_time_points,
                     sep = "_", times=as.character(seq(1,length(time_points))), idvar='id_boot')
    df_tv <- df_tv[order(df_tv$id_boot),]
    df_pred.data.a.y <- df_tv[,match(getvarnames(Y$call)$xvar,colnames(df_tv))]
    df_pred.data.a.y <- model.matrix(~.,data=df_pred.data.a.y)[,-1]
    
    # pred.data.astar.y
    pred.data.astar.y$id_boot <- seq(1:dim(pred.data.astar.y)[1])
    df_tv <- reshape(pred.data.astar.y, direction = "long", varying = vector_time_points,
                     sep = "_", times=as.character(seq(1,length(time_points))), idvar='id_boot')
    df_tv <- df_tv[order(df_tv$idno),]
    df_pred.data.astar.y <- df_tv[,match(getvarnames(Y$call)$xvar,colnames(df_tv))]
    df_pred.data.astar.y <- model.matrix(~.,data=df_pred.data.astar.y)[,-1]
    
    # pred.data.astar.a.astar.a.y
    pred.data.astar.a.astar.a.y$id_boot <- seq(1:dim(pred.data.astar.a.astar.a.y)[1])
    df_tv <- reshape(pred.data.astar.a.astar.a.y, direction = "long", varying = vector_time_points,
                     sep = "_", times=as.character(seq(1,length(time_points))), idvar='id_boot')
    df_tv <- df_tv[order(df_tv$idno),]
    df_pred.data.astar.a.astar.a.y <- df_tv[,match(getvarnames(Y$call)$xvar,colnames(df_tv))]
    df_pred.data.astar.a.astar.a.y <- model.matrix(~.,data=df_pred.data.astar.a.astar.a.y)[,-1]
    
    # pred.data.astar.a.y
    pred.data.astar.a.y$id_boot <- seq(1:dim(pred.data.astar.a.y)[1])
    df_tv <- reshape(pred.data.astar.a.y, direction = "long", varying = vector_time_points,
                     sep = "_", times=as.character(seq(1,length(time_points))), idvar='id_boot')
    df_tv <- df_tv[order(df_tv$idno),]
    df_pred.data.astar.a.y <- df_tv[,match(getvarnames(Y$call)$xvar,colnames(df_tv))]
    df_pred.data.astar.a.y <- model.matrix(~.,data=df_pred.data.astar.a.y)[,-1]
    
    #######################
  
    PredictY_DEIEM <- predict.aalen(object=Y, Z=df_pred.data.astar.a.y, times=t, n.sim=0)
    PredictY_TEDE_2 <- predict.aalen(object=Y, Z=df_pred.data.a.y, times=t, n.sim=0)
    PredictY_IEMIED <- predict.aalen(object=Y, Z=df_pred.data.astar.a.astar.a.y, times=t, n.sim=0)
    PredictY_IEDTE_1 <- predict.aalen(object=Y, Z=df_pred.data.astar.y, times=t, n.sim=0)
    
    DE <- mean(PredictY_DEIEM$S0) - mean(PredictY_TEDE_2$S0)
    IEM <- mean(PredictY_IEDTE_1$S0) - mean(PredictY_IEMIED$S0)
    IED <- mean(PredictY_IEMIED$S0) - mean(PredictY_DEIEM$S0)
    TE <- mean(PredictY_IEDTE_1$S0) - mean(PredictY_TEDE_2$S0)
    effects <- cbind(DE, IEM, IED, TE)
    effects
    }
