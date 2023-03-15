#' Gives test statistic value and its critical value for rank based test
#' @export
#' @param x1 is sample from Exp Population 1.
#' @param x2 is sample from Exp Population 2.
#' @param x3 is sample from Exp Population 3.
#' @param alpha is numeric value (level of significance)

Rank3EXP<-function(x1,x2,x3,alpha){

  k=3
  nsim<-rep(0,k)
  nsim[1]<-length(x1)
  nsim[2]<-length(x2)
  nsim[3]<-length(x3)

  N<-sum(nsim)


  bl=1000




  Q_boot<-rep(1,bl)






  c_function<-function(u){
    if(u<0)
      0
    else if(u==0)
      1/2
    else
      1
  }
  phi_n_function<-function(x1,x2,x3) {


    #####################creating upper indices ##################
    ui_11<-rank(c(x1))
    ui_12<-rank(c(x1,x2))
    ui_13<-rank(c(x1,x3))
    ui_21<-rank(c(x2,x1))
    ui_22<-rank(c(x2))
    ui_23<-rank(c(x2,x3))
    ui_31<-rank(c(x3,x1))
    ui_32<-rank(c(x3,x2))
    ui_33<-rank(c(x3))

    ######################## creating R_ij ^{li}##################
    R_1<-matrix(0,k,nsim[1])
    R_2<-matrix(0,k,nsim[2])
    R_3<-matrix(0,k,nsim[3])


    for (i in 1:nsim[1]) {
      R_1[1,i]<-ui_11[i]
    }
    for (i in 1:nsim[1]) {
      R_1[2,i]<-ui_21[nsim[2]+i]
    }
    for (i in 1:nsim[1]) {
      R_1[3,i]<-ui_31[i+nsim[3]]
    }




    for (i in 1:nsim[2]) {
      R_2[1,i]<-ui_12[nsim[1]+i]
    }
    for (i in 1:nsim[2]) {
      R_2[2,i]<-ui_22[i]
    }
    for (i in 1:nsim[2]) {
      R_2[3,i]<-ui_32[i+nsim[3]]
    }


    for (i in 1:nsim[3]) {
      R_3[1,i]<-ui_13[nsim[1]+i]
    }
    for (i in 1:nsim[3]) {
      R_3[2,i]<-ui_23[i+nsim[2]]
    }
    for (i in 1:nsim[3]) {
      R_3[3,i]<-ui_33[i]
    }

    ############## elements of R_i. ####################


    R_1_bar<-rep(0,k)
    for (i in 1:k) {
      R_1_bar[i]<-(rowSums(R_1))[i]/nsim[1]
    }

    R_2_bar<-rep(0,k)
    for (i in 1:k) {
      R_2_bar[i]<-(rowSums(R_2))[i]/nsim[2]
    }

    R_3_bar<-rep(0,k)
    for (i in 1:k) {
      R_3_bar[i]<-(rowSums(R_3))[i]/nsim[3]
    }

    ################################## w
    w<-function(a1,b1,c1){
      (1/a1)*(b1-((c1+1)/2))

    }


    w_11<-w(nsim[1],R_1_bar[1],nsim[1])
    w_21<-w(nsim[2],R_1_bar[2],nsim[1])
    w_31<-w(nsim[3],R_1_bar[3],nsim[1])
    w_12<-w(nsim[1],R_2_bar[1],nsim[2])
    w_22<-w(nsim[2],R_2_bar[2],nsim[2])
    w_32<-w(nsim[3],R_2_bar[3],nsim[2])
    w_13<-w(nsim[1],R_3_bar[1],nsim[3])
    w_23<-w(nsim[2],R_3_bar[2],nsim[3])
    w_33<-w(nsim[3],R_3_bar[3],nsim[3])

    wmatrix<-matrix(c(w_11,w_21, w_31,w_12, w_22, w_32,w_13, w_23, w_33), nrow = 3, ncol = 3)
    wmatrix
    ##################################### p_hat
    p_hat<-rep(0,k)
    for (i in 1:k) {
      p_hat[i]<- colSums(wmatrix)[i]/k
    }
    p_hat
    ################################# create D
    tau<-array(1,dim = c(k,k,k))
    D_1<-matrix(1,nrow = k,ncol = nsim[1])
    for (s in 1:k) {
      for (j in 1:nsim[1]) {
        D_1[s,j]<-(1/nsim[s])*((R_1[s,j]-R_1[1,j])-(R_1_bar[s]-((nsim[1]+1)/2)))
      }
    }
    D_2<-matrix(1,nrow = k,ncol = nsim[2])
    for (s in 1:k) {
      for (j in 1:nsim[2]) {
        D_2[s,j]<-(1/nsim[s])*((R_2[s,j]-R_2[2,j])-(R_2_bar[s]-((nsim[2]+1)/2)))
      }
    }
    D_3<-matrix(1,nrow = k,ncol = nsim[3])
    for (s in 1:k) {
      for (j in 1:nsim[3]) {
        D_3[s,j]<-(1/nsim[s])*((R_3[s,j]-R_3[3,j])-(R_3_bar[s]-((nsim[3]+1)/2)))
      }
    }



    ############################################## create tau

    tau1<-matrix(1,k,k)
    for (s in 1:k) {
      for (t in 1:k) {
        sum1<-0
        for (j in 1:nsim[1]) {
          sum1<-sum1+D_1[s,j]*D_1[t,j]
        }
        tau1[s,t]=sum1*sum(nsim)/(nsim[1]*(nsim[1]-1))
      }
    }



    tau2<-matrix(1,k,k)
    for (s in 1:k) {
      for (t in 1:k) {
        sum1<-0
        for (j in 1:nsim[2]) {
          sum1<-sum1+D_2[s,j]*D_2[t,j]
        }
        tau2[s,t]=sum1*sum(nsim)/(nsim[2]*(nsim[2]-1))
      }
    }


    tau3<-matrix(1,k,k)
    for (s in 1:k) {
      for (t in 1:k) {
        sum1<-0
        for (j in 1:nsim[3]) {
          sum1<-sum1+D_3[s,j]*D_3[t,j]
        }
        tau3[s,t]=sum1*sum(nsim)/(nsim[3]*(nsim[3]-1))
      }
    }





    tau_1<-as.vector(tau1)
    tau_2<-as.vector(tau2)
    tau_3<-as.vector(tau3)
    tau<-c(tau_1,tau_2,tau_3)

    tau_matrix<-array (c(tau), dim = c(k, k, k))




    ################################# creating matrix S


    S<-array(1,dim = c(k,k,k,k))




    for (i in 1:k) {
      for (j in 1:k) {
        for (m in 1:k) {
          for (o in 1:k) {
            if(m!=o){
              if(i!=j & i!=o & j==m)
                S[i,j,m,o]<- -tau_matrix[i,o,m]
              else if(i!=j & i==o & j==m)
                S[i,j,m,o]<- -tau_matrix[o,o,m]-tau_matrix[m,m,o]
              else if(i!=j & i==o & j!=m)
                S[i,j,m,o]<- -tau_matrix[m,j,i]
              else if(i==j & i!=o & j!=m)
                S[i,j,m,o]<- tau_matrix[m,o,i]
              else
                S[i,j,m,o]<-0
            }
            else{
              if(i==j & i!=m)
                S[i,j,m,o]<-tau_matrix[i,i,m]+tau_matrix[m,m,i]
              else if(i!=j & m!=i & m!=j)
                S[i,j,m,o]<-tau_matrix[i,j,m]
              else
                S[i,j,m,o]<-0
            }

          }
        }
      }

    }
    ######################################### Create V

    I<-rep(1,k)
    V<-matrix(1,nrow = k,ncol = k)
    for (i in 1:k) {
      for (j in 1:k) {
        V[i,j]<- (t(I)%*%S[,, i,j]%*%I)/(k^2)
      }
    }

    ##################################### Create T

    C_matrix<-matrix(0,nrow = k-1, ncol = k)
    for (i in 1:k-1) {
      for (j in 1:k) {
        if(i==j)
          C_matrix[i,j]<-1
        else if(j==i+1)
          C_matrix[i,j]<- -1
        else
          C_matrix[i,j]<-0
      }

    }
    C_inv <- inv( (C_matrix) %*% t(C_matrix))
    T_matrix<- (t(C_matrix)) %*% (C_inv) %*% (C_matrix)
    ##################################### Computing Q(T)


    TV<-T_matrix%*%V
    ev_TV<-eigen(TV)$values
    tr_TV<-tr(TV)

    Q<-(N/tr_TV)*( t(p_hat)%*%T_matrix%*%p_hat)



    #######################################Creating pseudo rank############################################
    G1_function<-function(x){
      a=0
      for (m in 1:nsim[1]) {
        a<-a+c_function(x-x1[m])
      }
      a/nsim[1]
    }

    G2_function<-function(x){
      a=0
      for (m in 1:nsim[2]) {
        a<-a+c_function(x-x2[m])
      }
      a/nsim[2]
    }

    G3_function<-function(x){
      a=0
      for (m in 1:nsim[3]) {
        a<-a+c_function(x-x3[m])
      }
      a/nsim[3]
    }
    G_function<-function(x){
      (G1_function(x)+G2_function(x)+G3_function(x))/3
    }

    Pseudo_R1<-rep(1,nsim[1])
    for (m in 1:nsim[1]) {
      Pseudo_R1[m]<-N*G_function(x1[m])+0.5
    }

    Pseudo_R2<-rep(1,nsim[2])
    for (m in 1:nsim[2]) {
      Pseudo_R2[m]<-N*G_function(x2[m])+0.5
    }

    Pseudo_R3<-rep(1,nsim[3])
    for (m in 1:nsim[3]) {
      Pseudo_R3[m]<-N*G_function(x3[m])+0.5
    }


    ############################################### S_i^2 ###########################
    S_1_sq=0
    for (m in 1:nsim[1]) {
      S_1_sq<-S_1_sq+(Pseudo_R1[m]-R_1[1,m]-mean(Pseudo_R1)+((nsim[1]+1)/2))^2
    }
    S_2_sq=0
    for (m in 1:nsim[2]) {
      S_2_sq<-S_2_sq+(Pseudo_R2[m]-R_2[2,m]-mean(Pseudo_R2)+((nsim[2]+1)/2))^2
    }
    S_3_sq=0
    for (m in 1:nsim[3]) {
      S_3_sq<-S_3_sq+(Pseudo_R3[m]-R_3[3,m]-mean(Pseudo_R3)+((nsim[3]+1)/2))^2
    }


    S_sq<-rep(0,k)
    S_sq[1]<-S_1_sq/(nsim[1]-1)
    S_sq[2]<-S_2_sq/(nsim[2]-1)
    S_sq[3]<-S_3_sq/(nsim[3]-1)

    f_1_hat<-(sum( (S_sq/(N-nsim)) ))^2/sum( ((S_sq/(N-nsim))^2)/(nsim-1) )
    ###################################################################################
    f_hat<-(tr(TV))^2 / (tr(TV %*% TV))

    df_hat<- qf(1-alpha,f_hat,f_1_hat,lower.tail = T)
    data.frame(statistic_a=Q,statistic_b= f_hat *  Q, statistic_c=Q, f_hat=f_hat, df_hat=df_hat)
  }



  phi_n<-phi_n_function(x1,x2,x3)$statistic_a
  N<-sum(nsim)
  theta_boot<-((nsim[1]*mean(x1))+(nsim[2]*mean(x2))+(nsim[3]*mean(x3)))/N
  phi_n_boot<-rep(0,bl)
  for (boot in 1:bl) {
    x1_boot<-rexp(nsim[1],1/theta_boot)
    x2_boot<-rexp(nsim[2],1/theta_boot)
    x3_boot<-rexp(nsim[3],1/theta_boot)
    phi_n_boot[boot]<-phi_n_function(x1_boot,x2_boot,x3_boot)$statistic_a
  }
  new_Q<- phi_n_boot[order(phi_n_boot)]
  C_star_Q<- new_Q[floor((1-alpha)*bl)]


  phintilde<-  phi_n_function(x1,x2,x3)$statistic_b

  phintildecv<-qchisq(1-alpha, df=phi_n_function(x1,x2,x3)$f_hat, lower.tail=TRUE)

  phin<-phi_n_function(x1,x2,x3)$statistic_c
  phincv<-phi_n_function(x1,x2,x3)$df_hat
  data.frame(test=c('PBQ','QB','QT'), test_statistic=c(phi_n, phintilde, phin), critical_value=c(C_star_Q,phintildecv,phincv))
  #data.frame(statistic_test_b=phintilde, critical_value_test_b=phintildecv)
}
