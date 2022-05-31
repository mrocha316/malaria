#####
#####
##### SINGLE STRAIN, N-TREATMENT GROUPS
#####
#####

fthreeN <- function(first, iterations = 50, n=0, conv = 1*exp(-16)){
  x = matrix(0, nrow=iterations, ncol=5+n)
  x[1,] = first
  if(n==0){
    colnames(x) <- c("Sh", "Ih", "Rh", "Sm", "Im")
    gam <- function(x){
      betaM*((x[2])/(sum(x[-c(4,5)])))
    }
  } else {
    colnames(x) <- c("Sh", "Ih", "Rh", "Sm", "Im", paste("T", 1:n, sep=""))
    gam <- function(x){
      betaM*((x[2]+ psi[1:n] %*% x[6:(5+n)])/(sum(x[-c(4,5)])))
    }
  }
  change <- matrix(c(LamH, 0, 0, LamM, 0, 
                     0, 0, 0, 0, 0, 
                     0, 1/step -MuH-delS-sigma, sigma, 0, 0, 
                     omega, 0, 1/step - MuH - omega, 0, 0, 
                     0, 0, 0, 0, 0, 
                     0, 0, 0, 0, 1/step -MuM), nrow=6, byrow=T)
  if(n != 0){
    change2 <- matrix(0, nrow=n, ncol=n)
    for(i in 1:n-1){
      change2[i,i] <- 1/step - MuH - delS - alpha[[i+1]]
      change2[i, i+1] <- alpha[[i+1]]
    }
    change2[n,n] <- 1/step - MuH - delS - alpha[[n+1]]
    change <- rbind(cbind(change, matrix(0, nrow=6, ncol=n)), cbind(matrix(0, nrow=n, ncol=5), change2))
    change[3,2] <- 1/step -MuH-delS-sigma - alpha[[1]]
    change[3,6] <- alpha[[1]]
    change[6+n,1] <- alpha[[n]]
  }
  for(i in 2:iterations){
    lam = betaH*(x[i-1,5]/(sum(x[i-1,-c(4,5)]))) 
    gamm = gam(x[i-1,])
    change[2,1] <- 1/step - MuH - lam
    change[2,2] <- lam
    change[5,4] <- 1/step - MuM - gamm
    change[5,5] <- gamm
    x[i,] <- c(1,x[i-1,]) %*% (change*step)
    if(sum(abs((x[i-1,] - x[i,]))) < conv){
      print(paste("Convergence Successful in ", i, "iterations."))
      break
    }
  }
  rownames(x) <- seq(1, nrow(x))
  #  return(change)
  return(x[1:i,])
}


poph <- 35*10^6
popm <- 3*poph

MuH <- 4.43*10^(-5)     ##Human Death Rate [0,1]
LamH <- MuH * poph      ##Human Birth Rate [0,1]  ##death rate * population
MuM <- 1/14     ##Mosquito Death Rate [0,1]
LamM <- MuM*popm     ##Mosquito Birth Rate [per capita]    ##death rate * population
omega <- 1/28    ##Rate at which humans lose natural immunity and return to susceptible [0,1]
sigma <- 1/33    ##Rate of natural recovery from infection [0,1]
betaH <- 0.02    ##Rate of mosquito to human infection [0,1]
betaM <- 0.01     ##Rate of human to mosquito infection [0,1]
delS <- 1.43*10^(-5)     ##additional disease-induced death rate [0,1]
alpha0 <- 0.05     #Treatment rate into first treated stage (Ih to T1) [0,1]
alpha1 <- 0.05     #Treatment rate into second treated stage (T1 to T2) [0,1]
alpha2 <- 0.05     #Treatment rate into third treated stage (T2 to T3) [0,1]
alpha3 <- 0.05     #Treatment rate into susceptible treated stage (T3 to Sh) [0,1]
psi1 <- 0.90       #Transmission reduction factor from human in treatment group 1 to mosquito
psi2 <- 0.75       #Transmission reduction factor from human in treatment group 2 to mosquito
psi3 <- 0.5

x0 <- c(.8*poph, .1*poph, .1*poph, .8*popm, .2*popm, 0, 0, 0)  ##Initial conditions

step <- 0.01 ##Precision of time step - to get exponentially smaller for precision visualization
time <- seq(0, 500, by=step)

alpha <- c(alpha0, alpha1, alpha2, alpha3)
psi <- c(psi1, psi2, psi3)


results2N <- fthreeN(x0, iterations=length(time), n=3, conv=1*exp(-16))

stbl <- min(length(time), nrow(results2N))

plot(x=time[1:stbl], results2N[,1], type="l", lwd=1, col="blue", xlab="Time", ylab="Humans", ylim=c(0,poph))
lines(x=time[1:stbl], results2N[,2], col="light blue")
lines(x=time[1:stbl], results2N[,3], col="dark blue")
lines(x=time[1:stbl], results2N[,6], col="mediumpurple")
lines(x=time[1:stbl], results2N[,7], col="purple")
lines(x=time[1:stbl], results2N[,8], col="midnightblue")
#legend("topright", c("Susceptible Humans", "Infected Humans", "Recovered Humans", "Treatment1", "Treatment2",
#                     "Treatment3"), lty=1, col=c("blue", "light blue", "dark blue", "mediumpurple", "purple", 
#                                                 "midnightblue"))
legend(225,0.8*poph, c("Susceptible Humans", "Infected Humans", "Recovered Humans", "Treatment1", "Treatment2",
                       "Treatment3"), lty=1, col=c("blue", "light blue", "dark blue", "mediumpurple", "purple", 
                                                   "midnightblue"))


plot(x=time[1:stbl], results2N[,4], type="l", lwd=1, col="red", xlab="Time", ylab="Mosquitos", ylim=c(0,popm))
lines(x=time[1:stbl], results2N[,5], col="dark red")
legend(200, 0.8*popm, c("Susceptible Mosquitos", "Infected Mosquitos"), lty=1, col=c("red",  "dark red"))

system.time(fthreeN(x0, iterations=length(time), n=3, conv=1*exp(-16)))[3]



head(as.data.frame(results2N))
test2 <- as.data.frame(results2N)

as.data.frame(results2N) %>% 
  ggplot(aes(x=seq(1:stbl)))+
  geom_line(aes(y=Sh), color="blue")+
  geom_line(aes(y=Ih), color="light blue")+
  geom_line(aes(y=Rh), color="dark blue")+
  geom_line(aes(y=T1), color="mediumpurple")+
  geom_line(aes(y=T2), color="purple")+
  geom_line(aes(y=T3), color="midnightblue")

colnames(results2N)[1]

as.data.frame(results2N) %>% 
  ggplot(aes(x=seq(1:stbl)))+
  geom_line(aes(y=get(colnames(results2N)[1])), color="blue")+
  geom_line(aes(y=get(colnames(results2N)[2])), color="light blue")+
  geom_line(aes(y=get(colnames(results2N)[3])), color="dark blue")+
  geom_line(aes(y=get(colnames(results2N)[6])), color="mediumpurple")+
  geom_line(aes(y=get(colnames(results2N)[7])), color="mediumpurple")+
  geom_line(aes(y=get(colnames(results2N)[8])), color="mediumpurple")
