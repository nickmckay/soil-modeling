library(SoilR)

t_start=0
t_end=10
t_step=1/12
t=seq(t_start,t_end,t_step)


t=seq(0,10,by=0.2)
k=0.8 # 1/time
C0=100 # mass
In = 30 # mass/time


Model1=OnepModel(t,k,C0,In)
Ct1=getC(Model1)

#my own version
cstock <- c()
for(ti in seq_along(t)){
  if(ti == 1){
    cstock[ti] <- C0
  }else{
    timestep <- t[ti] - t[ti-1]
    dCdT <- In - cstock[ti - 1]  * (k * timestep)
    cstock[ti] <- cstock[ti - 1] + dCdT
  }
}



plot(t, Ct1, type="l", ylab="Carbon stocks (kg)", xlab="Time (years)")
lines(t, cstock, col = "red")

