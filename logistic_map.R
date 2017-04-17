#l is half beat cycle (60s/heart/beat=1 cycle) a= amplitude
#x=x+p_r interval
#b=(2*l)/duration

#### step1: p wave
p_wave<-function(x,c)
{ l=c# half heartbeat cycle.eg.s/heartbeat=1 cycle
a=0.25#amplitude of p wave
x=x+(l/1.8)# x is the wave starting point shifted.
b=3# 2l/b is duration of p wave.
n=100 # fourier series levels, the bigher the more accurate
p1=0 # baseline
p2=0 # p wave
# fourier series to creat p wave
for (i in 1:n){
  harm1<-(((sin((pi/(2*b))*(b-(2*i))))/(b- (2*i))+(sin((pi/(2*b))*(b+(2*i))))
           /(b+(2*i)))*(2/pi))*cos((i*pi*x)/l)
  p2<-p2+harm1
}
pwav1=p1+p2
pwav=a*pwav1
return(pwav)
}
x<-seq(0,5,0.001) 
pwav<-p_wave(x, 0.5)
plot(x=x, y=pwav, type='l', xlim=c(0,5), ylim=c(0, 1), xlab="ECG wave", ylab="MV",main="ECG p wave simulation")

####step 2: simulation QRS wave:
####### q wave
q_wave<-function (x,c){
  l=c
  x=x+l/6
  a=0.03# amplitude 
  b=16#duration
  n=100
  q1=0 #(a/(2*b))*(2-b)
  q2=0
  for (i in 1:n){
    harm5=(((2*b*a)/(i*i*pi*pi))*(1-cos((i*pi)/b)))*cos((i*pi*x)/l);
    q2=q2+harm5}
  qwav=-1*(q1+q2)
  return (qwav)
}
x<-seq(0,5,0.001) 
qwav<-q_wave(x, 0.5)
plot(x=x, y=qwav, type='l', xlim=c(0,5), ylim=c(0, 1), xlab="ECG wave", ylab="MV",main="ECG q wave simulation")

####qrs wave
qrs_wave<-function(x,c){
  l=c
  a=1
  b=5
  n=100
  qrs1=0 #(a/(2*b))*(2-b)
  qrs2=0
  for (i in 1:n){
    harm=(((2*b*a)/(i*i*pi*pi))*(1-cos((i*pi)/b)))*cos((i*pi*x)/l);
    qrs2=qrs2+harm}
  qrswav=qrs1+qrs2
  return(qrswav)
}
x<-seq(0,5,0.001) 
qrswav<-qrs_wave(x, 0.5)
plot(x=x, y=qrswav, type='l', xlim=c(0,5), ylim=c(0, 1), xlab="ECG wave", ylab="MV",main="ECG qrs wave simulation")

##### s wave
s_wave<-function(x,c){
  l=c
  x=x-l/6
  a=0.25
  b=15
  n=100;
  s1=0 # (a/(2*b))*(2-b);
  s2=0;
  for (i in 1:n){
    harm3=(((2*b*a)/(i*i*pi*pi))*(1-cos((i*pi)/b)))*cos((i*pi*x)/l);
    s2=s2+harm3}
  swav=-1*(s1+s2)
  return(swav)}
x<-seq(0,2,0.001) 
swav<-s_wave(x, 0.5)
plot(x=x, y=swav, type='l', xlim=c(0,5), ylim=c(-2, 1), xlab="ECG wave", ylab="MV",main="ECG s wave simulation")

##### t wave

t_wave<-function(x,c){
  l=c
  a=0.35
  x=x-l/1.8
  b=7
  n=100
  t1=0 # 1/l
  t2=0
  for (i in 1:n){
    harm2=(((sin((pi/(2*b))*(b-(2*i))))/(b-(2*i))+(sin((pi/(2*b))*(b+(2*i))))/(b+(2*i)))*(2/pi))*cos((i*pi*x)/l);             
    t2=t2+harm2}
  twav1=t1+t2
  twav=a*twav1
  return (twav)}

x<-seq(0,2,0.001) 
twav<-t_wave(x, 0.5)
plot(x=x, y=twav, type='l', xlim=c(0,5), ylim=c(0, 1), xlab="ECG wave", ylab="MV",main="ECG t wave simulation")
# ###### adding waves
# ecg<-pwav+qrswav+twav+swav+qwav
# plot(x,ecg,xlab = 'x',ylab='Mv',type='l',main='ECG wave')

####### using logistic map to creat heart beat rhythm
T <- 200 # number of generations
#interval: 0.6-1s
p<-numeric(T)
length(p)
p<-rep(1.0,length.out=T)
p[1] <- 0.1  # starting population
r <- 3.9# reproduction rate; try different values, like 2.9, 3.1, 3.5, 3.9
for (i in 2:T) p[i] <- r * p[i - 1] * (1 - p[i - 1])
plot(p, type='l', ylim=c(0, 1), xlab="Generation", ylab="Population", main=sprintf("r = %0.2f", r))
#Q<-p[21:86]
#length(Q)
#plot(cumsum(Q)*60/sum(Q),rep(1,66))
abline(v=cumsum(Q))


############################### week 10
heart_rate = 60.0
Q<-p[21:(21 + heart_rate - 1)]# Q interval of heart beat
# Q/sum(Q) is percentage of total interval, *60 is within 60 seconds, the heart beat cycly.
Q<-Q/sum(Q) * 60.0
plot(cumsum(Q)*60/sum(Q),rep(1,60))
abline(v=cumsum(Q))
length(cumsum(Q)*60/sum(Q))
# Q is interval seconds.
x<-seq(0,60,0.008)
x_start = 0
print(x_start)
ecg<-numeric()
x_coordinate<-numeric()
i = 1
for (i in 1:length(Q)) {
  x_end = x_start + Q[i]
  if (i == length(Q)) {
    x_used = subset(x, x>=x_start)
  } else {
    x_used = subset(x, x>=x_start & x<x_end)
  }
  x_used = x_used - x_start
  pwav = p_wave(x_used, Q[i]/2)
  qrswav = qrs_wave(x_used, Q[i]/2)
  twav = t_wave(x_used, Q[i]/2)
  swav = s_wave(x_used, Q[i]/2)
  qwav = q_wave(x_used, Q[i]/2)
  combined_wav = pwav+qrswav+twav+swav+qwav
  ecg <- c(ecg, combined_wav)
  x_start = x_end
}
# Lift ecg graph to mostly above 0.
ecg = ecg + 0.2
plot(x,ecg,xlab = 'x',ylab='Mv', xlim=c(0,60), type='l',main='ECG wave')
plot(x,ecg,xlab = 'x',ylab='Mv', xlim=c(0,10), type='l',main='ECG wave')
abline(v=cumsum(Q),col='blue',lty=3)




# # Q is interval seconds.
# x<-seq(0,60,0.001)
# x_start = 0
# print(x_start)
# pwav<-numeric()
# qrswav<-numeric()
# twav<-numeric()
# swav<-numeric()
# qwav<-numeric()
# ecg<-numeric()
# i = 1
# for (i in 1:length(Q)) {
#   x_end = x_start + Q[i]
#   cat('x_start', x_start, 'x_end', x_end, '\n')
#   if (i == length(Q)) {
#     x_used = subset(x, x>=x_start)
#   } else {
#     x_used = subset(x, x>=x_start & x<x_end)
#   }
#   
#   pwav<-c(pwav, p_wave(x_used, Q[i]/2))
#   qrswav<-c(qrswav, qrs_wave(x_used, Q[i]/2))
#   twav<-c(twav, t_wave(x_used, Q[i]/2))
#   swav<-c(swav, s_wave(x_used, Q[i]/2))
#   qwav<-c(qwav, q_wave(x_used, Q[i]/2))
#  # pwav+qrswav+twav+swav+qwav
#   x_start = x_end
#   plot(x_used,pwav,xlab = 'x',ylab='Mv', xlim=c(0,3), type='l',main='ECG wave')
#   plot(x_used,qrswav,xlab = 'x',ylab='Mv', xlim=c(0,3), type='l',main='ECG wave')
#   plot(x_used,twav,xlab = 'x',ylab='Mv', xlim=c(0,3), type='l',main='ECG wave')
#   plot(x_used,swav,xlab = 'x',ylab='Mv', xlim=c(0,3), type='l',main='ECG wave')
#   plot(x_used,qwav,xlab = 'x',ylab='Mv', xlim=c(0,3), type='l',main='ECG wave')
# }
# ecg<-pwav+qrswav+twav+swav+qwav
# plot(x,ecg,xlab = 'x',ylab='Mv', xlim=c(0,60), type='l',main='ECG wave')
# plot(x,ecg,xlab = 'x',ylab='Mv', xlim=c(0,3), type='l',main='ECG wave')

#####plot 10
# ecg plot with in 50 heatbeats in 60s 
#shiny app
#write json file
#fit real ecg data
#writting more efficient codes, code optimzation. 
