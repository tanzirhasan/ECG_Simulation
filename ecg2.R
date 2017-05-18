

x<-seq(0,5,0.001)


generate_ecg <- function(x,bpm,a1,a2,a3,a4,a5){          #a1= amplitude of p, a2= amplitude of q, a3= amplitude of QR, a4= amplitude of S, a5= amplitude of T wave
  l= 60 / bpm 
  p_wav<-function(x){
  x=x+(1/1.8)# x is the wave starting point shifted .
  b=3# 2l/b is duration of p wave.
  n=100 # fourier series levels, the bigher the more accurate
  p1=1/l # baseline
  p2=0 # p wave
  # fourier series to creat p wave
  for (i in 1:n){
    harm1<-(((sin((pi/(2*b))*(b-(2*i))))/(b- (2*i))+(sin((pi/(2*b))*(b+(2*i))))
             /(b+(2*i)))*(2/pi))*cos((i*pi*x)/l)
    p2<-p2+harm1
  }
  pwav1=p1+p2
  pwav=a1*pwav1
  }
  
  pwav<-p_wav(x)
 
  ####### q wave
  q_wave<-function (x){
    #l=1
    x=x+1/6
    b=16#duration
    n=100
    q1=(a2/(2*b))*(2-b)
    q2=0
    for (i in 1:n){
      harm5=(((2*b*a2)/(i*i*pi*pi))*(1-cos((i*pi)/b)))*cos((i*pi*x)/l);
      q2=q2+harm5}
    qwav=-1*(q1+q2)
  }
  qwav<-q_wave(x)
  
  ####qrs wave
  qrs_wave<-function(x){
    b=5
    n=100
    qrs1=(a3/(2*b))*(2-b)
    qrs2=0
    for (i in 1:n){
      harm=(((2*b*a3)/(i*i*pi*pi))*(1-cos((i*pi)/b)))*cos((i*pi*x)/l);
      qrs2=qrs2+harm}
    qrswav=qrs1+qrs2
  }
  qrswav<-qrs_wave(x)
  
  ##### s wave
  s_wave<-function(x){
    x=x-1/6
    b=15
    n=100;
    s1=(a4/(2*b))*(2-b);
    s2=0;
    for (i in 1:n){
      harm3=(((2*b*a4)/(i*i*pi*pi))*(1-cos((i*pi)/b)))*cos((i*pi*x)/l);
      s2=s2+harm3}
    swav=-1*(s1+s2)
    }
  swav<-s_wave(x)
  
  ##### t wave
  
  t_wave<-function(x){
    x=x-1/1.8
    b=7
    n=20
    t1=1/l
    t2=0
    for (i in 1:n){
      harm2=(((sin((pi/(2*b))*(b-(2*i))))/(b-(2*i))+(sin((pi/(2*b))*(b+(2*i))))/(b+(2*i)))*(2/pi))*cos((i*pi*x)/l);             
      t2=t2+harm2}
    twav1=t1+t2
    twav=a5*twav1
    
    }
  twav<-t_wave(x)
  ###### adding waves
  ecg<-pwav+qrswav+twav+swav+qwav
  plot(x,ecg,type = 'l',xlab = 'x',ylab='Mv',main='ECG wave')}

generate_ecg(x,72,0.05,.025,1,.25,0.12)
   
