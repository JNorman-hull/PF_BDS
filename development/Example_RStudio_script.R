library(reshape)
library(MASS)
library(writexl)
library(BETS)
library(stringr) 
library(pracma)
#library(rstantools)
#library(ggstatsplot)
#library(WRS2)
#library(gridExtra)
library(readxl)
library(plotly)
library(dplyr)
#library(ggbiplot)
library(ggplot2)
library(tidyverse)
#library(plotrix)
#library(patchwork) # To display 2 charts together
#library(hrbrthemes)
#library(pwr)
#library(coin)
#library(MKpower)
#library(future.apply)
#library(effectsize)
library(gridExtra)
library(purrr)
####################################################################################################################################################################

#important
#loading all excel sheets
files = list.files(path="F:/Test_BDS_Josh/", pattern = ".csv")

files = str_c("F:/Test_BDS_Josh/", files)
names(files) = c("C360328110334", "C370328110457","C420328110536","C440328110507","U100328110504", "U130328110330","U170328110314","U240328110512")
data2= map_df(.x=files, .f=read.csv2, .id = "data_source",sep= ",", header = FALSE)

colnames(data2) <- c("data_source","Time", "PL","TL","PC" ,"TC","PR","TR",
                    "AccX", "AccY", "AccZ", "RotX",
                   "RotY", "RotZ")
#clean table
data2 = data2[-c(1), ]
data2 = data2[,-c(15:18)]
##################################################################################################################################################################
#calculating total pressure and acceleration magnitude and add to table - first sensor data as example
C36 <- subset(data2, data_source=="C360328110334")
C36 <-  read.csv2("C360328110334.csv",sep="," ,header = FALSE)
colnames(C36) <- c("Time", "PL","TL","PC" ,"TC","PR","TR",
                     "AccX", "AccY", "AccZ", "RotX",
                     "RotY", "RotZ")
C36 =C36[-c(1),]

C36  =  C36  %>% mutate_at(c('AccX', 'AccY','AccZ','PL','PC','PR'), as.numeric)

l=0
AccMag <-data.frame(AccMag=rep(NA,86289))
for (i in 1:86289){
  l=sqrt(C36[i,9]^2+C36[i,10]^2+C36[i,11]^2)
  AccMag[i,]=l
  i=i+1
}
C36 = cbind(C36,AccMag)

s=0
TotalPressure <-data.frame(TotalPressure=rep(NA,86289))
for (i in 1:86289){
  s=median(C36[i,3],C36[i,5],C36[i,7])
  TotalPressure[i,]=s
  i=i+1
}
C36 = cbind(C36,TotalPressure)

C36 =  C36 %>% mutate_at(c('Time', 'AccMag','TotalPressure'), as.numeric)

#######################################################################################################################################################
#calculating total pressure and acceleration magnitude and add to table for all
data2 =  data2 %>% mutate_at(c('AccX', 'AccY','AccZ','PL','PC','PR'), as.numeric)

l=0
AccMag <-data.frame(AccMag=rep(NA,690183))
for (i in 1:690183){
  l=sqrt(data2[i,9]^2+data2[i,10]^2+data2[i,11]^2)
  AccMag[i,]=l
  i=i+1
}
data2 = cbind(data2,AccMag)

s=0
TotalPressure <-data.frame(TotalPressure=rep(NA,690183))
for (i in 1:690183){
  s=median(C36[i,3],C36[i,5],C36[i,7])
  TotalPressure[i,]=s
  i=i+1
}
data2 = cbind(data2,TotalPressure)

data2 =  data2 %>% mutate_at(c('Time', 'AccMag','TotalPressure'), as.numeric)

###################################################################################################################################################
#ROI graph - first sensor data as example

b = ggplot()+
  geom_line(data=C36, aes(x=Time, y=TotalPressure, colour="blue"))
ggplotly(b) 


a = ggplot()+
  geom_line(data=C36, aes(x=Time, y=AccMag, colour="red"))
ggplotly(a) 



#find acc. peaks with min distance of 11 = 0.1sec
#series 15-17 12 steps = 0.1 sec
#series 6-14 11 steps = 0.1 sec
#1 sec = ~96 steps
t=C36[,2]  
d=findpeaks(x= C36$AccMag, minpeakheight=50, minpeakdistance= 12, sortstr=FALSE)
print(d)



###################################################################################################################################################
#creating differnece plot
#example the first three BDS data sets of data2

C36 <- subset(data2, data_source=="C360328110334")
C37 <- subset(data2, data_source=="C370328110457")
C37 = C37[-c(1), ]
C42 <- subset(data2, data_source=="C420328110536")
C42 = C42[-c(1), ]
C37  =  C37  %>% mutate_at(c('Time'), as.numeric)
C42  =  C42  %>% mutate_at(c('Time'), as.numeric)



C36_BN <- C36 %>% filter(Time < 162.386)
C36_AN <- C36 %>% filter(Time >= 162.386)
n1=normalize(C36_BN $Time, mode = "maxmin")
n2=normalize(C36_AN$Time, mode = "maxmin")
a=n1*0.5
b=(n2*0.5)+0.5
norm_timeC36=c(a,b)
d1 = data.frame(x=C36$AccMag, y=norm_timeC36)

C37_BN <- C37 %>% filter(Time < 364.940)
C37_AN <- C37 %>% filter(Time >= 364.940)
n3=normalize(C37_BN $Time, mode = "maxmin")
n4=normalize(C37_AN$Time, mode = "maxmin")
c=n3*0.5
d=(n4*0.5)+0.5
norm_timeC37=c(c,d)
d2 = data.frame(x=C37$AccMag, y=norm_timeC37)

C42_BN <- C42 %>% filter(Time < 431.051)
C42_AN <- C42 %>% filter(Time >= 431.051)
n5=normalize(C42_BN $Time, mode = "maxmin")
n6=normalize(C42_AN$Time, mode = "maxmin")
e=n5*0.5
f=(n6*0.5)+0.5
norm_timeC42=c(e,f)
d3 = data.frame(x=C42$AccMag, y=norm_timeC42)

min(nrow(d1),nrow(d2), nrow(d3))
#86289

data <-data.frame(NA_col=rep(NA, 86289))
for (i in 1:86289){
  s <- min( d1[i,1], d2[i,1],d3[i,1])
  m <- max( d1[i,1], d2[i,1],d3[i,1])
  dif= m-s
  data[i, ]<- dif
  i=i+1
}


data  =  data  %>% mutate_at(c('NA_col'), as.numeric)
bias <- mean(data$NA_col)
sd <- sd(data$NA_col)
lower <- bias - 1.96*sd
upper <- bias + 1.96*sd


d=ggplot()+
  geom_point(data=data, aes(x=norm_timeC36, y=NA_col, colour="red"))+
  geom_hline(yintercept = lower, color ="black", linetype= "dashed")+
  geom_hline(yintercept = upper, color ="black", linetype= "dashed")+
  geom_vline(xintercept = 0.51, color = "blue", linetype = "dashed")+
  geom_vline(xintercept = 0.49, colour = "blue", linetype = "dashed")+
  ylab("Difference in Acceleration Magnitute ")+
  xlab("normalised Time")
ggplotly(d)

#0.51 and 0.49 - time steps around Nadir (0.5s)
##############################################################################################################################################################
#creating nice plots with acceleration and pressure - first sensor data as example

C36 =C36 %>% filter(between(Time, 156.111, 220.636))# if you only want time between injection and tailwater
C36 =C36 %>% filter(between(Time, 162.386-0.5, 162.386+0.5)) #if you only want time between 0.5s before the nadir and 0.5s after nadir

c=ggplot()+
  geom_line(data=C36, aes(x=Time, y=AccMag+1000, color="Acceleration Magnitude"))+
  geom_line(data=C36, aes(x=Time, y=TotalPressure, color="Total Pressure"))+
  geom_vline(xintercept = 162.386+0.5,color = "black", linetype = "dashed", size=1)+
  geom_vline(xintercept =162.386-0.5,color = "black", linetype = "dashed", size=1)+
  #geom_vline(xintercept =0+3,color = "red", linetype = "dashed", size=1)+
  #geom_vline(xintercept = 156.111,color = "orange", linetype = "dashed", size=1)+
  geom_point(aes(x=162.386, y=934.78,color="Nadir"))+
  theme(legend.key.size = unit(2,'cm'),legend.title = element_text(size=18),
        legend.text = element_text(size=16), plot.title = element_text(size=22),
        axis.title.x = element_text(size =18),axis.title.y = element_text(size =18),
        axis.text.x = element_text(size = 16),axis.text.y = element_text(size = 16))+
  # scale_color_manual(breaks = c("Before", "Around_Nadir", "After", "Nadir"),
  #   values=c("red", "green", "blue","black"))+
  labs(title = "BDS C36")+
  scale_y_continuous(name ="Total Pressure [mbar]", breaks = seq(600,1300, 50),limits=c(600,1300) ) + scale_x_continuous(name ="Time [s]") 
#  xlab("Time [s]")+ ylab("Total Pressure [mbar]")
ggplotly(c)

#########################################################################################################################################################
#boxplot for pressure comparison

P =read_excel("BDS_Data_Denver_PRC.xlsx")


m = ggplot()+
  geom_boxplot(P, aes(x=operational_scenario, y=LRP))
ggplotly(m)  

l = ggplot()+
  #geom_boxplot(data=BDS_Data_Denver_PRC, aes(x=operational_scenario, y=max_pressure))
  geom_boxplot(data=BDS_Data_Denver_PRC, aes(x=operational_scenario, y=Nadir_Pressure))
ggplotly(l) 