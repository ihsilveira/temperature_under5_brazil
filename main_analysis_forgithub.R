
################################## DESCRIPTION ##########################################

# Working paper by Silveira et al. 
# Main analysis of temperature effects on under 5 mortality 
# We analyzed the effects by region, age group [neonatal (0-27 days), post-neonatal (28-364 days), childhood (1-4 years), total (under 5 years)] and causes of death
# we conducted space-time case-crossover analysis using conditional quasi-poisson regression models combined with dlnm models


################################## R PACKAGES REQUIRED ##########################################
if(!require(data.table)) install.packages("data.table"); library(data.table)
if(!require(dlnm)) install.packages("dlnm"); library(dlnm)
if(!require(splines)) install.packages("splines"); library(splines)
if (!require(gnm)) install.packages('gnm'); library(gnm)

################################## DATASET ##########################################

data_under5=fread("/home/ismael.silveira/mortality_under5/dataset/data_under5_20241206.csv") #esse aqui com causas externas
head(data_under5)
names(data_under5)

################################## GENERAL PARAMETERS FOR CROSS BASIS ##########################################
cb_base='ns'
lagn=7
lag_nk=2

result<-NULL # object to store the results

# before start the analysis we are creating a pdf for plots

pdf("/home/ismael.silveira/mortality_under5/main_and_additional_analysis/graph_result_teste_github.pdf",width=15,height=25)
layout(matrix(seq(5*3),nrow=5,byrow=T))
par(mex=0.8,mgp=c(2.5,1,0),las=1)

outcomes2<-c("0-27 days", "28-364 days", "1-4 years", "Total <5 years" )

################################## ANALYSIS WITH TOTAL UNDER FIVER MORTALITY BY REGION ##########################################

#-------------------------------------ANALISE BRAZIL--------------------------------------------------------

data_sub=data_under5

outcomes<-names(data_sub)[c(9,6,7,3)]

i=4


  
  data_sub <- data_under5
  
  data_sub <- data_sub[, count_cases := sum(get(outcomes[i])), by = stratum]
  
  data_sub=subset(data_sub, data_sub$count_cases>0)
  
  basis.temp <- crossbasis(data_sub$Daily_mean_temperature_XAVIER, lag= lagn, 
                           argvar=list(fun="ns",df=3),
                           arglag=list(fun="ns",knots=(logknots(lagn,lag_nk))))
  
  form<-as.formula(paste(outcomes[i],"~ basis.temp + ns(absumid_xavier,3) + feriado"))
  
  mod <- gnm(form, data=na.omit(data_sub), family = quasipoisson(), eliminate=factor(stratum))
  
  cp <- crosspred(basis.temp, mod)
  mmt=cp$predvar[which.min(cp$allRRfit)]
  mmp=sum(na.omit(cp$predvar) < mmt) / length(na.omit(cp$predvar))
  
  cp2 <- crosspred(basis.temp, mod,cen=mmt) 
  
  RR.eh <- unlist(crosspred(basis.temp, mod, at = quantile(cp2$predvar, 0.99,na.rm = T), cen = mmt)[c("allRRfit", "allRRlow", "allRRhigh")]) 
  RR.ec <- unlist(crosspred(basis.temp, mod, at = quantile(cp2$predvar, 0.01,na.rm = T), cen = mmt)[c("allRRfit", "allRRlow", "allRRhigh")])
  
  RR.table <- cbind(t(RR.eh),t(RR.ec),mmt,mmp)
  colnames(RR.table) <- c("RR.eh", "RR.eh.low", "RR.eh.high",
                          "RR.ec", "RR.ec.low", "RR.ec.high","mmt","mmp")
  
  RR.table<- round(RR.table,2)
  
  RR.table<-cbind("Brazil",outcomes[i],RR.table)
  
  result=rbind(result,RR.table)      
  
  y.lim.max<-round(max(RR.eh[1],RR.ec[1])+2,0)
  

  plot(cp2,"overall",ylim=c(0,y.lim.max),yaxt="n",lab=c(6,5,7),xlab="Temperature (°C)",ylab="RR",cex.lab=0.9,cex.axis = 0.9,main="Brazil", cex=1.0)
  
  axis(2)
  
  abline(v=mmt,lty=3)
  abline(v=c(quantile(cp2$predvar,c(0.01,0.99))),lty=2) ## We used predvar instead temperature variable
  

rm(list = ls()[!ls() %in% c("data_under5","result","cb_base","cb_per","lagn","lag_nk","outcomes2")])

#-------------------------------------CENTRAL-WEST REGION--------------------------------------------------------

data_region <- subset(data_under5, Sigla.da.UF %in% c("DF","GO","MS","MT"))

outcomes<-names(data_region)[c(9,6,7,3)] 

i=4 
  
  data_sub <- data_region 
  
  data_sub <- data_sub[, count_cases := sum(get(outcomes[i])), by = stratum]
  
  data_sub=subset(data_sub, data_sub$count_cases>0)
  
  basis.temp <- crossbasis(data_sub$Daily_mean_temperature_XAVIER, lag= lagn, 
                           argvar=list(fun="ns",df=3),
                           arglag=list(fun="ns",knots=(logknots(lagn,lag_nk))))
  
  form<-as.formula(paste(outcomes[i],"~ basis.temp + ns(absumid_xavier,3) + feriado"))
  
  mod <- gnm(form, data=data_sub, family = quasipoisson(), eliminate=factor(stratum))
  
  cp <- crosspred(basis.temp, mod)
  mmt=cp$predvar[which.min(cp$allRRfit)]
  mmp=sum(na.omit(cp$predvar) < mmt) / length(na.omit(cp$predvar))
  
  cp2 <- crosspred(basis.temp, mod,cen=mmt) #adicionei aqui pra tentar checar as predicoes
  
  RR.eh <- unlist(crosspred(basis.temp, mod, at = quantile(cp2$predvar, 0.99,na.rm = T), cen = mmt)[c("allRRfit", "allRRlow", "allRRhigh")]) 
  RR.ec <- unlist(crosspred(basis.temp, mod, at = quantile(cp2$predvar, 0.01,na.rm = T), cen = mmt)[c("allRRfit", "allRRlow", "allRRhigh")])
  
  RR.table <- cbind(t(RR.eh),t(RR.ec),mmt,mmp)
  colnames(RR.table) <- c("RR.eh", "RR.eh.low", "RR.eh.high",
                          "RR.ec", "RR.ec.low", "RR.ec.high","mmt","mmp")
  
  RR.table<- round(RR.table,2)
  
  RR.table<-cbind("Central-West",outcomes[i],RR.table)
  
  result=rbind(result,RR.table)      
  
  y.lim.max<-round(max(RR.eh[1],RR.ec[1])+2,0)
  
  plot(cp2,"overall",ylim=c(0,y.lim.max),yaxt="n",lab=c(6,5,7),xlab="Temperature (°C)",ylab="RR",cex.lab=0.9,cex.axis = 0.9,main="Central-West", cex=1.0)
  
  axis(2)
  
  abline(v=mmt,lty=3)
  abline(v=c(quantile(cp2$predvar,c(0.01,0.99))),lty=2) ## We used predvar instead temperature variable 
  
rm(list = ls()[!ls() %in% c("data_under5","result","cb_base","cb_per","lagn","lag_nk","outcomes2")])


#-------------------------------------NORTH REGION--------------------------------------------------------

data_region <- subset(data_under5, Sigla.da.UF %in% c("AC","AM","AP","PA","RO","RR","TO"))

outcomes<-names(data_region)[c(9,6,7,3)] 

i=4 
  
  data_sub <- data_region 
  
  data_sub <- data_sub[, count_cases := sum(get(outcomes[i])), by = stratum]
  
  data_sub=subset(data_sub, data_sub$count_cases>0)
  
  basis.temp <- crossbasis(data_sub$Daily_mean_temperature_XAVIER, lag= lagn, 
                           argvar=list(fun="ns",df=3),
                           arglag=list(fun="ns",knots=(logknots(lagn,lag_nk))))
  
  form<-as.formula(paste(outcomes[i],"~ basis.temp + ns(absumid_xavier,3) + feriado"))
  
  mod <- gnm(form, data=data_sub, family = quasipoisson(), eliminate=factor(stratum))
  
  cp <- crosspred(basis.temp, mod)
  mmt=cp$predvar[which.min(cp$allRRfit)]
  mmp=sum(na.omit(cp$predvar) < mmt) / length(na.omit(cp$predvar))
  
  cp2 <- crosspred(basis.temp, mod,cen=mmt) 
  
  RR.eh <- unlist(crosspred(basis.temp, mod, at = quantile(cp2$predvar, 0.99,na.rm = T), cen = mmt)[c("allRRfit", "allRRlow", "allRRhigh")]) 
  RR.ec <- unlist(crosspred(basis.temp, mod, at = quantile(cp2$predvar, 0.01,na.rm = T), cen = mmt)[c("allRRfit", "allRRlow", "allRRhigh")])
  
  RR.table <- cbind(t(RR.eh),t(RR.ec),mmt,mmp)
  colnames(RR.table) <- c("RR.eh", "RR.eh.low", "RR.eh.high",
                          "RR.ec", "RR.ec.low", "RR.ec.high","mmt","mmp")
  
  RR.table<- round(RR.table,2)
  
  RR.table<-cbind("North",outcomes[i],RR.table)
  
  result=rbind(result,RR.table)
  
  y.lim.max<-round(max(RR.eh[1],RR.ec[1])+2,0)
  
  plot(cp2,"overall",ylim=c(0,y.lim.max),yaxt="n",lab=c(6,5,7),xlab="Temperature (°C)",ylab="RR",cex.lab=0.9,cex.axis = 0.9,main="North", cex=1.0)
  
  axis(2)
  
  abline(v=mmt,lty=3)
  abline(v=c(quantile(cp2$predvar,c(0.01,0.99))),lty=2) ## We used predvar instead temperature variable 
  
rm(list = ls()[!ls() %in% c("data_under5","result","cb_base","cb_per","lagn","lag_nk","outcomes2")])

#-------------------------------------NORTHEAST REGION--------------------------------------------------------

data_region <- subset(data_under5, Sigla.da.UF %in% c("AL","BA","CE","MA","PB","PE","PI","RN","SE"))

outcomes<-names(data_region)[c(9,6,7,3)] 

i=4 
  
  data_sub <- data_region 
  
  data_sub <- data_sub[, count_cases := sum(get(outcomes[i])), by = stratum]
  
  data_sub=subset(data_sub, data_sub$count_cases>0)
  
  basis.temp <- crossbasis(data_sub$Daily_mean_temperature_XAVIER, lag= lagn, 
                           argvar=list(fun="ns",df=3),
                           arglag=list(fun="ns",knots=(logknots(lagn,lag_nk))))
  
  form<-as.formula(paste(outcomes[i],"~ basis.temp + ns(absumid_xavier,3) + feriado"))
  
  mod <- gnm(form, data=data_sub, family = quasipoisson(), eliminate=factor(stratum))
  
  cp <- crosspred(basis.temp, mod)
  mmt=cp$predvar[which.min(cp$allRRfit)]
  mmp=sum(na.omit(cp$predvar) < mmt) / length(na.omit(cp$predvar))
  
  cp2 <- crosspred(basis.temp, mod,cen=mmt) 
  
  RR.eh <- unlist(crosspred(basis.temp, mod, at = quantile(cp2$predvar, 0.99,na.rm = T), cen = mmt)[c("allRRfit", "allRRlow", "allRRhigh")]) 
  RR.ec <- unlist(crosspred(basis.temp, mod, at = quantile(cp2$predvar, 0.01,na.rm = T), cen = mmt)[c("allRRfit", "allRRlow", "allRRhigh")])
  
  RR.table <- cbind(t(RR.eh),t(RR.ec),mmt,mmp)
  colnames(RR.table) <- c("RR.eh", "RR.eh.low", "RR.eh.high",
                          "RR.ec", "RR.ec.low", "RR.ec.high","mmt","mmp")
  
  RR.table<- round(RR.table,2)
  
  RR.table<-cbind("Northeast",outcomes[i],RR.table)
  
  result=rbind(result,RR.table)
  
  y.lim.max<-round(max(RR.eh[1],RR.ec[1])+2,0)
  
  plot(cp2,"overall",ylim=c(0,y.lim.max),yaxt="n",lab=c(6,5,7),xlab="Temperature (°C)",ylab="RR",cex.lab=0.9,cex.axis = 0.9,main="Northeast", cex=1.0)
  
  axis(2)
  
  abline(v=mmt,lty=3)
  abline(v=c(quantile(cp2$predvar,c(0.01,0.99))),lty=2) ## We used predvar instead temperature variable 
  

rm(list = ls()[!ls() %in% c("data_under5","result","cb_base","cb_per","lagn","lag_nk","outcomes2")])

#-------------------------------------SOUTH REGION--------------------------------------------------------

data_region <- subset(data_under5, Sigla.da.UF %in% c("PR","RS","SC"))

outcomes<-names(data_region)[c(9,6,7,3)] 

i=4 
  
  data_sub <- data_region 
  
  data_sub <- data_sub[, count_cases := sum(get(outcomes[i])), by = stratum]
  
  data_sub=subset(data_sub, data_sub$count_cases>0)
  
  basis.temp <- crossbasis(data_sub$Daily_mean_temperature_XAVIER, lag= lagn, 
                           argvar=list(fun="ns",df=3),
                           arglag=list(fun="ns",knots=(logknots(lagn,lag_nk))))
  
  form<-as.formula(paste(outcomes[i],"~ basis.temp + ns(absumid_xavier,3) + feriado"))
  
  mod <- gnm(form, data=data_sub, family = quasipoisson(), eliminate=factor(stratum))
  
  cp <- crosspred(basis.temp, mod)
  mmt=cp$predvar[which.min(cp$allRRfit)]
  mmp=sum(na.omit(cp$predvar) < mmt) / length(na.omit(cp$predvar))
  
  cp2 <- crosspred(basis.temp, mod,cen=mmt) 
  
  RR.eh <- unlist(crosspred(basis.temp, mod, at = quantile(cp2$predvar, 0.99,na.rm = T), cen = mmt)[c("allRRfit", "allRRlow", "allRRhigh")]) 
  RR.ec <- unlist(crosspred(basis.temp, mod, at = quantile(cp2$predvar, 0.01,na.rm = T), cen = mmt)[c("allRRfit", "allRRlow", "allRRhigh")])
  
  RR.table <- cbind(t(RR.eh),t(RR.ec),mmt,mmp)
  colnames(RR.table) <- c("RR.eh", "RR.eh.low", "RR.eh.high",
                          "RR.ec", "RR.ec.low", "RR.ec.high","mmt","mmp")
  
  RR.table<- round(RR.table,2)
  
  RR.table<-cbind("South",outcomes[i],RR.table)
  
  result=rbind(result,RR.table)      
  
  y.lim.max<-round(max(RR.eh[1],RR.ec[1])+2,0)
  
  cp2 <- crosspred(basis.temp, mod,cen=mmt)
  
  plot(cp2,"overall",ylim=c(0,y.lim.max),yaxt="n",lab=c(6,5,7),xlab="Temperature (°C)",ylab="RR",cex.lab=0.9,cex.axis = 0.9,main="South", cex=1.0)
  
  axis(2)
  
  abline(v=mmt,lty=3)
  abline(v=c(quantile(cp2$predvar,c(0.01,0.99))),lty=2) ## We used predvar instead temperature variable
  

rm(list = ls()[!ls() %in% c("data_under5","result","cb_base","cb_per","lagn","lag_nk","outcomes2")])

#-------------------------------------SOUTHEAST REGION--------------------------------------------------------

data_region <- subset(data_under5, Sigla.da.UF %in% c("ES","MG","RJ","SP"))

outcomes<-names(data_region)[c(9,6,7,3)] 

i=4 
  
  data_sub <- data_region 
  
  data_sub <- data_sub[, count_cases := sum(get(outcomes[i])), by = stratum]
  
  data_sub=subset(data_sub, data_sub$count_cases>0)
  
  basis.temp <- crossbasis(data_sub$Daily_mean_temperature_XAVIER, lag= lagn, 
                           argvar=list(fun="ns",df=3),
                           arglag=list(fun="ns",knots=(logknots(lagn,lag_nk))))
  
  form<-as.formula(paste(outcomes[i],"~ basis.temp + ns(absumid_xavier,3) + feriado"))
  
  mod <- gnm(form, data=data_sub, family = quasipoisson(), eliminate=factor(stratum))
  
  cp <- crosspred(basis.temp, mod)
  mmt=cp$predvar[which.min(cp$allRRfit)]
  mmp=sum(na.omit(cp$predvar) < mmt) / length(na.omit(cp$predvar))
  
  cp2 <- crosspred(basis.temp, mod,cen=mmt) #adicionei aqui pra tentar checar as predicoes
  
  RR.eh <- unlist(crosspred(basis.temp, mod, at = quantile(cp2$predvar, 0.99,na.rm = T), cen = mmt)[c("allRRfit", "allRRlow", "allRRhigh")]) 
  RR.ec <- unlist(crosspred(basis.temp, mod, at = quantile(cp2$predvar, 0.01,na.rm = T), cen = mmt)[c("allRRfit", "allRRlow", "allRRhigh")])
  
  RR.table <- cbind(t(RR.eh),t(RR.ec),mmt,mmp)
  colnames(RR.table) <- c("RR.eh", "RR.eh.low", "RR.eh.high",
                          "RR.ec", "RR.ec.low", "RR.ec.high","mmt","mmp")
  
  RR.table<- round(RR.table,2)
  
  RR.table<-cbind("Southeast",outcomes[i],RR.table)
  
  result=rbind(result,RR.table)      
  
  y.lim.max<-round(max(RR.eh[1],RR.ec[1])+2,0)
  
  cp2 <- crosspred(basis.temp, mod,cen=mmt)
  
  plot(cp2,"overall",ylim=c(0,y.lim.max),yaxt="n",lab=c(6,5,7),xlab="Temperature (°C)",ylab="RR",cex.lab=0.9,cex.axis = 0.9,main="Southeast", cex=1.0)
  
  axis(2)
  
  abline(v=mmt,lty=3)
  abline(v=c(quantile(cp2$predvar,c(0.01,0.99))),lty=2) ## We used predvar instead temperature variable 
  

rm(list = ls()[!ls() %in% c("data_under5","result","cb_base","cb_per","lagn","lag_nk","outcomes2")])





#-------------------------------------ANALYSIS BRAZIL BY AGE GROUP--------------------------------------------------------

data_sub=data_under5

outcomes<-names(data_sub)[c(9,6,7,3)]


for (i in 1:3) {  
  
  data_sub <- data_under5
  
  data_sub <- data_sub[, count_cases := sum(get(outcomes[i])), by = stratum]
  
  data_sub=subset(data_sub, data_sub$count_cases>0)
  
  basis.temp <- crossbasis(data_sub$Daily_mean_temperature_XAVIER, lag= lagn, 
                           argvar=list(fun="ns",df=3),
                           arglag=list(fun="ns",knots=(logknots(lagn,lag_nk))))
  
  form<-as.formula(paste(outcomes[i],"~ basis.temp + ns(absumid_xavier,3) + feriado"))
  
  mod <- gnm(form, data=na.omit(data_sub), family = quasipoisson(), eliminate=factor(stratum))
  
  cp <- crosspred(basis.temp, mod)
  mmt=cp$predvar[which.min(cp$allRRfit)]
  mmp=sum(na.omit(cp$predvar) < mmt) / length(na.omit(cp$predvar))
  
  cp2 <- crosspred(basis.temp, mod,cen=mmt) 
  
  RR.eh <- unlist(crosspred(basis.temp, mod, at = quantile(cp2$predvar, 0.99,na.rm = T), cen = mmt)[c("allRRfit", "allRRlow", "allRRhigh")]) 
  RR.ec <- unlist(crosspred(basis.temp, mod, at = quantile(cp2$predvar, 0.01,na.rm = T), cen = mmt)[c("allRRfit", "allRRlow", "allRRhigh")])
  
  RR.table <- cbind(t(RR.eh),t(RR.ec),mmt,mmp)
  colnames(RR.table) <- c("RR.eh", "RR.eh.low", "RR.eh.high",
                          "RR.ec", "RR.ec.low", "RR.ec.high","mmt","mmp")
  
  RR.table<- round(RR.table,2)
  
  RR.table<-cbind("Brazil",outcomes[i],RR.table)
  
  result=rbind(result,RR.table)      
  
  y.lim.max<-round(max(RR.eh[1],RR.ec[1])+2,0)
  
  cp2 <- crosspred(basis.temp, mod,cen=mmt)
  
  plot(cp2,"overall",ylim=c(0,y.lim.max),yaxt="n",lab=c(6,5,7),xlab="Temperature (°C)",ylab="RR",cex.lab=0.9,cex.axis = 0.9,main=outcomes2[i], cex=1.0)
  
  axis(2)
  
  abline(v=mmt,lty=3)
  abline(v=c(quantile(cp2$predvar,c(0.01,0.99))),lty=2) ## We used predvar instead temperature variable 
  
  mtext("Brazil", cex=0.8)
  
}

rm(list = ls()[!ls() %in% c("data_under5","result","cb_base","cb_per","lagn","lag_nk","outcomes2")])


#-------------------------------------ANALYSIS BRAZIL BY CAUSE OF DEATH--------------------------------------------------------
data_under5=fread("/home/ismael.silveira/mortality_under5/dataset/data_under5_sub_causes_race.csv") 
head(data_under5)
names(data_under5)

# especificação da base
cb_base='ns'
lagn=7
lag_nk=2

outcomes2<-c("Infectious and parasitic diseases", "Diarrhea", "Undernutrition", "Malnutrition","Respiratory diseases" )

data_sub=data_under5

outcomes<-names(data_sub)[c(8:12)]


for (i in 1:length(outcomes)) {  
  
  data_sub <- data_under5
  
  data_sub <- data_sub[, count_cases := sum(get(outcomes[i])), by = stratum]
  
  data_sub=subset(data_sub, data_sub$count_cases>0)
  
  basis.temp <- crossbasis(data_sub$Daily_mean_temperature_XAVIER, lag= lagn, 
                           argvar=list(fun="ns",df=3),
                           arglag=list(fun="ns",knots=(logknots(lagn,lag_nk))))
  
  form<-as.formula(paste(outcomes[i],"~ basis.temp + ns(absumid_xavier,3) + feriado"))
  
  mod <- gnm(form, data=na.omit(data_sub), family = quasipoisson(), eliminate=factor(stratum))
  
  cp <- crosspred(basis.temp, mod)
  mmt=cp$predvar[which.min(cp$allRRfit)]
  mmp=sum(na.omit(cp$predvar) < mmt) / length(na.omit(cp$predvar))
  
  cp2 <- crosspred(basis.temp, mod,cen=mmt) #adicionei aqui pra tentar checar as predicoes
  
  RR.eh <- unlist(crosspred(basis.temp, mod, at = quantile(cp2$predvar, 0.99,na.rm = T), cen = mmt)[c("allRRfit", "allRRlow", "allRRhigh")]) 
  RR.ec <- unlist(crosspred(basis.temp, mod, at = quantile(cp2$predvar, 0.01,na.rm = T), cen = mmt)[c("allRRfit", "allRRlow", "allRRhigh")])
  
  RR.table <- cbind(t(RR.eh),t(RR.ec),mmt,mmp)
  colnames(RR.table) <- c("RR.eh", "RR.eh.low", "RR.eh.high",
                          "RR.ec", "RR.ec.low", "RR.ec.high","mmt","mmp")
  
  RR.table<- round(RR.table,2)
  
  RR.table<-cbind("Brazil",outcomes[i],RR.table)
  
  result=rbind(result,RR.table)      
  
  y.lim.max<-round(max(RR.eh[1],RR.ec[1])+2,0)
  
  cp2 <- crosspred(basis.temp, mod,cen=mmt)
  
  plot(cp2,"overall",ylim=c(0,y.lim.max),yaxt="n",lab=c(6,5,7),xlab="Temperature (°C)",ylab="RR",cex.lab=0.9,cex.axis = 0.9,main=outcomes2[i], cex=1.0)
  
  axis(2)
  
  abline(v=mmt,lty=3)
  abline(v=c(quantile(cp2$predvar,c(0.01,0.99))),lty=2) ## We used predvar instead temperature variable 
  
  #mtext("Brazil", cex=0.8)
  
}

dev.off()

fwrite(result,"/home/ismael.silveira/mortality_under5/main_and_additional_analysis/table_result_test_for_github.csv")



