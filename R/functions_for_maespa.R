#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
# function required to change parameters in MAESPA
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------

# function for create tree stand for MAESPA

make_stand<-function(){
  
  inter.row<-rep(seq(2.6,27,2.6),times=10) #create 100 trees in a row with 2.6m spacing between two trees
  intra.row<-rep(seq(3.85,40,3.85),each=10)#create 100 trees in a column with 3.85m spacing between two trees
  
  tree.stand<-data.frame(inter.row,intra.row) #get tree stand. 100 trees with 2.6*3.85m spacing
  names(tree.stand)<-NULL
  return(tree.stand)
}

#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------

# function to copy and paste original MAESPA input files in MAESPA directory for each simulation

startover <- function(){
  file.copy(dir("original_inputfiles", full.names=TRUE),"maespa",overwrite=TRUE)
}
cleanup <- function()unlink(Sys.glob("*.dat"))


#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------

# function to make structure file for MAESPA

make_str <- function(){
  
  replacePAR("maespa/str.dat", "elp", "lia", 1.0)  # spherical LIA
  replacePAR("maespa/str.dat", "cshape", "canopy", "CYL")  # canopy shape
  
}

#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------

# function to read and process WTC3 met data for MAESPA simulations

make_met <- function(ch){
  met <- subset(wtcflux, chamber == ch, select=c(DateTime, Tair_al, VPDair, PAR))
  met <- met[order(met$DateTime),-1]  # order by DateTime; remove DateTime
  
  # Fill NA with last non-NA
  met <- within(met, {
    Tair_al <- na.locf(Tair_al)
    VPDair <- na.locf(VPDair)
    PAR <- na.locf(PAR)
  })
  
  # CO2 must be added; data in wtcflux too noisy
  met$CA <- 400
  
  # Units VPD are Pa
  met$VPDair <- 1000 * met$VPDair
  
  replacemetdata(met, "maespa/met.dat", columns=c("TAIR","VPD","PAR","CA"), 
                 "maespa/met.dat", setdates=FALSE, khrs=24)
  replacePAR("maespa/met.dat", "khrsperday", "metformat", 24)
  replacePAR("maespa/met.dat", "startdate", "metformat", format(startDate, "%d/%m/%y"))
  replacePAR("maespa/met.dat", "enddate", "metformat", format(endDate, "%d/%m/%y"))
}


#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------

# finction to make tree stand for MAESPA. The trees file


make_trees <- function(ch){
  
  hd <- subset(treeh, chamber == ch)
  replacePAR("maespa/trees.dat","values","allhtcrown",round(hd$Plant_height / 100,2))
  replacePAR("maespa/trees.dat","nodates","allhtcrown",nrow(hd))
  replacePAR("maespa/trees.dat","dates","allhtcrown",format(hd$Date,"%d/%m/%y"))
  
  replacePAR("maespa/trees.dat","values","allradx",round(hd$crownwidth / 100 / 2, 2))
  replacePAR("maespa/trees.dat","nodates","allradx",nrow(hd))
  replacePAR("maespa/trees.dat","dates","allradx",format(hd$Date,"%d/%m/%y"))
  
  replacePAR("maespa/trees.dat","values","allrady",round(hd$crownwidth / 100 / 2, 2))
  replacePAR("maespa/trees.dat","nodates","allrady",nrow(hd))
  replacePAR("maespa/trees.dat","dates","allrady",format(hd$Date,"%d/%m/%y"))
  
  replacePAR("maespa/trees.dat", "values", "allhttrunk", 0.01)
  
  replacePAR("maespa/trees.dat", "values", "alldiam", round(hd$diameter/1000,2)) #convert mm to m
  replacePAR("maespa/trees.dat","nodates","alldiam",nrow(hd))
  replacePAR("maespa/trees.dat","dates","alldiam",format(hd$Date,"%d/%m/%y"))
  
  replacePAR("maespa/trees.dat", "values", "alllarea", round(hd$leafArea,4))
  replacePAR("maespa/trees.dat","nodates","alllarea",nrow(hd))
  replacePAR("maespa/trees.dat","dates","alllarea",format(hd$Date,"%d/%m/%y"))
  
  # replacePAR("maespa/trees.dat", "xycoords", "xy", make_stand())
}


#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------

# function to set simulation control parameters. The control file

make_con <- function(parlayers,nopoints){
  
  replacePAR("maespa/confile.dat", "startdate", "dates", format(startDate, "%d/%m/%y"))
  replacePAR("maespa/confile.dat", "enddate", "dates", format(endDate, "%d/%m/%y"))
  replacePAR("maespa/confile.dat","iohrly","control",1)
  # Use Medlyn 2011 model of gs.
  replacePAR("maespa/confile.dat", "modelgs", "model", 4)
  replacePAR("maespa/confile.dat", "modelrw", "model", 0)
  replacePAR("maespa/confile.dat", "nolay", "diffang", parlayers)
  replacePAR("maespa/confile.dat", "pplay", "diffang", nopoints)
  
  }


#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------

# function to set physiological parameters. The Physiology file

make_phy <- function(jmaxnodates,jmaxnolayers,jmax,jmaxdates,vcmaxnodates,vcmaxnolayers,vcmax,vcmaxdates,
                     theta,eavj,delsj,ajq,eavc,delsc,rdnodates,rday,rddates,dayresp=1,q10ndates,q10,q10dates,g1){
  
  replacePAR("maespa/phy.dat","g1","bbmgs",g1) 
  replacePAR("maespa/phy.dat","dates","bbmgs","14/09/13")
  # replacePAR("maespa/phy.dat","nolayers","bbmgs",g1layers)
  
  
  # Jmax parameters
  
  replacePAR("maespa/phy.dat","nodates","jmaxcon",jmaxnodates)
  replacePAR("maespa/phy.dat","nolayers","jmaxcon",jmaxnolayers)
  replacePAR("maespa/phy.dat","values","jmax",jmax) 
  replacePAR("maespa/phy.dat","dates","jmax",jmaxdates)                         
  
  
  # Vcmax parameters
  
  replacePAR("maespa/phy.dat","nodates","vcmaxcon",vcmaxnodates)
  replacePAR("maespa/phy.dat","nolayers","vcmaxcon",vcmaxnolayers) 
  replacePAR("maespa/phy.dat","values","vcmax",vcmax)
  replacePAR("maespa/phy.dat","dates","vcmax",vcmaxdates)
  
  
  # Temperature response of Vcmax and Jmax
  
  replacePAR("maespa/phy.dat","theta","jmaxpars",theta) #all parameters are for summer 2013
  replacePAR("maespa/phy.dat","eavj","jmaxpars",eavj)
  
  replacePAR("maespa/phy.dat","delsj","jmaxpars",delsj)
  replacePAR("maespa/phy.dat","ajq","jmaxpars",ajq)
  
  replacePAR("maespa/phy.dat","eavc","vcmaxpars",eavc)
  replacePAR("maespa/phy.dat","delsc","vcmaxpars",delsc)
  
  #Rday
  replacePAR("maespa/phy.dat","nodates","rdcon",rdnodates) 
  replacePAR("maespa/phy.dat","values","rd",rday) 
  replacePAR("maespa/phy.dat","dates","rd",rddates) 
  replacePAR("maespa/phy.dat","dayresp","rdpars",dayresp)
  
  # Tresponse of Rday
  replacePAR("maespa/phy.dat","rates","folq10",q10)
  replacePAR("maespa/phy.dat","ndates","folq10",q10ndates)
  replacePAR("maespa/phy.dat","dates","folq10",q10dates)
  


  
}

#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------

# function to run MAESPA for a given chamber
# change temperature response parameters using make_phy function for each run
# so the run_chamber function should update for each simulation 

run_chamber <- function(ch,parlayers=6,nopoints=12,jmaxnodates,jmaxnolayers,jmax,jmaxdates,vcmaxnodates,vcmaxnolayers,vcmax,vcmaxdates,
                        theta,eavj,delsj,ajq,eavc,delsc,rdnodates,rday,rddates,dayresp=dayresp,q10ndates,q10,q10dates,g1){
  o <- getwd()
  on.exit(setwd(o))
  
  startover()
  make_str()
  make_con(parlayers=parlayers,nopoints=nopoints)
  
  # The physiology file
  make_phy(jmaxnodates,jmaxnolayers,jmax,jmaxdates,vcmaxnodates,vcmaxnolayers,vcmax,vcmaxdates,
           theta,eavj,delsj,ajq,eavc,delsc,rdnodates,rday,rddates,dayresp=1,q10ndates,q10,q10dates,g1)
  
  # The trees file  
  make_trees(ch)
  
  # The met file
  make_met(ch)
  
  # run MAESPA
  setwd("maespa")
  shell("MAESPA64.exe", intern=TRUE)
  # h <- readlayflux()
  
  # merge with observed flux 
  # this section has to remove for Kashif's simulations
  
  #---------------------------------------------------------------------
  h <- readhrflux()
  message(sprintf("Chamber %s simulation complete.",ch))
  #---------------------------------------------------------------------
  
  h$chamber <- ch
  h$DateTime <- seq(from=as.POSIXct(as.character(startDate), tz="UTC"), length=nrow(h), by="1 hour")
  return(h)
}


#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------

# Function to read MAESPA layer output

run_chamber_layers <- function(ch,jmaxnodates,jmaxnolayers,jmax,jmaxdates,vcmaxnodates,vcmaxnolayers,vcmax,vcmaxdates,
                               theta,eavj,delsj,ajq,eavc,delsc,rdnodates,rday,rddates,dayresp=dayresp,q10ndates,q10,q10dates){
  o <- getwd()
  on.exit(setwd(o))
  
  startover()
  
  make_str()
  make_con()
  
  
  make_phy(jmaxnodates,jmaxnolayers,jmax,jmaxdates,vcmaxnodates,vcmaxnolayers,vcmax,vcmaxdates,
           theta,eavj,delsj,ajq,eavc,delsc,rdnodates,rday,rddates,dayresp=1,q10ndates,q10,q10dates)
  
  make_trees(ch)
  make_met(ch)
  
  setwd("maespa")
  shell("MAESPA64.exe", intern=TRUE)
  h <- readlayflux()
  h<-h[with(h, order(Date, Hour)),]
  
  
  # h <- readhrflux()
  message(sprintf("Chamber %s simulation complete.",ch))
  
  h$chamber <- ch
  h$DateTime <- seq(from=as.POSIXct(as.character(startDate), tz="UTC"), length=nrow(h), by="1 hour")
  return(h)
}


#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------


# Function to merge with measured chamber flux

merge_flux<-function(modelledflux){
  
  flux_control<-subset(wtcflux,wtcflux$chamber %in% c(insidechambers_con))
  df <- merge(flux_control, modelledflux, by=c("chamber","DateTime"), all=TRUE)
  # df <- cbind(flux_control, modelledflux)
  
  
  #get flux per leaf area
  
  df$hrPs<-df$hrPs-df$hrRmW #(substract wood respiration before convert flux per leaf area basis) 
  # df$hrPs<-df$hrPs-df$hrRf #(substract wood respiration assuming wood respiration=leaf respiration) 
  df$hrPs_la<-with(df,hrPs/leafArea) #units for MAESPA hourly flux data is mumol/tree/s. convert it to mumol/m2 (LA)/s
  
  df$hrPs_la_m<-with(df,FluxCO2*1000/leafArea) #units for wtc3 hourly flux data is mmol/s. convert it to mumol/m2 (LA)/s
  
  df$residuals<-with(df,hrPs_la_m-hrPs_la)
  
  return(df)
}

# df<-merge_flux(maes_con)
# # get flux at high PAR to plot the temperature response
# 
# dat.plot.a<-subset(df,df$PAR.x>1200)

#--------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------

# plot measured vs modelled flux

plot_flux<-function(Title){
  
  
  # windows(40,40);par(mfrow=c(1,1),mar=c(0,0,0,0),oma=c(5,7,2,2));par(las=1)
  with(df,plot(hrPs_la_m~hrPs_la,xlim=c(-5,25),ylim=c(-5,25),col=alpha("lightgray",0.1),xlab="",ylab="",axes=F))
  abline(0,1)
  
  magaxis(side=c(1:4),labels=c(1,1,0,0),frame.plot=T,las=1,cex.axis=1.3,ratio=0.4,tcl=0.5,majorn=3)
  
  lmod<-summary(lm(hrPs_la_m~hrPs_la,data=df))
  
  abline(lmod$coefficients[1],lmod$coefficients[2]
         ,col="black",lwd=3,lty=2)
  
  title(xlab=expression(A[modelled]~(mu*mol~m^2~s^-1)),outer=F,cex.lab=2)
  title(ylab=expression(A[measured]~(mu*mol~m^2~s^-1)),outer=F,cex.lab=2,line=1.7)
  
  
  mylabel_1 = bquote(italic(r)^2 == .(format(lmod$adj.r.squared, digits = 2)))
  mylabel_2 = bquote(italic(y) == .(format(lmod$coefficients[2], digits = 2))~italic(x) + .(format(lmod$coefficients[1], digits = 2)))
  
  # text(x = -1.0, y = 19, labels = mylabel_2)
  # text(x = -2.5, y = 18, labels = mylabel_1)
  
  text(x = 14, y = -2, labels = mylabel_2,cex=1.5)
  text(x = 11.5, y = -3.5, labels = mylabel_1,cex=1.5)
  
  
  # 
  title(main=Title,line=0.5,adj=.05,cex.main=1.6)
  
}

#----- ---------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------

# plot temperature response of canopy flux

plot_tresponse_wtc<-function(Legend=NULL){
  
  # windows(40,40);par(mfrow=c(1,1),mar=c(0,0,0,0),oma=c(5,7,2,2));par(las=1)
  smoothplot(TCAN,hrPs_la,pointcols=c(alpha("black",0.05)),
             linecol=c("black"),polycolor=c(alpha("gray",0.9))
             ,cex=1,main="",random="chamber",
             xlim=c(15,42),ylim=c(-5,25),xlab="",ylab="",
             data=dat.plot.a, kgam=4,axes=F)
  
  
  # box();axis(side=1,labels=T);axis(side=2,labels=T)
  magaxis(side=c(1:4),labels=c(1,1,0,0),frame.plot=T,las=1,cex.axis=1.3,ratio=0.4,tcl=0.5,majorn=3)
  
  
  title(xlab=expression(T[Canopy]~(degree*C)),outer=F,cex.lab=2)
  title(ylab=expression(A~(mu*mol~m^2~s^-1)),outer=F,cex.lab=2,line=1.7)
  
  if(is.null(Legend)==T){
    legend("topright",cex=1.5,pch=16,legend=c("measured","modelled"),col=c(alpha("red",0.5),alpha("black",0.5)),bty="n")
    
  }
  
  # else{legend("topright",legend=NULL,bty="n")}
  # title(main="AMBIENT_CHAMBERS",outer=T,cex.lab=2)
  
  par(new=T)
  
  #plot measured photosynthesis
  
  smoothplot(TCAN,hrPs_la_m,pointcols=c(alpha("red",0.05)),
             linecol=c("red"),polycolor=c(alpha("gray",0.5)),cex=1,main="",
             xlim=c(15,42),ylim=c(-5,25),xlab="",ylab="",
             data=dat.plot.a, kgam=4,axes=F,random="chamber")
  
}

# windows(100,60);par(mar=c(5.5,5.5,2,1),oma=c(0,0,0,0),mfrow=c(1,2));par(las=1)
# plot_flux()
# plot_tresponse_wtc()

#----------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------
