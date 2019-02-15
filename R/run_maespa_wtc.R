
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
# baseline MAESPA for wtc3

# parameters specified for one crown layer
# No seasonal acclimation
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------

# get functions and packages

source("R/data_processing.R")
source("R/functions_for_maespa.R")
source("R/loadLibraries.R")


# chamber codes. Only ambient-T treatment

insidechambers_con <- sprintf("C%02.f",(c(1,3,5,7,9,11)))


# start and end dates of simulation 

startDate <- as.Date("2013-09-14")
endDate <- as.Date("2014-05-26")


# Read WTC flux data
wtcflux <- make_or_read("cache/wtcflux.rds", make_wtcflux)
wtcflux$chamber<-factor(wtcflux$chamber)
wtcflux$DateTime<-as.POSIXct(wtcflux$DateTime,format="%d/%m/%Y %H:%M",tz="UTC")


wtcflux$Date<-as.Date(wtcflux$Date,format="%d/%m/%Y")
wtcflux$PAR<-ifelse(wtcflux$PAR<0,0,wtcflux$PAR) 

wtcflux<-subset(wtcflux,Date >= startDate & Date <= endDate)


# get tree height, crown width data 
treeh <-make_treeh()

# treeh <-readRDS(file="cache/treeh.rds")
treeh<-subset(treeh,Date >= startDate & Date <= endDate)


# maespa simulation
# parameter definitios
#' jmaxnodates: no of dates jmax specify
#' jmaxnolayers: no of crown layers jmax specify
#' jmaxdates: dates jmax specify
#' vcmaxnodates, vcmaxnolayers,vcmaxdates: same as jmax
#' rdnodates: no of dates day respiration rate specify
#' rddates: dates day respiration rate specify
#' theta: theta of electron transport
#' eavj:activation energy of jmax
#' delsj:  entropy of jmax
#' eavc: activation energy of vcmax
#' delsc:  entropy of vcmax
#' rday: respiration rate at 25C
#' dayresp: ratio of day respiration:dark respiration
#' q10: q10 of respiration
#' g1: g1 (Medlyn et al) 


maes_3 <-do.call(rbind,lapply(insidechambers_con, function(x)run_chamber(x,jmaxnodates=1,jmaxnolayers=1,
                                                                         jmax=c(178.3),
                                                                         jmaxdates=c("14/09/13"),
                                                                         vcmaxnodates=1,vcmaxnolayers=1,
                                                                         vcmax=c(103.6),
                                                                         vcmaxdates=c("14/09/13"),rdnodates=1,rddates=c("14/09/13"),
                                                                         theta=0.57,eavj=23800,delsj=627,ajq=0.26,eavc=59700,delsc=634,
                                                                         rday=c(1.4),dayresp=0.7,q10ndates=1,
                                                                         q10=c(0.07),q10dates=c("14/09/13"),g1=2.4)))

# # merge measured flux with modelled flux
# # convert units 
# # substract wood respiration

df<-merge_flux(maes_3)

# # get subset of high PAR (PAR>1200 for temperature response)

dat.plot.a<-subset(df,df$PAR.x>1200) # get subset of high PAR 

# plot measured vs modelled flux

plot_flux(Title=NULL)
legend("topleft",expression((bold(e))),bty="n",cex=1.5)

# plot temperature response 

plot_tresponse_wtc(Legend=F)
legend("topleft",expression((bold(f))),bty="n",cex=1.5)

legend("topright",cex=1.5,pch=16,legend=c("measured","modelled"),col=c(alpha("red",0.5),alpha("black",0.5)),bty="n")


#----------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------

#