
# Functions required for simulations


#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------

# function to install MAESPA

install_maespa <- function(topath=NULL){
  
  if(is.null(topath))stop("Specify installation path.")
  
  download.file("http://maespa.github.io/maespa/MAESPAV2.zip", 
                "tmp.zip",
                mode="wb")
  u <- unzip("tmp.zip", "MAESPA64.exe", exdir=topath)
  
}

#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------


# function for read and process WTC3 flux data

make_wtcflux <- function(){
  wtcflux <- downloadCSV("WTC_TEMP_CM_WTCFLUX_20130910-20140530_L2_V1.csv")
  
  # Problem : 5 extra rows in 2014-1-19 at midnight (timezone problem or something??)
  wtcflux <- summaryBy(. ~ DateTime + chamber, FUN=mean, data=wtcflux, keep.names=TRUE)
  
  wtcflux$Date <- as.Date(wtcflux$DateTime)
  
  # toss last day; has only one observation
  wtcflux <- subset(wtcflux, Date < as.Date("2014-5-27"))
  
  return(wtcflux)
}


#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------

# function for read a given dataset if available in WD or create/download

make_or_read <- function(fn, fun){
  
  if(!file.exists(fn)){
    out <- fun()
    saveRDS(out, fn)
  } else {
    out <- readRDS(fn)
  }
  
  return(out)
}


#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------

# function for read and process WTC3 tree data (height, diameter, crown width)
# 

make_treeh <- function(){
  # wtchd <- downloadCSV("WTC_TEMP_CM_TREE-HEIGHT-DIAMETER_20121129-20140527_L1.CSV")
  
  insidechambers <- sprintf("C%02.f",(c(1:12)))
  # wtchd <- downloadCSV("WTC_TEMP_CM_TREE-HEIGHT-DIAMETER_20121129-20140527_L1.CSV")
  wtchd <- read.csv("cache/WTC_TEMP_CM_TREE-HEIGHT-DIAMETER_20121129-20140527_L1.CSV")
  wtchd$Date <- as.Date(wtchd$DateTime)
  
  
  wtchd <- droplevels(subset(wtchd, Date >= startDate & Date <= endDate & chamber %in% insidechambers))
  # Geometric mean crown width
  wtchd$Crown_width <- with(wtchd, sqrt(Crown_width_EW * Crown_width_NS))
  
  # Multiple stem heights on some trees. Take height of maximum for every date/chamber.
  wtchd <- summaryBy(. ~ Date + chamber, FUN=max, keep.names=TRUE, data=wtchd)
  wtchd_means <- summaryBy(. ~ Date + chamber, FUN=mean, keep.names=TRUE, data=wtchd)
  
  # Use mean crown widths; sometimes in wrong field
  wtchd$Crown_width <- wtchd_means$Crown_width
  
  # Interpolated tree heights
  datesout <- seq(startDate, endDate, by="1 day")
  interp_h <- function(d, datesout){
    
    sp <- spline(x=d$Date, y=d$Plant_height, xout=datesout)
    return(data.frame(Date=datesout, chamber=unique(d$chamber), Plant_height=sp$y))
  }
  treeh <- do.call(rbind, lapply(split(wtchd, wtchd$chamber), interp_h, datesout=datesout))
  
  
  #get tree diameter (diameter measured at 130 cm height: assumed this height as DBH)
  
  interp_d <- function(d, datesout){
    
    sp_d <- spline(x=d$Date, y=d$X130, xout=datesout)
    return(data.frame(Date=datesout, chamber=unique(d$chamber), diameter=sp_d$y))
  }
  
  
  # interpollate tree diameter
  #------------------------------------------------------------------------------------------------
  #------------------------------------------------------------------------------------------------
  #re-read data. this is due to overcome dropping diameter values in C11 when getting averages above
  
  wtchd_dia <- read.csv("cache/WTC_TEMP_CM_TREE-HEIGHT-DIAMETER_20121129-20140527_L1.CSV")
  wtchd_dia$Date <- as.Date(wtchd_dia$DateTime,formate="%d-%m-%Y")
  wtchd_dia <- droplevels(subset(wtchd_dia, Date >= startDate & Date <= endDate & chamber %in% insidechambers))
  daiam<- do.call(rbind, lapply(split(wtchd_dia, wtchd_dia$chamber), interp_d, datesout=datesout))
  
  # wtchd_dia$X130<-ifelse(wtchd_dia$Date==as.Date("2012-11-29"),wtchd_dia$X15,wtchd_dia$X130)
  #
  #
  # # linear interpollation of initial plant diameter between measurements (from 29-11-2012 to 19-03-2013)
  # # initial diameter at 130 cm were assumed to similar as diameter at 15 cm
  #
  # diam_i<-subset(wtchd_dia,wtchd_dia$Date<=as.Date("2013-03-20"))
  # diam_i$X130<-ifelse(diam_i$X130==0,3.78,diam_i$X130) # chamber 5 diameter was zero. replaced with mean of ambient trees
  #
  # fff <- lm(X130 ~ chamber+Measurement+chamber:Measurement+chamber:I(Measurement), data=diam_i) #fit linear model with intercept and slope for each chamber
  # # summary(fff)
  # diam_i$X130<-suppressWarnings(predict(fff,diam_i))
  # diam_rest<-subset(wtchd_dia,wtchd_dia$Measurement>9)
  #
  #
  # wtchd_dia_all<-plyr::rbind.fill(diam_rest,diam_i)
  
  # daiam<- do.call(rbind, lapply(split(wtchd_dia_all, wtchd_dia_all$chamber), interp_d, datesout=datesout))
  
  #------------------------------------------------------------------------------------------------
  #------------------------------------------------------------------------------------------------
  
  
  
  # Fitted crown widths
  f <- lm(Crown_width ~ chamber + Plant_height + chamber:Plant_height + chamber:I(Plant_height^0.33), data=wtchd)
  treeh$crownwidth <- predict(f,treeh)
  
  
  # fitted crown width for the linear phase (data from first 25 measurements)
  
  dat<-subset(wtchd,wtchd$Measurement<25)
  ff <- lm(Crown_width ~ chamber + Plant_height + chamber:Plant_height, data=dat)
  
  
  treeh.1<-subset(treeh,Date>=as.Date("2013-03-19")) #  measured crown width data available. Use the allometric model based on all data
  treeh.2<-subset(treeh,Date<as.Date("2013-03-19")) # use allometric model fitted for first 25 measurements. Assume linear increase in CW with time
  treeh.2$crownwidth<-predict(ff,treeh.2)
  
  treeh<-rbind(treeh.2,treeh.1)
  
  # Subset of flux data with just tree leaf area
  treela <- summaryBy(leafArea ~ chamber + Date, FUN=mean, keep.names=TRUE, data=wtcflux)
  treeh <- merge(treeh, treela)
  treeh<-merge(treeh, daiam)
  
  treeh <- treeh[order(treeh$Date),]
  
  return(treeh)
}





#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------
# save a fprocessed data to cache
# treeh <-make_treeh()
# saveRDS(treeh, file = "cache/treeh.rds")
#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------



# with(treeh,plot(Date,leafArea,col=chamber))
#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------

# with(treeh,plot(Date,crownwidth))

#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------


# run_chamber <- function(ch,return=c("hrflux","layflux")){
#   o <- getwd()
#   on.exit(setwd(o))
#   
#   return<-match.arg(return)
#   
#   startover()
#   
#   make_str()
#   make_con()
#   make_phy()
#   
#   make_trees(ch)
#   make_met(ch)
#   
#   setwd("maespa")
#   shell("MAESPA64.exe", intern=TRUE)
#   h <- readlayflux()
#   h_h <- readhrflux()
#   message(sprintf("Chamber %s simulation complete.",ch))
#   
#   h$chamber <- ch
#   h_h$chamber <- ch
#   h_h$DateTime <- seq(from=as.POSIXct(as.character(startDate), tz="UTC"), length=nrow(h), by="1 hour")
#   
#   if(return=="hrflux"){
#     return(h_h)
#   }
#   
#   if(return=="layflux"){
#     return(h)
#   }
# }

#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------

#function to fit GAM and plot GAM

fitgam <- function(X,Y,dfr, k=-1, R=NULL){
  dfr$Y <- dfr[,Y]
  dfr$X <- dfr[,X]
  if(!is.null(R)){
    dfr$R <- dfr[,R]
    model <- 2
  } else model <- 1
  dfr <- droplevels(dfr)
  
  
  if(model ==1){
    g <- gam(Y ~ s(X, k=k), data=dfr)
  }
  if(model ==2){
    g <- gamm(Y ~ s(X, k=k), random = list(R=~1), data=dfr)
  }
  
  return(g)
}


#' Plot a generalized additive model
#' @param x Variable for X axis (unquoted)
#' @param y Variable for Y axis (unquoted)
#' @param data Dataframe containing x and y
#' @param kgam the \code{k} parameter for smooth terms in gam.
#' @param random An optional random effect (quoted)
#' @param log Whether to add log axes for x or y (but no transformations are done).
#' @param fitoneline Whether to fit only 
smoothplot <- function(x,y,g=NULL,data,
                       fittype=c("gam","lm"),
                       kgam=4,
                       random=NULL,
                       randommethod=c("lmer","aggregate"),
                       log="",
                       fitoneline=FALSE,
                       pointcols=NULL,
                       linecols=NULL, 
                       xlab=NULL, ylab=NULL,
                       polycolor=alpha("lightgrey",0.7),
                       axes=TRUE,PCH=16,
                       ...){
  
  fittype <- match.arg(fittype)
  randommethod <- match.arg(randommethod)
  
  if(!is.null(substitute(g))){
    data$G <- as.factor(eval(substitute(g),data))
  } else {
    fitoneline <- TRUE
    data$G <- 1
  }
  data$X <- eval(substitute(x),data)
  data$Y <- eval(substitute(y),data)
  data <- droplevels(data)
  
  data <- data[!is.na(data$X) & !is.na(data$Y) & !is.na(data$G),]
  
  if(is.null(pointcols))pointcols <- palette()
  if(is.null(linecols))linecols <- palette()
  
  if(is.null(xlab))xlab <- substitute(x)
  if(is.null(ylab))ylab <- substitute(y)
  
  # If randommethod = aggregate, average by group and fit simple gam.
  if(!is.null(random) && randommethod == "aggregate"){
    data$R <- data[,random]
    
    data <- summaryBy(. ~ R, FUN=mean, na.rm=TRUE, keep.names=TRUE, data=data,
                      id=~G)
    R <- NULL
  }
  
  
  if(!fitoneline){
    
    d <- split(data, data$G)
    
    if(fittype == "gam"){
      fits <- lapply(d, function(x)try(fitgam("X","Y",x, k=kgam, R=random)))
      if(!is.null(random))fits <- lapply(fits, "[[", "gam")
    } else {
      fits <- lapply(d, function(x)lm(Y ~ X, data=x))
    }
    hran <- lapply(d, function(x)range(x$X, na.rm=TRUE))
  } else {
    if(fittype == "gam"){
      fits <- list(fitgam("X","Y",data, k=kgam, R=random))
      if(!is.null(random))fits <- lapply(fits, "[[", "gam")
    } else {
      fits <- list(lm(Y ~ X, data=data))
    }
    hran <- list(range(data$X, na.rm=TRUE))
    
  }
  
  with(data, plot(X, Y, xaxt="n",yaxt="n", pch=PCH, col=pointcols[G],
                  xlab=xlab, ylab=ylab, ...))
  
  if(axes){
    if(log=="xy")magaxis(side=1:2, unlog=1:2)
    if(log=="x"){
      magaxis(side=1, unlog=1)
      axis(2)
      box()
    }
    if(log=="y"){
      magaxis(side=2, unlog=2)
      axis(1)
      box()
    }
    if(log==""){
      axis(1)
      axis(2)
      box()
    }
  }
  
  for(i in 1:length(fits)){
    
    if(fittype == "gam"){
      nd <- data.frame(X=seq(hran[[i]][1], hran[[i]][2], length=101))
      if(!inherits(fits[[i]], "try-error")){
        p <- predict(fits[[i]],nd,se.fit=TRUE)
        addpoly(nd$X, p$fit-2*p$se.fit, p$fit+2*p$se.fit, col=polycolor[i])
        lines(nd$X, p$fit, col=linecols[i], lwd=2)
      }
    }
    if(fittype == "lm"){
      pval <- summary(fits[[i]])$coefficients[2,4]
      LTY <- if(pval < 0.05)1 else 5
      predline(fits[[i]], col=linecols[i], lwd=2, lty=LTY)
    }
  }
  
  return(invisible(fits))
}


addpoly <- function(x,y1,y2,col=alpha("lightgrey",0.7),...){
  ii <- order(x)
  y1 <- y1[ii]
  y2 <- y2[ii]
  x <- x[ii]
  polygon(c(x,rev(x)), c(y1, rev(y2)), col=col, border=NA,...)
}

predline <- function(fit, from=NULL, to=NULL, col=alpha("lightgrey",0.7), ...){
  
  if(is.null(from))from <- min(fit$model[,2], na.rm=TRUE)
  if(is.null(to))to <- max(fit$model[,2], na.rm=TRUE)
  
  newdat <- data.frame(X = seq(from,to, length=101))
  names(newdat)[1] <- names(coef(fit))[2]
  
  pred <- as.data.frame(predict(fit, newdat, se.fit=TRUE, interval="confidence")$fit)
  
  addpoly(newdat[[1]], pred$lwr, pred$upr, col=col)
  
  #ablinepiece(fit, from=from, to=to, ...)
  lines(pred$fit~newdat[,1])
}

#'@title Add a line to a plot
#'@description As \code{abline}, but with \code{from} and \code{to} arguments. 
#'If a fitted linear regression model is used as asn argument, it uses the min and max values of the data used to fit the model.
#'@param a Intercept (optional)
#'@param b Slope (optional)
#'@param reg A fitted linear regression model (output of \code{\link{lm}}).
#'@param from Draw from this X value
#'@param to Draw to this x value
#'@param \dots Further parameters passed to \code{\link{segments}}
#'@export
ablinepiece <- function(a=NULL,b=NULL,reg=NULL,from=NULL,to=NULL,...){
  
  # Borrowed from abline
  if (!is.null(reg)) a <- reg
  
  if (!is.null(a) && is.list(a)) {
    temp <- as.vector(coefficients(a))
    from <- min(a$model[,2], na.rm=TRUE)
    to <- max(a$model[,2], na.rm=TRUE)
    
    if (length(temp) == 1) {
      a <- 0
      b <- temp
    }
    else {
      a <- temp[1]
      b <- temp[2]
    }
  }
  
  segments(x0=from,x1=to,
           y0=a+from*b,y1=a+to*b,...)
  
}


