# ----- Define a function for plotting a matrix ----- #
myImagePlot <- function(x, ...){
     min <- min(x)
     max <- max(x)
     yLabels <- rownames(x)
     xLabels <- colnames(x)
     title <-c()
  # check for additional function arguments
  if( length(list(...)) ){
    Lst <- list(...)
    if( !is.null(Lst$zlim) ){
       min <- Lst$zlim[1]
       max <- Lst$zlim[2]
    }
    if( !is.null(Lst$yLabels) ){
       yLabels <- c(Lst$yLabels)
    }
    if( !is.null(Lst$xLabels) ){
       xLabels <- c(Lst$xLabels)
    }
    if( !is.null(Lst$title) ){
       title <- Lst$title
    }
  }
# check for null values
if( is.null(xLabels) ){
   xLabels <- c(1:ncol(x))
}
if( is.null(yLabels) ){
   yLabels <- c(1:nrow(x))
}

layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(6,1), heights=c(1,1))

 # Red and green range from 0 to 1 while Blue ranges from 1 to 0
 ColorRamp <- rgb(seq(0,1,length=256),  
                  seq(0,1,length=256),  
                  seq(1,0,length=256))  
 ColorLevels <- seq(min, max, length=length(ColorRamp))

 # Reverse Y axis
 reverse <- nrow(x) : 1
 yLabels <- yLabels[reverse]
 x <- x[reverse,]

 # Data Map
#par(mar = c(2.5,25,0.5,1)) #for inter_BP
#par(mar=c(12,29,12,1)) #for high_BP
par(family="Times")
 par(mar=c(15,25,13,3)) #for high_kegg
#par(mar = c(2.5,25,2.5,1)) #for low_kegg
    # par(mar = c(5,25,5,1)) #for inter_kegg
#par(mar = c(2.5,25,2.5,1))

     image(1:length(xLabels), 1:length(yLabels), t(x), col=ColorRamp, xlab="",
 ylab="", axes=FALSE, zlim=c(min,max))
 if( !is.null(title) ){
    title(main=title)
 }
#     if (length(grep("KEGG",rownames(imagehigh)))>0){
 #        rownames(imagehigh) <- tolower(gsub("KEGG_","",rownames(imagehigh)))
  #   }else {
   #      rownames(imagehigh) <- tolower(rownames(imagehigh))
    # }
 #    yLabels <- rownames(imagehigh)
     axis(BELOW<-1, at=1:length(xLabels), labels=xLabels, cex.axis=0.6,)
     axis(LEFT <-2, at=1:length(yLabels), labels=yLabels, las= HORIZONTAL<-1,
          cex.axis=0.8)

 # Color Scale
# par(mar = c(3,3.2,2.5,3.5))
# image(1, ColorLevels,
     #matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
     #col=ColorRamp,
     #xlab="",ylab="",
     # xaxt="n")
 #layout(1)
}
# ----- END plot function ----- #
