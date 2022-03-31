library(ggplot2)
mydot <- read.table("dot.txt",header = FALSE) #read the coordinates for the dots
myline <- read.table("line.txt",header = FALSE) #read the coordinates for the lines
myxtics <- read.table("xticks.txt",header = TRUE) #read the x axis ticks
myytics <- read.table("yticks.txt",header = TRUE) #read the y-axis ticks
mycolor <- c("red","blue") #assign the colors to your own color scheme
names(mycolor) <- c("F","R") #name the colors with the forward and reverse codes. F and R should match the first and second color, respectively


pdf(file="dotplots.pdf", width=3, height=3, colormodel="cmyk")
myplot <- ggplot(mydot,aes(x = V1,y=V2,color=V3)) + geom_point(size=0.5) #create the ggplot object with the dots
#now create the plot. Parameters can be modified
myplot + geom_segment(data = myline, aes(x=V1,y=V2,xend=V3,yend=V4,color = V5)) + 
  scale_y_continuous(breaks = myytics$ypos,labels = myytics$yname, name = "New", expand=c(0,0)) +
  scale_x_continuous(breaks=myxtics$xpos, labels = NULL, name = 'Old', expand=c(0,0)) +
  theme_light() +
  theme(axis.ticks.x=element_blank(), 
        axis.ticks.y =element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  scale_color_manual(values = mycolor)
  

dev.off()

  