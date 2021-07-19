

png("linear_food_chain.png", height = 5, width = 5, units = "in", res = 300)

par(mai=c(0.8,0.8,0,0), omi=rep(0,4))
## boundaries
ss <- 5
nn <- 7
rr <- ss*3
cc <- ss*nn
## mid-points
xm <- ss/2 + seq(0,cc-ss,ss)
ymt <- rr - ss/2
ymb <- ss/2
## arrow locs
x0t <- seq(ss, by=2*ss, len=3)
x1t <- x0t + ss
## empty plot space
plot(c(0,rr), c(0,cc), type="n", xlab="", ylab="",
     xaxt="n", yaxt="n", bty="n")
## top row: state
arrows(y0 = x0t, y1 = x1t, x0 = ymt, code = 3,
       col = "black", lwd = 2, length = 0.1)
symbols(y=xm[c(1,3,5,7)], x=rep(ymt,4), circles=rep(ss/5,4),
        lty="solid",  fg=NA, bg=c("#339933", "#ff8100", "#488fdf", "#844870"),
        inches=FALSE, add=TRUE, lwd=3)

dev.off()



