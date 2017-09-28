.plotRCC <- function(rccMean, rccM2sd, srName, mErn, rcc)
{
  plot(rccMean, type="l", col="red", xlab="Rank", ylab="# genes recovered",
       ylim=c(0,(max(c(rccM2sd, rcc)))+1), lwd=1, main=srName)
  polygon(c(0,1:(length(rccM2sd)),length(rccM2sd)), c(0,rccM2sd,0),
          col="#aa660040", border=FALSE)
  lines(rccM2sd, type="l", col="green", lwd=1)

  lines(rcc, col="blue")
  points(mErn["x"], mErn["y"], type="h", lty=2, col="blue")
  lines(c(mErn["x"], mErn["x"]), c(rccM2sd[mErn["x"]], mErn["y"]),
        lwd=2, col="blue")
  lines(c(0, mErn["x"]), c(mErn["y"], mErn["y"]), lty=2, col="darkgrey")
  points(mErn["x"], mErn["y"], type="b", pch=16, col="blue")
}
