##### mess by B. Petitpierre, 2011 ####
### mess function calculates the MESS (i.e. extrapolation) as in Maxent ###
### arguments ###
### proj: projection datase t###
### cal: calibration dataset ###
### w: weight for each predictor (e.g. variables importance in SDM) ###
### return values ###
### MESS: mess as calculated in Maxent, i.e. the minimal extrapolation values ###
### MESSw: sum of negative MESS values corrected by the total number of predictors;
### if there is no negative values, MESSw is then the mean MESS ###
### MESSneg: number of predictors on which there is extrapolation ###


ecospat.mess <- function(proj, cal, w = "default") {
  if (!is.matrix(proj)) {
    proj <- as.matrix(proj)
    cal <- as.matrix(cal)
  }
  
  xy.proj <- proj[,1:2]
  xy.cal <- cal[,1:2] #Not used at the moment but could be to plot some additonal stuff
  proj <- proj[,-c(1:2)]
  cal <- cal[,-c(1:2)]
  
  if (w == "default") {
    w <- rep(1, ncol(proj))
  }

  minp <- apply(cal, 2, min)
  minp <- sapply(minp, rep, t = nrow(proj))
  maxp <- apply(cal, 2, max)
  maxp <- sapply(maxp, rep, t = nrow(proj))
  ecdf.cal <- apply(cal, 2, ecdf)

  fn <- function(k, seqdeseq, proj) {
    lapply(proj[, k], seqdeseq[[k]])
  }
  fi <- round(unlist(sapply(1:ncol(proj), fn, ecdf.cal, proj, simplify = TRUE)) *
    100)


  # x=fi,proj,minp,maxp
  messi <- function(x) {
    if (x[1] == 0) {
      MESSi <- (x[2] - x[3])/(x[4] - x[3]) * 100
    }
    if (x[1] > 0 & x[1] <= 50) {
      MESSi <- 2 * x[1]
    }
    if (x[1] > 50 & x[1] < 100) {
      MESSi <- 2 * (100 - x[1])
    }
    if (x[1] == 100) {
      MESSi <- (x[4] - x[2])/(x[4] - x[3]) * 100
    }
    return(MESSi)
  }

  count.neg <- function(x) {
    return(length(which(x < 0)))
  }
  total <- round(matrix(apply(cbind(fi, as.vector(proj), as.vector(minp), 
                                    as.vector(maxp)),1, messi), nrow = nrow(proj)))
  if (ncol(proj) > 1) {
    MESS <- apply(total, 1, min)
    MESSneg <- total
    MESSneg[MESSneg[] > 0] <- 0
    MESSneg[which(MESS >= 0), ] <- total[which(MESS >= 0), ]
    MESSw <- round(apply(MESSneg, 1, weighted.mean, w = w))
    MESSneg <- apply(total, 1, count.neg)
    return(cbind(xy.proj,MESS, MESSw, MESSneg))
  } else {
    return(cbind(xy.proj,total))
  }
}

##### plot.mess by B. Petitpierre, 2011 #### plot the MESS extrapolation index onto
##### the geographical space ### arguments ### xy: xy coordinates of the projection
##### dataset t### mess.object: dataframe returned by the mess() function ### return
##### values ### MESS: mess as calculated in Maxent, i.e. the minimal extrapolation
##### values (red= negative, blue= positive values) ### MESSw: sum of negative MESS
##### values corrected by the total number of predictors; if there is no negative
##### values, MESSw is then the mean MESS (red= negative, blue= positive values)###
##### MESSneg: number of predictors on which there is extrapolation ###

ecospat.plot.mess <- function (mess.object, cex = 1, pch = 15) 
{
  #Plot MESS
  col.mess.neg <- colorRampPalette(c("white", "red"))
  col.mess.pos <- colorRampPalette(c("white", "blue"))
  col.neg <- col.mess.neg(max(1 + abs(mess.object[, 3])))
  col.pos <- col.mess.pos(max(1 + abs(mess.object[, 3])))
  par(mfrow = c(2, 2))
  plot(mess.object[,1:2], cex = cex, pch = pch, col = 0, main = "MESS", xlab = paste("min=", 
                                                                                     min(mess.object[, 3]), " & max=", max(mess.object[, 3]), sep = ""), ylab = "")
  points(mess.object[,1:2][which(mess.object[, 3] < 0), ], cex = cex, pch = pch, 
         col = col.neg[mess.object[which(mess.object[, 3] < 0), 3]])
  points(mess.object[,1:2][which(mess.object[, 3] > 0), ], cex = cex, pch = pch, 
         col = col.pos[mess.object[which(mess.object[, 3] > 0), 3]])
  #Plot MESSw
  col.neg <- col.mess.neg(max(1 + abs(mess.object[, 4])))
  col.pos <- col.mess.pos(max(1 + abs(mess.object[, 4])))
  plot(mess.object[,1:2], cex = cex, pch = pch, col = 0, main = "MESSw", xlab = paste("min=", min(mess.object[, 4]), " & max=", max(mess.object[, 4]), sep = ""), ylab = "")
  points(mess.object[,1:2][which(mess.object[, 4] < 0), ], cex = cex, pch = pch, col = col.neg[mess.object[which(mess.object[, 4] < 0), 4]])
  points(mess.object[,1:2][which(mess.object[, 4] > 0), ], cex = cex, pch = pch, col = col.pos[mess.object[which(mess.object[, 4] > 0), 4]])
  #Plot MESSneg
  col.neg <- col.mess.neg(max(1 + abs(mess.object[, 5])))
  plot(mess.object[,1:2], cex = cex, pch = pch, col = col.neg[mess.object[, 5] + 1], main = "#MESSneg", xlab = paste("min=", min(mess.object[, 5]), " & max=", max(mess.object[, 5]), sep = ""), ylab = "")
}