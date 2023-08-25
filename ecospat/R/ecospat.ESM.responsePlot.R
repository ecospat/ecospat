ecospat.ESM.responsePlot <- 
  function(ESM.EnsembleModeling.output,
            ESM.modeling.output,
            fixed.var.metric = 'median'){
            #,parallel = FALSE){
    
    models <- ESM.modeling.output$models
    weights <- ESM.EnsembleModeling.output$weights
    weights.EF <- ESM.EnsembleModeling.output$weights.EF
    resp <- ESM.EnsembleModeling.output$ESM.fit$resp.var
    data <- ESM.modeling.output$data@data.env.var
    min.data <- apply(data,2,min)
    max.data <- apply(data,2,max)
    data.fixed.all <- as.data.frame(t(apply(data,2,fixed.var.metric)))
    data.fixed.all <- data.fixed.all[rep(1, each = 1000), ]
    var.names <- colnames(data)
    
    ### All except one variable are kept constant and projected for each bivariate model along the gradient of the variable under investigation
    #if(!parallel){ # parallel computation not possible because projections will be overwriten with parallel computation
    proj.fixed.list <- list()
    for(i in 1:ncol(data)){
      data.fixed <- data.fixed.all
      data.fixed[,i] <- seq(min.data[i],max.data[i],(max.data[i]-min.data[i])/999)
      proj.fixed <- ecospat.ESM.Projection(ESM.modeling.output,new.env=data.fixed)
      proj.fixed.list[[i]] <- ecospat.ESM.EnsembleProjection(ESM.prediction.output=proj.fixed,
                                                            ESM.EnsembleModeling.output=ESM.EnsembleModeling.output)
      proj.fixed.list[[i]] <- cbind(data.fixed[,i],proj.fixed.list[[i]])
    }
    # plot(proj.fixed.list[[i]]$EF~ proj.fixed.list[[i]][,1], xlab='',main= names(proj.fixed.list)[i], ylab='predicted value', ylim=c(0,1000), type='n',las = TRUE)
    # 
    # if(ncol(proj.fixed.list[[i]])>1){
    #   for(mod.i in models){
    #   points(proj.fixed.list[[i]][,mod.i]~ proj.fixed.list[[i]][,1],col='grey', lwd=2,type='l')
    # }}
    # points(proj.fixed.list[[i]]$EF~ proj.fixed.list[[i]][,1], xlab= names(proj.fixed.list)[i],col='red', lwd=2,  type='l')
    # rug(data[,i], col='black')
    # }}else{
    #   proj.fixed.list <- foreach(i = 1:ncol(data),
    #                       .packages = c("biomod2","raster",'ecospat')) %dopar% {
    #                                                       
    #     data.fixed <- data.fixed.all
    #     
    #     data.fixed[,i] <- seq(min.data[i],max.data[i],(max.data[i]-min.data[i])/999)
    #     proj.fixed <- ecospat.ESM.Projection(ESM.modeling.output,new.env=data.fixed)
    #     
    #     proj.fixed <- ecospat.ESM.EnsembleProjection(ESM.prediction.output=proj.fixed,
    #                                    ESM.EnsembleModeling.output=ESM.EnsembleModeling.output)
    #     
    #     proj.fixed <- cbind(data.fixed[,i],proj.fixed)
    #     
    #     names(proj.fixed)[1] <- colnames(data)[i]
    #     proj.fixed
    #                       }
    #   }
      
      names(proj.fixed.list) <- colnames(data)
      old.par <- graphics::par(no.readonly = TRUE) 
      on.exit(graphics::par(old.par))
      xs <- floor(sqrt(length(var.names)))
      ys <- ceiling(length(var.names) / xs)
      graphics::par(mfrow=c(xs, ys))
      
      ColModels <- c("#222E50","#007991","#439A86","#BCD8C1","#E9D985","#FF6978","#6D435A","#352D39","#6E8894","#FA7921","#FE9920")
    
      for (i in 1:ncol(data)) {
        if (length(models) == 1) {
          plot(proj.fixed.list[[i]][, 2] ~ proj.fixed.list[[i]][, 
                                                                1], xlab = "", main = names(proj.fixed.list)[i], 
               ylab = "predicted value", ylim = c(min(sapply(proj.fixed.list, 
                                                             function(x) {
                                                               min(x[, -1])
                                                             })), max(sapply(proj.fixed.list, function(x) {
                                                               max(x[, -1])
                                                             }))), type = "n", las = TRUE)
          points(proj.fixed.list[[i]][, 2] ~ proj.fixed.list[[i]][, 
                                                                  1], xlab = names(proj.fixed.list)[i], col = "red", 
                 lwd = 2, type = "l")
        }else{
          plot(proj.fixed.list[[i]]$EF ~ proj.fixed.list[[i]][, 
                                                              1], xlab = "", main = names(proj.fixed.list)[i], 
               ylab = "predicted value", ylim = c(min(sapply(proj.fixed.list, 
                                                             function(x) {
                                                               min(x[, -1])
                                                             })), max(sapply(proj.fixed.list, function(x) {
                                                               max(x[, -1])
                                                             }))), type = "n", las = TRUE)
          legend("topleft", legend = c("ensemble", models), 
                 fill = c("red", ColModels[1:length(models)]), 
                 box.lty = 0)
          points(proj.fixed.list[[i]]$EF ~ proj.fixed.list[[i]][, 
                                                                1], xlab = names(proj.fixed.list)[i], col = "red", 
                 lwd = 2, type = "l")
          for (mod.i in models) {
            points(proj.fixed.list[[i]][, mod.i] ~ proj.fixed.list[[i]][, 
                                                                        1], col = ColModels[which(models == mod.i)], 
                   lwd = 2, type = "l")
            rug(data[, i], col = "black")
          }
        }
      }
      
    ## delete temporary files whcih were produced to calculate response plots
    unlink(grep('proj_data.fixed.ESM.BIOMOD',
         list.files(paste0(getwd(),'/','ESM.BIOMOD.output_',ESM.EnsembleModeling.output$species),recursive  = TRUE,full.names = TRUE,include.dirs = TRUE))
         ,value=TRUE),
         recursive = TRUE)
    unlink(grep("ESM_Projections.data.fixed",
          list.files(paste0(getwd(),"/", "ESM.BIOMOD.output_", ESM.EnsembleModeling.output$species),full.names = TRUE)
          ,value = TRUE))
    
    return(proj.fixed.list = proj.fixed.list)
  }
