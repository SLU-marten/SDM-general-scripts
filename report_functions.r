# Create a bar chart of the number of trawl surveys per year
trawlPerYear <- function(data, source = NULL){
  if(source != "Total"){
    data <- dplyr::filter(data, source.x == source)
  }
  dataSum <- data |>
    group_by(Year, Quarter) |>
    summarise("Number of surveys" = n())
  theme_set(theme_classic())
  p <- ggplot(dataSum, aes(x = Year, y = `Number of surveys`,
                           fill = factor(Quarter, levels = c(4,3,2,1)))) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("#56b4e9", "#009e73", "#f0e442", "#0072b2")) +
    scale_y_continuous(expand = c(0, 0)) +
    guides(fill = guide_legend(title = "Quarter"))
  grid.newpage()
  grid.draw(ggplotGrob(p))
}
# Create two bar charts
# First, the number of catches per year, quarter and species
# Second, the yearly catches per trawl area 
plotYear <- function(data, spec){
  p1 <- ggplot(dplyr::filter(data, sp == spec), aes(x = Year, y = Number, fill = Quarter)) +
    geom_bar(stat = "identity") +
    guides(fill = guide_legend(title = "Quarter")) +
    theme(legend.position = c(0.1, 0.7))
  p2 <- ggplot(dplyr::filter(data, sp == spec), aes(x = Year, y = nrPerSweptArea)) +
    geom_bar(stat = "identity")
  grid.newpage()
  grid.draw(rbind(ggplotGrob(p1), ggplotGrob(p2), size = "last"))
}

## Correlation graph functions ####
# Smooth line for the lower part
lower_fn <- function(data, mapping, ...){
  ggplot(data, mapping) +
    geom_point() +
    geom_smooth(fill="red", color="red")
}
# Color scheme for the upper part
upper_fn <- function(data, mapping, method="p", use="pairwise", ...){
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  corr <- cor(x,y, method=method, use=use)
  colFn <- colorRampPalette(c("blue", "white", "red"), interpolate="spline")
  fill <- colFn(100)[findInterval(corr, seq(-1,1, length=100))]
  ggally_cor(data=data, mapping=mapping, ...) +
    theme_void() +
    theme(panel.background = element_rect(fill=fill))
}

## Biomod Visualization ####
visMod <- function(mod, ens){
  myEvalPlot <- bm_PlotEvalBoxplot(bm.out = mod, group.by = c('algo', 'algo'), do.plot = FALSE)$plot
  myVarImpPlot <- bm_PlotVarImpBoxplot(bm.out = mod, group.by = c('expl.var', 'algo', 'algo'), do.plot = FALSE)$plot
  
  cat("\n::: panel-tabset\n")
  cat("\n\n##### Modell output\n")
  cat("```{html}\n")
  cat(paste(capture.output(mod), collapse = "\n"))
  cat("\n```")
  cat("\n\n##### Validering \n")
  cat("\n::: panel-tabset \n")
  cat("\n\n##### variant 1 \n")
  bm_PlotEvalMean(mod)$plot
  cat("\n\n##### Variant 2 \n")
  plot(myEvalPlot)
  cat("\n:::\n")
  cat("\n\n##### Vaiable importance\n")
  plot(myVarImpPlot)
  cat("\n\n##### Responskurvor\n")
  if(ens){
    bm_PlotResponseCurves(mod, do.progress = F)$plot
  } else{
    cat("\n::: panel-tabset \n")
    for(modType in c("GLM", "GAM", "GBM", "RF", "XGBOOST")){
      cat("\n\n####", modType, " \n")
      plot(bm_PlotResponseCurves(mod, models.chosen = get_built_models(mod, algo = modType), do.progress = F)$plot)
    }
    cat("\n:::\n")
  }
  cat("\n:::\n")
}

modQuartoFuncBiomod_pre <- function(modelID, species, timePeriod, quarter){
  cat("\n::: panel-tabset\n")
  for(sp in species) {
    cat("\n\n#####", sp," \n")
    if(length(timePeriod) > 1 ) cat("\n::: panel-tabset\n")  # Deactivate if we only have one timeperiod
    for(tp in names(timePeriod)){
      if(length(timePeriod) > 1 ) cat("\n\n#####", tp," \n")  # Deactivate if we only have one timeperiod
      if(length(quarter) > 1 ) cat("\n::: panel-tabset\n")  # Deactivate if we only have one quarter
      for(q in names(quarter)){
        if(length(timePeriod) > 1 ) cat("\n\n#####", q," \n") # Deactivate if we only have one quarter
        # Check if model failed
        errorFiles <- list.files(paste0("model output/", modelID, "/Error log"))
        if(paste(sp, tp, q, "errorlog.txt", sep = "_") %in% errorFiles){
          errorMsg <- read.table(paste0("model output/", modelID, "/Error log/", paste(sp, tp, q, "errorlog.txt", sep = "_")))
          cat(paste0(unlist(errorMsg), collapse = "\n"))
          next
        }
        fileName <- paste0(sp, "_", tp, "_", q)
        folderPA <- paste0("model output/", modelID, "/PA/", sp)
        myBiomodModelOut <- get(load(paste0(folderPA, "/", sp, ".", fileName, ".models.out")))
        myBiomodEM <- get(load(paste0(folderPA, "/", sp, ".", fileName, ".ensemble.models.out")))
        myBiomodProj <- get(load(paste0(folderPA, "/proj_",  fileName, "/", sp, ".", fileName, ".projection.out")))
        myBiomodEMProj <- get(load(paste0(folderPA, "/proj_",  fileName, "/", sp, ".", fileName, ".ensemble.projection.out")))
        EMprojPA <- rast(myBiomodEMProj@proj.out@link[2]) # Selecting PA map from TSS threshold
        EMprojPA[EMprojPA == 0] <- NA
        
        cat("\n::: panel-tabset\n")
        cat("\n\n##### Prediktion\n")
        cat("\n::: panel-tabset\n")
        cat("\n\n##### Ensemble\n")
        plot(myBiomodEMProj)
        cat("\n\n##### Enskilda modeller\n")
        cat("\n::: panel-tabset\n")
        for(modType in c("GLM", "GAM", "GBM", "RF", "XGBOOST")){
          cat("\n\n#####", modType, " \n")
          plot(myBiomodProj, algo = modType)
        }
        cat("\n:::\n")
        cat("\n\n##### PA-mask \n")
        pp <- EMprojPA[[1]] |> 
          aqua::rasterPlot(sluCol = "green", rev=T, nClass = 1) |> 
          aqua::mapExtra()
        plot(pp)
        cat("\n:::\n")
        cat("\n\n##### Modellering\n")
        cat("\n::: panel-tabset\n")
        cat("\n\n##### Data\n")
        cat("\n\n##### Enskilda modeller\n")
        visMod(myBiomodModelOut, ens = F)
        cat("\n\n##### Ensemblemodeller\n")
        # Ensemble
        tryCatch({
          visMod(myBiomodEM, ens = T)
        },
        error = function(e){cat(capture.output(message(e)))})
        cat("\n:::\n")
        cat("\n:::\n")
      } 
      if(length(quarter) > 1 ) cat("\n:::\n")
    }
    if(length(timePeriod) > 1 ) cat("\n:::\n")
  }
  cat("\n:::\n") 
}

## SDM Visualization ####
insetMap <- function (map, rast){
  ins <- sf::st_as_sfc(sf::st_bbox(rast))
  europe <- rnaturalearth::ne_countries(continent = 'europe', returnclass = "sf", scale = "large")
  m2 <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = europe) +
    ggplot2::geom_sf(data = ins, color = "red", fill = NA, size = 2, linewidth = 1) +
    xlim(0, 30) +
    ylim(50, 70) +
    ggthemes::theme_map() +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white"),
                   panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth  = 1))
  
  m3 <- cowplot::ggdraw() +
    cowplot::draw_plot(map) +
    cowplot::draw_plot(m2, x = 0.47, y = 0.705, width = 0.3, height = 0.3)
  
  m3
}

modQuartoFuncSDM_pre <- function(modelID, species, timePeriod, quarter, data){
  cat("\n::: panel-tabset\n")
  for(sp in names(species)) {
    cat("\n\n#####", sp," \n")
    if(length(timePeriod) > 1 ) cat("\n::: panel-tabset\n") # Deactivate if we only have one timeperiod
    for(tp in 1:length(timePeriod)){
      if(length(timePeriod) > 1 ) cat("\n\n#####", names(timePeriod)[tp]," \n") # Deactivate if we only have one timeperiod
      if(length(quarter) > 1 ) cat("\n::: panel-tabset\n") # Deactivate if we only have full year
      for(q in 1:length(quarter)){
        if(length(quarter) > 1) cat("\n\n#####", names(quarter)[q]," \n") # Deactivate if we only have full year
        # Check if model failed
        errorFiles <- list.files(paste0("model output/", modelID, "/Error log"))
        if(paste(sp, names(timePeriod)[tp], names(quarter)[q], "errorlog.txt", sep = "_") %in% errorFiles){
          errorMsg <- read.table(paste0("model output/", modelID, "/Error log/", paste(sp, names(timePeriod)[tp], names(quarter)[q], "errorlog.txt", sep = "_")))
          cat(paste0(unlist(errorMsg), collapse = "\n"))
          next
        }
        
        # myBiomodEMProj <- get(load(paste0("model output/", modelID, "/PA/", sp, "/proj_Current/", sp, ".Current.ensemble.projection.out")))
        myBiomodEMProj <- get(load(
          paste0("model output/", modelID, "/PA/", sp, "/proj_",  sp, "_", names(timePeriod)[tp], "_", names(quarter)[q],
                 "/", sp, ".", sp, "_", names(timePeriod)[tp], "_", names(quarter)[q], ".ensemble.projection.out")
        ))
        EMprojPA <- rast(myBiomodEMProj@proj.out@link[2]) # Selecting PA map from TSS threshold
        EMprojPA[EMprojPA == 0] <- NA
        
        sdm_M <- readRDS(paste0("model output/", modelID, "/Abu/", sp, "/", sp, "_", names(timePeriod)[tp], "_", names(quarter)[q], "_sdm_M.rds"))
        sdm_P <- readRDS(paste0("model output/", modelID, "/Abu/", sp, "/", sp, "_", names(timePeriod)[tp], "_", names(quarter)[q], "_sdm_P.rds"))
        sdm_E <- rast(paste0("model output/", modelID, "/Abu/", sp, "/", sp, "_", names(timePeriod)[tp], "_", names(quarter)[q], "_Ensemble_model.tif"))
        
        trawlAreaRast <- rast("data/gis/trawl_area.tif")
        trawlAreaLine <- st_read("data/gis/trawl_area.shp", quiet = T)
        
        enDelta <- rast(paste0("model output/", modelID, "/Delta/", sp, "_", names(timePeriod)[tp], "_", names(quarter)[q],"_deltamodel.tif"))
        enDelta <- enDelta * trawlAreaRast
        
        dfAbu <- dplyr::filter(data,
                               Year >= as.numeric(timePeriod[tp][1, 1]),
                               Year <= as.numeric(timePeriod[tp][2, 1]),
                               Quarter %in% unlist(unique(quarter[q])))
        
        cat("\n::: panel-tabset\n")
        cat("\n\n##### Report figure \n")
        # ext(enDelta) <- c(9.5, 13.5, 55, 59.5)
        p <- aqua::rasterPlot(enDelta, sluCol = "yellow", rev = T, nClass = 5, classMethod = "quantile") |>
          aqua::mapExtra() +
          geom_sf(data = trawlAreaLine, fill = "transparent", linewidth = 0.5, col = "black") +
          coord_sf(xlim = c(9.5, 13.5), ylim = c(55, 59.5))
        p <- p |>  insetMap(enDelta)
        
        r <- rast(xmin = 9, xmax = 14, ymin=55, ymax = 60, ncols = 40, nrows = 80)
        sp_vect <- terra::vect(dfAbu, geom=c("ShootLong", "ShootLat"))
        sp_rast <- terra::rasterize(sp_vect, r, field = sp, fun = mean)
        p2 <- aqua::rasterPlot(sp_rast, sluCol = "red", rev = T, nClass = 5, classMethod = "jenks") |>
          aqua::mapExtra() +
          geom_sf(data = trawlAreaLine, fill = "transparent", linewidth = 0.5, col = "black") +
          coord_sf(xlim = c(9.5, 13.5), ylim = c(55, 59.5))
        p2 <- insetMap(p2, enDelta)     
        p3 <- ggpubr::ggarrange(p, p2, labels = c("A", "B"), ncol = 2)
        ggsave(paste0("model output/", modelID, "/Delta/", sp, "_", names(timePeriod)[tp], "_", names(quarter)[q],"_deltamodel.png"),
               p3, width = 22, height = 12, units = "cm")
        plot(p3)
        
        cat("\n\n##### Model\n")
        cat("```{html}\n")
        tryCatch(cat(paste(capture.output(sdm_M), collapse = "\n")), error = function(e){cat(capture.output(e))})
        cat("\n```")
        tryCatch(cat(capture.output(getModelInfo(sdm_M)), "```", sep="\n"), error = function(e){cat(capture.output(e))})
        
        cat("\n\n##### Partial response\n")
        cat("\n::: panel-tabset\n")
        cat("\n\n##### Ensemble\n")
        tryCatch(plot(rcurve(sdm_M, id = sdm_M@run.info$modelID)), error = function(e){cat(capture.output(message(e)))})
        for(l in 1:(length(sdm_M@setting@methods))){
          cat("\n\n#####", sdm_M@setting@methods[l], " \n")
          tryCatch(plot(rcurve(sdm_M, id = (1:sdm_M@setting@n.replicates) + (l - 1) * sdm_M@setting@n.replicates)), error = function(e){cat(capture.output(e))})
          if(sdm_M@setting@methods[l] == "rf"){
            tryCatch(plot(rcurve(sdm_M, id = (1:sdm_M@setting@n.replicates) + (l - 1) * sdm_M@setting@n.replicates), gg = F), error = function(e){cat(capture.output(e))})
          }
        }
        cat("\n:::\n")
        cat("\n\n##### Ind. predictions\n")
        cat("\n::: panel-tabset\n")
        cat("\n\n##### Ensemble\n")
        tryCatch(plot(sdm_E), error = function(e){cat(capture.output(e))})
        for(l in 1:length(sdm_M@setting@methods)){
          cat("\n\n#####", sdm_M@setting@methods[l], " \n")
          tryCatch(plot(sdm_P[[which(grepl(sdm_M@setting@methods[l], names(sdm_P)))]]), error = function(e){cat(capture.output(e))})
        }
        cat("\n:::\n")
        cat("\n:::\n")
      }
      if(length(quarter) > 1 ) cat("\n:::\n")
    }
    if(length(timePeriod) > 1 ) cat("\n:::\n")
  }
  cat("\n:::\n")
}

## Validation ####
reportTables <- function(modelID, species, timePeriod, quarter, df, dfAbu){
  errorFiles <- list.files(paste0("model output/", modelID, "/Error log"))
  # Table 1. Model results from presence/absence ensemble models
  PAensResult <- c()
  for(sp in names(species)){
    for(tp in 1:length(timePeriod)){
      for(q in 1:length(quarter)){
        # Check and skip if model failed
        if(paste(sp, names(timePeriod)[tp], names(quarter)[q], "errorlog.txt", sep = "_") %in% errorFiles) next
        myBiomodEM <- get(load(paste0("model output/", modelID, "/PA/", sp, "/", sp, ".", sp, "_", names(timePeriod)[tp], "_", names(quarter)[q], ".ensemble.models.out")))
        t1 <- data.frame(
          Species = sp,
          Timeperiod = paste0(timePeriod[1, tp], "-", timePeriod[2, tp]), 
          Season = names(quarter)[q], 
          AUC = myBiomodEM@models.evaluation@val$calibration[2],
          Threshold = myBiomodEM@models.evaluation@val$cutoff[1],
          Sensitivity = myBiomodEM@models.evaluation@val$sensitivity[1],
          Specificity = myBiomodEM@models.evaluation@val$specificity[1]
        )
        PAensResult <- rbind(PAensResult, t1)
      }
    }
  }
  
  # Table 2. Model results from abundance ensemble models
  spStats <- list()
  spSummary <- c()
  for(sp in names(species)){
    for(tp in 1:length(timePeriod)){
      for(q in 1:length(quarter)){
        # Check and skip if model failed
        if(paste(sp, names(timePeriod)[tp], names(quarter)[q], "errorlog.txt", sep = "_") %in% errorFiles) next
        
        sdm_M <- readRDS(paste0("model output/", modelID, "/Abu/", sp, "/", sp, "_", names(timePeriod)[tp], "_", names(quarter)[q], "_sdm_M.rds"))
        nIterations <- length(sdm_M@models[[1]][[1]])
        stats <- c()
        k = 0
        for (modType in names(sdm_M@models[[1]])){
          for (iter in 1:nIterations) {
            tryCatch({
              modInfo <- sdm_M@models[[sp]][[modType]][[as.character(k * nIterations + iter)]]@evaluation$test.dep
              mase <- Metrics::mase(modInfo@observed, modInfo@predicted)
              R2 <- cor(modInfo@observed, modInfo@predicted)^2
              st <- modInfo@statistics
              all <- data.frame(mod = modType, id = iter, cor = st$COR[1], pVal = st$COR[2], dev = st$Deviance, RMSE = st$RMSE,
                                MAE = st$MAE, MASE = mase, R2 = R2, obs = mean(modInfo@observed), pred = mean(modInfo@predicted))
              stats <- rbind(stats, all)
            },
            error = function(e) {
              cat(capture.output(message(e)))
            }
            )
          }
          k = k + 1
        }
        spStats[[paste(sp, names(timePeriod)[tp], names(quarter)[q], sep = ".")]] <- stats
        spSum <- data.frame(Species = sp,
                            Timeperiod = names(timePeriod)[tp],
                            Quarter = names(quarter)[q],
                            Abu = median(stats$obs), AbuSD = sd(stats$obs),
                            PreAbu = median(stats$pred), preSD = sd(stats$pred),
                            R2 = median(stats$R2), R2SD = sd(stats$R2),
                            MASE = median(stats$MASE), MASESD = sd(stats$MASE))
        spSummary <- rbind(spSummary, spSum)
      }
    }
  }
  
  # Tabell 3. Results for the delta model
  library(sp)
  tab3 <- c()
  for(sp in names(species)){
    for(tp in 1:length(timePeriod)){
      for(q in 1:length(quarter)){
        if(paste(sp, names(timePeriod)[tp], names(quarter)[q], "errorlog.txt", sep = "_") %in% errorFiles) next
        enDelta <- rast(paste0("model output/", modelID, "/Delta/", sp, "_", names(timePeriod)[tp], "_", names(quarter)[q], "_deltamodel.tif"))
        rawValsFilt <- dplyr::filter(df,
                                     Year >= as.numeric(timePeriod[tp][1, 1]),
                                     Year <= as.numeric(timePeriod[tp][2, 1]),
                                     Quarter %in% unlist(unique(quarter[q]))) |>
          dplyr::select(ShootLat, ShootLong, all_of(species)[sp])
        
        rawValsFiltShape <- rawValsFilt
        coordinates(rawValsFiltShape) <- ~ ShootLong + ShootLat
        terra::crs(rawValsFiltShape) <- terra::crs(enDelta)
        rasValue <- terra::extract(raster(enDelta), rawValsFiltShape)
        rasValue[is.na(rasValue)] <- 0
        rawValsFilt$mod <- rasValue
        linModVars <- cbind(rasValue, rawValsFilt[[3]]) |> na.omit()
        linMod <- lm(linModVars[,1] ~ linModVars[,2])
        cor <- cor(rasValue, rawValsFilt[[3]])
        
        dfAbu[sp] <- dfAbu[, sp] / dfAbu$sweptAreaDS_km2
        spMap <- dfAbu |> dplyr::filter(Year >= as.numeric(timePeriod[tp][1, 1]),
                                        Year <= as.numeric(timePeriod[tp][2, 1]),
                                        Quarter %in% unlist(unique(quarter[q])))
        r <- rast(xmin = 9, xmax = 14, ymin = 55, ymax = 60, ncols = 40, nrows = 80)
        sp_vect <- terra::vect(spMap, geom = c("ShootLong", "ShootLat"))
        sp_rast <- terra::rasterize(sp_vect, r, field = sp, fun = mean)
        rastExtent <- enDelta
        ext(rastExtent) <- c(9.5, 13.5, 55, 59.5)
        enDelta_resamp <- enDelta
        enDelta_resamp[is.na(enDelta_resamp)] <- 0
        crs(sp_rast) <- crs(enDelta_resamp)
        sp_rast_resamp <- terra::resample(sp_rast, enDelta_resamp)
        corMap <- data.frame(values(enDelta_resamp), values(sp_rast_resamp)) |> na.omit() |> cor() %>% .[1, 2] |> round(3)
        
        tab3 <- rbind(tab3, data.frame(Species = sp, Timeperiod = names(timePeriod[q]), Quarter = names(quarter[q]),
                                       pVal = summary(linMod)$coefficients[2, 4], cor, corMap = corMap, fStat = summary(linMod)$fstatistic[1]))
      }
    }
  }
  
  # Tabell 4. Variable importance PA-model
  tab4 <- c()
  for(sp in names(species)){
    for(tp in 1:length(timePeriod)){
      for(q in 1:length(quarter)){
        # Check and skip if model failed
        if(paste(sp, names(timePeriod)[tp], names(quarter)[q], "errorlog.txt", sep = "_") %in% errorFiles) next
        myBiomodEM <- get(load(paste0("model output/", modelID, "/PA/", sp, "/", sp, ".", sp, "_", names(timePeriod)[tp], "_", names(quarter)[q], ".ensemble.models.out")))
        varImp <- biomod2::get_variables_importance(myBiomodEM)
        varImp_mean <- tapply(varImp$var.imp, varImp$expl.var, mean) |> as.data.frame() |> t()
        tab4 <- rbind(tab4, data.frame(Species = sp, Timeperiod = names(timePeriod[q]), Quarter = names(quarter[q]),
                                       varImp_mean))
      }
    }
  }
  
  # Tabell 5. Variable importance abundansmodel
  tab5 <- c()
  for(sp in names(species)){
    for(tp in 1:length(timePeriod)){
      for(q in 1:length(quarter)){
        # Check and skip if model failed
        if(paste(sp, names(timePeriod)[tp], names(quarter)[q], "errorlog.txt", sep = "_") %in% errorFiles) next
        sdm_M <- readRDS(paste0("model output/", modelID, "/Abu/", sp, "/", sp, "_", names(timePeriod)[tp], "_", names(quarter)[q], "_sdm_M.rds"))
        tryCatch({
          varImp <- sdm::getVarImp(sdm_M)
          varImp_cor <<- varImp@varImportanceMean$corTest |> arrange(variables) %>% .[, 1:2] |> as_tibble() %>% 
            pivot_longer(-variables) %>% pivot_wider(names_from = variables , values_from = value) %>% .[, -1]
        }, error = function(e){
          vars <- sdm_M@data@sdmFormula@vars@numeric[,1:2] |> arrange(names) |> as_tibble() %>% 
            pivot_longer(-names) %>% pivot_wider(names_from = names , values_from = value) %>% .[, -1]
          vars[1,] <- NA
          varImp_cor <<- vars
        })
        tab5 <- rbind(tab5, data.frame(Species = sp, Timeperiod = names(timePeriod[q]), Quarter = names(quarter[q]), varImp_cor))
      }
    }
  }
  
  return(list(PAensResult, spSummary, tab3, tab4, tab5))
}