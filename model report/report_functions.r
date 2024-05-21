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
  cat("\n##### Modell output\n")
  cat("```{html}\n")
  cat(paste(capture.output(mod), collapse = "\n"))
  cat("\n```")
  cat("\n##### Validering \n")
  cat("\n::: panel-tabset \n")
  cat("\n##### variant 1 \n")
  bm_PlotEvalMean(mod)$plot
  cat("\n##### Variant 2 \n")
  plot(myEvalPlot)
  cat("\n:::\n")
  cat("\n##### Vaiable importance\n")
  plot(myVarImpPlot)
  cat("\n##### Responskurvor\n")
  if(ens){
    bm_PlotResponseCurves(mod, do.progress = F)$plot
  } else{
    cat("\n::: panel-tabset \n")
    for(modType in c("GLM", "GAM", "GBM", "RF", "XGBOOST")){
      cat("\n####", modType, " \n")
      plot(bm_PlotResponseCurves(mod, models.chosen = get_built_models(mod, algo = modType), do.progress = F)$plot)
    }
    cat("\n:::\n")
  }
  cat("\n:::\n")
}

modQuartoFuncBiomod_pre <- function(modelID, species, timePeriod, quarter){
  cat("\n::: panel-tabset\n")
  for(sp in species) {
    cat("\n#####", sp," \n")
    # cat("\n::: panel-tabset\n") # Deactivated since we only have one timeperiod
    for(tp in names(timePeriod)){
      # cat("\n#####", names(timePeriod)[j]," \n") # Deactivated since we only have one timeperiod
      # cat("\n::: panel-tabset\n") # Deactivated since we only have full year
      for(q in names(quarter)){
        # Check if model failed
        errorFiles <- list.files(paste0("../model output/", modelID, "/Error log"))
        if(paste(sp, tp, q, "errorlog.txt", sep = "_") %in% errorFiles){
          errorMsg <- read.table(paste0("../model output/", modelID, "/Error log/", paste(sp, tp, q, "errorlog.txt", sep = "_")))
          cat(paste0(unlist(errorMsg), collapse = "\n"))
          next
        }
        
        # cat("\n#####", names(quarter)[k]," \n") # Deactivated since we only have full year
        fileName <- paste0(sp, "_", tp, "_", q)
        folderPA <- paste0("../model output/", modelID, "/PA/", sp)
        myBiomodModelOut <- get(load(paste0(folderPA, "/", sp, ".", fileName, ".models.out")))
        myBiomodEM <- get(load(paste0(folderPA, "/", sp, ".", fileName, ".ensemble.models.out")))
        myBiomodProj <- get(load(paste0(folderPA, "/proj_",  fileName, "/", sp, ".", fileName, ".projection.out")))
        myBiomodEMProj <- get(load(paste0(folderPA, "/proj_",  fileName, "/", sp, ".", fileName, ".ensemble.projection.out")))
        EMprojPA <- rast(myBiomodEMProj@proj.out@link[2]) # Selecting PA map from TSS threshold
        EMprojPA[EMprojPA == 0] <- NA
        
        cat("\n::: panel-tabset\n")
        cat("\n##### Prediktion\n")
        cat("\n::: panel-tabset\n")
        cat("\n##### Ensemble\n")
        plot(myBiomodEMProj)
        cat("\n##### Enskilda modeller\n")
        cat("\n::: panel-tabset\n")
        for(modType in c("GLM", "GAM", "GBM", "RF", "XGBOOST")){
          cat("\n#####", modType, " \n")
          plot(myBiomodProj, algo = modType)
        }
        cat("\n:::\n")
        cat("\n##### PA-mask \n")
        pp <- EMprojPA[[1]] |> 
          aqua::rasterPlot(sluCol = "green", rev=T, nClass = 1) |> 
          aqua::mapExtra()
        plot(pp)
        cat("\n:::\n")
        cat("\n##### Modellering\n")
        cat("\n::: panel-tabset\n")
        cat("\n##### Data\n")
        cat("\n##### Enskilda modeller\n")
        visMod(myBiomodModelOut, ens = F)
        cat("\n##### Ensemblemodeller\n")
        # Ensemble
        tryCatch({
          visMod(myBiomodEM, ens = T)
        },
        error = function(e){cat(capture.output(message(e)))})
        cat("\n:::\n")
        cat("\n:::\n")
      } 
      # cat("\n:::\n")
    }
    # cat("\n:::\n")
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
  for(sp in species) {
    cat("\n#####", sp," \n")
    if(length(timePeriod) > 1 ) cat("\n::: panel-tabset\n") # Deactivate if we only have one timeperiod
    for(tp in 1:length(timePeriod)){
      if(length(timePeriod) > 1 ) cat("\n#####", names(timePeriod)[tp]," \n") # Deactivate if we only have one timeperiod
      if(length(quarter) > 1 ) cat("\n::: panel-tabset\n") # Deactivate if we only have full year
      for(q in 1:length(quarter)){
        if(length(quarter) > 1 ) cat("\n#####", names(quarter)[q]," \n") # Deactivate if we only have full year
        # Check if model failed
        errorFiles <- list.files(paste0("../model output/", modelID, "/Error log"))
        if(paste(sp, tp, q, "errorlog.txt", sep = "_") %in% errorFiles){
          errorMsg <- read.table(paste0("../model output/", modelID, "/Error log/", paste(sp, names(timePeriod)[tp], names(quarter)[q], "errorlog.txt", sep = "_")))
          cat(paste0(unlist(errorMsg), collapse = "\n"))
          next
        }
        
        # myBiomodEMProj <- get(load(paste0("../model output/", modelID, "/PA/", sp, "/proj_Current/", sp, ".Current.ensemble.projection.out")))
        myBiomodEMProj <- get(load(
          paste0("../model output/", modelID, "/PA/", sp, "/proj_",  sp, "_", names(timePeriod)[tp], "_", names(quarter)[q],
                 "/", sp, ".", sp, "_", names(timePeriod)[tp], "_", names(quarter)[q], ".ensemble.projection.out")
          ))
        EMprojPA <- rast(myBiomodEMProj@proj.out@link[2]) # Selecting PA map from TSS threshold
        EMprojPA[EMprojPA == 0] <- NA

        sdm_M <- readRDS(paste0("../model output/", modelID, "/Abu/", sp, "/", sp, "_", names(timePeriod)[tp], "_", names(quarter)[q], "_sdm_M.rds"))
        sdm_P <- readRDS(paste0("../model output/", modelID, "/Abu/", sp, "/", sp, "_", names(timePeriod)[tp], "_", names(quarter)[q], "_sdm_P.rds"))
        sdm_E <- rast(paste0("../model output/", modelID, "/Abu/", sp, "/", sp, "_", names(timePeriod)[tp], "_", names(quarter)[q], "_Ensemble_model.tif"))
        
        trawlAreaRast <- rast("../data/gis/trawl_area.tif")
        trawlAreaLine <- st_read("../data/gis/trawl_area.shp", quiet = T)
        
        enDelta <- rast(paste0("../model output/", modelID, "/Delta/", sp, "_", names(timePeriod)[tp], "_", names(quarter)[q],"_deltamodel.tif"))
        enDelta <- enDelta * trawlAreaRast

        dfAbu <- dplyr::filter(data,
                               Year >= as.numeric(timePeriod[tp][1, 1]),
                               Year <= as.numeric(timePeriod[tp][2, 1]),
                               Quarter %in% unlist(unique(quarter[q])))
        
        cat("\n::: panel-tabset\n")
        cat("\n##### Report figure \n")
        ext(enDelta) <- c(9.5, 13.5, 55, 59.5)
        p <- aqua::rasterPlot(enDelta, sluCol = "yellow", rev = T, nClass = 5, classMethod = "quantile") |>
          aqua::mapExtra() +
          geom_sf(data = trawlAreaLine, fill = "transparent", linewidth = 0.5, col = "black") +
          coord_sf(xlim = c(9.5, 13.5), ylim = c(55, 59.5))
        p <- p |>  insetMap(enDelta)
        plot(p)
        
        cat("\n##### Catches \n")
        r <- rast(xmin = 9, xmax = 14, ymin=55, ymax = 60, ncols = 40, nrows = 80)
        sp_vect <- terra::vect(dfAbu, geom=c("ShootLong", "ShootLat"))
        sp_rast <- terra::rasterize(sp_vect, r, field = sp, fun = mean)
        p2 <- aqua::rasterPlot(sp_rast, sluCol = "red", rev = T, nClass = 5, classMethod = "jenks") |>
          aqua::mapExtra() +
          geom_sf(data = trawlAreaLine, fill = "transparent", linewidth = 0.5, col = "black") +
          coord_sf(xlim = c(9.5, 13.5), ylim = c(55, 59.5))
        p2 <- insetMap(p2, enDelta)
        plot(p2)
        
        cat("\n##### Model\n")
        cat("```{html}\n")
        tryCatch(cat(paste(capture.output(sdm_M), collapse = "\n")), error = function(e){cat(capture.output(e))})
        cat("\n```")
        tryCatch(cat(capture.output(getModelInfo(sdm_M)), "```", sep="\n"), error = function(e){cat(capture.output(e))})
        
        cat("\n##### Partial response\n")
        cat("\n::: panel-tabset\n")
        cat("\n##### Ensemble\n")
        tryCatch(plot(rcurve(sdm_M, id = sdm_M@run.info$modelID)), error = function(e){cat(capture.output(message(e)))})
        for(l in 1:(length(sdm_M@setting@methods))){
          cat("\n#####", sdm_M@setting@methods[l], " \n")
          tryCatch(plot(rcurve(sdm_M, id = (1:sdm_M@setting@n.replicates) + (l - 1) * sdm_M@setting@n.replicates)), error = function(e){cat(capture.output(e))})
          if(sdm_M@setting@methods[l] == "rf"){
            tryCatch(plot(rcurve(sdm_M, id = (1:sdm_M@setting@n.replicates) + (l - 1) * sdm_M@setting@n.replicates), gg = F), error = function(e){cat(capture.output(e))})
          }
        }
        cat("\n:::\n")
        cat("\n##### Prediction\n")
        cat("\n::: panel-tabset\n")
        cat("\n##### Ensemble\n")
        tryCatch(plot(sdm_E), error = function(e){cat(capture.output(e))})
        for(l in 1:length(sdm_M@setting@methods)){
          cat("\n#####", sdm_M@setting@methods[l], " \n")
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
  # Table 1. Model results from presence/absence ensemble models
  PAensResult <- c()
  for(sp in names(species)){
    for(tp in 1:length(timePeriod)){
      for(q in 1:length(quarter)){
        myBiomodEM <- get(load(paste0("../model output/", modelID, "/PA/", sp, "/", sp, ".", sp, "_", names(timePeriod)[tp], "_", names(quarter)[q], ".ensemble.models.out")))
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
      sdm_M <- readRDS(paste0("../model output/", modelID, "/Abu/", sp, "/", sp, "_", names(timePeriod)[tp], "_", names(quarter)[q], "_sdm_M.rds"))
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
            df <- data.frame(mod = modType, id = iter, cor = st$COR[1], pVal = st$COR[2], dev = st$Deviance, RMSE = st$RMSE,
                             MAE = st$MAE, MASE = mase, R2 = R2, obs = mean(modInfo@observed), pred = mean(modInfo@predicted))
            stats <- rbind(stats, df)
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
        enDelta <- rast(paste0("../model output/", modelID, "/Delta/", sp, "_", names(timePeriod)[tp], "_", names(quarter)[q], "_deltamodel.tif"))
        
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
        sp_rast_resamp <- terra::resample(sp_rast, enDelta_resamp)
        r2 <- data.frame(values(enDelta_resamp), values(sp_rast_resamp)) |> na.omit() |> cor() %>% .[1, 2] |> round(3)
        
        tab3 <- rbind(tab3, data.frame(i, q, pVal = summary(lm1)$coefficients[2,4], cor, corMap = r2, summary(lm1)$fstatistic[1]))
      }
    }
  }
  tab3
  
  # Tabell 4. Variable importance PA-model
  tab4 <- c()
  for(i in names(species)){
    for(q in 1:length(quarter)){
      dat <- getData(i, k = q)
      t3 <- tapply(dat$myBiomodEM@variables.importance@val$var.imp, dat$myBiomodEM@variables.importance@val$expl.var, mean) |> as.data.frame() |> t()
      rownames(t3) <- i
      # print(paste(i, length(t3)))
      tab4 <- rbind(tab4, t3)
    }
  }
  tab4

  return(list(PAensResult, spSummary, spStats, tab3, tab4))
}


























# 
# 
# 
# 
# 
# 
# modQuartoFuncGAM <- function(species, sp2, predictors, PA = F){
#   cat("  \n::: panel-tabset  \n")
#   for (i in 1:length(species)) {
#     cat("\n#####", species[i]," \n")
#     cat("\n::: panel-tabset\n")
#     for(j in 1:4){
#       cat("\n##### Kvartal ", j," \n")
#       cat("\n::: panel-tabset\n")
#       cat("\n##### Summering\n")
#       mod <- gamFunc(as.character(sp2[i]), Quart = j)
#       formatv <- function(x, ...) {mapply(format, x, scientific = abs(x) < 0.001, ...)} 
#       a <- summary(mod)
#       x1 <- paste(format(a$formula), collapse = " ") |>
#         strsplit("(?<=[~+])", perl = TRUE) |>
#         knitr::kable("html", col.names = c("Formula:"))
#       x2table <- a$p.table |>  as.data.frame()
#       x2sign <- case_when(x2table[4] < 0.001 ~ "(***)",
#                           x2table[4]  < 0.01 ~ "(**)",
#                           x2table[4]  < 0.05 ~ "(*)",
#                           x2table[4]  < 0.1 ~ "(.)",
#                           x2table[4]  > 0.1 ~ "")
#       x2 <- x2table |> 
#         mutate(across(.fns = formatv, digit = 1, nsmall = 3),
#                " " = x2sign) |>
#         knitr::kable("html", align = "c") |>
#         kableExtra::kable_styling(position = "center")
#       
#       x3table <- a$s.table |>  as.data.frame()
#       x3sign <- case_when(x3table[4] < 0.001 ~ "(***)",
#                           x3table[4] < 0.01 ~ "(**)",
#                           x3table[4] < 0.05 ~ "(*)",
#                           x3table[4] < 0.1 ~ "(.)",
#                           x3table[4] > 0.1 ~ "")
#       x3 <- x3table |>
#         mutate(across(.fns = formatv, digit = 1, nsmall = 3),
#                " " = x3sign) |> 
#         knitr::kable("html", align = "c") |>
#         kableExtra::kable_styling(position = "center")
#       x4 <- data.frame(`R-sq` = round(a$r.sq, 4),
#                        "Deviance explained" = scales::label_percent(0.1)(a$dev.expl),
#                        "n" = as.integer(a$n),
#                        "AIC" = mod$aic) |>
#         t() |>
#         knitr::kable("html") |>
#         kableExtra::column_spec(1, bold = T)|>
#         kableExtra::kable_styling(position = "center")
#       cat(x1)
#       cat("**Parametric coefficients:**")
#       cat(x2)
#       cat("**Approximate significance of smooth terms:**")
#       cat(x3)
#       cat("**Model info:**")
#       cat(x4)
#       cat("\n##### GAM check\n")
#       plot(appraise(mod))
#       cat("\n##### Respons\n")
#       if(PA){
#         plot(mod, pages=4, shade = T, rug = T, residuals = F, shift = coef(mod)[1], trans = plogis)
#       } else {
#         plot(mod, pages=4, shade = T, rug = T, residuals = F)
#       }
#       # plot(mod, pages=2, shade = T, rug = T, residuals = F)
#       cat("\n##### Prediktion\n")
#       plot(raster::predict(predictors, mod, type="response"))
#       cat("\n:::\n")
#     }
#     cat("\n:::\n")
#   }
#   cat("\n:::\n")
# }
# 
# modQuartoFuncSDM <- function(data, species, timePeriod, quarter, predictors, log = F){
#   cat("  \n::: panel-tabset  \n")
#   for(i in 1:length(species)) {
#     cat("\n#####", names(species)[i]," \n")
#     cat("\n::: panel-tabset\n")
#     for(j in 1:length(timePeriod)){
#       cat("\n#####", names(timePeriod)[j]," \n")
#       cat("\n::: panel-tabset\n")
#       for(k in 1:length(quarter)){
#         cat("\n#####", names(quarter)[k]," \n")
#         master2 <- data |> filter(Year >= as.numeric(timePeriod[j][1,1]), 
#                                     Year <= as.numeric(timePeriod[j][2,1]),
#                                     Quarter %in% unlist(unique(quarter[k])))
#         if(nrow(master2) > 50 & sum(master2[, as.character(species[i])] != 0) > 20){
#           if(log){
#             form <- as.formula(paste0("log(", species[i], ")", "~ bTemp + bSal + bO2 + Depth + botCur + slope + fishInt"))
#           } else{
#             form <- as.formula(paste0(species[i], "~ bTemp + bSal + bO2 + Depth + botCur + slope + fishInt"))
#           }
#           
#           d1 <- sdmData(form, train = master2)
#           m1 <- sdm(form, 
#                     data = d1, 
#                     methods = c('glm', 'gam', 'brt','rf'), 
#                     replication = c("sub"), 
#                     test.percent = 25, # Lägg undan 25% av datat i varje replikation
#                     modelTypes = c('abundance'),
#                     n = 6, # Kör 10 itterationer med olika subsamplings av modellen
#                     # parallelsetting = list(ncore=4, method="parallel"),
#                     modelsettings = list(
#                       glm = list(family='gaussian'),
#                       brt = list(n.trees=500, distribution='gaussian'))
#                     )
#           p1Trunk <- predict(m1, predictors)
#           en1Trunk <- ensemble(m1, p1Trunk, setting = list(method='weighted', stat='cor'))
#           cat("\n::: panel-tabset\n")
#           cat("\n##### Modell\n")
#           cat("```{html}\n")
#           tryCatch(cat(paste(capture.output(m1), collapse = "\n")), 
#                    error = function(e){})
#           cat("\n```")
#           tryCatch(cat(capture.output(getModelInfo(m1)), "```", sep="\n"), 
#                    error = function(e){})
#           cat("\n##### Responskurvor\n")
#           cat("\n::: panel-tabset\n")
#           cat("\n##### Ensemble\n")
#           tryCatch(plot(rcurve(m1, id=m1@run.info$modelID)), 
#                    error = function(e){})
#           for(l in 1:(length(m1@setting@methods))){
#               cat("\n#####", m1@setting@methods[l], " \n")
#               tryCatch(plot(rcurve(m1, id = 1:m1@setting@n.replicates + (l - 1) * m1@setting@n.replicates)), 
#                        error = function(e){})
#           }
#           cat("\n:::\n")
#           cat("\n##### Prediktion\n")
#           cat("\n::: panel-tabset\n")
#           cat("\n##### Ensemble\n")
#           tryCatch(plot(en1Trunk), error = function(e){
#             cat(message(e))
#           })
#           for(l in 1:length(m1@setting@methods)){
#             cat("\n#####", m1@setting@methods[l], " \n")
#             tryCatch(plot(p1Trunk[[which(grepl(m1@setting@methods[l], names(p1Trunk)))]]), 
#                      error = function(e){
#                        cat(message(e))
#                      })
#           }
#           cat("\n:::\n")
#           cat("\n:::\n")
#         } else{
#           if(nrow(master2) <= 100){
#             cat("För få datapunkter\n")
#           } else {
#             cat("För få datapunkter med fångst\n")
#           }
#         }
#       }
#       cat("\n:::\n")
#     }
#     cat("\n:::\n")
#   }
#   cat("\n:::\n")
# }
# 
# modQuartoFuncBiomod <- function(myResp, myExpl, myRespXY, myRespName){
#   
#   invisible(capture.output(myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
#                                        expl.var = myExpl,
#                                        resp.xy = myRespXY,
#                                        resp.name = myRespName)))
#   cat("\n::: panel-tabset\n")
#   cat("\n##### Modellering\n")
#   cat("\n::: panel-tabset\n")
#   cat("\n##### Data\n")
#   plot(myBiomodData)
# 
#   invisible(capture.output(myBiomodOptions <- BIOMOD_ModelingOptions(
#     GLM = list(type = 'quadratic',
#                interaction.level = 0,
#                myFormula = NULL,
#                test = 'BIC',
#                family = 'binomial',
#                control = glm.control(epsilon = 1e-08,
#                                      maxit = 1000,
#                                      trace = FALSE)),
#     GAM = list(family = "quasipoisson"),
#     GBM = list(distribution = 'bernoulli',
#                n.trees = 7500,
#                interaction.depth = 9,
#                n.minobsinnode = 5,
#                shrinkage = 0.01,
#                bag.fraction = 0.5,
#                train.fraction = 1,
#                cv.folds = 3,
#                keep.data = FALSE,
#                verbose = FALSE,
#                perf.method = 'cv',
#                n.cores = 1),
#     RF = list(do.classif = TRUE,
#               ntree = 500,
#               mtry = 5,
#               sampsize = NULL,
#               nodesize = 5,
#               maxnodes = NULL))))
#   
#   cat("\n##### Enskilda modeller\n")
#   # Enskilda modeller med korsvalidering. CV kan bytas ut mot en klassisk datasplit
#   invisible(capture.output(myBiomodModelOut <- 
#     BIOMOD_Modeling(bm.format = myBiomodData,
#                     bm.options = myBiomodOptions,
#                     modeling.id = 'AllModels',
#                     models = c("GLM", "GBM", "GAM", "RF"),
#                     nb.rep = 2, # Antal repetitioner för varje kalibrering/Validering
#                     data.split.perc = 75, # Andel datapunkter i kalibrering
#                     # data.split.table = myBiomodCV,
#                     do.progress = FALSE,
#                     var.import = 2, #Hur många gånger ska varje variabel permuteras för att beräkna variable importance
#                     metric.eval = c('TSS','ROC'), # Valideringsmetoder
#                     do.full.models = FALSE #Anger om kalibrering ska ske över hela datasetet"
#                     # seed.val = 123, 
#                     # nb.cpu = 6 # Antal kärnor som ska användas
#     )))
#   
#   # get_evaluations(myBiomodModelOut)
#   
#   invisible(capture.output(myRespPlot <-
#     bm_PlotResponseCurves(bm.out = myBiomodModelOut,
#                           models.chosen = get_built_models(myBiomodModelOut),
#                           fixed.var = 'median', do.plot = FALSE)$plot))
#   invisible(capture.output(myEvalPlot <- bm_PlotEvalBoxplot(bm.out = myBiomodModelOut, group.by = c('algo', 'algo'), do.plot = FALSE)$plot))
#   invisible(capture.output(myVarImpPlot <- bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut, group.by = c('expl.var', 'algo', 'algo'), do.plot = FALSE)$plot))
#   
#   visMod <- function(mod, evalPlot, myVarImpPlot, resp){
#     cat("\n::: panel-tabset\n")
#     cat("\n##### Modell output\n")
#     cat("```{html}\n")
#     cat(paste(capture.output(mod), collapse = "\n"))
#     cat("\n```")
#     cat("\n##### Validering\n")
#     plot(evalPlot)
#     cat("\n##### Vaiable importance\n")
#     plot(myVarImpPlot)
#     cat("\n##### Responskurvor\n")
#     plot(resp)
#     cat("\n:::\n")
#   }
#   
#   visMod(myBiomodModelOut, myEvalPlot, myVarImpPlot, myRespPlot)
#   
#   cat("\n##### Ensemblemodeller\n")
#   # Ensemble
#   invisible(capture.output(myBiomodEM <- BIOMOD_EnsembleModeling(bm.mod = myBiomodModelOut,
#                                         models.chosen = 'all',
#                                         em.by = 'all',
#                                         em.algo = c('EMmean', 'EMmedian', 'EMwmean'),
#                                         metric.select = c('ROC'),
#                                         metric.select.thresh = c(0.7),
#                                         var.import = 2,
#                                         metric.eval = c('TSS', 'ROC'),
#                                         EMwmean.decay = 'proportional',
#                                         do.progress = FALSE
#                                         )))
#   
#   invisible(capture.output(myEMRespPlot <- bm_PlotResponseCurves(bm.out = myBiomodEM,
#                                         models.chosen = get_built_models(myBiomodEM),
#                                         fixed.var = 'median', do.plot = FALSE)$plot))
#   invisible(capture.output(myEMEvalPlot <- bm_PlotEvalBoxplot(bm.out = myBiomodEM, group.by = c('algo', 'algo'), do.plot = FALSE)$plot))
#   invisible(capture.output(myEMVarImpPlot <- bm_PlotVarImpBoxplot(bm.out = myBiomodEM, group.by = c('expl.var', 'algo', 'algo'), do.plot = FALSE)$plot))
# 
#   
#   
#   visMod(myBiomodEM, myEMEvalPlot, myEMVarImpPlot, myEMRespPlot)
#   
#   cat("\n:::\n")
#   
#   cat("\n##### Prediktion\n")
#   cat("\n::: panel-tabset\n")
#   cat("\n##### Enskilda modeller\n")
#   ### Prediktioner
#   
#   ###### Enskilda modeller
#   # Project single models
#   invisible(capture.output(myBiomodProj <- BIOMOD_Projection(bm.mod = myBiomodModelOut,
#                                     proj.name = 'Current',
#                                     new.env = myExpl,
#                                     models.chosen = 'all',
#                                     metric.binary = 'all',
#                                     metric.filter = 'all',
#                                     build.clamping.mask = TRUE)))
#   # myBiomodProj
#   plot(myBiomodProj)
#   
#   cat("\n##### Ensemble\n")
#   # Ensemble prediktion
#   invisible(capture.output(myBiomodEMProj <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM,
#                                                proj.name = 'CurrentEM',
#                                                new.env = myExpl,
#                                                models.chosen = 'all',
#                                                metric.binary = 'all',
#                                                metric.filter = 'all')))
#   # myBiomodEMProj
#   plot(myBiomodEMProj)
#   
#   cat("\n##### PA-mask \n")
#   EMproj <- rast(myBiomodEMProj@proj.out@link[1])[[3]]
#   cutoff <- myBiomodEM@models.evaluation@val |> 
#     filter(algo == "EMwmean", metric.eval == "TSS") |> 
#     dplyr::select(cutoff) |> 
#     as.numeric()
#   rclmat <- matrix(c(0, cutoff, NA, cutoff, 1000, 1), ncol=3, byrow=TRUE)
#   EMprojPA <- terra::classify(EMproj, rclmat)
#   predStackDelta <<- rast(predStack) * EMprojPA
#   plot(terra::classify(EMproj, matrix(c(0, cutoff, 0, cutoff, 1000, 1), ncol=3, byrow=TRUE)))
#   
#   cat("\n:::\n")
#   cat("\n:::\n")
#   
# }
# 
# modQuartoFuncBiomod_pre <- function(path, art, timePeriod, quarter){
#   cat("\n::: panel-tabset\n")
#   for(i in 1:length(art)) {
#     cat("\n#####", art[i]," \n")
#     # cat("\n::: panel-tabset\n")
#     for(j in 1:length(timePeriod)){
#       # cat("\n#####", names(timePeriod)[j]," \n")
#       # cat("\n::: panel-tabset\n")
#       for(k in 1:length(quarter)){
#         # cat("\n#####", names(quarter)[k]," \n")
#         myBiomodModelOut_path <- paste0(path, gsub(" ", "_", art[i]),"_PA_", names(timePeriod)[j], "_", names(quarter)[k], "_myBiomodModelOut.rds")
#         myRespPlot_path <- paste0(path, gsub(" ", "_", art[i]),"_PA_", names(timePeriod)[j], "_", names(quarter)[k], "_myRespPlot.rds")
#         myEvalPlot_path <- paste0(path, gsub(" ", "_", art[i]),"_PA_", names(timePeriod)[j], "_", names(quarter)[k], "_myEvalPlot.rds")
#         myVarImpPlot_path <- paste0(path, gsub(" ", "_", art[i]),"_PA_", names(timePeriod)[j], "_", names(quarter)[k], "_myVarImpPlot.rds")
#         myEMRespPlot_path <- paste0(path, gsub(" ", "_", art[i]),"_PA_", names(timePeriod)[j], "_", names(quarter)[k], "_myEMRespPlot.rds")
#         myEMEvalPlot_path <- paste0(path, gsub(" ", "_", art[i]),"_PA_", names(timePeriod)[j], "_", names(quarter)[k], "_myEMEvalPlot.rds")
#         myEMVarImpPlot_path <- paste0(path, gsub(" ", "_", art[i]),"_PA_", names(timePeriod)[j], "_", names(quarter)[k], "_myEMVarImpPlot.rds")
#         myBiomodProj_path <- paste0(path, gsub(" ", "_", art[i]),"_PA_", names(timePeriod)[j], "_", names(quarter)[k], "_myBiomodProj.rds")
#         myBiomodEM_path <- paste0(path, gsub(" ", "_", art[i]),"_PA_", names(timePeriod)[j], "_", names(quarter)[k], "_myBiomodEM.rds")
#         myBiomodEMProj_path <- paste0(path, gsub(" ", "_", art[i]),"_PA_", names(timePeriod)[j], "_", names(quarter)[k], "_myBiomodEMProj.rds")
#         EMprojPA_path <- paste0(path, gsub(" ", "_", art[i]),"_PA_", names(timePeriod)[j], "_", names(quarter)[k], "_EMprojPA.rds")
#         
#         lf <- list.files(path, full.names = T)
#         
#         if(all(c(myBiomodModelOut_path, myRespPlot_path, myEvalPlot_path, myVarImpPlot_path, myEMRespPlot_path, myEMEvalPlot_path,
#                  myEMVarImpPlot_path, myBiomodProj_path, myBiomodEM_path, myBiomodEMProj_path, EMprojPA_path) %in% lf)){
#           myBiomodModelOut <- readRDS(myBiomodModelOut_path)
#           myRespPlot<- readRDS(myRespPlot_path)
#           myEvalPlot <- readRDS(myEvalPlot_path)
#           myVarImpPlot <- readRDS(myVarImpPlot_path)
#           myEMRespPlot <- readRDS(myEMRespPlot_path)
#           myEMEvalPlot <- readRDS(myEMEvalPlot_path)
#           myEMVarImpPlot <- readRDS(myEMVarImpPlot_path)
#           myBiomodProj <- readRDS(myBiomodProj_path)
#           myBiomodEM <- readRDS(myBiomodEM_path)
#           myBiomodEMProj <- readRDS(myBiomodEMProj_path)
#           EMprojPA <- readRDS(EMprojPA_path)
#           
#           cat("\n::: panel-tabset\n")
#           cat("\n##### Prediktion\n")
#           cat("\n::: panel-tabset\n")
#           cat("\n##### Ensemble\n")
#           # invisible(capture.output(plot(myBiomodEMProj)))
#           plot(myBiomodEMProj)
#           cat("\n##### Enskilda modeller\n")
#           cat("\n::: panel-tabset\n")
#           for(modType in c("GLM", "GAM", "GBM", "RF", "XGBOOST")){
#             cat("\n#####", modType, " \n")
#             plot(myBiomodProj, algo = modType)
#           }
#           cat("\n:::\n")
#           
#           cat("\n##### PA-mask \n")
#           pp <- EMprojPA[[1]] |> aqua::rasterPlot(sluCol = "green", rev=T, nClass = 1) |> 
#             aqua::mapExtra()
#           plot(pp)
#           # invisible(capture.output(plot(EMproj[[1]])))
#           cat("\n:::\n")
#           
#           cat("\n##### Modellering\n")
#           cat("\n::: panel-tabset\n")
#           cat("\n##### Data\n")
#           cat("\n##### Enskilda modeller\n")
#           visMod <- function(mod, evalPlot, myVarImpPlot, resp, ens){
#             cat("\n::: panel-tabset\n")
#             cat("\n##### Modell output\n")
#             cat("```{html}\n")
#             cat(paste(capture.output(mod), collapse = "\n"))
#             cat("\n```")
#             cat("\n##### Validering \n")
#             cat("\n::: panel-tabset \n")
#             cat("\n##### variant 1 \n")
#             bm_PlotEvalMean(mod)$plot
#             cat("\n##### Variant 2 \n")
#             plot(evalPlot)
#             cat("\n:::\n")
#             cat("\n##### Vaiable importance\n")
#             plot(myVarImpPlot)
#             cat("\n##### Responskurvor\n")
#             if(ens){
#               bm_PlotResponseCurves(mod, do.progress = F)$plot
#             } else{
#               cat("\n::: panel-tabset \n")
#               for(modType in c("GLM", "GAM", "GBM", "RF", "XGBOOST")){
#                 cat("\n####", modType, " \n")
#                 plot(bm_PlotResponseCurves(mod, models.chosen = get_built_models(mod, algo = modType), do.progress = F)$plot)
#               }
#               cat("\n:::\n")
#             }
#             cat("\n:::\n")
#           }
#           visMod(myBiomodModelOut, myEvalPlot, myVarImpPlot, myRespPlot, ens = F)
#           
#           cat("\n##### Ensemblemodeller\n")
#           # Ensemble
#           tryCatch({
#             visMod(myBiomodEM, myEMEvalPlot, myEMVarImpPlot, myEMRespPlot, ens = T)
#           },
#           error = function(e){})
#           cat("\n:::\n")
#           cat("\n:::\n")
#         } else{
#           "Ingen fil kunde skapas"
#         }
#       } 
#       # cat("\n:::\n")
#     }
#     # cat("\n:::\n")
#   }
#   cat("\n:::\n") 
# }
# 
#   
# modQuartoFuncSDM_pre <- function(path, art, timePeriod, quarter, data){
#   cat("  \n::: panel-tabset  \n")
#   for(i in 1:length(art)) {
#     cat("\n#####", art[i]," \n")
#     # cat("\n::: panel-tabset\n")
#     for(j in 1:length(timePeriod)){
#       # cat("\n#####", names(timePeriod)[j]," \n")
#       # cat("\n::: panel-tabset\n")
#       for(k in 1:length(quarter)){
#         # cat("\n#####", names(quarter)[k]," \n")
#         EMprojPA_path <- paste0(path, gsub(" ", "_", art[i]),"_PA_", names(timePeriod)[j], "_", names(quarter)[k], "_EMprojPA.rds")
#         EMprojPA <- readRDS(EMprojPA_path)
#         trawlArea <- rast("C:/Github Projcects/NMK/WestCoast/indata/gis/trawl_area.tif")
#         data[gsub(" ", "", art[i])] <- data[,art[i]]/data$sweptAreaDS_km2
#         master2 <- data |> dplyr::filter(Year >= as.numeric(timePeriod[j][1,1]),
#                                          Year <= as.numeric(timePeriod[j][2,1]),
#                                          Quarter %in% unlist(unique(quarter[k])))
# 
#         if(nrow(master2) > 50 & sum(master2[, as.character(art[i])] != 0) > 20){
#           lf <- list.files(path, full.names = T)
#           d1p <- paste0(path, gsub(" ", "_", art[i]),"_Abu_", names(timePeriod)[j], "_", names(quarter)[k],"_d1_2.rds")
#           m1p <- paste0(path, gsub(" ", "_", art[i]),"_Abu_", names(timePeriod)[j], "_", names(quarter)[k],"_m1_2.rds")
#           p1Trunkp <- paste0(path, gsub(" ", "_", art[i]),"_Abu_", names(timePeriod)[j], "_", names(quarter)[k],"_p1Trunk_2.rds")
#           en1Trunkp <- paste0(path, gsub(" ", "_", art[i]),"_Abu_", names(timePeriod)[j], "_", names(quarter)[k],"_en1Trunk_2.rds")
#           if(all(c(d1p, m1p, p1Trunkp, en1Trunkp) %in% lf)){
#             d1 <-readRDS(d1p)
#             m1 <-readRDS(m1p)
#             p1Trunk <-readRDS(p1Trunkp)
#             en1Trunk <-readRDS(en1Trunkp)
# 
#             enDelta <- (en1Trunk * raster(EMprojPA[[1]]) * raster(trawlArea))
# 
#             cat("\n::: panel-tabset\n")
#             cat("\n##### Rapportfigur \n")
# 
#             spMap <- master2 |>
#               dplyr::select(Year, Month, Quarter, Depth, source.x, ShootLat, ShootLong, !!!syms(species))
#             r <- rast(xmin = 9, xmax = 14, ymin=55, ymax=60, ncols=40, nrows=80)
#             sp_vect <- terra::vect(spMap, geom=c("ShootLong", "ShootLat"))
#             sp_rast <- terra::rasterize(sp_vect, r, field = art[i], fun = mean)
#             # rast_df2 <- dplyr::arrange(na.omit(raster::as.data.frame(sp_rast, xy = TRUE)), desc(names(sp_rast)))
#             # rast_df <- dplyr::arrange(na.omit(raster::as.data.frame(enDelta, xy = TRUE)), desc(names(enDelta)))
#             rastExtent <- enDelta
#             # extent(rastExtent) <- c(min(rast_df$x), max(rast_df$x), min(rast_df$y), max(rast_df$y))
#             extent(rastExtent) <- c(9.5, 13.5, 55, 59.5)
#             insetMap2 <- function (map, rast){
#               ins <- sf::st_as_sfc(sf::st_bbox(rast))
#               europe <- rnaturalearth::ne_countries(continent = 'europe', returnclass = "sf", scale = "large")
#               m2 <- ggplot2::ggplot() +
#                 ggplot2::geom_sf(data = europe) +
#                 ggplot2::geom_sf(data = ins, color = "red", fill = NA, size = 2, linewidth = 1) +
#                 xlim(0, 30) +
#                 ylim(50, 70) +
#                 ggthemes::theme_map() +
#                 ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white"),
#                                panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth  = 1))
# 
#               m3 <- cowplot::ggdraw() +
#                 cowplot::draw_plot(map) +
#                 cowplot::draw_plot(m2, x = 0.47, y = 0.705, width = 0.3, height = 0.3)
# 
#               m3
#             }
#             p <- aqua::rasterPlot(enDelta, sluCol = "yellow", rev=T, nClass = 5, classMethod = "quantile") |>
#               aqua::mapExtra() +
#               geom_sf(data = trawl_area, fill = "transparent", linewidth = 0.5, col = "black") +
#               coord_sf(xlim = c(9.5,13.5), ylim = c(55,59.5))
#             p <- p |>  insetMap2(rastExtent)
#             plot(p)
#             cat("\n##### Antal per km2 \n")
#             tryCatch({
#               p2 <- aqua::rasterPlot(sp_rast, sluCol = "red", rev=T, nClass = 5, classMethod = "jenks") |>
#                 aqua::mapExtra() +
#                 geom_sf(data = trawl_area, fill = "transparent", linewidth = 0.5, col = "black") +
#                 coord_sf(xlim = c(9.5,13.5), ylim = c(55,59.5))
#               p2 <- insetMap2(p2, rastExtent)
#             },
#             error = function(cond){
#               p2 <- aqua::rasterPlot(sp_rast, sluCol = "red", rev=T, nClass = 3, classMethod = "jenks") |>
#                 aqua::mapExtra() +
#                 geom_sf(data = trawl_area, fill = "transparent", linewidth = 0.5, col = "black") +
#                 coord_sf(xlim = c(9.5,13.5), ylim = c(55,59.5))
#               p2 <- insetMap2(p2, rastExtent)
#             })
#             plot(p2)
#             cat("\n##### Modell\n")
#             cat("```{html}\n")
#             tryCatch(cat(paste(capture.output(m1), collapse = "\n")),
#                      error = function(e){})
#             cat("\n```")
#             tryCatch(cat(capture.output(getModelInfo(m1)), "```", sep="\n"),
#                      error = function(e){})
#             cat("\n##### Responskurvor\n")
#             cat("\n::: panel-tabset\n")
#             cat("\n##### Ensemble\n")
#             tryCatch(plot(rcurve(m1, id = m1@run.info$modelID)), error = function(e){})
#             for(l in 1:(length(m1@setting@methods))){
#               cat("\n#####", m1@setting@methods[l], " \n")
#               plot(rcurve(m1, id = (1:m1@setting@n.replicates) + (l - 1) * m1@setting@n.replicates))
#               if(m1@setting@methods[l] == "rf"){
#                 plot(rcurve(m1, id = (1:m1@setting@n.replicates) + (l - 1) * m1@setting@n.replicates), gg = F)
#               }
#             }
#             cat("\n:::\n")
#             cat("\n##### Prediktion\n")
#             cat("\n::: panel-tabset\n")
#             cat("\n##### Ensemble\n")
#             tryCatch(plot(en1Trunk), error = function(e){
#               cat(message(e))
#             })
#             for(l in 1:length(m1@setting@methods)){
#               cat("\n#####", m1@setting@methods[l], " \n")
#               tryCatch(plot(p1Trunk[[which(grepl(m1@setting@methods[l], names(p1Trunk)))]]),
#                        error = function(e){
#                          cat(message(e))
#                        })
#             }
#             cat("\n:::\n")
#             cat("\n:::\n")
#           }
#         } else{
#           if(nrow(master2) <= 100){
#             cat("För få datapunkter\n")
#           } else {
#             cat("För få datapunkter med fångst\n")
#           }
#         }
#       }
#       # cat("\n:::\n")
#     }
#     # cat("\n:::\n")
#   }
#   cat("\n:::\n")
# }
# 
# rapportfigurer <- function(path, art, timePeriod, quarter, data){
#     cat("  \n::: panel-tabset  \n")
#     for(i in 1:length(species)) {
#       cat("\n#####", art[i]," \n")
#       cat("\n::: panel-tabset\n")
#       for(j in 1:length(timePeriod)){
#         cat("\n#####", names(timePeriod)[j]," \n")
#         cat("\n::: panel-tabset\n")
#         for(k in 1:length(quarter)){
#           cat("\n#####", names(quarter)[k]," \n")
#           lf <- list.files(path, full.names = T)
#           
#           EMprojPA_path <- paste0(path, gsub(" ", "_", art[i]),"_PA_", names(timePeriod)[j], "_", names(quarter)[k], "_EMprojPA.rds")
#           EMprojPA <- readRDS(EMprojPA_path)
#           
#           d1p <- paste0(path, gsub(" ", "_", art[i]),"_Abu_", names(timePeriod)[j], "_", names(quarter)[k],"_d1_2.rds")
#           m1p <- paste0(path, gsub(" ", "_", art[i]),"_Abu_", names(timePeriod)[j], "_", names(quarter)[k],"_m1_2.rds")
#           p1Trunkp <- paste0(path, gsub(" ", "_", art[i]),"_Abu_", names(timePeriod)[j], "_", names(quarter)[k],"_p1Trunk_2.rds")
#           en1Trunkp <- paste0(path, gsub(" ", "_", art[i]),"_Abu_", names(timePeriod)[j], "_", names(quarter)[k],"_en1Trunk_2.rds")
#           if(all(c(d1p, m1p, p1Trunkp, en1Trunkp) %in% lf)){
#             d1 <- readRDS(d1p)
#             m1 <- readRDS(m1p)
#             p1Trunk <- readRDS(p1Trunkp)
#             en1Trunk <- readRDS(en1Trunkp)
#             enDelta <- (en1Trunk* raster(EMprojPA[[1]]))
#             
#             cat("\n::: panel-tabset\n")
#             cat("\n##### Rapportfigur \n")
#             data[gsub(" ", "", art[i])] <- data[,art[i]]/data$sweptAreaDS_km2
#             master2 <- data |> dplyr::filter(Year >= as.numeric(timePeriod[j][1,1]),
#                                              Year <= as.numeric(timePeriod[j][2,1]),
#                                              Quarter %in% unlist(unique(quarter[k])))
#             spMap <- master2 |>
#               dplyr::select(Year, Month, Quarter, Depth, source.x, ShootLat, ShootLong, !!!syms(species))
#             r <- rast(xmin = 9, xmax = 14, ymin=55, ymax=60, ncols=40, nrows=80)
#             sp_vect <- terra::vect(spMap, geom=c("ShootLong", "ShootLat"))
#             sp_rast <- terra::rasterize(sp_vect, r, field = art[i], fun = mean)
#             # rast_df2 <- dplyr::arrange(na.omit(raster::as.data.frame(sp_rast, xy = TRUE)), desc(names(sp_rast)))
#             
#             rast_df <- dplyr::arrange(na.omit(raster::as.data.frame(enDelta, xy = TRUE)), desc(names(enDelta)))
# 
#             rastExtent <- enDelta
#             # extent(rastExtent) <- c(min(rast_df$x), max(rast_df$x), min(rast_df$y), max(rast_df$y))
#             extent(rastExtent) <- c(9.5, 13.5, 55, 59.5)
#             
#             insetMap2 <- function (map, rast){
#               ins <- sf::st_as_sfc(sf::st_bbox(rast))
#               europe <- rnaturalearth::ne_countries(continent = 'europe', returnclass = "sf", scale = "large")
#               m2 <- ggplot2::ggplot() +
#                 ggplot2::geom_sf(data = europe) +
#                 ggplot2::geom_sf(data = ins, color = "red", fill = NA, size = 2, linewidth = 1) +
#                 xlim(0, 30) +
#                 ylim(50, 70) +
#                 ggthemes::theme_map() +
#                 ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white"),
#                                panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth  = 1))
#               
#               m3 <- cowplot::ggdraw() +
#                 cowplot::draw_plot(map) +
#                 cowplot::draw_plot(m2, x = 0.47, y = 0.705, width = 0.3, height = 0.3)
#               
#               m3
#             }
#             
#             p <- enDelta |> 
#               aqua::rasterPlot(sluCol = "yellow", rev=T, nClass = 5, classMethod = "quantile") |> 
#               aqua::mapExtra()  +
#               coord_sf(xlim = c(9.5,13.5), ylim = c(55,59.5))
#             
#             p <- p |>  insetMap2(rastExtent)
#             
#             plot(p)
#             
#             cat("\n##### Antal per km2 \n")
#             tryCatch({
#               p2 <- sp_rast |> 
#                 aqua::rasterPlot(sluCol = "red", rev=T, nClass = 5, classMethod = "jenks") |> 
#                 aqua::mapExtra() +
#                 geom_sf(data = trawl_area, fill = "transparent", linewidth = 0.5, col = "black") +
#                 coord_sf(xlim = c(9.5,13.5), ylim = c(55,59.5))
#               p2 <- p2 |> insetMap2(rastExtent)
#             }, 
#             error = function(cond){
#               p2 <- sp_rast |> 
#                 aqua::rasterPlot(sluCol = "red", rev=T, nClass = 3, classMethod = "jenks") |> 
#                 aqua::mapExtra() +
#                 geom_sf(data = trawl_area, fill = "transparent", linewidth = 0.5, col = "black") +
#                 coord_sf(xlim = c(9.5,13.5), ylim = c(55,59.5))
#               p2 <- p2 |> insetMap2(rastExtent)
#             })
#             
#             p2
#             cat("\n:::\n")
# 
#           }
#         }
#         cat("\n:::\n")
#       }
#       cat("\n:::\n")
#     }
#     cat("\n:::\n")
#   }
# 
# 
# rapportTabeller <- function(species, path, timePeriod, quarter, master, data){
#   
#   art <- names(species)
#   getData <- function(art, j = 1, k = 1){
#     return(list(
#       myBiomodModelOut = readRDS(paste0(path, gsub(" ", "_", art),"_PA_", names(timePeriod)[j], "_", names(quarter)[k], "_myBiomodModelOut.rds")),
#       myRespPlot = readRDS(paste0(path, gsub(" ", "_", art),"_PA_", names(timePeriod)[j], "_", names(quarter)[k], "_myRespPlot.rds")),
#       myEvalPlot = readRDS(paste0(path, gsub(" ", "_", art),"_PA_", names(timePeriod)[j], "_", names(quarter)[k], "_myEvalPlot.rds")),
#       myVarImpPlot = readRDS(paste0(path, gsub(" ", "_", art),"_PA_", names(timePeriod)[j], "_", names(quarter)[k], "_myVarImpPlot.rds")),
#       myEMRespPlot = readRDS(paste0(path, gsub(" ", "_", art),"_PA_", names(timePeriod)[j], "_", names(quarter)[k], "_myEMRespPlot.rds")),
#       myEMEvalPlot = readRDS(paste0(path, gsub(" ", "_", art),"_PA_", names(timePeriod)[j], "_", names(quarter)[k], "_myEMEvalPlot.rds")),
#       myEMVarImpPlot = readRDS(paste0(path, gsub(" ", "_", art),"_PA_", names(timePeriod)[j], "_", names(quarter)[k], "_myEMVarImpPlot.rds")),
#       myBiomodProj = readRDS(paste0(path, gsub(" ", "_", art),"_PA_", names(timePeriod)[j], "_", names(quarter)[k], "_myBiomodProj.rds")),
#       myBiomodEM = readRDS(paste0(path, gsub(" ", "_", art),"_PA_", names(timePeriod)[j], "_", names(quarter)[k], "_myBiomodEM.rds")),
#       myBiomodEMProj = readRDS(paste0(path, gsub(" ", "_", art),"_PA_", names(timePeriod)[j], "_", names(quarter)[k], "_myBiomodEMProj.rds")),
#       EMprojPA = readRDS(paste0(path, gsub(" ", "_", art),"_PA_", names(timePeriod)[j], "_", names(quarter)[k], "_EMprojPA.rds")),
#       d1p = readRDS(paste0(path, gsub(" ", "_", art),"_Abu_", names(timePeriod)[j], "_", names(quarter)[k],"_d1_2.rds")),
#       m1p = readRDS(paste0(path, gsub(" ", "_", art),"_Abu_", names(timePeriod)[j], "_", names(quarter)[k],"_m1_2.rds")),
#       p1Trunkp = readRDS(paste0(path, gsub(" ", "_", art),"_Abu_", names(timePeriod)[j], "_", names(quarter)[k],"_p1Trunk_2.rds")),
#       en1Trunkp = readRDS(paste0(path, gsub(" ", "_", art),"_Abu_", names(timePeriod)[j], "_", names(quarter)[k],"_en1Trunk_2.rds"))
#     ))
#   }
#   
#   # Tabell 1. Sammanställning av modellresultat för de ensemblemodeller 
#   tab1 <- c()
#   for(i in names(species)){
#     for(j in 1:length(quarter)){
#       dat <- getData(i, k=j)
#       t1 <- data.frame(
#         art = i,
#         säsong = names(quarter)[j],
#         AUC = dat$myBiomodEM@models.evaluation@val$calibration[2],
#         Tröskelvärde = dat$myBiomodEM@models.evaluation@val$cutoff[1],
#         Sensitivitet = dat$myBiomodEM@models.evaluation@val$sensitivity[1],
#         Specificitet = dat$myBiomodEM@models.evaluation@val$specificity[1]
#       )
#       tab1 <- rbind(tab1, t1)
#     }
#   }
#   tab1
#   
#   # Tabell 2. Resultat från test där ensemblemodellen för abundansdata testades genom att fem gånger träna modellen på 75% av data och testa den på resterande 25%. 
#   artStats <- list()
#   artSummary <- c()
#   for(art in names(species)){
#     for(q in 1:length(quarter)){
#       print(art)
#       dat <- getData(art, k = q)
#       antalItterationer <- length(dat$m1@models[[1]][[1]])
#       stats <- c()
#       k = 0
#       for (i in c("glm", "brt", "rf")){
#         for (j in 1:antalItterationer) {
#           print(paste(i, j))
#           tryCatch({
#             modInfo <- dat$m1@models[[gsub(" ", "", art)]][[i]][[as.character(k * antalItterationer + j)]]@evaluation$test.dep
#             mase <- Metrics::mase(modInfo@observed, modInfo@predicted)
#             R2 <- cor(modInfo@observed, modInfo@predicted)^2
#             st <- modInfo@statistics
#             df <- data.frame(mod = i, id = j, cor = st$COR[1], pVal = st$COR[2], dev = st$Deviance, RMSE = st$RMSE,
#                              MAE = st$MAE, MASE = mase, R2 = R2, obs = mean(modInfo@observed), pred = mean(modInfo@predicted))
#             stats <- rbind(stats, df)
#           },
#           error = function(e) { 
#             print("Error")
#           }
#           )
#         }
#         k = k + 1
#       }
#       artStats <- c(artStats, s <- list(stats))
#       artSum <- data.frame(Art = art,
#                            Säsong = names(quarter)[q],
#                            Abu = median(stats$obs), AbuSD = sd(stats$obs),
#                            PreAbu = median(stats$pred), preSD = sd(stats$pred),
#                            R2 = median(stats$R2), R2SD = sd(stats$R2),
#                            MASE = median(stats$MASE), MASESD = sd(stats$MASE))
#       artSummary <- rbind(artSummary, artSum)
#     }
#   }
#   names(artStats) <- names(species)
#   artStats
#   artSummary
#   
#   # Tabell 3. Resultat för den slutgiltiga kartprediktionen
#   # master <- readRDS("WestCoast/R/helpFiles/master_20230610.rds")
#   library(sp)
#   tab3 <- c()
#   for(i in names(species)){
#     # master2 <- master[master[species[[i]]]>0,]
#     for(q in 1:length(quarter)){
#       dat <- getData(i, k = q)
#       
#       ensMask <- dat$en1Trunkp * raster(dat$EMprojPA)
#       
#       rawVals_all <- master |>      
#         dplyr::filter(Year >= as.numeric(2005),
#                       Quarter %in% unlist(unique(quarter[q]))) |> 
#         dplyr::select(ShootLat, ShootLong, all_of(species[[i]]))
#       
#       rawVals_all2 <- rawVals_all
#       coordinates(rawVals_all2) <- ~ ShootLong + ShootLat
#       terra::crs(rawVals_all2) <- terra::crs(ensMask)
#       rasValue <-  terra::extract(ensMask, rawVals_all2)
#       rasValue[is.na(rasValue)] <- 0
#       rawVals_all$mod <- rasValue
#       # coordinates(rawVals_all) <- ~ ShootLong + ShootLat
#       # a <- log10(rawVals_all[[1]]) 
#       # a[a == -Inf] <- 0
#       a <- rawVals_all[[3]]
#       
#       aa <- cbind(rasValue, a) |> na.omit() 
#       
#       lm1 <- lm(aa[,1] ~ aa[,2])
#       cor <- cor(rasValue, a)
#       
#         trawlArea <- rast("C:/Github Projcects/NMK/WestCoast/indata/gis/trawl_area.tif")
#         data[gsub(" ", "", i)] <- data[, i]/data$sweptAreaDS_km2
#         master2 <- data |> dplyr::filter(Year >= as.numeric(timePeriod[1][1, 1]),
#                                          Year <= as.numeric(timePeriod[1][2, 1]),
#                                          Quarter %in% unlist(unique(quarter[1])))
#         enDelta <- (dat$en1Trunkp * raster(dat$EMprojPA[[1]]) * raster(trawlArea))
#         spMap <- master2 |>
#           dplyr::select(Year, Month, Quarter, Depth, source.x, ShootLat, ShootLong, !!!syms(species))
#         r <- rast(xmin = 9, xmax = 14, ymin=55, ymax=60, ncols=40, nrows=80)
#         sp_vect <- terra::vect(spMap, geom=c("ShootLong", "ShootLat"))
#         sp_rast <- terra::rasterize(sp_vect, r, field = i, fun = mean)
#         rastExtent <- enDelta
#         extent(rastExtent) <- c(9.5, 13.5, 55, 59.5)
#         
#         enDelta_resamp <- rast(enDelta)
#         enDelta_resamp[is.na(enDelta_resamp)] <- 0
#         sp_rast_resamp <- terra::resample(sp_rast, enDelta_resamp)
#         r2 <- data.frame(values(enDelta_resamp), values(sp_rast_resamp)) |> na.omit() |> cor() %>% .[1, 2] |> round(3)
#         
#       
#       tab3 <- rbind(tab3, data.frame(i, q, pVal = summary(lm1)$coefficients[2,4], cor, corMap = r2, summary(lm1)$fstatistic[1]))
#     }
#   }
#   tab3 
# 
#   # Tabell 4. Variable importance PA-model
#   tab4 <- c()
#   for(i in names(species)){
#     for(q in 1:length(quarter)){
#       dat <- getData(i, k = q)
#       t3 <- tapply(dat$myBiomodEM@variables.importance@val$var.imp, dat$myBiomodEM@variables.importance@val$expl.var, mean) |> as.data.frame() |> t()
#       rownames(t3) <- i
#       # print(paste(i, length(t3)))
#       tab4 <- rbind(tab4, t3)
#     }
#   }
#   tab4
#   
#   return(list(tab1, artSummary, tab3, tab4))
# }
# 
