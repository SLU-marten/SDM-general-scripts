# This file creates a presence/absence models and a abundance models for all selected species
# The file delta_model_report.qmd creates a html-report of the created models

# DATA PREPARATIONS ####
## Load packages ####
library(dplyr) # Data handling
library(terra) # Raster handling
library(sdm) # Abundance modeling
library(biomod2) # Presence/absence modeling
library(quarto)

# Create folder structure
modelID <- "20240910"
pathPA <- paste0(getwd(), "/model output/", modelID, "/PA")
dir.create(pathPA, recursive = TRUE)

## Select species from the loaded data file ####
species <- c("Plaice" = "Pleuronectes platessa", 
             "Cod" = "Gadus morhua"
             # "Herring" =  "Clupea harengus", 
             # "ThornySkate" = "Amblyraja radiata"
             )

## Select predictors from the loaded raster stack ####
predictors <-   c("bTemp", "bSal", "Depth", "botCur", "slope", "fishInt", "maxT")

## Load data ####
df <- readRDS("data/master_sample1K.rds") # The prepared data with catches and predictor values for each sample point
predStack <- readRDS("data/predStack_20240118.rds") # Raster stack of all predictor variables
predStack <- predStack[[predictors]]

# Creating a presence/absence data.frame of species occurence
dfPA <- dplyr::select(df, Year, Quarter, ShootLong, ShootLat, all_of(predictors), all_of(species)) |> 
  mutate(across(names(species), ~replace(., . > 0, 1)))

# Creating a data.frame for the abundance modeling with the variable n catches per swept area
dfAbu <- dplyr::select(df, Year, Quarter, sweptAreaDS_km2, all_of(predictors), all_of(species)) |> 
  mutate(across(names(species), ~ as.integer(.x / sweptAreaDS_km2)))

## Model parameters ####
timePeriod <- tibble("2005_2022" = c(2005, 2022)) # Which years should be included
quarter <- list("AllYear" = c(1,2,3,4)) # Which quarters should be included (Can be changed to months)
nRepPA_CV <- 2 # Number of repetitions of cross-validation points that will be drawn in the individual models
nRepPA_VI <- 2 # Number of permutations for to estimate variable importance in the individual models
nRepPAEns_VI <- 2 # Number of permutations for to estimate variable importance in the ensemble model
nRepAbu_SubSamp <- 2 # Number of subsampling in the abundance modeling
formula <- paste(" ~", paste(predictors, collapse = " + "))
biomodMethods <- c("GLM", "GBM", "GAM", "RF", "XGBOOST")
sdmMethods <- c('glm', 'brt', 'rf')
coordNames <- c("ShootLong", "ShootLat")

# Write file with selected model parameters
modPar <- c(paste("Timeperiod:", paste(timePeriod, collapse = ", ")),
            paste("Quarter:", paste(quarter, collapse = ", ")),
            paste("nRepPA_CV:", nRepPA_CV),
            paste("nRepPA_VI:", nRepPA_VI),
            paste("nRepPAEns_VI:", nRepPAEns_VI),
            paste("nRepAbu_SubSamp:", nRepAbu_SubSamp),
            paste("Formula:", formula),
            paste("Biomod modeltypes:", paste(biomodMethods, collapse = ", ")),
            paste("SDM modeltypes:", paste(sdmMethods, collapse = ", ")))

write.table(modPar, paste0(getwd(), "/model output/", modelID, "/modelParameters.txt"),
            col.names = FALSE, row.names = FALSE, quote = FALSE)

# MODELING ####
# Function for filtering out selected time period and quarter
dfFilter <- function(data, tp, q){
  dplyr::filter(data,
                Year >= as.numeric(timePeriod[tp][1, 1]), 
                Year <= as.numeric(timePeriod[tp][2, 1]),
                Quarter %in% unlist(unique(quarter[q]))
                )
}
# Loops and creates models for all quarters in all time periods for all species
for(tp in 1:length(timePeriod)){
  for(q in 1:length(quarter)){ 
    # Filtering out selected species, time period and quarter
    dfPAFilt <- dfFilter(dfPA, tp, q)
    dfAbuFilt <- dfFilter(dfAbu, tp, q)
    for(sp in names(species)){
      # Error handling. If one species/time period/quarter model fails it creates a text file with the error message and then the loop continues
      tryCatch({
        ## PRESENCE/ABSENCE ####
        ### Data preparations ####
        # Preparing data for biomod2. Here it's possible to set pseudo absence settings and predefined data for validation
        myBiomodData <- BIOMOD_FormatingData(resp.var = dfPAFilt[sp], # Data.frame of species presence absence data
                                             expl.var = dfPAFilt[, c(predictors)], # SpatRaster or data.frame of the predictors
                                             resp.xy = dfPAFilt[, coordNames], # Data.frame of the coordinates
                                             resp.name = sp, #Vector of species names
                                             dir.name = pathPA # Path of model output
                                             )
        ### Modeling ####
        #### Individual models ####
        # Individual modeling function
        # This example uses the strategy "bigboss", which is an optimized model tuning by the biomod team. You can also choose "default"
        # to use the default setting from each model. Other options are "user.defined" or "tuned" that lets you define your own model settings
        # and CV-settings for the individual models with the bm_ModelingOptions() function.
        myBiomodModelOut <- BIOMOD_Modeling(bm.format = myBiomodData,
                                            modeling.id = paste(sp, names(timePeriod)[tp], names(quarter)[q], sep = "_"),
                                            models = biomodMethods, # Included individual model techniques
                                            OPT.strategy = "bigboss",
                                            CV.nb.rep = nRepPA_CV,
                                            CV.perc = 0.75, # Cross validation split percentage
                                            CV.do.full.models = TRUE, # Defining whether models should be also calibrated and validated over the whole dataset
                                            var.import = nRepPA_VI,
                                            metric.eval = c("TSS","ROC"), # Validation methods
                                            do.progress = FALSE # Should a progress bar be displayed in the console
                                            # nb.cpu = 6 # Number of cpu-cores to be used
                                            )

        #### Ensemble model ####
        myBiomodEM <- BIOMOD_EnsembleModeling(bm.mod = myBiomodModelOut,
                                              models.chosen = 'all',
                                              em.by = 'all',
                                              em.algo = c('EMwmean'), # Individual model weight in the ensemble
                                              EMwmean.decay = 'proportional', # A value defining the relative importance of the weights
                                              metric.select = c('ROC'), # Use the ROC score to define which models to be included in the ensemble
                                              metric.select.thresh = c(0.7), # Only include models with ROC > 0.7
                                              var.import = nRepPAEns_VI,
                                              metric.eval = c('TSS', 'ROC'),
                                              do.progress = FALSE,
                                              # nb.cpu = 6 # Number of cpu-cores to be used
                                              )

        ### Predictions ####
        #### Individual models ####
        myBiomodProj <- BIOMOD_Projection(bm.mod = myBiomodModelOut,
                                          proj.name = paste(sp, names(timePeriod)[tp], names(quarter)[q], sep = "_"),
                                          new.env = rast(predStack),
                                          models.chosen = 'all',
                                          metric.binary = 'all',
                                          metric.filter = 'all',
                                          build.clamping.mask = T
                                          )

        #### Ensemble model ####
        myBiomodEMProj <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM,
                                                     bm.proj = myBiomodProj,
                                                     models.chosen = 'all',
                                                     metric.binary = 'all',
                                                     metric.filter = 'all'
                                                     )

        ## ABUNDANCE ####
        # Create folder structure
        pathAbu <- paste(getwd(), "model output", modelID, "Abu", sp, sep = "/")
        if(!dir.exists(pathAbu)) dir.create(pathAbu, recursive = TRUE)
        ### Data preparations ####
        form <- as.formula(paste0(sp, formula)) # Create model formula
        sdm_D <- sdmData(form,
                         train = dfAbuFilt,
                         filename = paste0(pathAbu, "/", sp, "_", names(timePeriod)[tp], "_", names(quarter)[q])
        )
        ### Modeling ####
        sdm_M <- sdm(form,
                     data = sdm_D,
                     methods = sdmMethods,
                     replication = c("sub"),
                     test.percent = 25,
                     n = nRepAbu_SubSamp,
                     modelTypes = c('abundance'),
                     parallelsetting = list(ncore = 4, method = "parallel")
                     # modelsettings = list(glm = list(family = "quasipoisson"),
                     #                      brt = list(n.trees = 500, distribution = 'quasipoisson'))
        )
        ### Prediction ####
        
        # Individual models
        sdm_P <- predict(sdm_M, stack(rast(predStack)))
        
        # Ensemble prediction weighted by the individual models correlation value
        sdm_E <- ensemble(sdm_M, sdm_P, setting = list(method = 'weighted', stat = 'cor', expr = 'cor > 0'))
        
        # Saving abundance model output
        saveRDS(sdm_D, paste0(pathAbu, "/", sp, "_", names(timePeriod)[tp], "_", names(quarter)[q], "_sdm_D.rds"))
        saveRDS(sdm_M, paste0(pathAbu, "/", sp, "_", names(timePeriod)[tp], "_", names(quarter)[q], "_sdm_M.rds"))
        saveRDS(sdm_P, paste0(pathAbu, "/", sp, "_", names(timePeriod)[tp], "_", names(quarter)[q], "_sdm_P.rds"))
        writeRaster(sdm_E, paste0(pathAbu, "/", sp, "_", names(timePeriod)[tp], "_", names(quarter)[q], "_Ensemble_model.tif"))

        ## DELTA ####
        
        # Create folder structure
        pathDelta <- paste(getwd(), "model output", modelID, "Delta/", sep = "/")
        if(!dir.exists(pathDelta)) dir.create(pathDelta, recursive = TRUE)
        
        # Importing PA mask from the PA model
        paMask <- rast(myBiomodEMProj@proj.out@link[2]) # Selecting PA map from TSS threshold
        paMask[paMask == 0] <- NA
        
        # Constrains the abundance model according to the PA model's areas of presence
        delta_M <- sdm_E * paMask
        
        # Save the delta model
        writeRaster(delta_M, paste0(pathDelta, sp, "_", names(timePeriod)[tp], "_", names(quarter)[q],"_deltamodel.tif"))
      },
      error = function(e){
        errorMsg <- capture.output(e)
        
        # Create folder structure
        pathError <- paste(getwd(), "model output", modelID, "Error log", sep = "/")
        if(!dir.exists(pathError)) dir.create(pathError, recursive = TRUE)
        
        # Write an error file
        write.table(errorMsg, paste0(pathError, "/", sp, "_", names(timePeriod)[tp], "_", names(quarter)[q],"_errorlog.txt"))
      })
    }
  }
}
