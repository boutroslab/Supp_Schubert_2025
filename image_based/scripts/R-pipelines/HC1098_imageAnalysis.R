# Image analysis script as used for HC1094 starting from February 2019 



## load dependecies 
library(stringr)
library(dplyr)
library(dbplyr)
library(EBImage)
library(parallel)
library(FNN)
library(RPostgreSQL)


## open database connection 
sc_db <- src_postgres(
  dbname = "HC1098",
  host = "b110-sc2sn01",
  user = "AntoniaS",
  password = "dvlror2"
)

av_db <- src_postgres(
  dbname = "HC1098",
  host = "b110-sc2sn01",
  user = "AntoniaS",
  password = "dvlror2"
)



## get mock names 
mockNames_cell <- tbl(sc_db,"average_cell") %>% colnames
mockNames_cell_av <- tbl(sc_db,"average_cell") %>% colnames # mockNames_cell[!mockNames_cell == "cell_area"] 
#mockNames_pop <- tbl(sc_db,"single_population") %>% colnames

## define variables for qc
emptyWell <-  F
noImage <-  F
noFeatures <- F

## image analysis and feature extraction are summarized in a function 'analyzeImages'
anlyzeImages <- function(x) {
  print(resDir)
  start.time <- Sys.time()
  images=unlist(strsplit(x=x,split="::"))
 
  ## grab image names     
  DAPI_image_name=images[1] #DAPI_image_name="/data2/raw_data/HC1094/HC10944t/HC10944t_1094_cl01_S01_FDA1_3_hc_2019.03.02.04.38.47/HC10944t_1094_cl01_S01_FDA1_3_hc_K16_1_DAPI.tif"#
  tub_image_name=images[3]
  actin_image_name=images[2] 
 
  ## convert physical barcode to database barcode
  identifier=sub(".*/(.+)_\\w+.tif","\\1",DAPI_image_name,perl=TRUE)


  barcodes=str_match(identifier,"HC10984f_\\w+?_(\\w+?)_(\\w+?)_(\\w+?)_(\\w+?)_(\\w+?)_(\\w+?)_(\\d+)$")
  print(identifier)



  ## unique identifiers for the field of interest are read from the physical barcodes
  field_id <- barcodes[8] 
  well_id <- barcodes[7] 
    
  screen_id <- barcodes[3] 
  cell_line_id <-barcodes[2]  
  cell_line_id <- paste0(substr(cell_line_id,1,2),
                        paste0("0",substr(cell_line_id,3,4)))
  library_id <-barcodes[4] 
  library_plate_id <- barcodes[5] 
  concentration <- barcodes[6] 
  concentration <- substr(concentration,1,1)
  
  
     
    if(file.exists(DAPI_image_name) && 
       file.exists(tub_image_name) && 
       file.exists(actin_image_name)){
        
        img_DAPI <- readImage(DAPI_image_name)
        img_Cy3 <- readImage(actin_image_name)
        img_FITC <- readImage(tub_image_name)
        
        ## test images for artifacts 
        strange_behavior <-  c("dark"=0,"fussel"=0)
        medi=median(img_DAPI)
        midi=mean(img_DAPI) 
        stat_diff <- abs(abs(medi)-abs(midi))
        
        if(medi<0.01){
            strange_behavior["dark"]=1 }
        
        max_sat_DAPI <- length(which(img_DAPI == 1)) / (2094 * 2094)
        max_sat_Cy3 <- length(which(img_Cy3 == 1)) / (2094 * 2094)
        
        if(max_sat_DAPI > 0.0025 | max_sat_Cy3 > 0.0025){
            strange_behavior["fussel"]=1}
        
        ## blur images with gaussian filter 
        img_DAPIsmooth = gblur(img_DAPI, radius = 51, sigma = 1)
        img_Cy3smooth = gblur(img_Cy3, radius  = 51, sigma = 4)
        img_FITCsmooth = gblur(img_FITC, radius  = 51, sigma = 4)
        
        
        
        ## segment nuclei using adaptive thresholding 
        ##  the window sized is chosen based on the size of the nuclei 
        ##      the first threshold is a bist smaller to sparate the nuclei
        ##      the second is larger to fill gaps 
        
        nucleusTresh = thresh(img_DAPIsmooth, 
                                  w = 20, h = 20, 
                                  offset = 0.004)
        
        nucleusTresh = fillHull(opening(nucleusTresh,
                                   kern=makeBrush(9, shape="disc")))
        
        nucleusFill = fillHull(thresh(img_DAPIsmooth, 
                                            w = 30, h = 30, 
                                            offset = 0.001))
        
        nucleusRegions = propagate(img_DAPIsmooth,
                                   seed=bwlabel(nucleusTresh),
                                   mask=nucleusFill)
        
        ## check wheter cells are present 
        if(max(nucleusRegions)>1) {
        
            ## segment cell bodies by adaptive and global thersholding 
            ##      the actin and tubulin channel are used   
            ##      regions with mising cytoplasm signal are filled with 
            ##          the nucleus binary mask and a fillHull command
            
            cytoplasmThresh = thresh(img_Cy3smooth,
                                     w = 28,
                                     h = 28,
                                     offset =  0.001)
            
            cytoplasmOpening = opening(cytoplasmThresh,
                                       kern=makeBrush(9,shape="disc"))
            
            # due to differences in staining intensities the global 
            #   thresholding for the cell bodies has to be somewhat adaptive
            cytoplasmOpening2 = opening(img_Cy3smooth > 0.03)
            
            
            nucleusRegions2 <- nucleusRegions
            nucleusRegions2[nucleusRegions != 0] <- 1
            

            cytoplasmCombined = cytoplasmOpening | cytoplasmOpening2  | nucleusRegions2  

            # obatoclax (in well N04 on the KIN plates) has an optical effect
            #       and has to be segmenetd differentially 
            if(well_id == "N04" && library_id %in% c("KIN1","KIN2","KIN3")){
                cytoplasmCombined = cytoplasmOpening | nucleusRegions2
            }
            
            storage.mode(cytoplasmCombined) = "integer"
    
            
            cytoplasmRegions = propagate(x = img_Cy3smooth,
                                         seeds = nucleusRegions,
                                         lambda=1e-04,
                                         mask=cytoplasmCombined)
            
            cytoplasmRegions = fillHull(cytoplasmRegions)
            
    
            
            
            ## merge channel and save as compressed image with segmentation borders  
            ImgColor = rgbImage(5*normalize(img_Cy3),
                                5*normalize(img_FITC),
                                5*normalize(img_DAPI))
            
            ImgOut = paintObjects(cytoplasmRegions,
                                  paintObjects(nucleusRegions,
                                               ImgColor,
                                               col='blue'),
                                  col='white')
            r <- 1
            attempt <- 1
            while( r == 1 && attempt <= 10 ) {
              attempt <- attempt + 1
              try(
                r <<- 
                  writeImage(ImgOut,
                             file=paste(resDir,"/",identifier,"_segmented",".png",sep=""),
                             type="png")
              )
              #Sys.sleep(3)
            }
	    
            stopifnot(file.exists(paste(resDir,"/",identifier,"_segmented",".png",sep="")))
            
            ## compute features
            F_DAPI = try(computeFeatures(nucleusRegions,
                                     ref=img_DAPI,
                                     xname="nuc",
                                     refnames="nuc",
                                     expandRef = NULL)
                        )

            F_tub = try(computeFeatures(cytoplasmRegions,
                                    ref=img_FITC,
                                    xname="cell",
                                    refnames="dsh",
                                    expandRef = NULL,
 				    methods.noref = c())
                        )
            
            F_actin = try(computeFeatures(cytoplasmRegions,
                                      ref=img_Cy3,
                                      xname="cell",
                                      refnames="act",
                                      expandRef = NULL)
                        )

            
            # remark: the subtraction of the mean is a normalization
            #         required for combining the images
            F_DAPI_tub = try(computeFeatures(cytoplasmRegions,
                                         ref=(img_DAPI - mean(img_DAPI)) * (img_FITC - mean(img_FITC)),
                                         xname="cell", refnames="nucdsh",
                                         expandRef = NULL,
 					 methods.noref = c())
                        )

	    ## stop image analysis if feature extraction failed 
            if(any(
                inherits(F_DAPI,"try-error",which=F),
                inherits(F_actin,"try-error",which=F),
                inherits(F_tub,"try-error",which=F),
                inherits(F_DAPI_tub,"try-error",which=F))
            ) {
                noFeatures <- T }

            ## get additional features
            ##      nuclear displacement
            ##      nearest nieghbour info

            nucleus.displacement = sqrt(
                rowSums(
                    (F_actin[,c("cell.act.m.cx", "cell.act.m.cy")] - F_DAPI[,c("nuc.nuc.m.cx", "nuc.nuc.m.cy")])^2
                )
            )

            names(nucleus.displacement) = c("nuclear.displacement")

            nearest.neighbours = try(get.knn(F_actin[,c("cell.act.m.cx", "cell.act.m.cy")],k=30)[["nn.dist"]][,c(10,20,30)],silent=T)
            if(inherits(nearest.neighbours,"try-error",which=F)){
                nearest.neighbours = cbind(dist.10.nn = rep(NA,nrow(F_DAPI)),
                                           dist.20.nn = rep(NA,nrow(F_DAPI)),
                                           dist.30.nn = rep(NA,nrow(F_DAPI)) )
            } else {
                colnames(nearest.neighbours) = c("dist.10.nn","dist.20.nn","dist.30.nn")
            }







            ## a data frame with all single cell features is created
            F_all <- cbind.data.frame(cell_line_id = rep(cell_line_id,nrow(F_DAPI)),
                                      screen_id = rep(screen_id,nrow(F_DAPI)),
                                      library_id = rep(library_id,nrow(F_DAPI)),
                                      library_plate_id = rep(library_plate_id,nrow(F_DAPI)),
                                      concentration = rep(concentration,nrow(F_DAPI)),
                                      well_id = rep(well_id,nrow(F_DAPI)),
                                      field_id = rep(field_id,nrow(F_DAPI)),
                                      F_DAPI,F_tub,F_actin,
                                      F_DAPI_tub,
                                      nucleus.displacement,nearest.neighbours,
                                      cell_area = F,
                                      stringsAsFactors = F
            )





    ################################################################################
    ################################################################################

            # summarize single cell features per well
            # Note: all objects wehre to area of the cell body and nucleaus are 
            #           identical will be excluded from the averages 
            featureNames <- colnames(F_all)[!colnames(F_all) %in% c("cell_line_id",
                                                                    "screen_id",
                                                                    "library_id",
                                                                    "library_plate_id",
                                                                    "concentration",
                                                                    "well_id",
                                                                    "field_id",
                                                                    "cell_area")]

            F_all_sum <- F_all[,featureNames]

            F_all_sum <- c(apply(F_all_sum[F_all$cell.0.s.area != F_all$nuc.0.s.area,],2,mean,trim=0.01, na.rm=TRUE),
                           apply(F_all_sum[F_all$cell.0.s.area != F_all$nuc.0.s.area,],2,sd,na.rm=TRUE),
                           nrow(F_all[F_all$cell.0.s.area != F_all$nuc.0.s.area,]),
                           nrow(F_all[F_all$cell.0.s.area == F_all$nuc.0.s.area,]))

            names(F_all_sum) <- c(paste("mean_",featureNames,sep=""),
                                  paste("sd_",featureNames,sep=""),
                                  "cell_count",
                                  "dead_cells")
            
            dead_cells <- F_all_sum["dead_cells"]
            
            F_all_sum <- F_all_sum[!names(F_all_sum) %in% "dead_cells"] 
            
            F_all_sum <- cbind.data.frame(cell_line_id = cell_line_id,
                                          screen_id = screen_id,
                                          library_id = library_id,
                                          library_plate_id = library_plate_id,
                                          concentration = concentration,
                                          well_id = well_id,
                                          field_id = field_id,
                                          as.list(F_all_sum),
                                          fussel = strange_behavior["fussel"],
                                          dark = strange_behavior["dark"],
                                          dead_cells,
                                          stringsAsFactors=F)

            
            # add the object count from the Dsh signal 
            puncTresh_global <- img_FITCsmooth > 0.025
            punct_count <- length(unique(bwlabel(puncTresh_global)[bwlabel(puncTresh_global)!=0]))
            F_all_sum$punct_count <- rep(punct_count,nrow(F_all_sum))






            ########################################################################
            r <- 1
            attempt <- 1
            while( r == 1 && attempt <= 10 ) {
              attempt <- attempt + 1
              try(
                r <<- write.table(F_all_sum,file=paste0(resDir,"/",identifier,".tab"),sep="\t",quote=FALSE,row.names = FALSE)
              )
             # Sys.sleep(3)
            } 
            stopifnot(file.exists(paste0(resDir,"/",identifier,".tab")))
	    
	    db_insert_into(con = av_db$con, table = "average_cell", values=F_all_sum)
   
          
	        
            ########################################################################

            remove(img_DAPI)
            remove(img_Cy3)
            remove(img_FITC)
            remove(img_DAPIsmooth)
            remove(img_Cy3smooth)
            remove(img_FITCsmooth)
            remove(ImgColor)
            remove(ImgOut)

            remove(F_DAPI)
            remove(F_tub)
            remove(F_actin)
            remove(F_DAPI_tub)
            remove(nucleus.displacement)
            remove(nearest.neighbours)
            remove(strange_behavior)
            remove(F_all)
            remove(F_all_sum)
            
        } else { 
            emptyWell <-  T}
            
    } else {
    noImage <- T
   }

        
###################################################################################################        
# write table with all NAs in the database if image analysis failed 
  


    #if(any(c(emptyWell,noImage,noFeatures))) {
     if(!file.exists(paste0(resDir,"/",identifier,".tab"))){       
                F_all <- data.frame(matrix(ncol = length(mockNames_cell))) 
                colnames(F_all) <- mockNames_cell
                F_all$cell_line_id <- cell_line_id
                F_all$screen_id <- screen_id
                F_all$library_id <- library_id
                F_all$library_plate_id <- library_plate_id
                F_all$concentration <- concentration
                F_all$well_id <-  well_id
                F_all$field_id <-  field_id
                
                featureNames <- colnames(F_all)
                
                # +2 for cell count and dead cells 
                F_all_sum <- data.frame(matrix(ncol = length(featureNames)))
                
                names(F_all_sum ) <- featureNames
                
                ####################################################################
                r <- 1
                attempt <- 1
                while( r == 1 && attempt <= 10 ) {
                  attempt <- attempt + 1
                  try(
                    r <<- write.table(F_all_sum,file=paste0(resDir,"/",identifier,".tab"),sep="\t",quote=FALSE,row.names = FALSE)
                  )
                #  Sys.sleep(3)
                } 

		stopifnot(file.exists(paste0(resDir,"/",identifier,".tab")))
		
          	db_insert_into(con = av_db$con, table = "average_cell", values=F_all_sum)

            	  
              
                
                ####################################################################
                
                
                remove(img_DAPI)
                remove(img_Cy3)
                remove(img_FITC)
                remove(img_DAPIsmooth)
                remove(img_Cy3smooth)
                remove(img_FITCsmooth)
                remove(ImgColor)
                remove(ImgOut)
                
                remove(F_DAPI)
                remove(F_tub)
                remove(F_actin)
                remove(F_DAPI_tub)
                remove(nucleus.displacement)
                remove(nearest.neighbours)
                remove(strange_behavior)
                remove(F_all)
                remove(F_all_sum)
        }
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)
}


emptyWell <-  F
noImage <-  F

args=commandArgs(trailingOnly = TRUE)
options("scipen"=100, "digits"=4,warn=-1)
resDir=args[length(args)];
args=as.list(args[-length(args)])
result=lapply(as.list(args),anlyzeImages)



