# explore and analyse Borneo Camera trap data ####

# clear environment ####

rm(list = ls())

# load libraries ####

library(ggplot2)
library(ade4)
library(data.table)
library(gridExtra)
library(AER)
library(MuMIn)
library(AICcmodavg)


# load data ####

detections <- rbind(read.csv("data/LEWS_CT_data.csv"), read.csv("data/SabahParks_CT_data.csv"), read.csv("data/Borneo_LEWS_ComForest_Jan20-Apr21.csv"))

muntjacs <- fread("data/Sabah_data_with_more_identified_Muntjacs.csv", )

meta_for_coordinates <- rbind(read.csv("data/p197_deployments.csv"), read.csv("data/p258_deployments.csv")) # to use coordinates because some are wrong in the detections - these files were downloaded on WI website

detection_distance <- read.csv("data/borneodetectiondistances.csv")

Sabah15_meta <- read.csv("data/Sabah_CT_metadata_2015.csv", fileEncoding = "UTF-8-BOM")
Sabah19_meta <- read.csv("data/Sabah CT metadata 2019.csv", fileEncoding = "UTF-8-BOM")

Lanjak_meta <- read.csv("data/Lanjak_Transect_metadata.csv")

Borneo_LEWS_meta <- read.csv("data/Borneo_LEWS_ComForest_Meta.csv")

species_grouping <- read.csv("data/Species List_Borneo_Rev12.21.22ogc.csv")

functional_trails <- read.csv("data/functionaltraits.csv")

endemic <- read.csv("data/endemic.csv")

# body_size <- read.csv("data/borneo_weights.csv")

# prepare and explore data ####

## open a pdf to keep all figures together

# pdf("results/all_figures.pdf", width = 10, height = 7)

# lanjak is a mess and I need to do generate a new meta file because one row in the meta corresponds to multiple deployment names...
Lanjak_meta_new <- detections
Lanjak_meta_new$elevation <- Lanjak_meta$Elevation[match( substr(detections$Deployment.Name, 1,3), Lanjak_meta$new_name)]
Lanjak_meta_new <- unique(Lanjak_meta_new[, c("Deployment.ID", "Deployment.Name", "elevation")])
Lanjak_meta_new <- Lanjak_meta_new[complete.cases(Lanjak_meta_new),]

# make one big metadata file
all_meta <- data.frame(data = c(rep("Lanjak", nrow(Lanjak_meta_new)),
                                rep("Sabah15", nrow(Sabah15_meta)),
                                rep("Sabah19", nrow(Sabah19_meta)),
                                rep("Borneo_LEWS", nrow(Borneo_LEWS_meta))),
                       Camera_ID = c(Lanjak_meta_new$Deployment.Name, Sabah15_meta$camera, Sabah19_meta$cam, Borneo_LEWS_meta$Camera_Location),
                       Deployment.ID = c(Lanjak_meta_new$Deployment.ID, Sabah15_meta$Deployment.ID, Sabah19_meta$Deployment.ID, Borneo_LEWS_meta$Deployment.ID),
                       elevation  = c(Lanjak_meta_new$elevation, Sabah15_meta$elev, Sabah19_meta$elev, Borneo_LEWS_meta$elev_m_asl))
all_meta$project <- factor(detections$Project[match(all_meta$Deployment.ID, detections$Deployment.ID)])

str(all_meta)


# check if there are duplicated Deployment.ID
all_meta$Deployment.ID[duplicated(all_meta$Deployment.ID)] # "" are okay because those are the ones with no data

all_meta$Camera_ID[duplicated(all_meta$Camera_ID)] # if only "CF37" "CF39" are showing up --> is fine


all_meta[all_meta$Deployment.ID %in% "d62215",] # these were originally an issue but we resolved it
all_meta[all_meta$Deployment.ID %in% "d62220",] # these were originally an issue but we resolved it

all_meta[all_meta$Camera_ID %in% "CF37",] # I thihk this is ok
all_meta[all_meta$Camera_ID %in% "CF39",] # I thihk this is ok





# prepare detection data ####

# add elevation 
detections$elevation <- all_meta$elev[match(detections$Deployment.ID, all_meta$Deployment.ID)]

all(table(unique(detections[, c("Project", "Actual.Lon", "Actual.Lat", "Deployment.ID", "elevation")])$Deployment.ID) == 1)# should be TRUE
unique(detections[is.na(detections$elevation), c("Project", "Subproject", "Deployment.Name", "Deployment.ID", "elevation")]) # should be EMTPY

# replace lat and lon (some are wrong in detections)
detections[, c("Actual.Lon", "Actual.Lat")] <- meta_for_coordinates[match(detections$Deployment.ID, gsub("^.*_", "", meta_for_coordinates$deployment_id)), c("longitude", "latitude")]

nrow(unique(detections[, c("Actual.Lon", "Actual.Lat")])) # 88 locations
table(unique(detections[, c("Actual.Lon", "Actual.Lat", "Project")])$Project) # 45 locations at Sabah, 43 at LEWS

nrow(detections) # 25402 detections
table(detections$Project) # 11424 detections at Sabah, 13978 at LEWS

tapply(detections$elevation, detections$Project, range)

ggplot(unique(detections[, c("Project", "Actual.Lon", "Actual.Lat", "Deployment.ID", "elevation", "Deployment.Name")]), aes(x = Actual.Lon, y = Actual.Lat, colour= log(elevation), shape = Project)) +
  geom_point() + facet_wrap( ~ Project, scales = "free") +
  scale_shape_manual(values=c(16,15)) + theme(legend.position="bottom", aspect.ratio=1)

ggsave("results/Sampling_locations.png", width = 10, height = 7, units = "in", dpi = 300)



# add detection distances 
detections$detection_distance <- detection_distance$Score[match(detections$Deployment.ID, detection_distance$Deployment.ID)]
if(any(is.na(detections$detection_distance))) stop("some deployments are missing detection distance")

prop.table(table(unique(detections[, c("Deployment.ID", "detection_distance")])$detection_distance))*100 # 67% of cameras were ranked as a "1"
table(unique(detections[, c("Deployment.ID", "detection_distance")])$detection_distance)

## drop deployments with detection distance = 4 
detections <- detections[detections$detection_distance < 4, ]

# format date and time ####
setDT(detections)

detections[, date:=as.Date(substr(Begin.Time, 1, 10)), ]
detections[, timestamp:=as.POSIXct(Begin.Time, "%Y-%m-%dT%H:%M:%S", tz = "UTC"),]

detections[is.na(timestamp),] # only 2 NAs... dropping them (that was the case before we removed the detection distance 4)
detections <- detections[!is.na(timestamp),]


detections[, range(timestamp), by = Project] # between January 2015 and June 2019 at Sabah and between January 2020 and November 2021 at LEWS 

# replace unknown Muntjacs species with what Andy further identified

unk_muntjacs <- detections[Common.Name %in% "Muntiacus species", ]


table(muntjacs$common_name) # a muntjac was changed to a sambar

intersect(unk_muntjacs$Sequence.ID, muntjacs$sequence_id) # sequences IDs have nothing to do!!! WI is annoying.... need to match by timestamp and location


nrow(unk_muntjacs)
nrow(muntjacs) #509 observations were changed by Andy

unk_muntjacs$timestamp

muntjacs[, timestamp := as.POSIXct(start_time, "%m/%d/%Y %H:%M:%S", tz = "UTC")]

table(muntjacs$timestamp %in% unk_muntjacs$timestamp ) # we want to see 510 --> Good!



unique(detections[match(paste(gsub("^.*_", "", muntjacs$deployment_id), muntjacs$timestamp), paste(detections$Deployment.ID, detections$timestamp) ), ]) # the chevrotain is because there was a muntjac and a chevrotain in the same picture --> need to address that


idx_Muntjacs_to_fix <- detections$Common.Name %in% "Muntiacus species" & paste(detections$Deployment.ID, detections$timestamp) %in% paste(gsub("^.*_", "", muntjacs$deployment_id), muntjacs$timestamp)

sum(idx_Muntjacs_to_fix) # 511 because one seqience has 2 muntjac (an adult and a young)

detections$Common.Name[idx_Muntjacs_to_fix] <- muntjacs$common_name[match(paste(detections$Deployment.ID, detections$timestamp)[idx_Muntjacs_to_fix], paste(gsub("^.*_", "", muntjacs$deployment_id), muntjacs$timestamp))]

detections$Species.Name[idx_Muntjacs_to_fix] <- paste(muntjacs$genus, muntjacs$species)[match(paste(detections$Deployment.ID, detections$timestamp)[idx_Muntjacs_to_fix], paste(gsub("^.*_", "", muntjacs$deployment_id), muntjacs$timestamp))]

unique(detections$Common.Name[match(paste(gsub("^.*_", "", muntjacs$deployment_id), muntjacs$timestamp), paste(detections$Deployment.ID, detections$timestamp) )])# look good (no unkown muntjacs) - and I chjecked further and the replacement worked


# change detection distance to 3 for d72105	after 2020-08-13 because it was knocked off (repositioneed for the October session - d74311)]

detections[Deployment.ID%in%"d72105" & date>=as.Date("2020-08-13"), detection_distance := 3]


# add species group
detections$SpGroup <- species_grouping$Group[match(detections$Species.Name, species_grouping$Species)]
# unique(detections[is.na(detections$SpGroup), c("Common.Name", "Species.Name")])


# add functional traits (that includes body weight)
setdiff(functional_trails$Scientific, detections$Species.Name) # list of not used in Bill's list
detections[, names(functional_trails)[-1]] <- functional_trails[match(detections$Species.Name, functional_trails$Scientific), -1]
# add body weight
# detections$body_weight <- body_size$weight[match(tolower(detections$Common.Name), body_size$species)]

# add endemic
setdiff(endemic$Common.Name, detections$Common.Name) # Should be empty
setdiff(detections$Common.Name, endemic$Common.Name) # should only be species we don't care abput

detections$endemic <- endemic$endemic[match(detections$Common.Name, endemic$Common.Name)]

# remove detection that are within 'min_cut' minutes of each other (for same species and same deployment ID)
min_cut = 30

detections <- detections[order(timestamp),,] # order by time stamp

detections[, timecut:=cut(timestamp, breaks = paste(min_cut, "min"), labels = F), by = list(Deployment.ID, Common.Name)]  # split timestamp into 30 minutes interval, for each deployment and species
detections <- detections[!duplicated(cbind(timecut, Deployment.ID, Common.Name)), ,] # keep only one "cut" for each deployment and species

any(duplicated(detections[, list(timecut, Deployment.ID, Common.Name),])) # should be FALSE


# cut deployments into 30 days sections

range(detections[, round(difftime(max(timestamp), min(timestamp), units = "days")), by = list(Deployment.ID)]$V1) # 0-783 days

day_cut = 30


detections[, daycut:=cut(timestamp, breaks = paste(day_cut, "day"), labels = F), by = list(Deployment.ID)]  # split timestamp into 30 minutes interval, for each deployment and species

detections <- merge(detections, 
                    detections[, .(startDate=min(timestamp),
                                   endDate=max(timestamp),
                                   Ndays=round(difftime(max(timestamp), min(timestamp), units = "days"))), by =list(Deployment.ID, daycut)],
                    by = c("Deployment.ID", "daycut"))

#  + remove any section < 10 days
min_days = 10

detections <- detections[Ndays>=min_days,,]

table(unique(detections[,.( Deployment.ID, daycut, Project)])$Project)


# get histogram of elevation per project
ggplot(unique(detections[,.(Project, Deployment.ID, elevation),])) + geom_violin(aes(Project, elevation))  + geom_point(aes(Project, elevation, col = elevation))

ggsave("results/Violin_elevation_per_project.png", width = 10, height = 7, units = "in", dpi = 300)


#list unknown species and those with too few detection locations, also list rat speciesm, and squirel species
min_location = 10 # min number of detection location for a species to be included (we are applying this for each site sepearetly)

x <- table(detections$Common.Name, detections$Deployment.ID)
x[x>0] = 1

min10dectect <- names(which(rowSums(x) <min_location ))

x <- table(detections$Common.Name, detections$Deployment.ID, detections$Project)
x[x>0] = 1
min10dectectEachSite <- names(which(rowSums(x[,,1]) <min_location | rowSums(x[,,2]) <min_location ))


unk_species_to_remove <- c("Animal Not on List",
                       "Calibration Photos",
                       "Camera Misfire",
                       "Camera Trapper",
                       "Human non-staff",
                       "No Animal",
                       "Other Bird species",
                       "Unknown Animal",
                       "Unknown Bird", "Unknown Cervid", "Unknown Chevrotain", "Unknown Civet_Genet_Linsang",
                       "Unknown Felid", "Unknown Flying Squirrel", "Unknown Ground-Dove",
                       "Unknown Ground Squirrel", "Unknown Langur", "Unknown Large Weasel_Badger_Otter",
                       "Unknown Macaque", "Unknown Mongoose", "Unknown Porcupine", "Unknown Small Rodent",
                       "Unknown Squirrel", 
                       "Muntiacus species", "Muntjac Species", "Rattus species", "Reptile species", "Treeshrew species"
                       
) 

rat_species <- c("Long-tailed Giant Rat",
                 "Moonrat",
                 "Mountain Treeshrew",
                 "Striped Treeshrew" ,
                 "Treeshrew species",
                 "Rattus species",
                 "Unknown Small Rodent")

squirrel_species <- c("Common Giant Flying Squirrel", 
                      "Horse-tailed Squirrel", 
                      "Pale Giant Squirrel", 
                      "Plantain Squirrel",
                      "Prevost's Squirrel",
                      "Three-striped Ground Squirrel",
                      "Tufted Ground Squirrel",
                      "Unknown Flying Squirrel",
                      "Unknown Squirrel")



# looking at species richness per elevation band (all known species) ####
detections_for_SRband <- detections[!detections$Common.Name %in% unk_species_to_remove,]
detections_for_SRband$elevationBand <- cut(detections_for_SRband$elevation, breaks = min(detections_for_SRband$elevation) + 200 * c(0:(max(detections_for_SRband$elevation)/200)), include.lowest = T, dig.lab = 5)

detections_for_SRband <- unique(detections_for_SRband[, .(Project, elevationBand, Common.Name, Deployment.ID, endemic)])

detections_for_SRband <- detections_for_SRband[, .(SpeciesRichness = .N, EndemismRate = sum(endemic)/.N), by = .(Project, elevationBand, Deployment.ID)]


ggplot(detections_for_SRband, aes(y = elevationBand, x = SpeciesRichness)) + 
  geom_boxplot(position = ) +
  xlab("Species Richness") +
  ylab("Elevation") +
  facet_wrap(~Project)
  
ggsave("results/Boxplot_SR_per_elevation.png", width = 10, height = 7, units = "in", dpi = 300)

ggplot(detections_for_SRband, aes(y = elevationBand, x = EndemismRate)) + 
  geom_boxplot(position = ) +
  xlab("Endemism rate") +
  ylab("Elevation") +
  facet_wrap(~Project)

ggsave("results/Boxplot_Endemism_per_elevation.png", width = 10, height = 7, units = "in", dpi = 300)

# looking at species richness ~ elevation (all known species) ####
detections_for_SR <- detections[!detections$Common.Name %in% unk_species_to_remove,]

detections_for_SR <- unique(detections_for_SR[, .(Project, elevation, Common.Name, Deployment.ID, endemic)])

# detections_for_SR <- detections_for_SR[, .SD[unique(detections[, .(Deployment.ID, Ndays)])[, .(Ndays = sum(Ndays)), by = Deployment.ID], on=c("Deployment.ID")]]


detections_for_SR <- detections_for_SR[, .(SpeciesRichness = .N, EndemismRate = sum(endemic)/.N), by = .(Project, elevation, Deployment.ID)]


detections_for_SR$Ndays <- with(unique(detections[,.(Deployment.ID, daycut, Ndays)]), tapply(Ndays, Deployment.ID, sum))[detections_for_SR$Deployment.ID]


p <- ggplot(detections_for_SR, aes( x = elevation, y = SpeciesRichness)) +
  geom_point(aes(fill = Project, colour = Project, shape = Project), alpha = .5) + 
  labs(y = "Species richness") + 
  scale_shape_manual(values=c(21,22), labels = c("Sabah", "LEWS")) + 
  scale_colour_discrete(labels = c("Sabah", "LEWS")) + 
  scale_fill_discrete(labels = c("Sabah", "LEWS"))

max_elev = max(detections_for_SR$elevation) # 1500m elevation cuttoff or no cutoff max(elevation) 

mSRP <- glm(SpeciesRichness ~  Project*poly(elevation,2) - poly(elevation,2) - Project, offset = log(Ndays), data = detections_for_SR[elevation < max_elev,], family = poisson())

mSR <- glm(SpeciesRichness ~ poly(elevation,2), offset = log(Ndays), data = detections_for_SR[elevation < max_elev,], family = poisson())

(overdispersion = mSRP$deviance/mSRP$df.residual) ; AER::dispersiontest(mSRP)# a little bit over 1... not too bad but let's use it when we calculated uncertainty or compare models

# Compare with smaller models
mSRP1 <- glm(SpeciesRichness ~  Project*poly(elevation,1) - poly(elevation,1) - Project, offset = log(Ndays), data = detections_for_SR[elevation < max_elev,], family = poisson())
mSRP2<- glm(SpeciesRichness ~  Project, offset = log(Ndays), data = detections_for_SR[elevation < max_elev,], family = poisson())
mSRP3<- glm(SpeciesRichness ~  Project+poly(elevation,2), offset = log(Ndays), data = detections_for_SR[elevation < max_elev,], family = poisson())
mSRP4<- glm(SpeciesRichness ~  Project+poly(elevation,1), offset = log(Ndays), data = detections_for_SR[elevation < max_elev,], family = poisson())

mSR1 <- glm(SpeciesRichness ~ poly(elevation,1), offset = log(Ndays), data = detections_for_SR[elevation < max_elev,], family = poisson())
mSR2 <- glm(SpeciesRichness ~ 1, offset = log(Ndays), data = detections_for_SR[elevation < max_elev,], family = poisson())


MuMIn::QAICc(mSRP, chat = overdispersion)
AICcmodavg::c_hat(mSRP)
modcomp <- AICcmodavg::aictab(list('SR ~ elevation2' = mSR, 
                                   'SR ~ elevation' = mSR1,
                                   'SR ~ 1' = mSR2,
                                   'SR ~ elevation2*Project' = mSRP,
                                   'SR ~ elevation*Project' = mSRP1,
                                   'SR ~ Project' = mSRP2,
                                   'SR ~ elevation2+Project' = mSRP3,
                                   'SR ~ elevation+Project' = mSRP4
                                   
                                   ), second.ord = T, c.hat = overdispersion)

format(modcomp, digits = 2)

write.csv(format(modcomp, digits = 2), "results/SR_model_comparison.csv", row.names = F)

coef_m <- summary(mSRP3, dispersion = overdispersion)$coefficients


newdata = detections_for_SR[elevation < max_elev, .( Project, elevation, Ndays = 100)] # detections_for_SR[elevation < max_elev, .(Project, elevation)]
predictions <- predict(mSRP3, newdata = newdata, se.fit = T)
newdata$PoissonFit  <-  exp(predictions$fit)
newdata$SpeciesRichness <- newdata$PoissonFit # to be able to plot the ribbons, apparently
newdata$LCI  <-  exp(predictions$fit - 1.96 * predictions$se.fit* sqrt(overdispersion))
newdata$UCI  <-  exp(predictions$fit + 1.96 * predictions$se.fit* sqrt(overdispersion))
# nights$PoissonP  <- coef_m[match(paste0("elevation:Common.Name", nights$Common.Name, "Project", nights$Project), rownames(coef_m)),4]
newdata$pElev1<- coef_m["poly(elevation, 2)1",4] #newdata$pElev1<- coef_m[match(paste0("Project", newdata$Project, ":poly(elevation, 2)1"), rownames(coef_m)),4] #
newdata$pElev2<- coef_m["poly(elevation, 2)2",4 ] #newdata$pElev2<-coef_m[match(paste0("Project", newdata$Project, ":poly(elevation, 2)2"), rownames(coef_m)),4] #

newdata[, PoissonLinetype := ifelse(pElev2<0.05, "significant quadratic",
                                    ifelse(pElev1<0.05, "significant linear",
                                           " "))]

# p$data <- detections_for_SR[elevation < max_elev,]
g <- p + geom_ribbon(data = newdata[!PoissonLinetype %in% " ",], aes(ymin  = LCI, ymax = UCI, fill = Project), alpha = 0.3) + 
  geom_line(data = newdata, aes(y = PoissonFit, color = Project), linewidth = 1) + # , linetype = PoissonLinetype,
  scale_linetype_manual(values = c("significant quadratic" = "solid", "significant linear" = "dotted",  " " = "blank")) + 
  # labs(title = "Poisson (<max_elev m, no outliers)", linetype = "Elevation relationship")  +#   labs(title = "Poisson (<max_elev m, no outliers, locations where zero-inflated binomial <0.9)")
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.text.x.top = element_text(hjust = 0)) +
  labs(x = "elevation [m]")



# detections_for_SR <- detections_for_SR[elevation < max_elev,]
detections_for_SR$elevationBand <- cut(detections_for_SR$elevation, breaks = min(detections_for_SR$elevation) + 200 * c(0:(max(detections_for_SR$elevation)/200)), include.lowest = T, dig.lab = 5)


SampleEffortPerElevationBand <- detections_for_SR[, .(Ndays = sum(Ndays)), by = elevationBand ]
SampleEffortPerElevationBand[, xmin := as.numeric(regmatches(elevationBand, regexpr("\\d*(?=,)",elevationBand, perl = T)))] 
SampleEffortPerElevationBand[, xmax := as.numeric(regmatches(elevationBand, regexpr("(?<=,)\\d*",elevationBand, perl = T)))] 
SampleEffortPerElevationBand[, xmean := (xmin + xmax)/2] 

g +  geom_rect(data = SampleEffortPerElevationBand, inherit.aes = FALSE, mapping = aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = 1), color = "black", fill = "transparent") +
  annotate("text", label = SampleEffortPerElevationBand$Ndays, x = SampleEffortPerElevationBand$xmean, y = .5) +
  annotate("text", label = "N camera days", x = 400, y = -.5) 

ggsave("results/SpeciesRichness_against_elevation.tiff", width = 9, height = 7, units = "in", dpi = 600)

# SampleEffortPerElevationBand <-  data.frame(Ndays = tapply(detections_for_SR$Ndays, detections_for_SR$elevationBand, sum))

cowplot::plot_grid(g, ggplot(detections_for_SR[elevation < max_elev,], aes(elevationBand,y  = Ndays)) + geom_bar(stat  = "sum", show.legend = F)+
  labs(y = "N camera days", x = "elevation [m]"), nrow = 2) 


ggsave("results/SpeciesRichness_against_elevation_hidsEffort.png", width = 8, height = 7, units = "in", dpi = 300)


# look at endemism rate

pe <- ggplot(detections_for_SR, aes( x = elevation, y = EndemismRate)) +
  geom_point(aes(fill = Project, colour = Project, shape = Project), alpha = .5) + 
  labs(y = "Endemism rate") + 
  scale_shape_manual(values=c(21,22), labels = c("Sabah", "LEWS")) + 
  scale_colour_discrete(labels = c("Sabah", "LEWS")) + 
  scale_fill_discrete(labels = c("Sabah", "LEWS"))

max_elev = max(detections_for_SR$elevation) # 1500m elevation cuttoff or no cutoff max(elevation) 




# Compare with smaller models
mSRP <- glm(EndemismRate ~  Project*poly(elevation,2) - poly(elevation,2) - Project, data = detections_for_SR[elevation < max_elev,], family = binomial())
mSRP1 <- glm(EndemismRate ~  Project*poly(elevation,1) - poly(elevation,1) - Project, data = detections_for_SR[elevation < max_elev,], family = binomial())
mSRP2<- glm(EndemismRate ~  Project, data = detections_for_SR[elevation < max_elev,], family = binomial())
mSRP3<- glm(EndemismRate ~  Project+poly(elevation,2), data = detections_for_SR[elevation < max_elev,], family = binomial())
mSRP4<- glm(EndemismRate ~  Project+poly(elevation,1), data = detections_for_SR[elevation < max_elev,], family = binomial())

mSR <- glm(EndemismRate ~ poly(elevation,2), data = detections_for_SR[elevation < max_elev,], family = binomial())
mSR1 <- glm(EndemismRate ~ poly(elevation,1), data = detections_for_SR[elevation < max_elev,], family = binomial())
mSR2 <- glm(EndemismRate ~ 1, data = detections_for_SR[elevation < max_elev,], family = binomial())


modcomp <- AICcmodavg::aictab(list('endemism rate ~ elevation2' = mSR, 
                                   'endemism rate ~ elevation' = mSR1,
                                   'endemism rate ~ 1' = mSR2,
                                   'endemism rate ~ elevation2*Project' = mSRP,
                                   'endemism rate ~ elevation*Project' = mSRP1,
                                   'endemism rate ~ Project' = mSRP2,
                                   'endemism rate ~ elevation2+Project' = mSRP3,
                                   'endemism rate ~ elevation+Project' = mSRP4
                                   
), second.ord = T)

format(modcomp, digits = 2)

write.csv(format(modcomp, digits = 2), "results/EndemismRate_model_comparison.csv", row.names = F)




# Summary plots and tables for known species ####
detections <- detections[!detections$Common.Name %in% c(unk_species_to_remove),]

table(detections$Project); nrow(detections)

## plot
stdzed_det <- detections[, .(mean_det_rate = mean(.N/as.numeric(Ndays))), by= list(Common.Name, SpGroup, elevation, Project)] # mean detection rate at each elevation
stdzed_det[, det_rate_per_100_night := round(mean_det_rate*100)] # round n detection you would have for 100 days

stdzed_det <- stdzed_det[, .(elevation = rep(elevation, det_rate_per_100_night)), by = list(Common.Name, Project, SpGroup)] # repeat elevation value as many times as n detections per 100 days so we can use it in a violin plot


stdzed_det <- droplevels(stdzed_det[!Common.Name %in% min10dectect, ]) # only keep those that are dected in a lieast 10 locations for this figure

stdzed_det$Common.Name <- reorder(stdzed_det$Common.Name, stdzed_det$elevation, median)


p <- ggplot(stdzed_det, aes(x = elevation, y = Common.Name, fill = SpGroup)) + 
  geom_boxplot()  + 
  geom_rug(color = "black", sides = "b") + coord_cartesian(ylim = c(-1, nlevels(stdzed_det$Common.Name)), clip = "off") + 
  labs(fill = "Guild",
       x = "elevation [m]",
       y = "")

p + labs(title = "Both sites")

ggsave("results/Boxplot_elevation_per_species.png", width = 10, height = 7, units = "in", dpi = 600)


p + facet_grid(~Project, labeller = labeller(Project = c(
  "Sabah Parks Elevational Mammal Survey Project" = "Sabah",
  "Smithsonian Borneo Mammal Survey at LEWS Project" = "LEWS")))

ggsave("results/Boxplot_elevation_per_species_and_site.png", width = 10, height = 7, units = "in", dpi = 600)

## table
summaryTable <- detections[, .(deployments = length(unique(Deployment.ID)),
                               detections = .N,
                               elevation_mean = round(mean(elevation), 2),
                               elevation_min = min(elevation),
                               elevation_max = max(elevation)), by= list(Species.Name, Common.Name, SpGroup, Project)] 
summaryTable[, Project := ifelse(grepl("LEWS", Project), "L", "S")]
summaryTable[, elevation_range := paste(elevation_min, "-", elevation_max)]
summaryTable[, c("elevation_min", "elevation_max"):=NULL]

summaryTable <- reshape(summaryTable, direction = "wide", idvar = c("Species.Name", "Common.Name", "SpGroup"), timevar = "Project")


setDF(summaryTable)
names(summaryTable) <- gsub("SpGroup", "Guild", names(summaryTable))
summaryTable[, c("deployments.S", "deployments.L", "detections.S", "detections.L")][is.na(summaryTable[, c("deployments.S", "deployments.L", "detections.S", "detections.L")])] <- 0

setDT(summaryTable)

summaryTable <- summaryTable[, .(Species.Name, Common.Name, Guild, deployments.S, deployments.L, detections.S, detections.L, elevation_mean.S, elevation_mean.L, elevation_range.S, elevation_range.L)]


write.csv(summaryTable, "results/Summary_table.csv", row.names = F)

# Analyse only for species that are detected at least 10 locations at each sites and 100 detections minimum####

# calculate camera nights per species and per elevation
nights <- unique(detections[ !Common.Name %in% c(min10dectectEachSite, unk_species_to_remove), .(.N), by = .(Project, SpGroup, BodyMass, Common.Name, Deployment.ID, detection_distance, daycut, elevation, Ndays)])
nights[, Ndays := as.numeric(Ndays)]

# keep only when at least 100 detections
nights <- droplevels(nights[, if(sum(N)>=100) .SD, by= Common.Name])

# detect outliers (before we fill all the zeroes)
nights[, outlier := ifelse(N/Ndays > mean(N/Ndays)+(3*sd(N/Ndays)), TRUE, FALSE), by = Common.Name]


# add zeroes to species not seen
nights <- nights[, .SD[unique(nights[, .(Common.Name, SpGroup, BodyMass)]), on=c("Common.Name", "SpGroup", "BodyMass")], by = .(Project, elevation, Deployment.ID, daycut, detection_distance, Ndays)]

unique(table(nights$Common.Name)) # should be 783

nights[is.na(N), N := 0]
nights[is.na(outlier), outlier := FALSE]

# remove outliers
nights <- nights[outlier %in% FALSE, ]
nights[, outlier := NULL]

# plot

p <- ggplot(nights, aes( x = elevation, y = N/Ndays)) +
  geom_point(aes(fill = Project, colour = Project, shape = Project), alpha = .5) + 
  # geom_point(data = nights[outlier %in% TRUE,], aes(colour = Project, shape = Project), alpha = .5) +
  facet_wrap(~Common.Name, scales = "free_y",
             labeller = labeller(Common.Name = label_wrap_gen(width = 20))) +
  labs(y = "N detections / day") + 
  scale_shape_manual(values=c(21,22), labels = c("Sabah", "LEWS")) + 
  scale_colour_discrete(labels = c("Sabah", "LEWS")) + 
  scale_fill_discrete(labels = c("Sabah", "LEWS"))

## now do poisson regression with offset only using data with "true zeroes"
max_elev =  1500# 1500m elevation cuttoff or no cutoff max(elevation) 
mP <- glm(N ~  detection_distance + Project*poly(elevation,2)*Common.Name - poly(elevation,2)*Project - poly(elevation,2)*Common.Name - Project - poly(elevation,2) - Common.Name, offset = as.numeric(Ndays), data = nights[elevation < max_elev,], family = poisson())

m <- glm(N ~  detection_distance + poly(elevation,2)*Common.Name-poly(elevation,2) , offset = as.numeric(Ndays), data = nights[elevation < max_elev,], family = poisson())
# m <- glm(N ~  detection_distance + poly(elevation,2)*Common.Name-poly(elevation,2) + offset(as.numeric(Ndays)), data = nights[pStrZero < 0.9 & elevation < max_elev], family = poisson())

(overdispersion = mP$deviance/mP$df.residual) ; AER::dispersiontest(mP)# a little bit over 1... not too bad but let's use it when we calculated uncertainty or compare models



# Compare with smaller models
m1P <- glm(N ~  poly(elevation,2)*Project*Common.Name - poly(elevation,2)*Project - poly(elevation,2)*Common.Name - Project - poly(elevation,2) - Common.Name, offset = as.numeric(Ndays), data = nights[elevation < max_elev,], family = poisson())
m2P<- glm(N ~  detection_distance + poly(elevation,1)*Project*Common.Name - poly(elevation,1)*Project - poly(elevation,1)*Common.Name - poly(elevation,1)  - Common.Name, offset = as.numeric(Ndays), data = nights[elevation < max_elev,], family = poisson())
m3P<- glm(N ~  detection_distance + poly(elevation,2)*Project, offset = as.numeric(Ndays), data = nights[elevation < max_elev,], family = poisson())
m4P<- glm(N ~  detection_distance + Project, offset = as.numeric(Ndays), data = nights[elevation < max_elev,], family = poisson())
m5P <- glm(N ~ poly(elevation,2)*Project, offset = as.numeric(Ndays), data = nights[elevation < max_elev,], family = poisson())
m6P <- glm(N ~ elevation*Project, offset = as.numeric(Ndays), data = nights[elevation < max_elev,], family = poisson())
m7P<- glm(N ~  Project, offset = as.numeric(Ndays), data = nights[elevation < max_elev,], family = poisson())

m1 <- glm(N ~  poly(elevation,2)*Common.Name-poly(elevation,2), offset = as.numeric(Ndays), data = nights[elevation < max_elev,], family = poisson())
m2<- glm(N ~  detection_distance + poly(elevation,1)*Common.Name-poly(elevation,1), offset = as.numeric(Ndays), data = nights[elevation < max_elev,], family = poisson())
m3<- glm(N ~  detection_distance + poly(elevation,2), offset = as.numeric(Ndays), data = nights[elevation < max_elev,], family = poisson())
m4<- glm(N ~  detection_distance, offset = as.numeric(Ndays), data = nights[elevation < max_elev,], family = poisson())
m5 <- glm(N ~ poly(elevation,2), offset = as.numeric(Ndays), data = nights[elevation < max_elev,], family = poisson())
m6 <- glm(N ~ elevation, offset = as.numeric(Ndays), data = nights[elevation < max_elev,], family = poisson())
m7<- glm(N ~  1, offset = as.numeric(Ndays), data = nights[elevation < max_elev,], family = poisson())



MuMIn::QAICc(mP, chat = overdispersion)
AICcmodavg::c_hat(mP)
modcomp <- AICcmodavg::aictab(list('N ~ detectiondist + elevation2*species' = m, 
                        'N ~ elevation2*species' = m1,
                        'N ~ detectiondist + elevation*species' = m2,
                        'N ~ detectiondist + elevation2' = m3,
                        'N ~ detectiondist' = m4,
                        'N ~ elevation2' = m5,
                        'N ~ elevation' = m6,
                        'N ~ 1' =m7,
                        'N ~ detectiondist + elevation2*species*Project' = mP, 
                        'N ~ elevation2*species*Project' = m1P,
                        'N ~ detectiondist + elevation*species*Project' = m2P,
                        'N ~ detectiondist + elevation2*Project' = m3P,
                        'N ~ detectiondist + Project' = m4P,
                        'N ~ elevation2*Project' = m5P,
                        'N ~ elevation*Project' = m6P,
                        'N ~ Project' =m7P), second.ord = T, c.hat = overdispersion)

format(modcomp, digits = 2)

write.csv(format(modcomp, digits = 2), "results/model_comparison.csv", row.names = F)

coef_m <- summary(mP, dispersion = overdispersion)$coefficients


newdata = nights[N>0 & elevation < max_elev, .(Project, elevation, Common.Name, Ndays = mean(Ndays), detection_distance = 1)]
predictions <- predict(mP, newdata = newdata, se.fit = T)
newdata$PoissonFit  <-  exp(predictions$fit)
newdata$N <- newdata$PoissonFit # to be able to plot the ribbons, apparently
newdata$LCI  <-  exp(predictions$fit - 1.96 * predictions$se.fit* sqrt(overdispersion))
newdata$UCI  <-  exp(predictions$fit + 1.96 * predictions$se.fit* sqrt(overdispersion))
# nights$PoissonP  <- coef_m[match(paste0("elevation:Common.Name", nights$Common.Name, "Project", nights$Project), rownames(coef_m)),4]
newdata$pElev1<- coef_m[match(paste0("Project", newdata$Project, ":poly(elevation, 2)1:Common.Name", newdata$Common.Name), rownames(coef_m)),4]
newdata$pElev2<- coef_m[match(paste0("Project", newdata$Project, ":poly(elevation, 2)2:Common.Name", newdata$Common.Name), rownames(coef_m)),4]

newdata[, PoissonLinetype := ifelse(pElev2<0.05, "significant quadratic",
                                   ifelse(pElev1<0.05, "significant linear",
                                          " "))]

p$data <- nights[elevation < max_elev,]
p + geom_ribbon(data = newdata[!PoissonLinetype %in% " ",], aes(ymin  = LCI, ymax = UCI, fill = Project), alpha = 0.3) + 
  geom_line(data = newdata, aes(y = PoissonFit, linetype = PoissonLinetype, color = Project), linewidth = 1) +
  scale_linetype_manual(values = c("significant quadratic" = "solid", "significant linear" = "dotted",  " " = "blank")) + 
  # labs(title = "Poisson (<max_elev m, no outliers)", linetype = "Elevation relationship")  +#   labs(title = "Poisson (<max_elev m, no outliers, locations where zero-inflated binomial <0.9)")
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.text.x.top = element_text(hjust = 0))


ggsave("results/Model_predictions.png", width = 10, height = 7, units = "in", dpi = 300)



###################################################


### multivariate analysis with mixed quantitative variables and factors.

All_ddToPlot <- NULL

for(what in c("Both sites", "Sabah", "LEWS")) {
  
  if(what %in% "Both sites") x <- stdzed_det[, .(elevation = median(elevation)), by = list(Common.Name, SpGroup)] else x <- stdzed_det[grepl(what, Project),.(elevation = median(elevation)), by = list(Common.Name, SpGroup)]
 
  # y <-  detections[, .SD[1, .(Herbivore_score, Nocturnal,CrepuscularDiurnal, BodyMass)], by = Common.Name] # with activity pattern
  
  y <-  detections[, .SD[1, .(Herbivore_score, BodyMass)], by = list(Common.Name)]
  
  
  x <- merge(x, y)
  
  x[, BodyMass := log10(BodyMass)]
  
   x <- x[complete.cases(x), ]
  
  dd <- dudi.mix(x[,setdiff(names(x), c("Common.Name", "SpGroup")), with = F], scannf = F)
  
  if(what %in% "Both sites") dd_both <- dd
  
  
  ddToPlot <- as.data.table(mapply(c, dd$co[, c(1,2)]*3,dd$li[, c(1,2)]))
  
  ddToPlot[, label := c(rownames(dd$co),  as.character(x$Common.Name))]
  
  # ddToPlot[, Project := c(rep("", nrow(dd$co)), x$Project)]
  
  # ddToPlot[, color := c(rep("#000000", nrow(dd$co)), spColors[x$Common.Name])]
  ddToPlot$guild <- detections$SpGroup[match(ddToPlot$label, detections$Common.Name)]
  
  ddToPlot[, size := c(rep(3.5, nrow(dd$co)), rep(3, nrow(dd$li)))]
  
  ddToPlot[, c("xend", "yend") := c(rep(0,  nrow(dd$co)), rep(NA, nrow(dd$li)))]
  
  ddToPlot$what = what
  
  All_ddToPlot <- rbind(All_ddToPlot,ddToPlot )
}

summary(dd_both)

Projected_inertia <- round(dd_both$eig*100/sum(dd_both$eig))

p <- ggplot(data = All_ddToPlot[what %in% "Both sites", ], aes(x = Comp1, y = Comp2, label =label, color = guild)) +
  # geom_vline(aes(xintercept = 0)) +
  # geom_hline(aes(yintercept = 0)) +
  geom_segment(data = All_ddToPlot[what %in% "Both sites" & is.na(guild), ], aes(xend=xend, yend=yend)) +
  geom_point() +
  ggrepel::geom_text_repel(aes(size = size), nudge_x = 0, nudge_y = .1) +
  # scale_color_identity() +
  scale_size_identity() +
  labs(x = paste0("Axis 1 (", Projected_inertia[1], "%)"), y = paste0("Axis 2 (", Projected_inertia[2], "%)")) +
  theme(aspect.ratio = 1) +
  xlim(-3, 3) + ylim(-3,3)


p

ggsave("results/Multivariate_results.tiff", width = 8, height = 7, units = "in", dpi = 600)




write.csv(round(dd_both$c1,2), "results/Multivariate_column_normed_scores.csv")

p$data <- All_ddToPlot[what %in% c("Sabah", "LEWS"), ]
p + facet_grid(~what)


dev.off()


###### LOOK AT RATS ####
rats <- detections[ Common.Name %in% rat_species, ]
rats <- unique(rats[ , .(.N), by = .(Project, SpGroup,BodyMass,Deployment.ID, detection_distance, daycut, elevation, Ndays)])
rats[, Ndays := as.numeric(Ndays)]


# add zeroes
rats <- rats[, .SD[unique(nights[, .(SpGroup, BodyMass)]), on=c("SpGroup", "BodyMass")], by = .(Project, elevation, Deployment.ID, daycut, detection_distance, Ndays)]


rats[is.na(N), N := 0]

# keep only when at least 100 detections
rats <- droplevels(rats[, if(sum(N)>=100) .SD])

# detect outliers (before we fill all the zeroes)
rats[, outlier := ifelse(N/Ndays > mean(N/Ndays)+(3*sd(N/Ndays)), TRUE, FALSE)]

# remove outliers
rats <- rats[outlier %in% FALSE, ]
rats[, outlier := NULL]

# plot

p <- ggplot(rats, aes( x = elevation, y = N/Ndays)) +
  geom_point(aes(fill = Project, colour = Project, shape = Project), alpha = .5) + 
  labs(y = "N detections / day") + 
  scale_shape_manual(values=c(21,22), labels = c("Sabah", "LEWS")) + 
  scale_colour_discrete(labels = c("Sabah", "LEWS")) + 
  scale_fill_discrete(labels = c("Sabah", "LEWS"))

p

## now do poisson regression with offset only using data with "true zeroes"

# Compare models
m1P <- glm(N ~  detection_distance + poly(elevation,2)*Project - poly(elevation,2) - Project - poly(elevation,2), offset = as.numeric(Ndays), data = rats, family = poisson())
m2P<- glm(N ~  poly(elevation,2)*Project - poly(elevation,2) - Project , offset = as.numeric(Ndays), data = rats, family = poisson())
m3P<- glm(N ~  detection_distance + poly(elevation,1)*Project - poly(elevation,1) - Project, offset = as.numeric(Ndays), data = rats, family = poisson())
m4P<- glm(N ~  poly(elevation,1)*Project - poly(elevation,1) - Project , offset = as.numeric(Ndays), data = rats, family = poisson())
m5P <- glm(N ~  detection_distance + poly(elevation,2) + Project, offset = as.numeric(Ndays), data = rats, family = poisson())
m6P <-  glm(N ~  detection_distance + poly(elevation,1) + Project, offset = as.numeric(Ndays), data = rats, family = poisson())
m7P<- glm(N ~  poly(elevation,2) + Project , offset = as.numeric(Ndays), data = rats, family = poisson())
m8P<- glm(N ~  poly(elevation,1) + Project, offset = as.numeric(Ndays), data = rats, family = poisson())
m9P<- glm(N ~  Project , offset = as.numeric(Ndays), data = rats, family = poisson())

m1 <- glm(N ~  detection_distance + poly(elevation,2), offset = as.numeric(Ndays), data = rats, family = poisson())
m2 <- glm(N ~  poly(elevation,2), offset = as.numeric(Ndays), data = rats, family = poisson())
m3 <- glm(N ~  detection_distance + poly(elevation,1), offset = as.numeric(Ndays), data = rats, family = poisson())
m4 <- glm(N ~  poly(elevation,1), offset = as.numeric(Ndays), data = rats, family = poisson())
m5 <- glm(N ~  detection_distance, offset = as.numeric(Ndays), data = rats, family = poisson())
m6 <-  glm(N ~ 1, offset = as.numeric(Ndays), data = rats, family = poisson())

overdispersion = AICcmodavg::c_hat(mP)
MuMIn::QAICc(mP, chat = overdispersion)
AICcmodavg::c_hat(mP)
modcomp <- AICcmodavg::aictab(list('N ~ detectiondist + elevation2*Project' = m1P,
                                   'N ~ detectiondist + elevation2' = m1,
                                   'N ~ elevation2*Project' = m2P,
                                   'N ~ elevation2' = m2,
                                   'N ~ detectiondist + elevation*Project' = m3P,
                                   'N ~ detectiondist + elevation' = m3,
                                   'N ~ elevation*Project' = m4P,
                                   'N ~ elevation' =m4,
                                   'N ~ ddetectiondist + elevation2+Project' = m5P,
                                   'N ~ ddetectiondist' = m5,
                                   'N ~ detectiondist + elevation + Project' = m6P,
                                   'N ~ 1' = m6,
                                   'N ~ elevation2+Project' = m7P,
                                   'N ~ elevation+Project' = m8P,
                                   'N ~ Project' = m9P), second.ord = T)

format(modcomp, digits = 2)

# # write.csv(format(modcomp, digits = 2), "results/model_comparison.csv", row.names = F)
# 
coef_m <- summary(m1P, dispersion = overdispersion)$coefficients


newdata = rats[N>0, .(Project, elevation, Ndays = mean(Ndays), detection_distance = 1)]
predictions <- predict(m1P, newdata = newdata, se.fit = T)
newdata$PoissonFit  <-  exp(predictions$fit)
newdata$N <- newdata$PoissonFit # to be able to plot the ribbons, apparently
newdata$LCI  <-  exp(predictions$fit - 1.96 * predictions$se.fit* sqrt(overdispersion))
newdata$UCI  <-  exp(predictions$fit + 1.96 * predictions$se.fit* sqrt(overdispersion))
# # nights$PoissonP  <- coef_m[match(paste0("elevation:Common.Name", nights$Common.Name, "Project", nights$Project), rownames(coef_m)),4]
newdata$pElev1<- coef_m[match(paste0("poly(elevation, 2)1:Project", newdata$Project), rownames(coef_m)),4]
newdata$pElev2<- coef_m[match(paste0("poly(elevation, 2)2:Project", newdata$Project), rownames(coef_m)),4]

newdata[, PoissonLinetype := ifelse(pElev2<0.05, "significant quadratic",
                                    ifelse(pElev1<0.05, "significant linear",
                                           " "))]
# 
p$data <- rats
p + geom_ribbon(data = newdata[!PoissonLinetype %in% " ",], aes(ymin  = LCI, ymax = UCI, fill = Project), alpha = 0.3) +
  geom_line(data = newdata, aes(y = PoissonFit, linetype = PoissonLinetype, color = Project), linewidth = 1) +
  scale_linetype_manual(values = c("significant quadratic" = "solid", "significant linear" = "dotted",  " " = "blank")) +
  # labs(title = "Poisson (<max_elev m, no outliers)", linetype = "Elevation relationship")  +#   labs(title = "Poisson (<max_elev m, no outliers, locations where zero-inflated binomial <0.9)")
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.text.x.top = element_text(hjust = 0))


# # ggsave("results/Model_predictions.png", width = 10, height = 7, units = "in", dpi = 300)
# 
