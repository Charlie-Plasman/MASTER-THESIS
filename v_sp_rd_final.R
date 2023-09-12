### VIRTUAL SPECIES GRTS

library(sp)
library(raster)
library(SDMPlay)
library(spsurvey)
library(SOmap)
library(caret)
library(dplyr)

n_repeats = 10
sample_size <- 50
full_n_pres_list = list()
full_kappa_v_list = list()
full_predictions_list = list()
full_sampling_sites_list = list()
basemap = SOmap(bathy_legend = F, ice = T, trim = -45, fronts = "Park", border = F) 
wap = extent(-105, -50, -90, -60)
stack_linear_cropped = crop(stack_linear, wap)

# Étape 1 : Génération de l'espèce virtuelle
setwd("C:/Users/charl/Desktop")
virtual_species = raster("MYCTO_PRED.tif")
# Calcul des quantiles pour définir les seuils
threshold <- quantile(virtual_species, probs = 0.99)
# Classification binaire en utilisant les seuils
binary_virtual_species <- reclassify(virtual_species, c(0,threshold, 0))
binary_virtual_species <- reclassify(binary_virtual_species, c(threshold,1, 1))
binary_virtual_species_cropped = crop(binary_virtual_species, wap)
# Convertir les valeurs logiques en valeurs numériques (0 pour FALSE, 1 pour TRUE)
numeric_values <- as.numeric(values(binary_virtual_species) == 1)
# Compter le nombre total de pixels
total_pixels <- ncell(binary_virtual_species) 
# Compter le nombre de pixels ayant une valeur de 1
pixels_with_value_1 <- sum(!is.na(values(binary_virtual_species)) & values(binary_virtual_species) == 1, na.rm = TRUE)
# Calculer le pourcentage
percentage <- (pixels_with_value_1 / total_pixels) * 100
# Afficher le pourcentage
print(paste("Pourcentage de pixels avec valeur 1:", round(percentage, 2), "%"))

# Add the new raster to the existing stack
envi_data = stack_linear_cropped
envi_data <- addLayer(envi_data, binary_virtual_species_cropped)

for(i in 1:n_repeats){
  # Étape 2 : Récolte données initiales
  true_distrib_points <- rasterToPoints(binary_virtual_species_cropped)
  true_distrib_points <- data.frame(true_distrib_points)
  # Visualiser la distribution des valeurs de profondeur
  hist(true_distrib_points$MYCTO_PRED)
  initial_data = 25
  n_strata <- c("0" = initial_data,"1" = initial_data)
  # Convert the data frame 'distrib_points' to an 'sf' object
  true_distrib_points_sf <- st_as_sf(true_distrib_points, coords = c("x", "y"))
  print(st_crs(true_distrib_points_sf))
  # Définir le CRS géographique approprié (WGS84)
  geographic_crs <- st_crs("+proj=longlat +datum=WGS84 +no_defs")
  # Assigne le CRS géographique à ton objet
  st_crs(true_distrib_points_sf) <- geographic_crs
  # Vérifier que le CRS a été correctement défini
  print(st_crs(true_distrib_points_sf))
  # Définir le CRS projeté approprié (exemple : UTM)
  target_crs <- "+proj=utm +zone=30 +datum=WGS84 +units=m +no_defs"
  # Convertir les coordonnées géographiques en coordonnées projetées
  true_distrib_points_sf <- st_transform(true_distrib_points_sf, crs = target_crs)
  # Vérifier que le CRS a été correctement défini
  st_crs(true_distrib_points_sf)
  
  eqprob_strat <- grts(true_distrib_points_sf, n_base = n_strata, stratum_var = "MYCTO_PRED")
  # Afficher les résultats
  print(eqprob_strat)
  plot(eqprob_strat, siteuse = "Base")
  ###
  coord = data.frame(latitude = eqprob_strat$sites_base$lat_WGS84, longitude = eqprob_strat$sites_base$lon_WGS84)
  coord_sp = SpatialPointsDataFrame(coords = coord[,c("longitude","latitude")],data =coord, proj4string =  CRS("+proj=longlat +datum=WGS84"))
  plot(binary_virtual_species_cropped)
  plot(coord_sp, add = T)
  ###
  true_sites = data.frame(latitude = eqprob_strat$sites_base$lat_WGS84, 
                          longitude = eqprob_strat$sites_base$lon_WGS84)
  true_sites = cbind.data.frame(true_sites$longitude,true_sites$latitude)
  colnames(true_sites) = c("longitude", "latitude")
  
  # Étape 5 : Ajout des nouvelles données
  realdistrib = stack(stack_linear_cropped, binary_virtual_species_cropped)
  true_sites_df = SDMPlay:::SDMtab(xydata=true_sites, predictors=realdistrib, unique.data=F, same=F, background.nb = 0)
  count <- sum(true_sites_df$MYCTO_PRED == 1)
  # Supprimer la première colonne
  true_sites_df <- true_sites_df[, -1]
  # Placer la dernière colonne en première position
  last_column <- ncol(true_sites_df)
  true_sites_df <- true_sites_df[, c(last_column, 1:(last_column-1))]
  # Changer le nom de la première colonne (qui était la dernière)
  new_first_column_name <- "id"
  colnames(true_sites_df)[1] <- new_first_column_name
  # Calculer la taille du data frame
  #taille_df <- nrow(true_sites_df)
  # Calculer l'indice à partir duquel vous voulez changer les valeurs en 0 (deuxième moitié)
  #indice_debut_modification <- ceiling(taille_df / 2)
  # Utiliser dplyr pour modifier les valeurs
  #true_sites_df <- true_sites_df %>%
  #mutate(id = ifelse(row_number() > indice_debut_modification, 0, id))
  training_data = true_sites_df
  # Étape 2 : Modèle de base
  # Entraînement d'un modèle de base (p.ex. Random Forest)
  initial_model <- SDMPlay:::compute.brt(x = training_data, proj.predictors = stack_linear_cropped,
                                         tc = 2, lr = 0.0001, bf = 0.75, n.trees = 50, 
                                         step.size = 50, n.folds = 5)
  initial_prediction = initial_model$raster.prediction
  # Calcul des quantiles pour définir les seuils
  threshold <- quantile(initial_prediction, probs = 0.90)
  # Classification binaire en utilisant les seuils
  binary_initial_prediction <- reclassify(initial_prediction, c(0,threshold, 0))
  binary_initial_prediction <- reclassify(binary_initial_prediction, c(threshold,1, 1))
  binary_initial_prediction_cropped = crop(binary_initial_prediction, wap)
  plot(binary_initial_prediction_cropped)
  #plot(basemap)
  #SOplot(first_prediction, col = hcl.colors(80, "viridis"), legend = F)
  
  # Itérations
  iterations = 10
  new_presences_list <- list()
  kappa_values_list <- list()
  predictions_list = list()
  sampling_sites_list = list()
  predictions_list[[1]] = binary_initial_prediction_cropped
  sampling_sites_list[[1]] = training_data
  for(j in 1:iterations){
    # Étape 4 : Simulation de l'échantillonnage
    # Charger ta zone d'étude 
    study_area <- crop(envi_data$gb_depth, wap)
    # Générer des points aléatoires dans la zone d'étude
    random_points <- sampleRandom(study_area, size = sample_size, xy = TRUE)
    random_points=data.frame(random_points[,-3])
    # Afficher les premiers points générés
    head(random_points)
    sampling_sites = cbind.data.frame(random_points$x,random_points$y)
    colnames(sampling_sites) = c("longitude", "latitude")
    # Étape 5 : Ajout des nouvelles données
    sampled_sites_df = SDMPlay:::SDMtab(xydata=sampling_sites, predictors=envi_data, unique.data=F, same=F, background.nb = 0)
    # Remplacer les NA par zéros dans une colonne spécifique
    sampled_sites_df$MYCTO_PRED[is.na(sampled_sites_df$MYCTO_PRED)] <- 0
    count <- sum(sampled_sites_df$MYCTO_PRED == 1)
    new_presences_list[[j]] <- count
    # Supprimer la première colonne
    sampled_sites_df <- sampled_sites_df[, -1]
    # Placer la dernière colonne en première position
    last_column <- ncol(sampled_sites_df)
    sampled_sites_df <- sampled_sites_df[, c(last_column, 1:(last_column-1))]
    # Changer le nom de la première colonne (qui était la dernière)
    new_first_column_name <- "id"
    colnames(sampled_sites_df)[1] <- new_first_column_name
    sampling_sites_list = append(sampling_sites_list,list(sampled_sites_df))
    new_training_data = rbind(training_data,sampled_sites_df)
    
    # Étape 6 : Construction du modèle itératif
    iterative_model <- SDMPlay:::compute.brt(x = new_training_data, proj.predictors = stack_linear_cropped,
                                             tc = 2, lr = 0.0001, bf = 0.75, n.trees = 50, 
                                             step.size = 50, n.folds = 5)
    iterative_prediction = iterative_model$raster.prediction
    plot(iterative_prediction)
    # Calcul des quantiles pour définir les seuils
    threshold <- quantile(iterative_prediction, probs = 0.90)
    # Classification binaire en utilisant les seuils
    binary_iterative_prediction <- reclassify(iterative_prediction, c(0,threshold, 0))
    binary_iterative_prediction <- reclassify(binary_iterative_prediction, c(threshold,1, 1))
    plot(binary_iterative_prediction)
    #plot(basemap)
    #SOplot(binary_iterative_prediction, col = hcl.colors(80, "viridis"), legend = F)
    predictions_list = append(predictions_list,binary_iterative_prediction)
    
    # Réajuster les données d'entraînement pour la prochaine itération
    binary_initial_prediction_cropped = binary_iterative_prediction
    training_data = new_training_data
    
    # Étape 7 : Calcul des métriques
    # Load your binary raster predictions and known distribution raster
    prediction_raster <- binary_iterative_prediction
    known_distribution_raster <- binary_virtual_species_cropped
    
    prediction_matrix <- as.matrix(prediction_raster)
    known_distribution_matrix <- as.matrix(known_distribution_raster)
    prediction_vector <- as.vector(prediction_matrix)
    known_distribution_vector <- as.vector(known_distribution_matrix)
    
    # Convert vectors to factors with levels "0" and "1"
    prediction_factor <- factor(prediction_vector, levels = c("0", "1"))
    known_distribution_factor <- factor(known_distribution_vector, levels = c("0", "1"))
    confusion_matrix <- confusionMatrix(prediction_factor, known_distribution_factor)
    
    # Calculate Cohen's Kappa
    kappa_value <- confusion_matrix$overall['Kappa']
    # Display the Cohen's Kappa coefficient
    print(paste("Cohen's Kappa Coefficient:", kappa_value))
    kappa_values_list[[j]] = kappa_value
    
  }
  full_n_pres_list[[i]] = new_presences_list 
  full_kappa_v_list[[i]] = kappa_values_list
  full_predictions_list[[i]] = predictions_list    
  full_sampling_sites_list[[i]] = sampling_sites_list
}

first_kappas_list = list()
for(i in seq_along(full_predictions_list)){
  raster = full_predictions_list[[i]][[1]]
  prediction_raster <- raster
  known_distribution_raster <- binary_virtual_species_cropped
  
  
  prediction_matrix <- as.matrix(prediction_raster)    
  known_distribution_matrix <- as.matrix(known_distribution_raster)
  prediction_vector <- as.vector(prediction_matrix)
  known_distribution_vector <- as.vector(known_distribution_matrix)
  
  # Convert vectors to factors with levels "0" and "1"
  prediction_factor <- factor(prediction_vector, levels = c("0", "1"))
  known_distribution_factor <- factor(known_distribution_vector, levels = c("0", "1"))
  confusion_matrix <- confusionMatrix(prediction_factor, known_distribution_factor)
  
  # Calculate Cohen's Kappa
  kappa_value <- confusion_matrix$overall['Kappa']
  # Display the Cohen's Kappa coefficient
  print(paste("Cohen's Kappa Coefficient:", kappa_value))
  first_kappas_list[[i]] = kappa_value
}

full_kappa_v_list_save = full_kappa_v_list
for(i in seq_along(full_kappa_v_list)){
  full_kappa_v_list[[i]] <- c(first_kappas_list[[i]][[1]], full_kappa_v_list[[i]])
}

# Transposer les données
transposed_data <- sapply(full_kappa_v_list, unlist)
# Convertir en dataframe
kappa_dataframe <- as.data.frame(transposed_data)
# Transposer le dataframe pour avoir chaque réplicat en colonne
transposed_dataframe <- as.data.frame(t(kappa_dataframe))
# Changer les noms des colonnes en "1", "2", "3", ..., "11"
colnames(transposed_dataframe) <- 1:11
# Créer un boxplot à partir du dataframe transposé
boxplot(transposed_dataframe, main = "Random - 25 échantillons",
        xlab = "Nombre d'itérations", ylab = "Accord du modèle avec la 'vraie distribution' virtuelle")
transposed_pres <- t(sapply(full_n_pres_list, unlist))
boxplot(transposed_pres, main = "Random - 25 échantillons",
        xlab = "Nombre d'itérations", ylab = "Nombres de nouvelles présences trouvées")


#save.image("C:/Users/charl/Desktop/FINAL_RANDOM.RData")

plot(binary_virtual_species_cropped)
plot(binary_initial_prediction)
plot(binary_iterative_prediction)
plot(iterative_prediction)
