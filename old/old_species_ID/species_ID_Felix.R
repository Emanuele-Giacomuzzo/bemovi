#Define the first species (to use for training data for species identification)
data_traits$species <- ifelse(data_traits$Species_1 %in% c("Pau_0g", "Pau_4g"), "Pau", "Col") #Tell monocultures 

#Perform the species identification
#Get the training data (80 % of the monoculture data), verification data (rest of
#monoculture data) and prediction data (communities)
#Perform separately for every environment, to account for plastic effects
{
  #Perform identification for populations without stress
  #Get the training data
  {
    #Filter data to keep only monoculture data for training
    data.filt <- filter(data_traits, Culture=="monoculture") #Select monoculture data
    #Set a seed
    set.seed(294) 
    
    #Create a random order of the rows, and get 80 and 20 % of these numbers
    rows <- sample(nrow(data.filt)) #Random numbers to the rows
    training <- rows[1:round(length(rows)*0.8)] #Take first 80% OF ROWS
    verification <- rows[1+round(length(rows)*0.8):length(rows)] #20% FOR VERIFICATION 
    
    #Select the training data
    training_data <- filter(data.filt[training, ], Culture=="monoculture") %>% #SELECT THE VARIABLES TO TRAIN THE MODEL 
      dplyr::select(species, mean_area, mean_major, mean_minor, mean_turning, net_speed, gross_speed, sd_turning, sd_gross_speed, Stress_treatment)
  }
  
  #Get the verification data
  verification_data <- filter(data.filt[verification, ], Culture=="monoculture") %>% #SAME BUT FOR THE VERIFICATION DATA 
    dplyr::select(species, mean_area, mean_major, mean_minor, mean_turning, net_speed, gross_speed, sd_turning, sd_gross_speed, species, Stress_treatment)
  
  #Perform svm (Support Vector Machines (SVM) with a linear kernel) algorithm to training data
  #Method was tested for svm, random forest, linear discriminant analysis, k-nearest neighbour and classification and regression trees
  #svm yielded best overall result, and highest correct prediction for both species
  #Simpler methods (e.g. lda) tended to have only slightly lower overall success, but resulted in clearly higher false predictions
  #for Paramecium (15% false versus 6 % false)
  {
    # Run algorithms using 10-fold cross validation
    control <- trainControl(method="cv", number=10) #DEFINE THE CONTROL VARIABLES 
    metric <- "Accuracy" #ACCURACY TO SELECT BEST MODEL 
    
    # Fit the model
    model <- train(species~., data=training_data, method="svmRadial", metric=metric, trControl=control) #FIT THE MODEL USING SVM (BUT YOU MIGHT USE ALSO ANOTHER METHOD)
    
    #Verify the data using the rest of the monoculture data
    test.data <- predict(model, verification_data) #PREDICT THE OUTCOME OF THE VERIFICATION DATASET 
    confusionMatrix(test.data, as.factor(verification_data$species)) #SHOW THE ACCURACY OF THE METHOD 
    }
}

#Now perform prediction on community data
{
  # Make predictions
  predictions <- model %>% predict(filter(data_traits, Culture!="monoculture")) #PREDICT THE COMMUNITY DATA
  
  #Add the species information to the community data
  predicted <- filter(data_traits, Culture!="monoculture") #ADD PREDICTIONS TO THE DATAFRAME 
  predicted$species <- as.character(predictions) #ADD PREDICTIONS 
}

#Summarize the data to get densities (order by treatments, and calculate population
#density), by scaling for emasured volume and number of frames
data_communities <- predicted %>% group_by(ID, Species_1, Species_2, Stress_treatment, experiment_time, species) %>% #CALCULATE POPULATION DENSITY IN EACH VIDEO 
  summarize(indiv_per_mL=sum(N_frames)/250*1000/34.4, mean_area=mean(mean_area))
data_communities <- data_communities %>% arrange(ID, experiment_time) #ORDER IT
