pierceBCA <- function(stdAbs, stdConcs, std_N = 2, sampleNames, sampleAbs, sample_N = 2, wellVol, DF = 5, LBX = 5, target_ug = NULL) {
  
  # Author: Eric Prince
  # Date: 2017-07-17
  
  # Western blot assays are commonplace tools in molecular biology research used to observe changes in protein expression.  
  # In order to correctly complete this assay, we need to ensure that we load a consistent amount of protein in each lane
  # of our SDS-PAGE or bis-tris gels.  To get an idea of our protein concentration, we use an external calibration curve that
  # relates the concentration of protein in µg/mL to an absorbance value at 490 nm wavelength.  These values must then be collected,
  # plotted, and analyzed to determine what the proper volume of lysate is for each given sample.  Often times, we will tinker with
  # these values and determine what the highest amount of protein is, given that our lowest concentrated sample is our limiting
  # reagent.  Also, we find missing standard at times and therefor change our x-axis values, making an excel style template some-
  # what ill-suited compared to alternative options.  Here, I present a function to streamline this process.  After feeding in the
  # required values, the function will optimize for the highest lysate mass and deliver a statistical model and results table for
  # final sample preparation.
  
  # This function results in a list.  It is best to save the output to the list and extract the plot as the 1st item in the list,
  # and the statistical table as the 2nd item in the list.  The linear model statistics will be printed, as well as the optimal
  # protein mass.
  
  ###################################################################
  #*****************************************************************#
  #                     FUNCTION ARGS                               #
  #*****************************************************************#
  # stdAbs: a vector list of absorbance values obtained for standards
  # stdConcs: a vector list of unique standard concentrations (given in the same order as stdAbs), typically given in µg/mL
  # std_N: an integer number of replicates for standards, typically N = 2
  # sampleNames: a vector list of unique sample names to be quantified
  # sampleAbs: a vector list of absorbance values obtained for samples
  # sample_N: an integer number of replicates for samples
  # wellVol: Volume of the lane well in µL
  # DF: dilution factor for BCA sample prep.  This is usually 5X, as I add 2 µL of lysate into 8 µL of RIPA for analysis
  # LBX: dilution factor for loading buffer; typically a 5X loading buffer is used
  # target_ug: the target mass of lysate used in the sample.  If NULL, the target mass will be optimized for the maximum ug in
  #            the lowest concentration sample.
  ###################################################################
  ###################################################################
  
  library(tidyverse)
  library(magrittr)
  library(cowplot)
  
  # iterate standard concentration values alongside absorbance data
  std_concentrations <- c()
  for (i in stdConcs) {
    std_concentrations <- c(std_concentrations, rep(i, std_N))
  }
  
  # iterate sample names for absorbance value list
  sample_names <- c()
  for (i in sampleNames) {
    sample_names <- c(sample_names, rep(i, sample_N))
  }
  
  # Error handling
  if (length(std_concentrations) != length(stdAbs)) {
    print("The length of the standard concentrations list does not the length of the standard absorbance values")
  }
  
  # Merge information about the standards together and calculate mean values
  standard_df <- data_frame(std_concentrations,
                   stdAbs)
  
  summarized_df <- standard_df %>%
    group_by(std_concentrations) %>%
    summarise(mean_value = mean(stdAbs))
  
  # plot the summarized data for reference
  plt <- summarized_df %>%
    ggplot(aes(std_concentrations, mean_value)) +
    geom_point() +
    geom_smooth(method = lm) +
    labs(title = paste("BCA Standard Curve - ", Sys.Date()),
         y = "Absorbance, 490 nm",
         x = "Lysate Concentration, µg/mL") +
    theme_cowplot()
    
  # Determine linear model statistics
  lm_result <- lm(summarized_df$mean_value ~ summarized_df$std_concentrations)
  print(summary(lm_result))
  
  # merge and summarize data regarding samples
  sample_df <- data_frame(sample_names, sampleAbs) %>%
    group_by(sample_names) %>%
    summarise(mean_value = mean(sampleAbs))
  
  # Use the inverse equation of the linear model to calculate the sample concentration based on absorbance values
  sample_df$sample_concs <- (sample_df$mean_value - lm_result$coefficients[1]) / lm_result$coefficients[2]
  # Adjust concentrations to µg/µL
  sample_df$sample_concs_df.corr <- sample_df$sample_concs * DF / 1000
  # Add column for loading buffer volumes
  sample_df$loading_buffer <- rep((wellVol/LBX), nrow(sample_df))
  
  # Check if target concentration has been given
  if (!is.null(target_ug)) {
    sample_df$sample_uL <- target_ug / sample_df$sample_concs_df.corr
  }
  
  # If target concentration was not given, optimize the mass of protein available
  else {
    # determine sample with lowest concentration
    min_abs_sample <- sample_df[which.min(sample_df$sample_concs_df.corr),]
    # calculate the available volume in the well by subtracting loading buffer volume from the well volume
    vol_lysate <- wellVol - (wellVol/LBX)
    # multiply the volume available by the lysate concentration to get µg protein, and round down to the nearest int
    optim_ug <- trunc(vol_lysate * min_abs_sample$sample_concs_df.corr)
    # caluclate the volume required to obtain the optimal µg of lysate
    sample_df$sample_uL <- optim_ug / sample_df$sample_concs_df.corr
    # Print what the optimal µg value is
    print(paste("The optimal mass of protein is", optim_ug, "µg"))
  }
  
  # Calculate how much ripa will be required to have an even volume for each sample
  sample_df$ripa_uL <- wellVol - sample_df$loading_buffer - sample_df$sample_uL
  
  # Return each object together in a list; Only way to return both from the function..
  return(list(plt, sample_df))
}
