# LunchinatoRs April 6, 2018
# Erica Baken

setwd("~/Documents/School/Presentations/Lunchinators April 6 2018")
library(phytools)
library(stringr)

# Read in the data
BB <- read.tree("BBPruned.tre")
Habitats <- read.csv("HabitatsPruned.csv")
MajRule <- read.csv("MajRuleAveragesClass3.csv") # This is calculated over the 1000 posterior trees, which is why they aren't whole numbers

# Modify data to the right format 
BB$tip.label <- str_replace_all(BB$tip.label, "[[:punct:]]"," ") 
Class3 <- Habitats$Class3
  names(Class3) <- Habitats$Species    

Info <- summary(Class3) 
      
rownames(MajRule) <- MajRule[,1] 
MajRule <- MajRule[,-1]
      
# Calculating Null Expected
        # This method is taking all known transitions and evenly distributing them
        # according the microhabitat type distribution across rows. More populated microhabitats
        # have an inherently higher chance of transitioning out of that microhabitat just 
        # because there are more of them

DeltaTot <- sum(na.omit(unlist(MajRule))) # Adding up all transitions according to MajRule = 90.331
n <- 332 # Number of species
m <- 5 # Number of microhabitat types
nx <- Info
    
      Deltaxy <- matrix(NA, nrow = 5, ncol = 5)
      for (i in 1:5){ Deltaxy[i,] <- (DeltaTot*nx[i]/n)/(m-1) }
      diag(Deltaxy) <- 0
      rownames(Deltaxy) <- colnames(Deltaxy) <- names(nx)
      Deltaxy
    
# Chi Square Tests
    # Overall Differences
    chi2 <- sum(na.omit(as.vector((MajRules- Deltaxy)^2/Deltaxy))) # for each cell, (Expected - Observed)^2/Expected 
    pchisq(chi2, df = 16, lower.tail = F) # df = (r-1)(c-1); p = 1.258 X 10^-7
    
    # Transitions TO Arboreality
    chi2 <- sum((MajRules[-1,1] - Deltaxy[-1,1])^2/ Deltaxy[-1,1])
    pchisq(chi2, df = 3, lower.tail = F) # p = 8.759 X 10-8
    
    # Transitions FROM Arboreality
    chi2 <- sum((MajRules[1,-1] - Deltaxy[1,-1])^2/ Deltaxy[1,-1])
    pchisq(chi2, df = 3, lower.tail = F) # p = 0.009
