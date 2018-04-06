# LunchinatoRs April 6, 2018
# Erica Baken

setwd("~/Documents/School/Presentations/2018 LunchinatoR")
library(phytools)
library(stringr)

# Read in the data
BB <- read.tree("BBPruned.tre")
Habitats <- read.csv("HabitatsPruned.csv")
MajRule <- read.csv("MajRuleAveragesClass3.csv") # This is calculated over the 1000 posterior trees, which is why they aren't whole numbers

# Modify data to the right format 
BB$tip.label <- str_replace_all(BB$tip.label, "[[:punct:]]"," ") 
rownames(MajRule) <- MajRule[,1] 
  MajRule <- MajRule[,-1]
Class3 <- Habitats$Class3
  names(Class3) <- Habitats$Species    
      
# Calculating Null Expected
        # This method is taking all known transitions and evenly distributing them
        # according the microhabitat type distribution across rows. More populated microhabitats
        # have an inherently higher chance of transitioning out of that microhabitat just 
        # because there are more of them

sums <- apply(MajRule, 1, function(x){
  sum(na.omit(x))
})


Deltaxy <- matrix(NA, nrow = 5, ncol = 4)

for (i in 1:length(sums)){
x<- summary(Class3)[-i]
Deltaxy[i,] <- sums[i]*(x/sum(x))
}

MajRule2 <- t(apply(MajRule,1,function(x){na.omit(x)}))

# Chi Square Tests
    # Overall Differences
    round((MajRule2- Deltaxy)/sqrt(Deltaxy), 2)
    chi2 <- sum(na.omit(as.vector((MajRule2- Deltaxy)^2/Deltaxy))) # for each cell, (Expected - Observed)^2/Expected 
    pchisq(chi2, df = 15, lower.tail = F) # df = 20 - 5; p = 1.258 X 10^-7
    
    # Transitions TO Arboreality 
    ToA <- MajRule[-1,1]
    ToANull <- Deltaxy[-1,1]
    chi2 <- sum((ToA - ToANull)^2/ToANull)
    pchisq(chi2, df = 3, lower.tail = F) # p = 0.43
    
    # Transitions FROM Arboreality
    FromA <- MajRule2[1,]
    FromANull <- Deltaxy[1,]
    chi2 <- sum((FromA - FromANull)^2/ FromANull)
    pchisq(chi2, df = 3, lower.tail = F) # p = 0.35
