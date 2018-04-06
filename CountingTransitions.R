# Counting Transitions of Discrete Traits
# Erica Baken

setwd("~/Documents/School/Presentations/Lunchinators April 6 2018")
library(phytools)
library(stringr)

# Reading in and manipulating data
BB <- read.tree("BBPruned.tre")
BB$tip.label <- str_replace_all(BB$tip.label, "[[:punct:]]"," ")    
Habitats <- read.csv("HabitatsPruned.csv")
Class3 <- Habitats$Class3
names(Class3) <- Habitats$Species    

# MAJORITY RULES
    
MajorityRules <- function(from, to, data, phy) {

      ASR <- ace(x = data, phy = phy, type = "discrete", method = "ML") # Ancestral Character Estimation
      x <- which(ASR$lik.anc[ ,from] > 0.5)  # Which internal nodes have higher than 50% chance of being terrestrial
      FromOver50 <- x + Ntip(phy) # Assigning the correct label numbers to each node
      Branches1 <- match(phy$edge[ ,1], FromOver50)  
      BranchNum1 <- which(Branches1 != "NA") # All the branch numbers where the beginning node is over 50% terrestrial
      
      y <- which(ASR$lik.anc[ ,to] > 0.5) # Which internal nodes have higher than 50% chance of being arboreal
      ToOver50 <- y + Ntip(phy) # Assigning the correct label numbers to each node
      ToSubset <- names(data[data==to]) # Adding tips of arboreal species to list of ASR nodes that are over 50% arboreal. 
      # Necessary for counting transitions that happen at the tips of the tree
      To <- match(phy$tip.label, ToSubset) # Matching the Tree's tip names to the Arboreal Subset species names
      BranchNumTo <- which(To != "NA")
      ToOver50Plus <- c(ToOver50, BranchNumTo) # Combining tips with internal nodes
      Branches2 <- match(phy$edge[ ,2], ToOver50Plus) # Matching the tree's tip names to all nodes or species that are arboreal
      BranchNum2 <- which(Branches2 != "NA") # All the edge numbers where the ending nodes are over 50% arboreal
      
      Transitions <- intersect(BranchNum1, BranchNum2) # Collects the edges (row number) where the edge beginning is >50% terrestrial and the edge end is >50% arboreal
      
      # A to T
      Branches3 <- match(phy$edge[ ,1], ToOver50)
      BranchNum3 <- which(Branches3 != "NA") # All the branch numbers where the beginning node is over 50% arboreal
      FromSubset <- names(data[data==from]) # Adding tips of terrestrial species to list of ASR nodes that are over 50% terrestrial 
      # Necessary for counting transitions that happen at the tips of the tree
      From <- match(phy$tip.label, FromSubset) # Matching the Tree's tip names to the Terrestrial Subset species names
      BranchNumFrom <- which(From != "NA")
      FromOver50Plus <- c(FromOver50, BranchNumFrom) # Combining tips with internal nodes
      Branches4 <- match(phy$edge[ ,2], FromOver50Plus) # Matching the Tree's tip names to all nodes or species that are arboreal
      BranchNum4 <- which(Branches4 != "NA") # All the edge numbers where the ending nodes are over 50% arboreal
      
      Transitions.Back <- intersect(BranchNum3, BranchNum4) # Collects the edges (row number) where the edge beginning is >50% arboreal and the branch end is >50% terrestrial
      
      Out <- list(length(Transitions), length(Transitions.Back), Transitions, Transitions.Back)
      names(Out) <- c('TotalTransitions', 'TotalBackTransitions', 'TransitionBranches', 'BackTransitionBranches')
      Out
}

    Class3TtoA <- MajorityRules(from = "T",to = "A", data = Class3, phy = BB)
    Class3FtoA <- MajorityRules(from = "F",to = "A", data = Class3, phy = BB)
    Class3RtoA <- MajorityRules(from = "R",to = "A", data = Class3, phy = BB)
    Class3WtoA <- MajorityRules(from = "W",to = "A", data = Class3, phy = BB)
    Class3TtoF <- MajorityRules(from = "T",to = "F", data = Class3, phy = BB)
    Class3TtoR <- MajorityRules(from = "T",to = "R", data = Class3, phy = BB)
    Class3TtoW <- MajorityRules(from = "T",to = "W", data = Class3, phy = BB)
    Class3FtoR <- MajorityRules(from = "F",to = "R", data = Class3, phy = BB)
    Class3FtoW <- MajorityRules(from = "F",to = "W", data = Class3, phy = BB)
    Class3RtoW <- MajorityRules(from = "R",to = "W", data = Class3, phy = BB)
          
      Class3TtoA # 31 forward; 6 back
      Class3FtoA # 0 forward; 0 back
      Class3RtoA # 0 forward; 2 back
      Class3WtoA # 0 forward; 0 back
      Class3TtoF # 19 forward; 1 back
      Class3TtoR # 16 forward; 1 back
      Class3TtoW # 7 forward; 2 back
      Class3FtoR # 0 forward; 2 back
      Class3FtoW # 0 forward; 1 back
      Class3RtoW # 0 forward; 2 back


# STOCHASTIC MAPPING
Cols2 <- c("Arboreal" = "#7fbf7b", "Saxicolous" = "#fee090", 
           "Terrestrial" = "#af8dc3", "Aquatic" = "#91bfdb", 
           "Fossorial" = "#762a83")

names(Cols2) <- c("A", "R", "T", "W", "F") # Naming the colors for matching Class3 classifications
pi <- c(0,0,1,0,0)
names(pi) <- c("A", "R","T", "W", "F") 
levels(Class3) # Look at the order in which Class3 factors appear and this should be the order of pi 
pi2 <- pi[c(1,5,2,3,4)]
              
# FITTING THESE MODELS TAKES A LONG TIME
fit <- fitMk(BB, Class3, model="ARD", pi = pi2) # around 10 minutes; fitting an MK model; "ARD" stands for All Rates Different model.  
fit$logLik # -320.6549
    
fit2 <- fitMk(BB, Class3, model="SYM", pi = pi2)  # around 3 minutes
fit2$logLik # -356.8753
                  
fit3 <- fitMk(BB, Class3, model="ER", pi = pi2)  # around 1 min
fit3$logLik # -383.0709

extRemes::lr.test(x = fit$logLik, y = fit2$logLik, df = 10) # df is number of free parameters. ARD is better than SYM
extRemes::lr.test(x = fit$logLik, y = fit3$logLik, df = 19) # ARD is better than ER
                                
treesARD <- make.simmap(BB, Class3, pi = pi2, Q = "empirical", model = "ARD", nsim = 100) # around 9 minutes
treesER <- make.simmap(BB, Class3, pi = pi2, Q = "empirical", model = "ER", nsim = 100)
treesSYM <- make.simmap(BB, Class3, pi = pi2, Q = "empirical", model = "SYM", nsim = 100)
 
sumARD <- summary(treesARD)                
Pies <- sumARD$ace # can use this for plotting Ancestral State Estimations on phylogeny
                  
                  
# Can do fun plotting with binary traits                  
Class1 <- Habitats$Class1
names(Class1) <- Habitats$Species  
pi <- c(0,1)
names(pi) <- c("A", "T")
treesER.Class1 <- make.simmap(BB, Class1, pi = pi, Q = "empirical", model = "ER", nsim = 100)
TreeSummary <- describe.simmap(treesER.Class1)
TreeSummary$count
TreeSummary$ace
plot(TreeSummary, type = "fan")

XX <- densityMap(treesER.Class1, res = 100, plot=F)  
n <- length(XX$cols) # Making a vector so you can color it in gradient
XX$cols[1:n] <- colorRampPalette(Cols2, space="Lab")(n) # Making color gradient
plotSimmap(XX$tree, colors = XX$cols, type = "fan", fsize = .01, lwd = 3)