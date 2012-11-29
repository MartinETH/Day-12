# Assignment 12 #
# Martin Gubler #



## Event sequence object (copy/paste)

library(TraMineR)
data(biofam)
biofam$cohort <- cut(biofam$birthyr, c(1900,1930,1940,1950,1960),
                     labels=c("1900-1929", "1930-1939", "1940-1949", "1950-1959"),
                     right=FALSE)
bf.states <- c("Parent", "Left", "Married", "Left/Married", "Child",
               "Left/Child", "Left/Married/Child", "Divorced")
bf.shortlab <- c("P","L","M","LM","C","LC", "LMC", "D")
biofam.seq <- seqdef(biofam[,10:25], states=bf.shortlab,
                     labels=bf.states, weights=biofam$wp00tbgs)
weight <- attr(biofam.seq, "weight")

tm <- c(
  ## "P",  "L",   "M",   "LM",  "C",     "LC",    "LMC",   "D"
  ##------------------------------------------------------------
  "P",   "L",   "M",   "L,M", "C",     "L,C",   "L,M,C", "M,D",
  "P",   "L",   "P,M", "M",   "P,C",   "C",     "M,C",   "M,D",
  "",    "D,L", "M",   "L",   "D,C",   "D,L,C", "L,C",   "D",
  "P,D", "D",   "P",   "L,M", "D,P,C", "D,C",   "C",     "D",
  "P",   "L",   "M",   "L,M", "C",     "L",     "L,M",   "M,D",
  "P",   "",    "M",   "M",   "P",     "L,C",   "M",     "D",
  "P,D", "L",   "",    "",    "D,P",   "D",     "L,M,C", "D",
  "P",   "L",   "P,M", "M",   "P,C",   "C",     "M,C",   "D"
)
tm <- matrix(tm, nrow=8, ncol=8, byrow=TRUE)
dimnames(tm) <- dimnames(seqetm(biofam.seq, method="transition",
                                use.labels=FALSE))
tm

bf.seqe <- seqecreate(biofam.seq, tevent=tm, use.labels=FALSE)


# 1. Compare the 10 most frequent subsequences between ages 20 and 29 to the 10
# most frequent between 26 and 35, considering only subsequences with at least two events. 
# (Tips: Minimum and maximum ages should be specified relatively to the first position;
# e.g., since sequence start is at 20 years old, to specify age 29, you should give 9 (= 29 􀀀 20.)

constraint.2029 <- seqeconstraint(ageMin=0,ageMaxEnd=9)
constraint.2635 <- seqeconstraint(ageMin=6, ageMaxEnd=15)
biofam.subseqe.2029 <- seqefsub(bf.seqe, constraint=constraint.2029, pMinSupport=0.01)
biofam.subseqe.2635 <- seqefsub(bf.seqe, constraint=constraint.2635, pMinSupport=0.01)
biofam.subseqe.2029[1:10,]
biofam.subseqe.2635[1:10,]

# 2. Compare the counts of the three patterns (M) --> (C), (C) --> (M).

subseq <- c("(M)-(C)","(C)-(M)")
biofam.subseq <- seqefsub (bf.seqe, strsubseq = subseq)
biofam.subseq

# 3. Find thefive association rules exhibiting the highest lift value among those with a
# 15% minimum support. Are these rules interesting?.

biofam.subseq <- seqefsub(bf.seqe, pMinSupport = 0.15)
biofam.rules <- TraMineR:::seqerules(biofam.subseq)
biofam.rules[order(biofam.rules$Lift, decreasing = TRUE)[1:10], 1:4]

## Hardly any "lift" visible (not much increased likelihood of experiencing one of these states after starting sequence with "P")
## Contentwise, these links do not seem very surprising (apart from the fact that some people seem to be staying at home despite having married (and even becoming a parent))

# 4. Define a vector of event weights inversely proportional to the event frequencies

## Extracting frequencies of involved events
print (event.frequency <- seqefsub(bf.seqe, maxK = 1, minSupport = 0))

## Computing weights and sort weigths according to alphabet

freq <- event.frequency$data$Count

event.alphabet <- levels(bf.seqe)
event.weight <- matrix(sum(freq)/freq, length(freq), 1)

event.names <- as.character(event.frequency$subseq)
event.names <- sub("(", "", event.names, fixed = TRUE)
event.names <- sub(")", "", event.names, fixed = TRUE)
pm <- pmatch(event.alphabet, event.names)
rownames(event.weight) <- event.names
event.weight <- event.weight[pm, ]
round(event.weight, 2)

# 5. Using the previously computed event weights together with a time shift cost equal
# to 0.1, compute the matrix of normalized OME dissimilarities.

library(TraMineRextras)
biofam.indel.costs <- sqrt(event.weight)
dOME <- seqedist(bf.seqe, idcost = biofam.indel.costs, vparam = 0.1, norm = "YujianBo")
summary(as.vector(dOME))
dOME[1:10,1:10]


# 6. Compute the four-cluster solution using the weighted PAM method and label the
# clusters with their respective medoid sequence. 
# (Tips: Retrieve the medoids with function disscenter together with argument medoids.index="first".)


library(WeightedCluster)
set.seed(1)
biofam.seqe.pam <- wcKMedoids(dOME, k = 4, weight = seqeweight(bf.seqe))
event.cluster4 <- biofam.seqe.pam$clustering
medoids <- disscenter(dOME, group=event.cluster4, medoids.index="first")
event.cluster4.label <- as.character(bf.seqe[medoids])
event.cluster4.factor <- factor(event.cluster4, labels = event.cluster4.label)

#### I was able to build my own code following the sample solution, but I could not produce lines 107/108 myself.

# 7. Render the clusters using (a) the parallel coordinate plot, and (b) the survival
# probability curves until first occurrence of each event except `P'.


## PCP
alphabet <- c("P","L","M","C","D")
seqpcplot(bf.seqe, group=event.cluster4.factor, alphabet=alphabet,
          filter = list(type = "function",
                        value = "cumfreq",
                        level = 0.5),
          order.align = "first",
          ltype = "non-embeddable",
          cex = 1.5, lwd = .9,)


## Survival of events
seqedplot(bf.seqe, group = event.cluster4.factor, lwd = 4, ignore = c("P"))