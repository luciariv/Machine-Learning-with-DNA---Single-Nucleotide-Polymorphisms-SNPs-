# Machine-Learning-with-DNA---Single-Nucleotide-Polymorphisms-SNPs-
install.packages("vcfR")
library(vcfR)
library(ade4)
library(adegenet)

VCF <- read.vcfR("filogenia_milano.vcf")
# Since VCF data does not typically include any sort of population information, I have to load this data
# separately from the text-delimited file (my popmap file)
pop.data <- read.table("filogenia_milano.tsv", sep = "\t", header = TRUE)


# We can now check that all the samples in the VCF and the population data frame are included:
all(colnames(VCF@gt)[-1] == pop.data$V1)

#genetic scores per individual
C.genind<- vcfR2genind(VCF)

#and check the population information is there
pop(C.genind) <- pop.data$Population

library("poppr")
data("C.genind")

library(ggplot2)

#PCA
data(C.genind)
x.hakes <- tab(C.genind, freq=TRUE, NA.method="mean")
pca.hakes <- dudi.pca(x.hakes, center=TRUE, scale=FALSE)
3

dudi.pca(df = x.hakes, center = TRUE, scale = FALSE, scannf = FALSE, nf = 3)
pca.hakes

#Visualization of genotypes coordinates:
s.label(pca.hakes$li)


# Distribution of the individuals from different groups.
s.class(pca.hakes$li, fac=pop(C.genind),
        col= transp(funky(9),.9),
        axesel=FALSE, cstar=0, cpoint=3, grid = FALSE, sub="PC1,PC2")

add.scatter.eig(pca.hakes$eig[1:50],3,1,2, ratio=.2 , posi=c("bottomright"))

#funky originlaly 15



#Let us examine the second and third axes:
#What is the major factor of genetic differentiation in hake stocks?
#What is the second one? What is the third one?

s.class(pca.hakes$li, fac=pop(C.genind),
        xax=2, yax=3, col=transp(funky(9),.9),
        axesel=FALSE, cstar=0, cpoint=3, sub = "PC2,PC3", grid = FALSE)
add.scatter.eig(pca.hakes$eig[1:50],3,2,3, ratio=.2,  posi=c("bottomright"))


# PC1 & PC3

s.class(pca.hakes$li, fac=pop(C.genind),
        xax=1, yax=3, col=transp(funky(9),.9),
        axesel=FALSE, cstar=0, cpoint=3, sub = "PC1,PC3", grid = FALSE)
add.scatter.eig(pca.hakes$eig[1:50],3,2,3, ratio=.2,  posi=c("bottomright"))

#Variance of first and second  principal components:
pca.hakes$eig[1]

pc1 <- pca.hakes$li[,1]
var(pc1)

var(pc1)*703/704

mean(pc1^2) 
n <- length(pc1)
0.5*mean(dist(pc1)^2)*((n-1)/n)

#Variances (abs values to percentages)

eig.perc <- 100*pca.hakes$eig/sum(pca.hakes$eig)
head(eig.perc)

#Allele contributions (frequencies) / 1 axis at a time:

loadingplot(pca.hakes$c1^2)


## How many genotypes possess one copy of ONE SPECIFIC allele?

X <- tab(C.genind) 
class(X)
dim(X)
SNP3520_fpt_61.0 <- X[, "3520_fpt_61.0"]
table(SNP3520_fpt_61.0)

#Which individuals are they?

rownames(X)[SNP3520_fpt_61.0 == 0.5]


#Contribution of alleles 1st Axis

myLoadings <- pca.hakes$c1[, 1]^2
names(myLoadings) <- rownames(pca.hakes$c1)
loadingplot(myLoadings)
loadingplot(myLoadings, xlab = "Alleles")
loadingplot(myLoadings, xlab = "Alleles", ylab = "Weight of the alleles",
            main = "Contribution of alleles \n to the first sPCA axis")

#Contribution of alleles 2nd Axis
loadingplot(pca.hakes$c1^2, axis=2, main = "Allele contributions to the PC2")

#Contribution of alleles 3rd Axis
loadingplot(pca.hakes$c1^3, axis=3, main = "Allele contributions to the PC3")

#5% alleles contributing most to showing the diversity within hake stocks:

temp <- loadingplot(myLoadings, threshold = quantile(myLoadings, + 0.95), xlab = "Alleles", ylab = "Weight of the alleles",
                    main = "Greatest contribution of alleles in hakes diversity (5%) ")

temp

# Assess the average contribution of each marker

boxplot(myLoadings ~ C.genind$loc.fac, las=3, main="Contributions by markers \nto the first global score")


#allele frequencies per population

hakes.genpop <- adegenet::genind2genpop(C.genind)
Freq <- adegenet::makefreq(hakes.genpop)
round(Freq[1:4, 1:26])

#Checking
apply(Freq, MARGIN = 1, FUN = sum)    # Just checking

write.csv(Freq,file="Freqperpop.txt",quote = F)

#allele frequencies per population (ok)
genpop <- genind2genpop(C.genind)
genpop
adegenet::makefreq(genpop)
