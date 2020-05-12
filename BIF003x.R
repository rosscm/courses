# edX BIF003x
# Statistical Analysis in Bioinformatics

########
## Week 2 - Basic DNA Analysis
########

library(seqinr)
library(ape)
library(rentrez)
library(Peptides)

# Read in FASTA file
cox1 <- read.fasta("cox1.fasta", seqtype = "AA")

# Retrieve GENBANK sequence as binary object
AB003468 <- read.GenBank("AB003468", as.character = "TRUE")

# Save sequence in FASTA
write.dna(AB003468, file ="AB003468.fasta", format = "fasta", append = FALSE, nbcol = 6, colsep = " ", colw = 10)

# Search for sequences using rentrez_search
entrez_search(db = "nucleotide", term = "human superoxide dismutase")

# Convert sequence to string for analysis
CloningVector <- AB003468[[1]]

# Get total count of each nucleotide
count <- count(CloningVector, 1)

# See all dinucleotide or trinucleotide combos
count(CloningVector, 2)
count(CloningVector, 3)

# Retrieve GC content
## Tells us the ratio of G/C residues compared to A/T residues
## important because coding regions tend to be higher in GC
## also affects the "melting" temperature of DNA
GC <- GC(CloningVector)

# Calculate GC content by windows of the sequence
## the above value is global; would be more interesting to see how GC content
## varies over the course of a sequence (eg entire genome)
slidingwindowGCplot <- function(windowsize, inputseq) {
  GCwindow <- seq(1, length(inputseq)-windowsize, by = windowsize)

  # Get number of sequence chunks
  n <- length(GCwindow)
  Chunks <- numeric(n)

  # GC content per chunk
  for (i in 1:n) {
    chunk <- inputseq[GCwindow[i]:(GCwindow[i]+(windowsize-1))]
    chunkGC <- GC(chunk)
    print(chunkGC)
    Chunks[i] <- chunkGC
  }

  # Plot
  plot(GCwindow, Chunks, type = "b", xlab = "Nucleotide start position", ylab = "GC content",
       main = paste("GC plot with windowsize", windowsize))
}

# Run function
slidingwindowGCplot(200, CloningVector)

# Calculate basic protein sequence statistics
# Determie amino acid composition
aaComp(cox1[1])
aaComp(cox1)

# Get the aliphatic index (an indicator of thermostability of globular proteins)
aIndex(cox1)

# Predict net charge of the protein
# Ouputs 1 net charge for each sequence
charge(cox1)

# Calculate hydrophobicity
hydrophobicity(cox1)

########
## Week 3 - Sequence Alignment
########

library(Biostrings)
library(seqinr)
library(msa)
library(ape)
library(phangorn)

# Read in sequence file
prokaryotes <- read.fasta(file = "prok.fasta", seqtype = "DNA")

# Split into individual sequences for first pairwise alignment
seq1 <- as.character(prokaryotes[[1]])
seq1 <- paste(seq1, collapse = "")
seq2 <- as.character(prokaryotes[[2]])
seq2 <- paste(seq2, collapse = "")

# Align seq1 and seq2
pairalign <- pairwiseAlignment(pattern = seq2, subject = seq1)

# Get alignment summary stats
summary(pairalign)

# Convert alignment to strings to export FASTA
pairalignString <- BStringSet(c(toString(subject(pairalign)),toString(pattern(pairalign))))
writeXStringSet(pairalignString, "aligned.txt", format = "FASTA")

# Create dotplot of two aligned sequences
coxgenes <- read.fasta(file = "cox1multi.fasta", seqtype = "AA")
cox1 <- as.character(coxgenes[[1]])
cox2 <- as.character(coxgenes[[2]])
dotPlot(cox1, cox2, main = "Human vs Mouse Cox1 Dotplot")
dotPlot(cox1, cox2, wsize = 3, wstep = 3, nmatch = 3, main = "Human vs Mouse Cox1 Dotplot\nwsize = 3, wstep = 3, nmatch = 3")

# Only use residues 1-100
dotPlot(cox1[1:100], cox2[1:100], wsize = 3, wstep = 3, nmatch = 3, main = "Human vs Mouse Cox1 first 100 AA Dotplot\nwsize = 3, wstep = 3, nmatch = 3")

# Read in multiple sequence files using msa
coxAA <- readAAStringSet("cox1multi.fasta")
prokDNA <- readDNAStringSet("prok.fasta")

# Align sequences using CLUSTALW (default)
coxAligned <- msa(coxAA)
prokAligned <- msa(prokDNA)

# Print untruncated alignment
print(prokAligned, show = "complete")

# Genertae color-coded printout
msaPrettyPrint(prokAligned)

# Cast msa alignment to StringSet to export as FASTA
prokAlignStr <- as(prokAligned, "DNAStringSet")
writeXStringSet(prokAlignStr, file = "prokAligned.fasta")

coxAlignStr <- as(coxAligned, "AAStringSet")
writeXStringSet(coxAlignStr, file = "coxAligned.fasta")

# Export alignment in PHYLIP format
write.phylip(coxAligned, "coxAligned.phylip")

# Convert prokaryotic alignment to seqinr format for phylogeny analysis
prokAligned2 <- msaConvert(prokAligned, type = "seqinr::alignment")

# Generate distance matrix
# Distances calculated based on whether or not the positions are identical
prokdist <- dist.alignment(prokAligned2, "identity")

# Use neighbour-joining method to construct basic distance phylogenetic tree
# Uses pairwise distances to cluster the closest sequence fist, then
# iterively adds sequences to that pair
prokTree <- nj(prokdist)
plot(prokTree)

# A more evolutionary-centric approach is maximum parsimony, which actually looks at shared
# characters between sequences and maximizes similarity while minimizing the number of steps
# required to construct a the phylogenetic history between the "cades", or species.
# NOTE very computationally expensive (<30 sequences)

# Genrate phylogenetic data, phyDat, for parsimony tree
prokAligned3 <- msaConvert(prokAligned, type = "phangorn::phyDat")

# Pratchet automatically does branch rearrangements to find "better"
# (more parsimonious) trees
ParsTree <- pratchet(prokAligned3)
plot(ParsTree)

# One variation on the distance approach is the maximum likelihood method.
# Maximum likelihood lets you use a distance matrix (as previously) but also includes the full sequence
# alignment and calculates the statistical likelihood of a tree given the full data.
# As such maximum likelihood is far more computationally intense, but it may also provide more
# accurate reconstructions than "conventional" distance approaches.
fit <- pml(prokTree, prokAligned3)

# Now that we have evaluated the statistical likelihood, or fit, of the tree using pml, we can optimize
# this tree using optim.pml.

# optim.pml requires us to provide the "model" of evolution that we think will best fit our data.
# The model we will use for this exercise is the Jukes-Cantor model, which assumes that all bases
# are found in equal numbers, and that all mutations are equally likely (e.g. A to T or C or G, etc.).
# Another "popular" model is K80, or the Kimura two-parameter model, which states that
# transition mutations (A <-> G)  will have a different rate from transversions (C<->T).

# Rearrangement indicates what kind of branch swapping the algorithm will try
# when attempting to improve the three – the options are stochastic, NNI (nearest neighbor),
# and ratchet.
fitJC <- optim.pml(fit, model = "JC", rearrangement = "stochastic")
plot(fitJC)
fitJC <- optim.pml(fit, model = "K80", rearrangement = "stochastic")
plot(fitJC)

# Bootstrap for statistical power
# Use fit ML tree with 100 bootstrap replicates and nearest-neighbour
# interchange for branch swapping
bootstrapped <- bootstrap.pml(fitJC, bs=100, optNni=TRUE, multicore=TRUE, control = pml.control(trace=0))

# Use the phangorn bootstrap plotting routine and takes all of the fitJC trees generated in the
# bootstrapped dataset and displays results where the support (number of trees that agree out of the
# 100) is greater than 50; the type controls the tree display (here we use p for phylogram).
plotBS(midpoint(fitJC$tree), bootstrapped, p = 50, type = "p")
