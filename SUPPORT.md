**InRelate** is an R package that implements the likelihood method of Sethuraman (2018) to estimate pairwise genetic relatedness while accounting for the presence of population structure, inbreeding, and admixture. Follow the steps below to run InRelate.

**Step 1**: Use your program of choice to estimate population structure under the admixture model - for e.g. STRUCTURE (Pritchard et al., 2000), MULTICLUST (Sethuraman 2013), ADMIXTURE (Alexander et al., 2009); pick the optimal number of subpopulations (often denoted as K) using your method of choice - e.g. Evanno et al., 2005, Dent and VonHoldt 2012, Alexander and Lange 2013 (cross validation errors), or Sethuraman 2013 (bootstrapping followed by model choice via AIC/BIC). Then use the outputs of the above programs to format the subpopulation allele frequencies (here called the "pkla" file) and admixture proportions (here called the "indivq" file). We will also require the STRUCTURE formatted multilocus genotype file (here called the "str" file) as input.

**Step 2**: Install the InRelate package - you can do this by downloading the .tar.gz file and sourcing it. Or if you're running it in Rstudio, just click on Tools->Install Packages. Under Install From, select Package Archive File (.tar.gz), click Install. This will install InRelate onto your machine. 

**Step 3**: Read in your datasets and calculate relatedness using “calculateRelate” by providing the .str file, .indivq file, and .pkla file as inputs. 

.deltas and .relat files will then be output, containing the IBD proportions (.deltas) and pairwise relatedness (.relat).

For a detailed vignette, please see the R vignette on the GitHub page. If you use the InRelate package, please cite us:

Nguyen, Khuyen, Morrissey, Christopher, and Arun Sethuraman. "InRelate - An R Package for Estimating Pairwise Genetic Relatedness in Structured Populations." (2024)

**References**:

Sethuraman, Arun. "Estimating genetic relatedness in admixed populations." G3: Genes, Genomes, Genetics 8.10 (2018): 3203-3220.

Sethuraman, Arun. "On inferring and interpreting genetic population structure-applications to conservation, and the estimation of pairwise genetic relatedness." (2013).

Evanno, Guillaume, Sebastien Regnaut, and Jérôme Goudet. "Detecting the number of clusters of individuals using the software STRUCTURE: a simulation study." Molecular ecology 14.8 (2005): 2611-2620.

Earl, Dent A., and Bridgett M. VonHoldt. "STRUCTURE HARVESTER: a website and program for visualizing STRUCTURE output and implementing the Evanno method." Conservation genetics resources 4 (2012): 359-361.

Alexander, David H., John Novembre, and Kenneth Lange. "Fast model-based estimation of ancestry in unrelated individuals." Genome research 19.9 (2009): 1655-1664.

Alexander, David H., and Kenneth Lange. "Enhancements to the ADMIXTURE algorithm for individual ancestry estimation." BMC bioinformatics 12 (2011): 1-6.

Pritchard, Jonathan K., Matthew Stephens, and Peter Donnelly. "Inference of population structure using multilocus genotype data." Genetics 155.2 (2000): 945-959.
