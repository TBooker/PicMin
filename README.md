PicMin
======

In this repository you'll find the an R package for applying PicMin, a method for identifying repeated patterns of evidence for adaptation from genome scans.

If you are interested in the method, take a look at the vignette detailing an analysis of Arabidopsis data from Bohutínská et al (2021). It's an RMarkdown script and walks the user through the process of applying PicMin to data. Several of the steps are specific to the Arabidopsis data, but this is 

The method is detailed in our manuscript available on bioRXiv: 
[https://www.biorxiv.org/content/10.1101/2022.03.24.485690v1](https://www.biorxiv.org/content/10.1101/2022.03.24.485690v1).

________________________
### Analysis windows or SNPs

We designed PicMin with analyses that compare the results of window-based genome scans in many species.  Window-based analyses lend themselves more readily to comparisons across lineages than individual genetic markers (e.g. using genes as windows, where orthology can be more readily identified). Depending on the phylogenetic distance between one's samples, individual SNPs may or may not be present in all lineages. Additionally, restricting analyses to individual SNPs may cause some signal of repeated adaptation to be lost if, for example, species A exhibited evidence at a particular SNP, but species B only showed evidence at a tightly linked SNP. In this case, one may want to conclude that there has been repeated adaptation acting on the same genomic element in the different lineages, but this would go undetected if one were looking only at individual SNPs. There may be ways to do this, but I have not thought them through. Get in touch if you would like to discuss this though. 
 
The use of analysis windows entails a difficult decision - how wide or narrow should analysis windows be? As in any analysis of a single lineage, comparing the results of genome scans across lineages will be sensitive to certain analysis choices. If windows are narrow relative to the decay of linkage disequilibrium, multiple windows may be associated with a given gene, necessitating stringent multiple comparisons correction. Alternatively, if windows are too wide the signal of selection or adaptation within a window may be diluted by linked sites. Ultimately the choice of analysis window size should be informed by local recombination rates or the extent of linkage disequilibrium in genome, the number of genes or genomic units that can be compared across lineages and considerations of multiple comparison corrections.


