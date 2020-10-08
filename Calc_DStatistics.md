# Performing various ABBA-BABA tests to the infer the extent of introgression between *Zosterops tenuirostris* and *Zosterops lateralis* on Norfolk Island.
A. Sendell-Price, Oct 2020.

## STEP 1: Calculate genome-wide estimate of introgression using DSuite
In our first assessment of hybridization across the genome we will calculate a single D-statistic for the genome. This will give us an idea of the overall extent of introgression between the two species. For this we will use the program [DSuite](https://github.com/millanek/Dsuite) which allows for the fast calculation of the D-statistic from SNP data in VCF format.

First, set path to DSuite and specify the VCF file
```
DSUITE=~/Dropbox/DPhil/BIN/Dsuite/Build/Dsuite
VCF=../VCFs/ZFified_Norfolk_Hybridization_indv73_pos7019400.vcf.gz
```
To run the analysis we will need to create a tab delimited file ("Pops.txt") that assigns samples to populations. That file looks like this, where "Outgroup" is used to assigns individuals to the outgroup (in this case we use the Reunion grey white-eye a.k.a *Z.borbonicus*) and remaining samples are assigned to their corresponding populations:

```
head -n 3 Pops.txt
```
```
15-179	Outgroup
PN1	NNZ_Zlat
PN10	NNZ_Zlat
PN11	NNZ_Zlat
```

Calculate the D-statistic (ABBA/BABA) and f4-ratio (f_G) statistics for all trios of species/pops in the dataset (the outgroup being fixed). The results are as definded in Patterson et al. 2012 (equivalent to Durand et al. 2011 when the Outgroup is fixed for the ancestral allele). Run Dsuite using Dtrios to calculate a single genome-wide estimate of Patterson’s D, aka the ABBA-BABA statistic:

```
$DSUITE Dtrios $VCF Pops.txt
```

This will take a little while to run, following which DStuite will write several files. The file we are interested in is "Pops__Dmin.txt", which contains the lowest estimate of D and a p-value estimated via DSuite's jacknifing procedure. The absence of introgression can be rejected if the p-value is below the significance level. In our case significant genome-wide levels of divergence HAVE been detected. As this value is positive, we can infer that introgression from *Z.tenuirostris* into *Z.lateralis* has taken place.

```
cat Pops__Dmin.txt
```

```
P1  P2  P3  Dstatistic  p-value f_G
NNZ_Zlat  NI_Zlat NI_Zten 0.0397274 5.66069e-12 0.0181002
```


## STEP 2: Determining how patterns of introgression vary across the genome

So far we have determined that at the genome level there is significant levels of introgression between *Z.tenuirostris* and *Z.lateralis* on Norfolk Island, the next step will be to calculate the D-statistic in non-overlapping windows so that we can assess how patterns of introgression vary across the genome.

This could be conducted using DSuite's "Dinvestigate" command, however I think this method has a severe limitation in that "Dinvestigate" uses the number of SNPs rather than SNP coordinates when determining window ranges. As a result windowed D-statistics calculated with DSuite are calculated for windows of varying size - which is not ideal.

So, instead we will make use of Simon H. Martin's "genomics_general" script-set which includes several useful python scripts for running windowed ABBA-BABA tests. See: https://github.com/simonhmartin/genomics_general.

Using bioconda create two python environments, one running python2 and the other running python 3. We will need to switch between these as Simon's scripts assume different python versions.

```
conda create --name py2 python=2.7
conda create --name py3 python=3.8
```

Activate python2 environment
```
conda activate py2
```

Use parseVCF.py to convert VCF file into ".geno" format (format required for downstream analyses).

Set path to parseVCF.py and output file:
```
parseVCF=~/Dropbox/DPhil/BIN/genomics_general/VCF_processing/parseVCF.py
GENO=../VCFs/ZFified_Norfolk_Hybridization_indv73_pos7019400.geno.gz
```

Convert vcf file:
```
python $parseVCF -i $VCF --skipIndels | bgzip > $GENO
```

The outputted file will be formatted like this:
CHROM	POS	    15-179	N101	N102	N103
1      107144	C/C	    C/C   N/N	  C/C
1	    107145	G/G	    G/G	  N/N   G/G

Missing data is denoted as N, and phased and unphased genotypes are shown conventionally with | and /. For our data all sites are unphased.

Calculate D-statistic using ABBABABAwindows.py
First set path to python script, and define necessary variables:
```
ABBABABAwindows=~/Dropbox/DPhil/BIN/genomics_general/ABBABABAwindows.py
WINDOW_SIZE=100000
MIN_SNPs=100  # <- Previous studies suggest 100 SNPs is the absolute minumum
P1=NNZ_Zlat
P2=NI_Zlat
P3=NI_Zten
```

Run ABBABABAwindows.py
```
python $ABBABABAwindows \
-g $GENO \
-f phased \
-o D-stats_windowsize${WINDOW_SIZE}_min${MIN_SNPs}.csv \
-w $WINDOW_SIZE \
-m $MIN_SNPs \
-P1 $P1 -P2 $P2 -P3 $P3 -O Outgroup \
--minData 0.5 --popsFile Pops.txt
```

################################################################################
STEP 3: Calculate DFS
################################################################################

By calculating genome-wide and windowed D-statistics we now know the extent to which the Z.tenuirostris and Z.lateralis genomes are introgressed and how patterns of introgression vary across the genome. However, the D-statistic itself is not overly informative of the timing or direction of this introgression. So, now we will calculate D frequency spectrum (DFS) - an extension of the D-statistic in which D is partitioned according to the frequencies of derived allles. The DFS is strongly altered by different ages and directions of introgression and is also sensitive to demographic events such as population bottlenecks (all explained by Simon Martin here: https://doi.org/10.1093/molbev/msaa239).

This is a three step process. First run freq.py to compute allele frequencies at each site in each population:

Set path to freq.py:
```
FREQ=~/Dropbox/DPhil/BIN/genomics_general/freq.py
```

Compute allele frequencies:
```
python $FREQ \
-g $GENO \
-p $P1 -p $P2 -p $P3 -p Outgroup \
--popsFile Pops.txt | gzip > basecounts.tsv.gz
```

Next run sfs.py, which produces Site Frequency Spectra. By including the --outgroup flag, we are telling it to use the outgroup population to polarise allele frequencies for the ingroups. We then tell it to output just one 3D SFS, representing populations P1 (NNZ Zlat), P2 (NI Zlat) and P3 (NI Zten).

Activate python3 environment:
```
conda activate py3
```

Set path to sfs.py:
SFS=../genomics_general-0.4/sfs.py

Compute SFS, subsample the data down to 22 haplotypes for P1 and P2 and 14 for P3. Note: subsampling is in haploid.

```
python $SFS \
-i basecounts.tsv.gz \
--inputType baseCounts \
--outgroup Outgroup \
--FSpops NNZ_Zlat NI_Zlat NI_Zten \
--subsample 22 22 14 \
--pref NorfolkIsHybrid. \
--suff .subsample22_22_14.sfs
```

Finally, plot DFS using R Script plot_DFS_from_SFS.R, Which looks like this:

```{r}
dfs_data <- get.DFS(base_counts=FS[,-4], #base counts are the first three columns (i.e everything minus column 4)
                    site_counts=FS[,4], # site counts are column
                    Ns = c(22,22,14)) #Ns provide the haploid sample sizes of each population (1 and 2 must always be equal)
```

### plot

pdf("NorfolkIsHybrid.NNZ_Zlat_NI_Zlat_NI_Zten.subsample22_22_14.pdf")
par(mar=c(1,4,1,1))
plotDFS(dfs_data$DFS, dfs_data$weights, method="lines", col_D="red", no_xlab=F)
dev.off()

### code for exporting a table of plotted values
write.table(data.frame(D=round(dfs_data$DFS,4),
                        weight=round(dfs_data$weights,4)),
             file="NorfolkIsHybrid.NNZ_Zlat_NI_Zlat_NI_Zten.subsample22_22_14.csv",
             sep=",",row.names=FALSE,quote=FALSE)




#STEP 3: Calculate windowed divergence/diversity stats
python ../genomics_general-0.4/popgenWindows.py \
-g ../VCFs/ZFified_Norfolk_Hybridization_indv73_pos7019400.geno.gz \
-f phased \
-w 10000 \
-m 100 \
-o Fst_dxy_pi_10kb.csv \
-p NNZ_Zlat -p NI_Zlat -p NI_Zten -p Zbor \
--popsFile Pops.txt
