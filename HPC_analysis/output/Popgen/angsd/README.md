
[Chris H.'s github link](https://github.com/chollenbeck/king_scallop_popgen_2022/tree/main/R)

# Workflow using vcftools 

**About** 
 
* start with the vcf.gz file produced from angsd raw_snps.vcf.gz 
	
	* review the angsd workflow for detials on parameters used
	
	* *note*: angsd outputs a bcf file, this was converted to vcf to run below
	
* with ```vcftools``` run stepwise filtering and omissions, along the way sanity check plotted in R
	
* labeling scheme is numeric for each step

## Preparation

* if running in on the HPC, ssh into a high perfomance node 

* prompted for your password
	
```
ssh himem01
```

* load tools, we will inly require ```vcftools```


```
module load <route to vcf tools>	
```

## Step pre 

* check the mean depth for each site and each indiv, common depth threshold of 1x to be considered acceptable  - should report the mean depth for all the samples as well

* export the depth output file to gather these data in R

*..in R*
```
vcftools --vcf raw_snps.vcf.gz -- site-mean-depth --out mean_depth.txt
```

*OR*

*..in bash*
```
vcftools --vcf raw_snps.vcf.gz -- site-mean-depth --out mean_depth.txt
```

## Step 1: Filter by depth & quality

* **objective**: use ```--minDP``` and ```--minGQ``` in ```vcftools``` to filter genotypes with depth < 10 and genotype quality < 20, respectively

*..in R*
```
system("vcftools --gzvcf raw_snps.vcf.gz --out out.1 --minDP 10 --minGQ 20 --recode --recode-INFO-all")
```

*OR*

*..in bash*
```
vcftools --gzvcf raw_snps.vcf.gz --out out.1 --minDP 10 --minGQ 20 --recode --recode-INFO-all
```

## Step 2: Filter monomorphic sites

* **objective**: use ```--maf``` in ```vcftools``` to filter out sites that were made monomorphic by the previous filter

*..in R*
```
system("vcftools --vcf out.1.recode.vcf --maf 0.001 --out out.2 --recode --recode-INFO-all")
```

*OR*

*..in bash*
```
vcftools --vcf out.1.recode.vcf --maf 0.001 --out out.2 --recode --recode-INFO-all
```


## Step 3: Identify individuals with missing data 

* **objective**: use ```--max-missing``` in ```vcftools``` to remove individuals with more than 50% missing data

* think of % missingness as, "we are allowing a max missingess of XX %" therefore the smaller the number the more strict

*..in R*
```
system("vcftools --vcf out.2.recode.vcf --out out.3 --max-missing 0.5 --recode --recode-INFO-all")
```

*OR*

*..in bash*
```
vcftools --vcf out.2.recode.vcf --out out.3 --max-missing 0.5 --recode --recode-INFO-all
```


* **objective**: use ```--missing-indv``` in ```vcftools``` to output an imiss file

*..in R*
```
system("vcftools --vcf out.3.recode.vcf --out out.3 --missing-indv")
```

*OR*

*..in bash*
```
vcftools --vcf out.3.recode.vcf --out out.3 --missing-indv
```

* **objective**: import the imiss file into R session, diagnose which bam files (samples) have more than 70% missing data and omit them

* **note** - also important here to run differnt thresholds to determine a reasonable cut-off

	* First get a list of all of the duplicate individuals so that these can be saved from filtering in this step. 
	All of the duplicates have a 'b' at the end of the individual name
	
	* next plot a histogram of the data, visualize and determine cut-off, call above and below as well for sanity checking
	
	* filter individuals and call the IDs
	
*..in R*
```
library(readr)

out_3_imiss <- read_tsv("out.3.imiss") # load missingness file

dups_a <- str_subset(out_3_imiss$INDV, "b$") # nothing was output here 
dups_b <- gsub(pattern = "b", replacement = "", x = dups_a) # nothing output here as well
dups <- c(dups_a, dups_b)

qplot(out_3_imiss$F_MISS) # plot histogram, view these data determine cut-off 

miss_20 <- filter(out_3_imiss, F_MISS > 0.2, ! INDV %in% dups) %>% select(INDV) # filter and call
miss_50 <- filter(out_3_imiss, F_MISS > 0.5, ! INDV %in% dups) %>% select(INDV) # filter and call
miss_70 <- filter(out_3_imiss, F_MISS > 0.7, ! INDV %in% dups) %>% select(INDV) # filter and call

write_delim(miss_70, "remove.3.inds", col_names = FALSE) # write the indivualds to remove
```

**note** - above we only export those with 70% missigness and this removed 18 individuals; removing
 50% and 20% each removed the same 19 individuals


## Step 4: Filter individuals with missing data

* **objective**: use ```--remove``` in ```vcftools``` to remove individuals with >70% missing data

* navigate back into the bash sessionwith vcftools, create out.4 as the vcf witht he individuals removed; 
**note**, 18 individuals were removed

*..in R*
```
system("vcftools --vcf out.3.recode.vcf --out out.4 --remove remove.3.inds --recode --recode-INFO-all")
```

*OR*

*..in bash*
```
vcftools --vcf out.3.recode.vcf --out out.4 --remove remove.3.inds --recode --recode-INFO-all
```

## Step 4.2: Investigate high depth filter

* **objective** - in R, use a mean depth filter of > 300 to omit high-depth sites

*..in R*
```
library(readr)
library(tidyverse)

site_depth_4 <- read_tsv("out.4.ldepth") %>%  # Read in the site depth file,
				mutate(MEAN_DEPTH = SUM_DEPTH / 394) # note we have 394 individuals! Calc depth

qplot(site_depth_4$MEAN_DEPTH) # plot histogram to visualize mean site depth per indiv

mean_site_depth_4 <- mean(site_depth_4$MEAN_DEPTH) # what is te mean? = 16.8

to_keep_4 <- filter(site_depth_4, MEAN_DEPTH < 300) # keep those iwh < 300


to_filter_4 <- filter(site_depth_4, MEAN_DEPTH >= 300) %>%
                  select(CHROM, POS) # make a list of sites to filter
				  
nrow(to_filter_4) # there are NONE! Good.

write_delim(to_filter_4, "remove.4.sites", col_names = FALSE) # did not output, I saw there was non to omit
```

* if there were sites to remove and you output 'remove.4.sites' (**there were none!**) 
then remove sites using ```vcftools``` in the following

* shown below output as out.5, however we did not run this!

*..in R*
```
system("vcftools --vcf out.4.vcf --out out.5 --exclude-positions remove.4.sites --recode --recode-INFO-all")
```

*OR*

*..in bash*
```
vcftools --vcf out.4.vcf --out out.5 --exclude-positions remove.4.sites --recode --recode-INFO-all
```

## Step 5: Investigate sites with missing data

* **objective**: use ```--max-missing``` in ```vcftools``` to remove sites with more than 75% missing data

*..in R*
```
system("vcftools --vcf out.4.recode.vcf --out out.5 --max-missing 0.75 --recode --recode-INFO-all")

```

*OR*

*..in bash*
```
vcftools --vcf out.4.recode.vcf --out out.5 --max-missing 0.75 --recode --recode-INFO-all

```

* **outputs**: 
	- OUTPUT: After filtering, kept 3640 out of a possible 3897 Sites
	- OUTPUT: After filtering, kept 394 out of 394 Individuals


## Step 5.2: Build file in R to omit based on missingness 

* **objective**: use ```--missing-indv``` in ```vcftools``` to calculate individual missingness - this is done again since the pool of sites is lower

* **note** -  this creates imiss and .kig files, does not overwrite the 5.recode

*..in R*
```
system("vcftools --vcf out.5.recode.vcf --out out.5 --missing-indv") 
```

*OR*

*..in bash*
```
vcftools --vcf out.5.recode.vcf --out out.5 --missing-indv
```

* **outputs**: 
	- OUTPUT:After filtering, kept 394 out of 394 Individuals
	- OUTPUT:Outputting Individual Missingness
	- OUTPUT: After filtering, kept 3640 out of a possible 3640 Sites
       
	   
* **objective** - in R, open the imiss file to call and output individuals with missingness 

*..in R*
```
library(readr)
library(tidyverse)

out_5_imiss <- read_tsv("out.5.imiss") %>%
                  extract(INDV, "pop", "(\\w+)_", remove = FALSE)

out_5_imiss %>% count(pop)

ggplot(out_5_imiss, aes(x = F_MISS, fill = pop)) + geom_histogram() #Plot a quick histogram of the data from the plot we see that 60% is the target filter for omission

miss_60 <- filter(out_5_imiss, F_MISS > 0.60) %>%  select(INDV) # Select individuals with more than 60% missing data

write_delim(miss_60, "remove.5.inds", col_names = FALSE) # write out
```

* **outputs**: 
	- remove.5.inds listed three individuals to omit!

## Step 6: Filter individuals with missing data

* **objective**: use ```--remove``` in ```vcftools``` to remove the individuals with more than 60% missing genotype

*..in R*
```
system("vcftools --vcf out.5.recode.vcf --out out.6 --remove remove.5.inds --recode --recode-INFO-all")
```

*OR*

*..in bash*
```
vcftools --vcf out.5.recode.vcf --out out.6 --remove remove.5.inds --recode --recode-INFO-all
```

* **outputs**: 
	- out.6 now at 391 indivuals from 394! 


## Step 6.2: Check for duplicates and relatedness filter

* **objective**: use ```--relatedness``` in ```vcftools``` to check for duplicate individuals using relatedness


*..in R*
```
system("vcftools --vcf out.6.recode.vcf --relatedness --out out.6")
```

*OR*

*..in bash*
```
vcftools --vcf out.6.recode.vcf --relatedness --out out.6
```

* hop into R and use a AJK > 0.7 (abnormally high relatedness) to check for relatedness, output file  of ids to omit

*..in R*
```
library(readr)
library(tidyverse)

out_6_relate <- read_tsv("out.6.relatedness") # load the data

relatives <- filter(out_6_relate, INDV1 != INDV2) %>% 
              filter(RELATEDNESS_AJK > 0.7) # # there were no relatived following this criteria of > 0.7 on this scale!
```

## Step 7: Convert to vcf.gz and rename


* **objective** - convert to .gz and rename

*..in R*
```
system("bcftools sort our.6.recode.vcf  -Oz -o all_final.vcf.gz")
```

*OR*

*..in bash*
```
bcftools sort our.6.recode.vcf  -Oz -o all_final.vcf.gz
```

# Summary: 
	- 18 individuals removed based on site missingness >70%
	- 3 individuals removed based on based on >60% missing genotype
	- no individulas shown to be abnormally related AJK > 0.8
