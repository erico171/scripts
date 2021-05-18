# scripts

Here I'll put relatively simple, miscellaneous scripts.

## runTreemix.R

### Description

runTreemix.R automates the execution and plotting of Treemix analyses for multiple migration events, from zero to the maximum allowed. It also automatically converts VCF files before running the analysis, if necessary.

### Dependencies

- Treemix: https://bitbucket.org/nygcresearch/treemix/wiki/Home

- R environment: https://cran.r-project.org/ (or at a Linux repository near you ;)

- R package 'optparse': if you don't already have it, the script will automatically try to install it once you run it; if for any reason it doesn't work, you can do it manually with the command `install.packages("optparse")` inside R.

### Running

Most of runTreemix options are actually Treemix arguments passed directly to it, with a few exceptions. Run it with `Rscript runTreemix.R -h` to see all the options. If you are going to use the convertion from VCF option, keep in mind that it is a fairly simple converter and just read the SNPs in a straight forward way, so it assumes you want all the SNPs present there in the very same order they are. So if you need to filter / reorganize your SNPs, do it before running runTreemix.


## gphocs2tree.pl

This is a very simple Perl script that generates a nexus tree with divergence times (and 95% confidence intervals) from a G-PhoCS run. It takes the information about the phylogenetic relationship between populations from the G-PhoCS control file and tau estimates from the trace file.

Usage: `perl gphocs2tree.pl [G-PhoCS control file] [G-PhoCS trace file] [root name] [burn-in (decimal)] [output file]`

G-PhoCS: https://github.com/gphocs-dev/G-PhoCS
