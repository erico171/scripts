#!/usr/bin/env Rscript

suppressWarnings(
  if(!require(optparse,quietly = T)){
    cat(" Dependency 'optparse' not found, trying to install...\n")
    install.packages("optparse")
    library(optparse)
  }
)


optionsList = list(
  make_option(opt_str = c("-i","--input"),type = "character",default = NULL,help = "Treemix input or VCF to be converted (see -y and -a options)"),
  make_option(opt_str = c("-y","--type"),type = "character",default = "treemix",help = "Type of the input file.\n\t\tCurrently supported types are 'treemix' (default) and 'vcf' (Variant Call Format)"),
  make_option(opt_str = c("-r","--root"),type = "character",default = NULL,help = "Outgroup pop name"),
  make_option(opt_str = c("-m","--migs"),type = "character",default = NULL,help = "Maximum migration events to be added"),
  make_option(opt_str = c("-p","--poporder"),type = "character",default = NULL,help = "File with poporder to plot residual fit"),
  make_option(opt_str = c("-o","--basename"),type = "character",default = "out",help = "Base name for the outputs (default: out)"),
  make_option(opt_str = c("-x","--plusxaxis"),type = "character",default = "0.1",help = "Extra \"space\" for the X axis on plot_tree (default: 0.1)"),
  make_option(opt_str = c("-t","--treemix"),type = "character",default = NULL,help = "Full path for the treemix binary (default: 'treemix' from the environment)"),
  make_option(opt_str = c("-f","--plotf"),type = "character",default = NULL,help = "Full path for the \"plotting_funcs.R\" script (default: [treemix_path]/plotting_funcs.R)"),
  make_option(opt_str = c("-v","--vcf2treemix"),type = "logical",action = "store_true",default = F,help = "Just converts VCF to Treemix input and quit"),
  make_option(opt_str = c("-b","--bootstrap"),type = "logical",action = "store_true",default = F,help = "Generate one bootstrap (default: disabled)"),
  make_option(opt_str = c("-s","--stderror"),type = "logical",action = "store_true",default = F,help = "Calculate standard errors of migration weights (default: disabled)"),
  make_option(opt_str = c("-a","--assignment"),type = "character",default = NULL,help = "Needed only with the '-y vcf' option.\n\t\tA text file containing the names of the populations to which each of the samples belong,\n\t\tin the same order as they appear in the VCF columns"),
  make_option(opt_str = c("-k","--kmers"),type = "integer",default = 1,help = "Number of SNPs per block for estimation of covariance matrix (default: 1)"),
  make_option(opt_str = c("-l","--global"),type = "logical",action = "store_true",default = F,help = "Treemix's -global: do a round of global rearrangements after adding all populations (default: disabled)"),
  make_option(opt_str = c("-c","--tf"),type = "character",default = NULL,help = "Treemix's -tf [file name]: read the tree topology from a file, rather than estimating it"),
  make_option(opt_str = c("-g","--graph"),type = "character",default = NULL,help = "Treemix's -g [vertices file name] [edges file name]: read the graph from a previous TreeMix run"),
  make_option(opt_str = c("-e","--micro"),type = "logical",action = "store_true",default = F,help = "Treemix's -micro: microsatellite data"),
  make_option(opt_str = c("-w","--cormig"),type = "character",default = NULL,help = "Treemix's -cor_mig [file]: list of known migration event to include (also use -climb)"),
  make_option(opt_str = c("-n","--noss"),type = "logical",action = "store_true",default = F,help = "Treemix's -noss: Turn off sample size correction"),
  make_option(opt_str = c("-d","--seed"),type = "integer",default = NULL,help = "Treemix's -seed [int]: set the seed for random number generation"),
  make_option(opt_str = c("-j","--nwarn"),type = "integer",default = NULL,help = "Treemix's -n_warn [int]: display first N warnings")
  
)

  
opt = (OptionParser(option_list = optionsList,add_help_option = T,))
arguments = parse_args(opt)

if(is.null(arguments$input)){
  print_help(opt)
  stop("Inform the path for Treemix input (should be a GZIP file)")
}
if(!arguments$vcf2treemix){
  if(is.null(arguments$root)){
    print_help(opt)
    stop("Inform the name of the outgroup population")
  }
  if(is.null(arguments$migs)){
    print_help(opt)
    stop("Inform the maximum number of migrations")
  }
  if(is.null(arguments$poporder)){
    print_help(opt)
    stop("Inform the name of the poporder file")
  }
  if(is.null(arguments$treemix)){
    arguments$treemix = Sys.readlink(Sys.which("treemix"))
  }
  if(is.null(arguments$plotf)){
    arguments$plotf = paste(dirname(arguments$treemix),"plotting_funcs.R",sep = "/")
  }  
}

if(arguments$type == "vcf"){
  if(is.null(arguments$assignment)){
    print_help(opt)
    stop("You must provide the name of the assingment file with the -a option")
  }else{
    
    cat(" Converting VCF input to Treemix...\n")
    
    pops = as.character(read.table(file = arguments$assignment)[,1])
    myVCF = as.matrix(read.table(file = arguments$input,comment.char = "#"))
    
    firstCol = grep(pattern = "^[01]/[01]",x = myVCF[1,],perl = T)[1]
    myVCF = myVCF[,firstCol:ncol(myVCF)]
    
    nSamples = ncol(myVCF)
    nPops = length(pops)
    if(nSamples != nPops){
      stop(paste("ERROR: The population assignment file provided contains",nPops,"population names, but your VCF contains",nSamples,"samples."))
    }
    
    pops.unique = unique(pops)
    
    matriz = matrix(data = NA,nrow = nrow(myVCF),ncol = length(pops.unique),dimnames = list(NULL,pops.unique))
    
    for(snp in 1:nrow(myVCF)){
      alelos1 = sub(pattern = "^([01]).+",replacement = "\\1",x = myVCF[snp,])
      alelos2 = sub(pattern = "^[01]/([01]).+",replacement = "\\1",x = myVCF[snp,])
      for(pop in pops.unique){
        indexes = grep(pattern = paste0("^",pop,"$"),x = pops)
        count0 = length(grep(pattern = "0",x = alelos1[indexes],fixed = T)) + length(grep(pattern = "0",x = alelos2[indexes],fixed = T))
        count1 = length(grep(pattern = "1",x = alelos1[indexes],fixed = T)) + length(grep(pattern = "1",x = alelos2[indexes],fixed = T))
        
        matriz[snp,pop] = paste(count0,count1,sep = ",")
      }
    }
    
    newInputName = paste0(arguments$input,".converted")
    write.table(x = matriz,file = newInputName,quote = F,row.names = F)
    R.utils::gzip(filename = newInputName,overwrite=T,remove=F)
    arguments$input = paste0(newInputName,".gz")
    
    if(arguments$vcf2treemix){
      stop("Input VCF file converted to Treemix. Stoping..")
    }
    
  }
}else if(arguments$type == "treemix"){
  if(!R.utils::isGzipped(arguments$input)){
    R.utils::gzip(filename = arguments$input,overwrite=T,remove=F)
    arguments$input = paste0(arguments$input,".gz")
  }
}

source(arguments$plotf)

pdf(file = paste0(arguments$basename,".pdf"))
lliks = numeric()

for(m in 0:arguments$migs){

  nome = paste0(arguments$basename,"_m",m)
  comando = paste(arguments$treemix,"-i",arguments$input,"-o",nome,"-root",arguments$root,"-m",m,"-k",arguments$kmers)
  if(arguments$bootstrap){
    cat(" IT'S A BOOTSTRAP!\n")
    comando = paste(comando,"-bootstrap")
  }
  if(arguments$stderror){
    comando = paste(comando,"-se")
  }
  if(arguments$global){
    comando = paste(comando,"-global")
  }
  if(!is.null(arguments$tf)){
    comando = paste(comando,"-tf",arguments$tf)
  }
  if(!is.null(arguments$graph)){
    comando = paste(comando,"-g",arguments$graph)
  }
  if(arguments$micro){
    comando = paste(comando,"-micro")
  }
  if(!is.null(arguments$cormig)){
    comando = paste(comando,"-cor_mig",arguments$cormig)
  }
  if(arguments$noss){
    comando = paste(comando,"-noss")
  }
  if(!is.null(arguments$seed)){
    comando = paste(comando,"-seed",arguments$seed)
  }
  if(!is.null(arguments$nwarn)){
    comando = paste(comando,"-n_warn",arguments$nwarn)
  }
  # print(comando)
  system(command = comando)
  plot_tree(stem = nome,plus = as.numeric(arguments$plusxaxis)); title(main = nome)
  plot_resid(stem = nome,pop_order = arguments$poporder); title(main = nome)
  llikFile = paste0(nome,".llik")
  llikDF = readLines(con = file(llikFile))
  llik = sub(pattern = ".+: ([\\d\\.]+) $",replacement = "\\1",x = llikDF[2],perl = T)
  lliks[as.character(m)] = llik

}

plot(x = as.numeric(labels(lliks)),y = lliks,xlab = "# migration events",ylab = "Log(likelihood)",yaxt="n")
minVal = as.numeric(min(lliks))
maxVal = as.numeric(max(lliks))
axis(side = 2,at = seq(from=floor(minVal),to=ceiling(maxVal),by=round((maxVal - minVal)/10)))

# print(max(lliks))
dev.off()