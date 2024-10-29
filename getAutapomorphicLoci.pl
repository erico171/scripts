#!/usr/bin/perl

use strict;
use Bio::SeqIO;
use List::Util qw(shuffle);

$| = 1;

my $usage = " Usage: getAutopomorphic.pl [fastas_directory] [pop_regex OR lookup_file] [min_samples_per_pop] [shuffle_reps]\n";

die $usage unless $#ARGV == 3;

my $dir = $ARGV[0];
my $rePops = $ARGV[1];
my $minSamples = $ARGV[2];
my $shuffle = $ARGV[3];

my $regexMode = 1;
my %lookup;
my @pops = ();
if(-e $rePops){
    $regexMode = 0;
    open(IF,$rePops) or die " ERROR opening lookup file $rePops\n";
    while(my $line = <IF>){
        chomp $line;
        my @splitLine = split(/\s/,$line);
        $lookup{$splitLine[0]} = $splitLine[1];
        push(@pops,$splitLine[1]);
    }
    close IF;
}

my $outFile = "autapomorphies.tsv";
# foreach my $l (keys %lookup){
#     print "$l $lookup{$l}\n";
# }

#1 ler todos os fastas para pegar todas as populacoes presentes, anotando o numero de populacoes em cada um para depois descartar os que nao possuem todas

opendir(DIR,$dir) or die "\n Could not open directory: $dir\n$usage";
my @files = readdir DIR;
closedir DIR;
my @fastas = ();
foreach my $f (@files){
	if($f =~ /^.+\.fa[sta]*$/){
		push(@fastas,$f);
	}
}

if($regexMode){

    my @goodFastas = ();
    foreach my $f (@fastas){
        print STDERR "\r Reading fastas... $f                                       ";
        my $seqio_object = Bio::SeqIO->new(-file => "$dir/$f", -format => "fasta", -alphabet => "dna");
        my @thisPops = ();
        while(my $seq_object = $seqio_object->next_seq){
            my $name = $seq_object->display_id();
            if($regexMode){
                if($name =~ /$rePops/){
                    push(@thisPops,"$1$2");
                }
                else{
                    die " ERROR finding pop name in seq $name, file $f. Please check your regex.\n";
                }
            }
            else{
                push(@thisPops,$lookup{$name});
    #             print "$lookup{$name}\n";
            }
        }
    #     @pops = getUniques(@thisPops);
        my $status = 1;
        my @thisPopsCount = getUniquesCount(@thisPops);
        foreach my $t (@thisPopsCount){
            if($t < $minSamples){
                $status = 0;
                last;
            }
        }
        push(@pops,getUniques(@thisPops));
        push(@goodFastas,$f) if($status);
    }
    print STDERR "\n";
    @fastas = @goodFastas;
}

@pops = getUniques(@pops);


if($shuffle > 0){
    
    my @nAutapomorphicLoci = ();
    for(my $i=1;$i<=$shuffle;$i++){
        $nAutapomorphicLoci[$i] = 0;
    }
    foreach my $f (@fastas){
        my $found = 0;
        my $soma = 0;
        foreach my $a (@nAutapomorphicLoci){
            $soma += $a;
        }
        my $seqio_object = Bio::SeqIO->new(-file => "$dir/$f", -format => "fasta", -alphabet => "dna");
        my @thisPops = ();
        my @names = ();
        my @seqs = ();
        while(my $seq_object = $seqio_object->next_seq){
            my $name = $seq_object->display_id();
            if($regexMode){
                $name =~ /$rePops/;
                push(@thisPops,"$1$2");
                
            }
            else{
                push(@thisPops,$lookup{$name});
            }
            push(@names,$name);
            push(@seqs,uc($seq_object->seq()));
                
        }
        
        for(my $s=0;$s<length($seqs[0]);$s++){
        
            my @basesSitio = ();
            my %countValidBases;
            foreach my $p (@pops){
                $countValidBases{$p} = 0;
            }
            for(my $n=0;$n<=$#names;$n++){
                my $base = substr($seqs[$n],$s,1);
                if($base =~ /[ATCG]/){
                    push(@basesSitio,$base);
                    $countValidBases{$thisPops[$n]}++;
                }
            }
            my @uniqueBasesSitio = getUniques(@basesSitio);
            if($#uniqueBasesSitio > 0){
            
                for(my $i=1;$i<=$shuffle;$i++){
                    
                    print STDERR "\r File: $f    Site: $s    Shuffle: $i    Sum: $soma    ";
            
                
                
                    my $goodSite = 1;
                    foreach my $v (keys %countValidBases){
                        if($countValidBases{$v} < $minSamples){
                            $goodSite = 0;
                            last;
                        }
                    }
                    if($goodSite){
                        
                        @thisPops = shuffle(@thisPops);
                        foreach my $base (@uniqueBasesSitio){
                            my @thisSitePops = ();
                            for(my $n=0;$n<=$#names;$n++){
                                if($basesSitio[$n] eq $base){
                                    push(@thisSitePops,$thisPops[$n]);
                                }
                            }
                            my @uniqueThisSitePops = getUniques(@thisSitePops);
                            if($#uniqueThisSitePops==0){
                                my $pop = $uniqueThisSitePops[0];
                                my @thisPopBases = ();
                                for(my $n=0;$n<=$#names;$n++){
                                    if($thisPops[$n] eq $pop and $basesSitio[$n] =~ /[ATCG]/){
                                        push(@thisPopBases,$basesSitio[$n]);
                                    }
                                }
                                my @uniqueThisPopBases = getUniques(@thisPopBases);
                                if($#uniqueThisPopBases == 0){
                                    $nAutapomorphicLoci[$i]++;
                                    $found = 1;
                                }
                            }
                        }
                    }
                }
            }
            if($found){
                last;
            }
        }
    }
    print STDERR "\n";
    print "@nAutapomorphicLoci\n";
}
else{
    #2 criar um hash de loci, que contera uma sequencia de numeros, um para cada pop, indicando quantas autapomorfias (SNPs) aquele locus possui para cada pop
            
#     my %autapomorphies;

    #3 abrir cada fasta e, para cada sitio, verificar se é variavel e se contem o minimo de amostras para cada pop.
    #4 em caso positivo, e caso haja pops com BASES FIXADAS E EXCLUSIVAS, somar uma autapomorfia para elas nesse locus.

    print "Locus";
    foreach my $p (@pops){
        print "\tcount\_$p\tpos\_$p";
    }
    print "\n";

    foreach my $f (@fastas){
        print STDERR "\r Searching for autapomorphies... $f           ";
        my $seqio_object = Bio::SeqIO->new(-file => "$dir/$f", -format => "fasta", -alphabet => "dna");
        my @thisPops = ();
        my @names = ();
        my @seqs = ();
        my $printFlag = 0;
        while(my $seq_object = $seqio_object->next_seq){
            my $name = $seq_object->display_id();
            if($regexMode){
                $name =~ /$rePops/;
                push(@thisPops,"$1$2");
                
            }
            else{
                push(@thisPops,$lookup{$name});
            }
            push(@names,$name);
            push(@seqs,uc($seq_object->seq()));
                
        }
        my %countAutopomorphies;
        my %autapomorphicSites;
        foreach my $p (@pops){
            $countAutopomorphies{$p} = 0;
            $autapomorphicSites{$p} = "";
        }
        for(my $s=0;$s<length($seqs[0]);$s++){
            my @basesSitio = ();
            my @popsSitio = ();
            my %countValidBases;
            foreach my $p (@pops){
                $countValidBases{$p} = 0;
            }
            for(my $n=0;$n<=$#names;$n++){
                my $base = substr($seqs[$n],$s,1);
                if($base =~ /[ATCG]/){
                    push(@basesSitio,$base);
                    push(@popsSitio,$thisPops[$n]);
                    $countValidBases{$thisPops[$n]}++;
                }
            }
            my $goodSite = 1;
            foreach my $v (keys %countValidBases){
                if($countValidBases{$v} < $minSamples){
                    $goodSite = 0;
                    last;
                }
            }
            if($goodSite){
                my @uniqueBasesSitio = getUniques(@basesSitio);
                if($#uniqueBasesSitio > 0){
                    foreach my $base (@uniqueBasesSitio){
                        my @thisSitePops = ();
                        for(my $n=0;$n<=$#names;$n++){
                            if($basesSitio[$n] eq $base){
                                push(@thisSitePops,$thisPops[$n]);
                            }
                        }
                        my @uniqueThisSitePops = getUniques(@thisSitePops);
                        if($#uniqueThisSitePops==0){
                            my $pop = $uniqueThisSitePops[0];
                            my @thisPopBases = ();
                            for(my $n=0;$n<=$#names;$n++){
                                if($thisPops[$n] eq $pop and $basesSitio[$n] =~ /[ATCG]/){
                                    push(@thisPopBases,$basesSitio[$n]);
                                }
                            }
                            my @uniqueThisPopBases = getUniques(@thisPopBases);
                            if($#uniqueThisPopBases == 0){
                                $countAutopomorphies{$pop}++;
                                $autapomorphicSites{$pop} .= ($s+1)."_";
                                $printFlag = 1;
                            }
                        }
                    }
                }
            }
        }
        
        if($printFlag){
            print "$f";
            foreach my $p (@pops){
        #       $countsText .= "$countAutopomorphies{$p}\t";
                chop $autapomorphicSites{$p};
                print "\t$countAutopomorphies{$p}\t$autapomorphicSites{$p}";
            }
            print "\n";
        }
    }

    print STDERR "\n";
}

sub getUniquesCount { #(@array)

    my @array = sort(@_);
    my @counts = (1);
    for(my $i=1;$i<=$#array;$i++){
        if($array[$i] ne $array[$i-1]){
            push(@counts,1);
        }
        else{
            $counts[$#counts]++;
        }
    }
    
    return @counts;
}

sub getUniques { #(@array)

    my @array = sort(@_);
    my @uniques = ($array[0]);
    for(my $i=1;$i<=$#array;$i++){
        if($array[$i] ne $array[$i-1]){
            push(@uniques,$array[$i]);
        }
    }
    
    return @uniques;
}
