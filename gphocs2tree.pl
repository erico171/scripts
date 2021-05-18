#!/usr/bin/perl

use strict;

validateArgs();

my $ctrlFile = $ARGV[0];
my $traceFile = $ARGV[1];
my $rootName = $ARGV[2];
my $burnin = $ARGV[3];
my $outfile = $ARGV[4];

my $nLines = `wc -l $traceFile`;
$nLines--;
my $skip = int($nLines * $burnin);

open(TRACE,$traceFile) or die;
my $header = <TRACE>;
$header =~ s/ //g;
my @header = split(/\t/,$header);
my @tauCols = ();
for(my $i=0;$i<=$#header;$i++){
	if($header[$i] =~ /^tau_(.+)/){
		$header[$i] = $1;
		push(@tauCols,$i);
	}
}

my %taus;
my %medias;
for(my $i=0;$i<=$#tauCols;$i++){
	$medias{$header[$tauCols[$i]]} = 0;
	@{$taus{$header[$tauCols[$i]]}} = ();
}

my $count = 0;
while(my $line = <TRACE>){
	if($skip > 0){
		$skip--;
	}
	else{
		$count++;
		$line =~ s/ //g;
		my @splitLine = split(/\t/,$line);
		for(my $i=0;$i<=$#tauCols;$i++){
			$medias{$header[$tauCols[$i]]} += $splitLine[$tauCols[$i]];
			push(@{$taus{$header[$tauCols[$i]]}},$splitLine[$tauCols[$i]]+0);
		}
	}
}
close TRACE;


my %ic95;
foreach my $tau (keys %taus){
	my @sorted = sort { $a <=> $b } @{$taus{$tau}};

	my $tam = scalar(@sorted);
	my $inf = int(0.025*$tam);
	my $sup = int(0.975*$tam);
	
	$ic95{$tau} = "[&CI95={$sorted[$inf],$sorted[$sup]}]";
}

foreach my $media (keys %medias){
	$medias{$media} /= $count;
}

open(CTRL,$ctrlFile) or die;

my $ancFlag = 0;
my $curFlag = 0;

my %ancestors;
my @currents = ();
my $lastAnc = "";

while(my $line = <CTRL>){
	
	$curFlag = 1 if($line =~ /CURRENT-POPS-START/);
	$curFlag = 0 if($line =~ /CURRENT-POPS-END/);
	$ancFlag = 1 if($line =~ /ANCESTRAL-POPS-START/);
	$ancFlag = 0 if($line =~ /ANCESTRAL-POPS-END/);
	
	if($line =~ /name\s+(\S+)/){
		push(@currents,$1) if($curFlag);
		if($ancFlag){
			$lastAnc = $1;
		}
	}
	
	if($line =~ /children\s+(\S+)\s+(\S+)/ and $ancFlag){
		$ancestors{$lastAnc} = [($1,$2)];
	}
}
close CTRL;

my $newick = "";
buildTree($rootName);

open(OF,">$outfile") or die;
print OF "#NEXUS\nBegin trees;\ntree TREE1 = $newick\nEnd;\n";
close OF;

sub buildTree { #(pop)
	my $anc = shift;
	my $c0 = ${$ancestors{$anc}}[0];
	my $c1 = ${$ancestors{$anc}}[1];
	if($anc eq $rootName){
		$newick = "($c0:$medias{$anc},$c1:$medias{$anc})$ic95{$rootName};";
		
		if(isAncestor($c0)){
			$newick =~ s/$c0:$medias{$anc}/$c0$ic95{$c0}:$medias{$anc}/;
		}
		if(isAncestor($c1)){
			$newick =~ s/$c1:$medias{$anc}/$c1$ic95{$c1}:$medias{$anc}/;
		}
	}
	else{
		die "regex error: $newick\n" unless $newick =~ /[\(,]$anc\[[^\]]+\]:([\d\.\-eE]+)/;
		my $oldtime = $1;
		my $newtime = $oldtime - $medias{$anc};
		$newick =~ s/([\(,])$anc(\[[^\]]+\]):$oldtime/$1($c0:$medias{$anc},$c1:$medias{$anc})$2:$newtime/;
		if(isAncestor($c0)){
			$newick =~ s/$c0:$medias{$anc}/$c0$ic95{$c0}:$medias{$anc}/;
		}
		if(isAncestor($c1)){
			$newick =~ s/$c1:$medias{$anc}/$c1$ic95{$c1}:$medias{$anc}/;
		}
	}
	buildTree($c0) if(isAncestor($c0));
	buildTree($c1) if(isAncestor($c1));
}

sub isAncestor {
	my $pop = shift;
	foreach my $cur (@currents){
		if($pop eq $cur){
			return 0;
		}
	}
	return 1;
}

sub validateArgs {

	my $usage = " Usage: gphocs2tree.pl [G-PhoCS control file] [G-PhoCS trace file] [root name] [burn-in (decimal)] [output file]\n";
	
	die $usage if($#ARGV != 4);
	
	open(CTRL,$ARGV[0]) or die "$usage\n ERROR: Could not open G-PhoCS control file $ARGV[0]\n";
	local $/ = undef;
	my $ctrlText = <CTRL>;
	close CTRL;
	die "$usage\n ERROR: File $ARGV[0] does not seem to be a G-PhoCS control file\n" if($ctrlText !~ /CURRENT\-POPS\-START/);
	$/ = "\n";
	
	open(TRACE,$ARGV[1]) or die "$usage\n ERROR: Could not open G-PhoCS trace file $ARGV[1]\n";
	my $firstLine = <TRACE>;
	close TRACE;
	die "$usage\n ERROR: File $ARGV[1] does not seem to be a G-PhoCS trace file\n" if($firstLine !~ /^Sample/);
	
	die "$usage\n ERROR: root name '$ARGV[2]' not found inside G-PhoCS control file $ARGV[0]\n" if($ctrlText !~ /name\s+$ARGV[2]/);
	
	if($ARGV[3] =~ /^[\d\.\-eE]+$/){
		my $num = $ARGV[3]*1;
		die "$usage\n ERROR: burn-in must be a numeric value between 0 and 1\n" if($num >= 1 or $num < 0);
	}
	else{
		die "$usage\n ERROR: burn-in must be a numeric value\n";
	}
	
	if(-e $ARGV[4]){
		print " File $ARGV[4] already exists. Overwrite? [Y/n] ";
		my $overwrite = <STDIN>;
		exit "\n" if($overwrite =~ /[nN]/);
	}
}
