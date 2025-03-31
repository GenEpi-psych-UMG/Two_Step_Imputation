#!/usr/bin/perl

#######################
#
# Extract imputation info from Michigan VCF HRC imputation files (.gz)
#
# (c) Alexander Teumer 2021 ateumer@uni-greifswald.de
#
#######################

$sep = "\t";   # separator in output file

# get filenames
$dirname = $ARGV[0];
if ($dirname eq "") {
  print STDERR "Input directory required!\n";
  exit 1;
}

$prefix = $ARGV[1];
if ($prefix eq "") {
  print STDERR "Prefix required (e.g. SHIP-0_R4)!\n";
  exit 1;
}

$filemask = $prefix.".*.vcf.gz\$";
opendir(DIR, "$dirname");
@files = grep(/$filemask/,readdir(DIR));
closedir(DIR);

# Get output file name
$outname = $dirname."/".$prefix.".all.info";
print "output name: ".$outname."\n";
open (OUT,">$outname.txt");
print OUT "SNPID".$sep."REF".$sep."ALT".$sep."Genotyped".$sep."MAF".$sep."R2".$sep."AVG_CS".$sep."ER2".$sep."ID"."\n";

$col_chr  = 0;
$col_position  = 1;
$col_id = 2;
$col_alleleA  = 3;
$col_alleleB  = 4;
$col_format  = 7; # INFO

foreach $filename (@files) {

print "Processing: $filename";

$cnt=-1;
open (F,"zcat $dirname/$filename|");
while(<F>) {
  
  my $line = $_;
  #$line=~s/[\n\r]//g;
  
  if ((($line) =~ /^#CHROM\t/) || ($cnt>=0)) {
  	$cnt++;
  }
  if ($cnt <= 0) { 
	next;
  }	

  my @vals = split("\t",$line,9);

  # generate modout values
  my $chr = $vals[$col_chr];
  $chr=~s/^0+//g; # remove leading zeros
  my $pos = $vals[$col_position];
  my $format = $vals[$col_format];
  my $id = $vals[$col_id];
  my $allele_a = $vals[$col_alleleA];
  my $allele_b = $vals[$col_alleleB];


  my $type = "SNP"; # HRC includes SNPs only so far

  my $snpid = $chr.':'.$pos.':'.$allele_a.':'.$allele_b;
  
  # effallele, freq, n
  $ref_all = $vals[$col_alleleA];
  $alt_all = $vals[$col_alleleB];

  # imputation info
  $oevar_imp = "NA";
  my @buf = split(";",$format);
  for ($i=scalar @buf - 1; $i>=0; $i--) {
  	($key,$val) = split("=",$buf[$i],2);
	if ($key eq "R2") {
		$oevar_imp = $val;
		last;
	}


  }
 
# imputation infoe
  $oevar_impe = "NA";
  my @buf = split(";",$format);
  for ($i=scalar @buf - 1; $i>=0; $i--) {
  	($key,$val) = split("=",$buf[$i],2);
	if ($key eq "ER2") {
		$oevar_impe = $val;
		last;
	}


  }

# Average call score
  $oevar_cs = "NA";
  my @buf = split(";",$format);
  for ($i=scalar @buf - 1; $i>1; $i--) {
  	($key,$val) = split("=",$buf[$i],2);
	if ($key eq "AVG_CS") {
		$oevar_cs = $val;
		last;
	}


  }

# Minor allele frequency
  $oevar_MAF = "NA";
  my @buf = split(";",$format);
  for ($i=scalar @buf - 1; $i>1; $i--) {
  	($key,$val) = split("=",$buf[$i],2);
	if ($key eq "MAF") {
		$oevar_MAF = $val;
		last;
	}


  }

# Genotyped
  $genotype = "Genotyped";
  my @buf = split(";",$format);
  for ($i=scalar @buf; $i>1; $i--) {
  	$key = $buf[1];
	if ($key eq "IMPUTED") {

		$genotype = "Imputed";
		last;
	}


  }


 
  # output
  print OUT $snpid.$sep.$ref_all.$sep.$alt_all.$sep.$genotype.$sep.$oevar_MAF.$sep.$oevar_imp.$sep.$oevar_cs.$sep.$oevar_impe.$sep.$id."\n";

}
close(F);
print " $cnt result lines processed.\n";

}
close(OUT);

print "Done.\n";