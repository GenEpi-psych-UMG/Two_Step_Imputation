#!/usr/bin/perl

#######################
#
# Extract imputation info from Beagle5 VCF imputation files (.gz)
#
# Adapted for Beagle5 specific INFO field format (DR2=0.66;AF=0.3624;IMP)
# April 2025
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
# Only include fields that are relevant for Beagle5 output
print OUT "SNPID".$sep."REF".$sep."ALT".$sep."Genotyped".$sep."AF".$sep."R2".$sep."ID"."\n";

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

  my $snpid = $chr.':'.$pos.':'.$allele_a.':'.$allele_b;
  
  # Reference and alternate alleles
  $ref_all = $vals[$col_alleleA];
  $alt_all = $vals[$col_alleleB];

  # Parse INFO field for Beagle5 format (DR2=0.66;AF=0.3624;IMP)
  my $dr2 = "NA";    # Will be mapped to R2 in output
  my $af = "NA";     # Allele frequency
  my $genotype = "Genotyped";  # Default is genotyped
  
  # Check INFO field entries
  my @buf = split(";",$format);
  foreach my $entry (@buf) {
    if ($entry =~ /^DR2=(.+)$/) {
      $dr2 = $1;
    }
    elsif ($entry =~ /^AF=(.+)$/) {
      $af = $1;
    }
    elsif ($entry eq "IMP") {
      $genotype = "Imputed";
    }
  }
  
  # Output only the required fields
  print OUT $snpid.$sep.$ref_all.$sep.$alt_all.$sep.$genotype.$sep.$af.$sep.$dr2.$sep.$id."\n";
}
close(F);
print " $cnt result lines processed.\n";

}
close(OUT);

print "Done.\n";