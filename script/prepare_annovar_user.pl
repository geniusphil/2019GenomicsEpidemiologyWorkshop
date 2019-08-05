#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;

our $REVISION = '$Revision: fcb810b8464da61f2ecb40d529b9fafcdfa13742 $';
our $DATE =	'$Date: 2018-07-08 23:15:35 -0400 (Sun,  8 Jul 2018) $';  
our $AUTHOR =	'$Author: Kai Wang <wangk@biocluster.wglab.org> $';

our ($verbose, $help, $man);
our ($dbtype, $dbfile, $outfile, $vcffile, $dbsnpfile, $buildver, $transcriptid, $cdotfile, $afstring, $acstring, $infostring, $twopos, $dbnsfpver);

our %iupac = (R=>'AG', Y=>'CT', S=>'GC', W=>'AT', K=>'GT', M=>'AC', A=>'AA', C=>'CC', G=>'GG', T=>'TT', B=>'CGT', D=>'AGT', H=>'ACT', V=>'ACG', N=>'ACGT', '.'=>'-', '-'=>'-');

GetOptions('verbose'=>\$verbose, 'help|h'=>\$help, 'man'=>\$man, 'outfile=s'=>\$outfile, 'dbtype=s'=>\$dbtype, 'vcffile=s'=>\$vcffile,
	'dbsnpfile=s'=>\$dbsnpfile, 'buildver=s'=>\$buildver, 'transcriptid=s'=>\$transcriptid, 'cdotfile=s'=>\$cdotfile, 'afstring=s'=>\$afstring,
	'acstring=s'=>\$acstring, 'infostring=s'=>\$infostring, 'twopos'=>\$twopos, 'dbnsfpver=s'=>\$dbnsfpver) or pod2usage ();

$help and pod2usage (-verbose=>1, -exitval=>1, -output=>\*STDOUT);
$man and pod2usage (-verbose=>2, -exitval=>1, -output=>\*STDOUT);
@ARGV or pod2usage (-verbose=>0, -exitval=>1, -output=>\*STDOUT);
@ARGV == 1 or pod2usage ("Syntax error");

($dbfile) = @ARGV;

$dbtype or pod2usage ("Error in argument: please specify --dbtype argument");
$dbtype eq 'cosmic' and $vcffile || pod2usage ("Error in argument: please specify --vcffile argument");

if ($outfile) {
	open (STDOUT, ">$outfile") or die "Error: cannot write to output file $outfile: $!\n";
}



if ($dbtype eq 'clinvar') {
	prepareClinVar ($dbfile);
} elsif ($dbtype eq 'clinvar2') {
	prepareClinVar2 ($dbfile);
} elsif ($dbtype eq 'clinvar_preprocess') {
	ClinvarPreprocess ($dbfile);
} elsif ($dbtype eq 'clinvar_preprocess2') {
	ClinvarPreprocess2 ($dbfile);
} elsif ($dbtype eq 'cosmic') {
	prepareCosmic ($dbfile, $vcffile);
} else {
	pod2usage ("Error: the -dbtype of $dbtype is not supported yet");
}





sub prepareClinVarOld {
	my ($dbfile) = @_;
	open (FH, "convert2annovar.pl -format vcf4 -include $dbfile |") or die "Error: cannot read from dbfile $dbfile: $!\n";
	
	if ($outfile) {
		open (STDOUT, ">$outfile") or die "Error: cannot write to output file: $!\n";
	}
	
	my %sig = (0 => 'unknown', 1 => 'untested', 2 => 'non-pathogenic', 3 => 'probable-non-pathogenic', 4 => 'probable-pathogenic', 5 => 'pathogenic', 6 => 'drug-response', 7 => 'histocompatibility', 255 => 'other');
	while (<FH>) {
		m/^#/ and next; 
		s/[\r\n]+$//;
		my @field  = split(/\t/, $_);
		
		$field[12]=~m/(CLNDBN=([^;]+);.*CLNACC=([^;]+))/ or die "Error: invalid record found in avinputfile: <$field[12]>";
		my $clndbnacc = $1;
		
		$field[12]=~m/(CLNDSDB=([^;]+);.*CLNDSDBID=([^;]+))/ or die "Error: invalid record found in avinputfile: <$field[12]>";
		my $clndsdb = $1;
		
		$field[12]=~m/CLNSIG=([\d\|\,]+)/ or die "Error: invalid record found in avinput file: <$field[12]>";	#CLNSIG=5|5|5,0|0; in RCV000013623.23 
		my @sigall = split (/\|/, $1);		#CLNSIG=5|2|2;
		my (@clnsig, $clnsig);
		for my $i (0 .. @sigall-1) {
			if ($sigall[$i] =~ m/,/) {
				my @temp = split (/,/, $sigall[$i]);
				@temp = map {$sig{$_} || 'unknown'} @temp;
				push @clnsig, join(",", @temp);
			} else {
				push @clnsig, $sig{$sigall[$i]} || 'unknown';
			}
		}
		$clnsig = join ('|', @clnsig);
		
		print join ("\t", @field[0..4]), "\t", "CLINSIG=$clnsig;$clndbnacc;$clndsdb\n";

	}
}




sub ClinvarPreprocess {
#the goal is to pre-process clinvar VCF (after splitting by VT) and delete unrelevant mutations from the file, before doing left-normalization
#prepare_annovar.pl   -dbtype clinvar_preprocess ~/temp8 -out ~/temp9
	my ($dbfile) = @_;
	open (FH, $dbfile) or die "Error: cannot read from inputfile $dbfile: $!\n";
	
	if ($outfile) {
		open (STDOUT, ">$outfile") or die "Error: cannot write to output file: $!\n";
	}
	
	my ($posindex, $prestring) = (0, '');
	
	while (<FH>) {
		
		if (m/^#/) {		#header lines are printed identically
			print;
			next;
		}
		
		if (m/(OLD_MULTIALLELIC=.+)$/) {	#only if we use vt for splitting multi-allelic variants in VCF files
			if ($prestring eq $1) {
				$posindex++;
			} else {
				$posindex = 0;		#found a new multiallelic variant
				$prestring = $1;
			}
		} else {
			($posindex, $prestring) = (0, '');
		}
		
		my @field  = split(/\t/, $_);
		$field[7]=~m/(CLNDBN=([^;]+);.*CLNACC=([^;]+))/ or die "Error: invalid record found in avinputfile: <$field[7]>";
		my ($clndbn, $clnacc) = ($2, $3);
		
		$field[7]=~m/(CLNDSDB=([^;]+);.*CLNDSDBID=([^;]+))/ or die "Error: invalid record found in avinputfile: <$field[7]>";
		my ($clndsdb, $clndsdbid) = ($2, $3);
		
		$field[7]=~m/CLNSIG=([^;]+)/ or die "Error: invalid record found in avinput file: <$field[7]>";
		my ($clnsig) = ($1);
	
		#the key here is to delete mutations are that are not really in the record itself
		#an example below that the C allele is actually not in the annotation
		#7       150648198       rs1137617       A       C,G,T   .       .       RS=1137617;RSPOS=150648198;RV;dbSNPBuildID=86;SSR=0;SAO=1;VP=0x05017800030515053e110100;GENEINFO=KCNH2:3757;WGT=1;VC=SNV;PM;TPA;PMC;SLO;REF;SYN;ASP;VLD;G5;HD;GNO;KGPhase1;KGPhase3;LSD;OM;CLNALLE=2,3;CLNHGVS=NC_000007.13:g.150648198A>G,NC_000007.13:g.150648198A>T;CLNSRC=.,.;CLNORIGIN=1,1;CLNSRCID=.,.;CLNSIG=2|2|3,5;CLNDSDB=MedGen|MedGen|MeSH:MedGen:SNOMED_CT,MedGen;CLNDSDBID=CN169374|CN230736|D008133:C0023976:9651007,CN221809;CLNDBN=not_specified|Cardiovascular_phenotype|Long_QT_syndrome,not_provided;CLNREVSTAT=mult|single|single,single;CLNACC=RCV000181727.2|RCV000253499.1|RCV000283094.1,RCV000181832.1;CAF=0.2278,.,0.7722,.;COMMON=1
		if ($prestring) {	#if there is a multi-allelic variant
			$field[7] =~ m/CLNHGVS=([^;]+);/ or die "Error: invalid record found in avinput file: <$field[7]>";
			my @hgvs = split (/,/, $1);
			#defined $hgvs[$posindex] or warn "posindex ($posindex) not found in $field[12] (@field[0..4])" and next;	#something is wrong
			defined $hgvs[$posindex] or next;
			
			if ($hgvs[$posindex] =~ m/(?:g|m)\.\d+(\w+)>(\w+)/) {
				my ($nt1, $nt2) = ($1, $2);
				if ($nt2 ne $field[4]) {
					$posindex--;
					next;
				}
			} elsif ($hgvs[$posindex] =~ m/g\.\d+(\w+)\\x3d/) {	#RS=1042714 CLNHGVS=NC_000005.9:g.148206473G\x3d
				1;	#nothing needs to be done, assume that the first occurence is the one that we are interested in
			
			} elsif ($hgvs[$posindex] =~ m/(?:g|m)\.\d+dup(\w)$/) {	#RS=142323886 CLNHGVS=NC_000012.11:g.25359046_25359047dupAA,NC_000012.11:g.25359046dupA
				if (not $field[4] =~ m/\w$1$/) {
					$posindex--;
					next;
				}
			} elsif ($hgvs[$posindex] =~ m/(?:g|m)\.\d+_\d+dup(\w+)$/) {	#RS=113564356 CLNHGVS=NC_000003.11:g.150645419_150645420dupAC,NC_000003.11:g.150645419_150645422dupACAC
				if (not $field[4] =~ m/^\w$1$/) {
					$posindex--;
					next;
				}
			} elsif ($hgvs[$posindex] =~ m/g\.(\d+)_(\d+)dup$/) {	#RS=113564356 CLNHGVS=NC_000003.11:g.150645419_150645420dupAC,NC_000003.11:g.150645419_150645422dupACAC
				if (length($field[4])-length($field[3]) != $2-$1+1) {
					$posindex--;
					next;
				}
			} elsif ($hgvs[$posindex] =~ m/g\.\d+dup$/) {	#RS=113564356 CLNHGVS=NC_000003.11:g.150645419_150645420dupAC,NC_000003.11:g.150645419_150645422dupACAC
				if (length($field[4])-length($field[3]) != 1) {
					$posindex--;
					next;
				}
			} elsif ($hgvs[$posindex] =~ m/g\.\d+_\d+ins(\w+)$/) {	#RS=779985493 CLNHGVS=NC_000001.10:g.241663898_241663899insGAGAGA
				if ($field[3].$1 ne $field[4]) {
					$posindex--;
					next;
				}
			} elsif ($hgvs[$posindex] =~ m/g\.\d+_\d+del(\d+)$/) {	#RS=267608120 CLNHGVS=NC_000002.11:g.48033292_48033308del17
				if (length($field[3]) - length($field[4]) != $1) {
					$posindex--;
					next;
				}
			} elsif ($hgvs[$posindex] =~ m/g\.\d+_\d+del(\w+)$/) {	#RS=267608120 CLNHGVS=NC_000002.11:g.48033748_48033751delCAAG
				if (not ($field[3] =~ m/\w$1$/ and $field[4] =~ m/^\w$/)) {
					$posindex--;
					next;
				}
			
			} elsif ($hgvs[$posindex] =~ m/g\.\d+_\d+del(\w+)ins(\w+)$/) {	#RS=794728490 CLNHGVS=NC_000007.13:g.150648198_150648200delATAinsGTG
				if (not ($field[3] =~ m/\w$1$/ and $field[4] =~ m/^\w$2$/)) {
					$posindex--;
					next;
				}
			} elsif ($hgvs[$posindex] =~ m/g\.\d+del(\w)/) {	#RS=267608120 CLNHGVS=NC_000002.11:g.48033748_48033751delCAAG
				if (not ($field[3] =~ m/\w$1$/ and $field[4] =~ m/^\w$/)) {
					$posindex--;
					next;
				}
			} elsif ($hgvs[$posindex] =~ m/g\.\d+_\d+inv(\w+)/) {	#RS=879353161 CLNHGVS=NC_000010.10:g.73115941_73115942invTG
				if (not $field[3] =~ m/\w$1$/) {
					$posindex--;
					next;
				}
			} else {
				warn "WARNING: non-canonical pattern found in HGVS: $field[7]>";
			}
			
			my @temp;
			@temp = split (/,/, $clnsig);
			$clnsig = $temp[$posindex];
			@temp = split (/,/, $clndbn);
			$clndbn = $temp[$posindex];
			@temp = split (/,/, $clndsdb);
			$clndsdb = $temp[$posindex];
			@temp = split (/,/, $clndsdbid);
			$clndsdbid = $temp[$posindex];
			@temp = split (/,/, $clnacc);
			$clnacc = $temp[$posindex];
			
			
			$field[7] =~ s/CLNDBN=([^;]+)/CLNDBN=$clndbn/;
			$field[7] =~ s/CLNACC=([^;]+)/CLNACC=$clnacc/;
			$field[7] =~ s/CLNDSDB=([^;]+)/CLNDSDB=$clndsdb/;
			$field[7] =~ s/CLNDSDBID=([^;]+)/CLNDSDBID=$clndsdbid/;
			$field[7] =~ s/CLNSIG=([^;]+)/CLNSIG=$clnsig/;
			
			$field[7] =~ s/CLNHGVS=([^;]+)/CLNHGVS=$hgvs[$posindex]/;
			
		}
		
		print STDOUT join("\t", @field);
	}
}

sub ClinvarPreprocess2 {
#the goal is to pre-process clinvar VCF (after splitting by VT) and delete unrelevant mutations from the file, before doing left-normalization
#prepare_annovar.pl   -dbtype clinvar_preprocess ~/temp8 -out ~/temp9
	my ($dbfile) = @_;
	open (FH, $dbfile) or die "Error: cannot read from inputfile $dbfile: $!\n";
	
	if ($outfile) {
		open (STDOUT, ">$outfile") or die "Error: cannot write to output file: $!\n";
	}
	
	my ($posindex, $prestring) = (0, '');
	
	while (<FH>) {
		
		if (m/^#/) {		#header lines are printed identically
			print;
			next;
		}
		
		if (m/(OLD_MULTIALLELIC=.+)$/) {	#only if we use vt for splitting multi-allelic variants in VCF files
			if ($prestring eq $1) {
				$posindex++;
			} else {
				$posindex = 0;		#found a new multiallelic variant
				$prestring = $1;
			}
		} else {
			($posindex, $prestring) = (0, '');
		}
		

##INFO=<ID=ALLELEID,Number=1,Type=Integer,Description="the ClinVar Allele ID">
##INFO=<ID=CLNDN,Number=.,Type=String,Description="ClinVar's preferred disease name for the concept specified by disease identifiers in CLNDISDB">
##INFO=<ID=CLNDNINCL,Number=.,Type=String,Description="For included Variant : ClinVar's preferred disease name for the concept specified by disease identifiers in CLNDISDB">
##INFO=<ID=CLNDISDB,Number=.,Type=String,Description="Tag-value pairs of disease database name and identifier, e.g. OMIM:NNNNNN">
##INFO=<ID=CLNDISDBINCL,Number=.,Type=String,Description="For included Variant: Tag-value pairs of disease database name and identifier, e.g. OMIM:NNNNNN">
##INFO=<ID=CLNHGVS,Number=.,Type=String,Description="Top-level (primary assembly, alt, or patch) HGVS expression.">
##INFO=<ID=CLNREVSTAT,Number=.,Type=String,Description="ClinVar review status for the Variation ID">
##INFO=<ID=CLNSIG,Number=.,Type=String,Description="Clinical significance for this single variant">
##INFO=<ID=CLNSIGCONF,Number=.,Type=String,Description="Conflicting clinical significance for this single variant">
##INFO=<ID=CLNSIGINCL,Number=.,Type=String,Description="Clinical significance for a haplotype or genotype that includes this variant. Reported as pairs of VariationID:clinical significance.">
##INFO=<ID=CLNVC,Number=1,Type=String,Description="Variant type">
##INFO=<ID=CLNVCSO,Number=1,Type=String,Description="Sequence Ontology id for variant type">
##INFO=<ID=CLNVI,Number=.,Type=String,Description="the variant's clinical sources reported as tag-value pairs of database and variant identifier">
##INFO=<ID=DBVARID,Number=.,Type=String,Description="nsv accessions from dbVar for the variant">
##INFO=<ID=GENEINFO,Number=1,Type=String,Description="Gene(s) for the variant reported as gene symbol:gene id. The gene symbol and id are delimited by a colon (:) and each pair is delimited by a vertical bar (|)">
##INFO=<ID=MC,Number=.,Type=String,Description="comma separated list of molecular consequence in the form of Sequence Ontology ID|molecular_consequence">
##INFO=<ID=ORIGIN,Number=.,Type=String,Description="Allele origin. One or more of the following values may be added: 0 - unknown; 1 - germline; 2 - somatic; 4 - inherited; 8 - paternal; 16 - maternal; 32 - de-novo; 64 - biparental; 128 - uniparental; 256 - not-tested; 512 - tested-inconclusive; 1073741824 - other">
##INFO=<ID=RS,Number=.,Type=String,Description="dbSNP ID (i.e. rs number)">
##INFO=<ID=SSR,Number=1,Type=Integer,Description="Variant Suspect Reason Codes. One or more of the following values may be added: 0 - unspecified, 1 - Paralog, 2 - byEST, 4 - oldAlign, 8 - Para_EST, 16 - 1kg_failed, 1024 - other">

		my @field  = split(/\t/, $_);
		
		$field[7]=~m/CLNDNINCL/ and next;	#this is for included variant
			#Interpretations may be made on a single variant or a set of variants, such as a haplotype. Variants
	 		#that have only been interpreted as part of a set of variants (i.e. no direct interpretation for the
			#variant itself) are considered "included" variants. The VCF files include both variants with a direct
			#interpretation and included variants. Included variants do not have an associated disease (CLNDN, 
			#CLNDISDB) or a clinical significance (CLNSIG). Instead there are three tags are specific to the 
			#included variants - CLNDNINCL, CLNDISDBINCL, and CLNSIGINCL (see below).

		
		$field[7]=~m/ALLELEID=([^;]+)/ or die "Error: invalid record found in avinputfile: <$field[7]>";
		my $alleleid = $1;
		
		$field[7]=~m/CLNDN=([^;]+)/ or die "Error: invalid record found in avinputfile: <$field[7]>";
		my $clndn = $1;
		
		$field[7]=~m/CLNDISDB=([^;]+)/ or die "Error: invalid record found in avinputfile: <$field[7]>";
		my $clndisdb = $1;
		
		$field[7]=~m/CLNREVSTAT=([^;]+)/ or die "Error: invalid record found in avinputfile: <$field[7]>";
		my $clnrevstat = $1;
		
		$field[7]=~m/CLNSIG=([^;]+)/ or die "Error: invalid record found in avinputfile: <$field[7]>";
		my $clnsig = $1;
		
	
		#the key here is to delete mutations are that are not really in the record itself
		#an example below that the C allele is actually not in the annotation
		#7       150648198       rs1137617       A       C,G,T   .       .       RS=1137617;RSPOS=150648198;RV;dbSNPBuildID=86;SSR=0;SAO=1;VP=0x05017800030515053e110100;GENEINFO=KCNH2:3757;WGT=1;VC=SNV;PM;TPA;PMC;SLO;REF;SYN;ASP;VLD;G5;HD;GNO;KGPhase1;KGPhase3;LSD;OM;CLNALLE=2,3;CLNHGVS=NC_000007.13:g.150648198A>G,NC_000007.13:g.150648198A>T;CLNSRC=.,.;CLNORIGIN=1,1;CLNSRCID=.,.;CLNSIG=2|2|3,5;CLNDSDB=MedGen|MedGen|MeSH:MedGen:SNOMED_CT,MedGen;CLNDSDBID=CN169374|CN230736|D008133:C0023976:9651007,CN221809;CLNDBN=not_specified|Cardiovascular_phenotype|Long_QT_syndrome,not_provided;CLNREVSTAT=mult|single|single,single;CLNACC=RCV000181727.2|RCV000253499.1|RCV000283094.1,RCV000181832.1;CAF=0.2278,.,0.7722,.;COMMON=1
		if ($prestring) {	#if there is a multi-allelic variant
			$field[7] =~ m/CLNHGVS=([^;]+);/ or die "Error: invalid record found in avinput file: <$field[7]>";
			my @hgvs = split (/,/, $1);
			#defined $hgvs[$posindex] or warn "posindex ($posindex) not found in $field[12] (@field[0..4])" and next;	#something is wrong
			defined $hgvs[$posindex] or next;
			
			if ($hgvs[$posindex] =~ m/(?:g|m)\.\d+(\w+)>(\w+)/) {
				my ($nt1, $nt2) = ($1, $2);
				if ($nt2 ne $field[4]) {
					$posindex--;
					next;
				}
			} elsif ($hgvs[$posindex] =~ m/g\.\d+(\w+)\\x3d/) {	#RS=1042714 CLNHGVS=NC_000005.9:g.148206473G\x3d
				1;	#nothing needs to be done, assume that the first occurence is the one that we are interested in
			
			} elsif ($hgvs[$posindex] =~ m/(?:g|m)\.\d+dup(\w)$/) {	#RS=142323886 CLNHGVS=NC_000012.11:g.25359046_25359047dupAA,NC_000012.11:g.25359046dupA
				if (not $field[4] =~ m/\w$1$/) {
					$posindex--;
					next;
				}
			} elsif ($hgvs[$posindex] =~ m/(?:g|m)\.\d+_\d+dup(\w+)$/) {	#RS=113564356 CLNHGVS=NC_000003.11:g.150645419_150645420dupAC,NC_000003.11:g.150645419_150645422dupACAC
				if (not $field[4] =~ m/^\w$1$/) {
					$posindex--;
					next;
				}
			} elsif ($hgvs[$posindex] =~ m/g\.(\d+)_(\d+)dup$/) {	#RS=113564356 CLNHGVS=NC_000003.11:g.150645419_150645420dupAC,NC_000003.11:g.150645419_150645422dupACAC
				if (length($field[4])-length($field[3]) != $2-$1+1) {
					$posindex--;
					next;
				}
			} elsif ($hgvs[$posindex] =~ m/g\.\d+dup$/) {	#RS=113564356 CLNHGVS=NC_000003.11:g.150645419_150645420dupAC,NC_000003.11:g.150645419_150645422dupACAC
				if (length($field[4])-length($field[3]) != 1) {
					$posindex--;
					next;
				}
			} elsif ($hgvs[$posindex] =~ m/g\.\d+_\d+ins(\w+)$/) {	#RS=779985493 CLNHGVS=NC_000001.10:g.241663898_241663899insGAGAGA
				if ($field[3].$1 ne $field[4]) {
					$posindex--;
					next;
				}
			} elsif ($hgvs[$posindex] =~ m/g\.\d+_\d+del(\d+)$/) {	#RS=267608120 CLNHGVS=NC_000002.11:g.48033292_48033308del17
				if (length($field[3]) - length($field[4]) != $1) {
					$posindex--;
					next;
				}
			} elsif ($hgvs[$posindex] =~ m/g\.\d+_\d+del(\w+)$/) {	#RS=267608120 CLNHGVS=NC_000002.11:g.48033748_48033751delCAAG
				if (not ($field[3] =~ m/\w$1$/ and $field[4] =~ m/^\w$/)) {
					$posindex--;
					next;
				}
			
			} elsif ($hgvs[$posindex] =~ m/g\.\d+_\d+del(\w+)ins(\w+)$/) {	#RS=794728490 CLNHGVS=NC_000007.13:g.150648198_150648200delATAinsGTG
				if (not ($field[3] =~ m/\w$1$/ and $field[4] =~ m/^\w$2$/)) {
					$posindex--;
					next;
				}
			} elsif ($hgvs[$posindex] =~ m/g\.\d+del(\w)/) {	#RS=267608120 CLNHGVS=NC_000002.11:g.48033748_48033751delCAAG
				if (not ($field[3] =~ m/\w$1$/ and $field[4] =~ m/^\w$/)) {
					$posindex--;
					next;
				}
			} elsif ($hgvs[$posindex] =~ m/g\.\d+_\d+inv(\w+)/) {	#RS=879353161 CLNHGVS=NC_000010.10:g.73115941_73115942invTG
				if (not $field[3] =~ m/\w$1$/) {
					$posindex--;
					next;
				}
			} else {
				warn "WARNING: non-canonical pattern found in HGVS: $field[7]>";
			}
			
			my @temp;
			@temp = split (/,/, $alleleid);
			$alleleid = $temp[$posindex];
			
			@temp = split (/,/, $clndn);
			$clndn = $temp[$posindex];
			
			@temp = split (/,/, $clndisdb);
			$clndisdb = $temp[$posindex];
			
			@temp = split (/,/, $clnrevstat);
			$clnrevstat = $temp[$posindex];
			
			@temp = split (/,/, $clnsig);
			$clnsig = $temp[$posindex];
					
			
			$field[7] =~ s/ALLELEID=([^;]+)/ALLELEID=$alleleid/;
			$field[7] =~ s/CLNDN=([^;]+)/CLNDN=$clndn/;
			$field[7] =~ s/CLNDISDB=([^;]+)/CLNDISDB=$clndisdb/;
			$field[7] =~ s/CLNREVSTAT=([^;]+)/CLNREVSTAT=$clnrevstat/;
			$field[7] =~ s/CLNSIG=([^;]+)/CLNSIG=$clnsig/;
			
			$field[7] =~ s/CLNHGVS=([^;]+)/CLNHGVS=$hgvs[$posindex]/;
			
		}
		
		print STDOUT join("\t", @field);
	}
}

sub prepareClinVar {
	my ($dbfile) = @_;
	open (FH, "convert2annovar.pl -format vcf4 -include $dbfile |") or die "Error: cannot read from dbfile $dbfile: $!\n";
	
	if ($outfile) {
		open (STDOUT, ">$outfile") or die "Error: cannot write to output file: $!\n";
	}
	
	#my %sig = (0 => 'unknown', 1 => 'untested', 2 => 'non-pathogenic', 3 => 'probable-non-pathogenic', 4 => 'probable-pathogenic', 5 => 'pathogenic', 6 => 'drug-response', 7 => 'histocompatibility', 255 => 'other');
	my %sig = (0 => 'Uncertain significance', 1 => 'not provided', 2 => 'Benign', 3 => 'Likely benign', 4 => 'Likely pathogenic', 5 => 'Pathogenic', 6 => 'drug response', 7 => 'histocompatibility', 255 => 'other');
	while (<FH>) {
		m/^#/ and next; 
		s/[\r\n]+$//;
		my @field  = split(/\t/, $_);
		
		$field[12]=~m/(CLNDBN=([^;]+);.*CLNACC=([^;]+))/ or die "Error: invalid record found in avinputfile: <$field[12]>";
		my ($clndbn, $clnacc) = ($2, $3);
		
		$field[12]=~m/(CLNDSDB=([^;]+);.*CLNDSDBID=([^;]+))/ or die "Error: invalid record found in avinputfile: <$field[12]>";
		my ($clndsdb, $clndsdbid) = ($2, $3);
		
		#the goal below is to split both "|" and ",". I should first split on "," before splitting on "|", but did not in an opposite way. It does not affect results though
		$field[12]=~m/CLNSIG=([\d\|\,]+)/ or die "Error: invalid record found in avinput file: <$field[12]>";	#CLNSIG=5|5|5,0|0; in RCV000013623.23 
		my @sigall = split (/\|/, $1);		#CLNSIG=5|2|2;
		my (@clnsig, $clnsig);
		for my $i (0 .. @sigall-1) {
			if ($sigall[$i] =~ m/,/) {
				my @temp = split (/,/, $sigall[$i]);
				@temp = map {$sig{$_} || 'unknown'} @temp;
				push @clnsig, join(",", @temp);
			} else {
				push @clnsig, $sig{$sigall[$i]} || 'unknown';
			}
		}
		$clnsig = join ('|', @clnsig);
		
		$clnsig =~ s/,/\\x2c/g;
		$clndbn =~ s/,/\\x2c/g;
		$clnacc =~ s/,/\\x2c/g;
		$clndsdb =~ s/,/\\x2c/g;
		$clndsdbid =~ s/,/\\x2c/g;
		
		print join ("\t", @field[0..4], $clnsig, $clndbn, $clnacc, $clndsdb, $clndsdbid), "\n";

	}
}

sub prepareClinVar2 {
	my ($dbfile) = @_;
	open (FH, "convert2annovar.pl -format vcf4 -include $dbfile |") or die "Error: cannot read from dbfile $dbfile: $!\n";
	
	if ($outfile) {
		open (STDOUT, ">$outfile") or die "Error: cannot write to output file: $!\n";
	}
	
	while (<FH>) {
		m/^#/ and next; 
		s/[\r\n]+$//;
		my @field  = split(/\t/, $_);
		
		$field[12]=~m/ALLELEID=([^;]+)/ or die "Error: invalid record found in avinputfile: <$field[7]>";
		my $alleleid = $1;
		
		$field[12]=~m/CLNDN=([^;]+)/ or die "Error: invalid record found in avinputfile: <$field[7]>";
		my $clndn = $1;
		
		$field[12]=~m/CLNDISDB=([^;]+)/ or die "Error: invalid record found in avinputfile: <$field[7]>";
		my $clndisdb = $1;
		
		$field[12]=~m/CLNREVSTAT=([^;]+)/ or die "Error: invalid record found in avinputfile: <$field[7]>";
		my $clnrevstat = $1;
		
		$field[12]=~m/CLNSIG=([^;]+)/ or die "Error: invalid record found in avinputfile: <$field[7]>";
		my $clnsig = $1;
		
				
		$alleleid =~ s/,/\\x2c/g;
		$clndn =~ s/,/\\x2c/g;
		$clndisdb =~ s/,/\\x2c/g;
		$clnrevstat =~ s/,/\\x2c/g;
		$clnsig =~ s/,/\\x2c/g;
		
		print join ("\t", @field[0..4], $alleleid, $clndn, $clndisdb, $clnrevstat, $clnsig), "\n";

	}
}
sub prepareCosmic {
	my ($dbfile, $vcffile) = @_;
	my %mut;		#key=mutation ID value=chr+pos+ref+alt
	my (%cosmic, %cosmicid, %nonfound);
	my ($idprefix);
	
	if ($vcffile =~ m/\.vcf\.gz$/) {
		open (VCF, "gunzip -c $vcffile |") or die "Error: cannot read from gunzip: $!\n";
	} else {
		open (VCF, $vcffile) or die "Error: cannot read from VCF file: $!\n";
	}
	$_ = <VCF>;
	m/^##fileformat=VCFv4/ or die "Error: the supplied VCF file does not have valid version 4 header\n";
	while (<VCF>) {
		m/^#/ and next;
		my @field = split (/\t/, $_);
		
		#following paragraph is no longer necessary with new version of COSMIC
		#if ($field[2] =~ s/^COSM//) {		#for coding variants
		#	$idprefix = "COSM";
		#} elsif ($field[2] =~ s/^COSN//) {		#for non-coding variants
		#	$idprefix = "COSN";
		#} else {
		#	die "Error: invalid record found in VCF file: ID does not start with COSM or COSN: <$_>\n";
		#}
		
		
		my ($chr, $start, $id, $ref_allele, $mut_allele) = @field[0..4];
		my $end = $start + length($ref_allele) - 1;
		
		$ref_allele =~ m/[^ACGTacgt]/ and next;
		$mut_allele =~ m/[^ACGTacgt]/ and next;
		
		if (length ($ref_allele) > 1 or length ($mut_allele) > 1) {
			if(length($ref_allele) > length ($mut_allele)) { 		# deletion or block substitution
				my $head = substr($ref_allele, 0, length ($mut_allele));
				if ($head eq $mut_allele) {
					$start = $start+length($head);			#then change start position
					
					$ref_allele = substr ($ref_allele, length ($mut_allele));

					$mut_allele = "-";
				} elsif (substr ($ref_allele, 0, 1) eq substr ($mut_allele, 0, 1)) {	#first base is identical (this is used in many VCF files)
					$start++;
					$ref_allele = substr ($ref_allele, 1);
					$mut_allele = substr ($mut_allele, 1);
				}
					
			} elsif(length($mut_allele) >= length ($ref_allele)) { 		# insertion or block substitution
				my $head = substr ($mut_allele, 0, length ($ref_allele));
				if ($head eq $ref_allele) {
					$start = $start+length($ref_allele)-1;
					
					$mut_allele = substr ($mut_allele, length ($ref_allele));
					$ref_allele = '-';
				} elsif (substr ($ref_allele, 0, 1) eq substr ($mut_allele, 0, 1)) {	#first base is identical (this is used in many VCF files)
					$start++;
					$ref_allele = substr ($ref_allele, 1);
					$mut_allele = substr ($mut_allele, 1);
				}
			}
		}
		

		$mut{$id} = join ("\t", $chr, $start, $end, $ref_allele, $mut_allele);
	}
	print STDERR "NOTICE: Finished reading ", scalar (keys %mut), " mutation ID from the VCF file $vcffile\n";
	
	open (DB, $dbfile) or die "Error: cannot read from DB file: $!\n";
	$_ = <DB>;
	my @field = split (/\t/, $_);
	my ($mutid, $id_tumor, $primary_site, $filetype);
	
	if ($field[16] eq 'Mutation ID') {
		$mutid = $field[16];
		$filetype = 'coding';
	} elsif ($field[11] eq 'ID_NCV') {
		$mutid = $field[11];
		$filetype = 'noncoding';
	} else {
		die "Error: COSMIC MutantExport format error: column 17 or 12 should be 'Mutation ID' or 'ID_NCV'\n";
	}
	while (<DB>) {
		@field = split (/\t/, $_);
		
		if ($filetype eq 'coding') {
			$mutid = $field[16];
			$id_tumor = $field[6];
			$primary_site = $field[7];
		} elsif ($filetype eq 'noncoding') {
			$mutid = $field[11];
			$id_tumor = $field[2];
			$primary_site = $field[3];
		}
		
		my $chrstring = $mut{$mutid};
		if (not defined $chrstring) {
			$nonfound{$mutid}++;
			next;
		}
		$cosmic{$chrstring} .= ";$id_tumor,$primary_site";		#these two columns are "ID_tumour" and "Primary site"
		#$cosmicid{$chrstring} .= ",$idprefix$field[16]";	#this is no longer necessary with new version of COSMIC
		$cosmicid{$chrstring} .= ",$mutid";
	}
	print STDERR "NOTICE: Finished reading ", scalar (keys %cosmic), " COSMIC records in DB file $dbfile\n";
	print STDERR "WARNING: ", scalar (keys %nonfound), " COSMIC ID from MutantExport file cannot be found in VCF file (this may be normal if the VCF file only contains coding or noncoding variants\n";
	
	#the following is to eliminate duplicate entries. For example, COSM256593 and COSM256594 refers to the same mutation in the same sample, but are annotated twice, once with CDK11B/NM_001787 and once with CDC2L2/ENST00000357760
	
	for my $key (keys %cosmic) {
		my @id_site = split (/;/, $cosmic{$key});
		shift @id_site;
		my (%found_tumorid, %found_site);
		for my $i (0 .. @id_site-1) {
			my ($tumorid, $site) = split (/,/, $id_site[$i]);
			$found_tumorid{$tumorid} and next;		#the idea is that for duplicate entries, the tumor ID must be identical
			$found_tumorid{$tumorid}++;
			$found_site{$site}++;
		}
		
		my @cosmicid;
		@cosmicid = split (/,/, $cosmicid{$key});
		shift @cosmicid;
		my %cosmicid = map {$_, 1} @cosmicid;
		
		print $key, "\t", "ID=", join (",", keys %cosmicid), ";", "OCCURENCE=";
		my $occurence;
		for my $site (keys %found_site) {
			$occurence .= "," . $found_site{$site} . "(" . $site . ")";
		}
		$occurence =~ s/^,//;
		print "$occurence\n";
	}			
}


sub sum {
	my $sum = 0;
	for (@_) {
		$sum += $_;
	}
	return $sum;
}


sub revcom {
	my ($seq) = @_;
	$seq = reverse $seq;
	$seq =~ tr/acgtACGT/tgcaTGCA/;
	return ($seq);
}


=head1 SYNOPSIS

 prepare_annovar_user.pl [arguments] <dbfile | stdin>

 Optional arguments:
        -h, --help                      print help message
        -m, --man                       print complete documentation
        -v, --verbose                   use verbose output
            --outfile <file>		output file name
            --dbtype <string>		specify database type (cosmic, clinvar, clinvar2)

 Function: reformat and prepare ANNOVAR annotation database
        
 Example: # For preparing COSMIC database for annotation in ANNOVAR
          prepare_annovar_user.pl -dbtype cosmic CosmicMutantExport.tsv  -vcf CosmicCodingMuts.vcf > hg38_cosmic76.txt
          prepare_annovar_user.pl -dbtype cosmic CosmicMutantExport.tsv  -vcf CosmicNonCodingVariants.vcf >> hg38_cosmic76.txt
          
          # For preparing CLINVAR database for annotation in ANNOVAR
          wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20180603.vcf.gz
          wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20180603.vcf.gz.tbi
          vt decompose clinvar_20180603.vcf.gz -o temp.split.vcf
          prepare_annovar_user.pl   -dbtype clinvar_preprocess2 temp.split.vcf -out temp.split2.vcf
          vt normalize temp.split2.vcf -r ~/project/seqlib/GRCh38/old/GRCh38.fa -o temp.norm.vcf -w 2000000
          prepare_annovar_user.pl -dbtype clinvar2 temp.norm.vcf -out hg38_clinvar_20180603_raw.txt
          index_annovar.pl hg38_clinvar_20180603_raw.txt -out hg38_clinvar_20180603.txt -comment comment_20180708.txt
 
 Version: $Date: 2018-07-08 23:15:35 -0400 (Sun,  8 Jul 2018) $

=head1 OPTIONS

=over 8

=item B<--help>

print a brief usage message and detailed explanation of options.

=item B<--man>

print the complete manual of the program.

=item B<--verbose>

use verbose output.

=item B<--outfile>

specify the output file name. By default, output is written to STDOUT.

=item B<--dbtype>

specify database type, such as esp, cosmicgff, hgmdgff, etc

=back

=head1 DESCRIPTION

This program is used to convert various databases to a format that can be used 
by ANNOVAR.

For questions or comments, please contact kai@openbioinformatics.org.

=cut