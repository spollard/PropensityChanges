# Takes all substitutions from a suboutfile and gives data frame of substitions
# including probability. Can include a probability threshold. Choose 0 for no
# threshold. 

use strict; 
use warnings;
use 5.010.001;


local $\ = "\n";


my $suboutfile = "suboutfile";
my $subsFile = "all_subs";

my $gens = 500;
my $burnin = 100;
my $pThreshold = 0;

my $threshold = $gens * $pThreshold;

my $fh = initSubs($suboutfile, $burnin);
open my $fh_out, ">$subsFile"or die $!;

print $fh_out "parent\tchild\tsite\tAA_parent\tAA_child\tprobability";

my $nodeSubs = {};
my $gen = 0;
while ($gen < $gens) {
	my $line = <$fh>;
	if ($line =~ /\/\//) { $gen++; print $gen; next } 
	if ($line =~ "parent") { next }
	chomp $line;
	
	my $fields = [split " ", $line];
	my $tag = shift @$fields; # tag contains parent..child e.g. 132..142
	
	foreach my $sub (@$fields) {
		$nodeSubs->{$tag . $sub}++;
	}
}

foreach my $nodeSub ( keys %$nodeSubs) {
	if ($nodeSubs->{$nodeSub} >= $threshold) { 
		$nodeSub =~ /(\d+)..(\d+)([ACDEFGHIKLMNPQRSTVWY])(\d+)([ACDEFGHIKLMNPQRSTVWY])/;
		print $fh_out "$1\t$2\t$4\t$3\t$5\t" . $nodeSubs->{$nodeSub} / $gens;
	} 
}


sub initSubs { my ($file, $burnin) = @_;
#	print $file;exit;
	open my $fh, "<", $file or die "cannot open substitution file"; 
	my $g = 0; my $line;
	<$fh>; # burn header line
	while ($g != $burnin and defined( $line = <$fh> )) {
		if ($line =~ /\/\//) {
#			print "found a generation $g" if $control->{debug}; 
			$g++;
		}

	}
	die "End of burnin not reached; check burnin" if $g != $burnin;
	return $fh;
}
