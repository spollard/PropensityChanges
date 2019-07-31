# Convert an old suboutfile with branch labels (e.g. ##..##) to new suboutfiles
# with node labels (e.g. Node_136 or Canis_lupus_8234). 

use strict;
use warnings;

package Suboutfile;

unless (caller) {
	my $suboutfile_old = $ARGV[0] // "../suboutfile_highprob_ASRV_updatebls7";
	my $nodes_file = $ARGV[1] // "nodes";
	my $suboutfile_new = $ARGV[2] // ($suboutfile_old . "_new");  
	
	OldToNew($suboutfile_old, $nodes_file, $suboutfile_new);
}

sub OldToNew {
	my ($suboutfile_old, $nodes_file, $suboutfile_new) = @_;
	
	my $nodes = ReadNodes($nodes_file);
	ConvertSuboutfile($suboutfile_old, $nodes, $suboutfile_new);
}

sub ReadNodes {
	my ($nodes_file) = @_;
	open my $nodes_in, $nodes_file or die $!;
	<$nodes_in>; # Burn header line;
	
	my $nodes = {};
	while (<$nodes_in>) {
		my ($id, $name) = split;
		$nodes->{$id} = $name;
	}
	return $nodes;
}

sub ConvertSuboutfile {
	my ($suboutfile_old, $nodes, $suboutfile_new) = @_;
	open my $fh, $suboutfile_old or die $!;
	
	<$fh>; # Burn header line;
	open my $fh_out, ">", $suboutfile_new; 
	
	print $fh_out "Node\tSubstitutions\n";
	
	while (my $line = <$fh>) {
		$line =~ s/\d+\.\.(\d+)/$nodes->{$1}/;
		print $fh_out $line;
	}
}