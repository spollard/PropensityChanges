# Read a tree file, prints the internal id assigned to each node. 
# Only works on old tree out files that have branch ids.
# For example 
# (((Tarentola_mauritanica_8569 11..1:0.4049072,(Coleonyx_variegatus_52435 10..2:0.2211825


use strict; 
use warnings;


package Treeoutfile; 


unless (caller) {
	my $tree_file = $ARGV[0] // "../Data/pipeline/treeoutfile_head";
	my $nodes_file = $ARGV[1] // "nodes";
	
	ToNodes($tree_file, $nodes_file);
}

sub ToNodes {
	my ($tree_file, $nodes_file) = @_;
	my $nodes = ReadNodes($tree_file);
	
	PrintNodes($nodes, $nodes_file);
}


sub ReadNodes { 
	my ($tree_file) = @_;
	open my $tree_in, '<', $tree_file or die $!;
	
	my $line = scalar <$tree_in>;
	
	my $matches = [$line =~ /\w* \d*..\d*/g];
	
	my $nodes = {0 => "Node_0"};
	
	for my $match (@$matches) {
		$match =~ /(\w*) \d*..(\d*)/;
		my $name = $1;
		my $id = $2;
		$nodes->{$2} = ($1 || "Node_$2");
	}
	
	return $nodes
}

sub PrintNodes {
	my ($nodes, $nodes_file) = @_;
	
	open my $nodes_out, '>', $nodes_file;
	print $nodes_out "ID\tName\n";
	
	while (my ($id, $name) = each %$nodes) {
		print $nodes_out "$id\t$name\n";	
	}
}
