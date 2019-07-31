# Find the amino acid pairs from common substitution sites in
# the branch pairs

use strict;
use warnings;


use Node;
use NodeHash;

local $\ = "\n";


my $treeFile =  $ARGV[0] // "../Data/treeoutfile_old";
my $subFile =  $ARGV[1] // "../Data/suboutfile_randomized_keep_anc_old";
my $nodepairs = $ARGV[2] // "../Data/nodepairs";
my $NodePairSubs_file = $ARGV[3] // "../Data/NodePairSubs_randomized_keep_anc_old";

open my $fhNPSout, '>', $NodePairSubs_file;
open my $fhTree, "<", $treeFile;
open my $fhSub, "<", $subFile;

(my $nodes, $fhTree) = readNextTree($fhTree);

my $rootId = findRootID($nodes);
$nodes->findAllAncestors($rootId);


open my $fhNP, $nodepairs or die $!;
my $nodePairs = readNodePairs($fhNP);
initNodePairs($nodePairs);

# Reset the file positions
$fhTree = initTree($treeFile, 0);
$fhSub = initSubs($subFile, 0);

$fhTree = updateDistances($nodes, $fhTree);
$nodes->propAllDist($rootId);

$fhSub = updateSubs($nodes, $fhSub);

updateNodePairs($nodes, $nodePairs);


print "Finishing";
print $fhNPSout writeNodePairs($nodePairs);

# Return 1 (true)
1;

##		End Main		 ##

#####
# Top level subroutines
####

sub initNodePairs {
	print STDOUT "Initializing Node Pairs";
	my ($nodePairs) = @_;

	foreach my $nodePair (@$nodePairs) {
		$nodePair->{distance} = 0;
		$nodePair->{subs} = "";
		$nodePair->{sumBLs} = 0;
	}
	return $nodePairs;
}

sub readNodePairs {
	print STDOUT "Reading Node Pairs";
	my ($fh) = @_;

	my $nodePairs = [];

	<$fh>; # Burn header line

	while (my $line = <$fh>) {
		chomp $line;
		next if not $line;

		my $fields = [split " ", $line];

		my $nodePair = {
			firstID => $fields->[0],
			secondID => $fields->[1],
			MRCA => $fields->[2],
		};
		push @$nodePairs, $nodePair;
	}
	return $nodePairs;
}

sub updateNodePairs {
	print "Updating Node Pairs";
	my ($nodes, $nodePairs, $g) = @_;

	my $np=0;
	foreach my $branchPair (@$nodePairs) {
		$branchPair->{distance} +=
		  $nodes->{$branchPair->{firstID}}->{distanceToRoot} +
		  $nodes->{$branchPair->{secondID}}->{distanceToRoot} -
		  2 * $nodes->{$branchPair->{MRCA}}->{distanceToRoot};
		$branchPair->{sumBLs} +=
		  $nodes->{$branchPair->{firstID}}->{distance} +
		  $nodes->{$branchPair->{secondID}}->{distance};
	
		my $subs1 = $nodes->{$branchPair->{firstID}}->{subs};
		my $subs2 = $nodes->{$branchPair->{secondID}}->{subs};

		my $commonSites = [intersection([keys %$subs1], [keys %$subs2])];

		#print "Common sites@$commonSites";

		foreach my $site (@$commonSites) {
			$branchPair->{subs} .=
			    $subs1->{$site}->{from}
			  . $subs2->{$site}->{from}
			  . $site
			  . $subs1->{$site}->{to}
			  . $subs2->{$site}->{to} . ',';
		}
		$branchPair->{subs} .= ';';
	}
}

sub intersection {
	my ($a, $b) = @_;
	my %a = map { $_ => 1} @$a;
	return grep { $a{$_} } @$b;
}

sub updateSubs {
	my ($branches, $fh) = @_;
	return $fh if not $fh;
	while (my $line = <$fh>) {
		if ($line =~ /\/\//) { return $fh }
		chomp $line;
		my $branch = readBranch($line);
		$branches->{$branch->{to}}->{subs} = $branch->{subs};
	}
	die "End of substitution file reached. Burnin might be too high.";
}

sub readBranch {

	#	print STDOUT "reading branch";
	my ($string) = @_;

	#	print $string;exit;
	my $fields = [split " ", $string];

	#print "f:@$fields"; exit;
	(shift @$fields) =~ /(\d*)..(\d*)/ or die $!;

	my $branch = {from => $1, to => $2};
	foreach my $field (@$fields) {
		my $sub = readSub($field);
		$branch->{subs}->{$sub->{site}} = $sub;
	}

	return $branch;
}

sub writeBranch {
	my ($branch) = @_;
	my $string = "$branch->{from}..$branch->{to}";

	if (exists $branch->{subs}) {
		foreach my $site (sort { $a <=> $b } keys %{$branch->{subs}}) {
			$string .= "\n" . writeSub(
				$branch->{subs}->{$site}
			  );
		}
	}
	return $string;
}

sub readSub {

	#	print STDOUT "Reading Substitution";
	my ($string) = @_;

	$string =~ /([ACDEFGHIKLMNPQRSTVWY])(\d*)([ACDEFGHIKLMNPQRSTVWY])/;

	return newSub($2, $1, $3);
}

sub writeSub {
	my ($sub) = @_;
	return "$sub->{site}\t$sub->{from}\t$sub->{to}";
}

sub newSub {
	my ($site, $from, $to) = @_;

	# Check for good amino acid
	die "unknown aa found: $from" if not($from =~ /[ACDEFGHIKLMNPQRSTVWY]/);
	die "unknown aa found: $to" if not($to =~ /[ACDEFGHIKLMNPQRSTVWY]/);

	my $sub = {
		site => $site,
		from => $from,
		to => $to,
	};

	return $sub;
}

sub updateDistances {
	my ($nodes, $fh) = @_;
	my $line = <$fh>;
	if (not defined $line or $line =~ /^\s*$/){

		#		print "found last tree...";
		return $fh

		  # Don't do more work than I have to
		  #		seek $fh, 0, 0;
		  #		$line = <$fh>;
	}

	#	print $line;
	chomp $line;

	my $ids = [$line =~ /..(\d*):/g];
	my $dists = [$line =~ /:(\d*\.?\d*)/g];

	foreach (0 .. $#$ids) {
		$nodes->{$ids->[$_]}->{distance} = $dists->[$_];
	}

	return $fh;
}

sub initSubs {
	my ($file, $burnin) = @_;
	print STDOUT "Initializing Substitutions File";

	#	print $file;exit;
	open my $fh, "<", $file or die "cannot open substitution file";
	my $g = 0;
	my $line;
	<$fh>; # burn header line
	while ($g != $burnin and defined( $line = <$fh> )) {
		if ($line =~ /\/\//) {
			$g++;
		}

	}
	die "End of burnin not reached; check burnin" if $g != $burnin;
	return $fh;
}

# Reads one tree at a time
sub readNextTree {
	print STDOUT "Reading next tree";
	my ($fh) = @_;
	my $line = <$fh>;
	if (not defined $line or $line =~ /^\s*$/){

		#		print "found last tree...";
		seek $fh, 0, 0;
		$line = <$fh>;
	}
	chomp $line;
	my $nodes = parseTree($line);

	return $nodes, $fh;
}

sub initSeq {
	my ($file, $burnin) = @_;
	open(my $fh, "<$file") or die "cannot open sequence file";
	my $g = 0;
	my $line;
	while ($g != $burnin and $line = <$fh>) {
		if ($line =~ /\/\//) {
			$g++;
		}

	}
	die "End of burnin not reached; check burnin" if $g != $burnin;
	return $fh;
}

sub initTree {
	my ($file, $burnin) = @_;
	print STDOUT "Initializing Tree";
	open(my $fh, "<$file") or die "cannot open tree file";
	<$fh>; # burn first line
	my $line2 = <$fh>;
	if (not defined $line2 or $line2 =~ /^\s*$/) {

		#		print "found only one tree";
		seek $fh, 0, 0;
	}
	else {

		#		print "found multiple trees";
		seek $fh, 0, 0;
		for (my $t = 0; $t < $burnin; $t++) {
			<$fh>;
		}
	}
	return $fh;
}

sub removeMaskedSubs {
	my ($nodes, $mask) = @_;

	my ($s, $subs);
	foreach my $node (values %$nodes) {
		$subs = $node->{subs};
		$s = 0;
		while ($s <= $#$subs){
			if (not $mask->[$subs->[$s]]) { splice(@$subs,$s,1) }
			else { $s++; }
		}
	}
}

sub dbpause {
	if (my $message = shift){
		print $message;
	}
	else { print "Waiting for user to press Enter";}
	<STDIN>;
}

sub setupInfo {
	my ($info, $nodes) = @_;
	nameIdHashes($nodes, $info);
	$info->{rootId} = findRootID($nodes);
}

sub writeNodePairs {
	my ($nodePairs) = @_;

	#my $string = "FirstID\tSecondID\tSumBLs\tDistance\tSubs\n";
	my $string = "NodeID1\tNodeID2\tDistance\tSubstitutions\n";
	foreach my $nodePair (@$nodePairs) {
		$string .= $nodePair->{firstID} . "\t" . $nodePair->{secondID} . "\t"

		  #. $nodePair->{sumBLs} . "\t"
		  . $nodePair->{distance} . "\t" . $nodePair->{subs} . "\n";
	}
	return $string;
}

sub printNodePairs {
	my ($fh, $nodePairs) = @_;
	print $fh
"firstID\tsecondID\tMRCA\tdivergence\tconvergence\tpseudoconverge\tdistance";
	foreach my $nodePair (@$nodePairs) {
		print $fh $nodePair->{firstID}, "\t", $nodePair->{secondID}, "\t",
		  $nodePair->{MRCA}, "\t", $nodePair->{diverge}, "\t",
		  $nodePair->{converge},
		  "\t", $nodePair->{pseudoConverge}, "\t",  $nodePair->{distance};
	}
}

sub printSites {
	my ($fh, $entropies, $siteConv) = @_;
	print $fh "entropy\tconvergence";
	for (my $s = 0; $s <= $#$entropies; $s++) {
		print $fh $entropies->[$s], "\t", $siteConv->[$s];
	}
}

sub findRootID {
	my ($nodes) = @_;
	foreach my $node (values %$nodes) {
		if ($node->{isRoot}) {
			return $node->{id};
		}
	}
}

sub nameIdHashes {
	my ($nodes, $info) = @_;
	hashifempty($info, "idToName");
	hashifempty($info, "nameToId");
	foreach my $id ( keys %$nodes ) {
		$info->{idToName}->{$id} = $nodes->{$id}->{name};
		$info->{nameToId}->{ $nodes->{$id}->{name} } = $id;
	}
	makekeys($info, "nodeNames", $info->{nameToId});
	makekeys2($info, "nodeIds", $info->{idToName});
}

# Calculate total distance of tree
sub totalTreeDist {
	my ($nodes) = @_;
	my $dist = 0;
	foreach my $node (values %$nodes) {
		next if $node->{isRoot};
		$dist += $node->{distance};
	}
	return $dist;
}

sub normalize {
	my ($reference, $divisor) = @_;

	# 'Values' works on both hash references and array references
	my $values;
	if (ref $reference eq "ARRAY") { $values = values @$reference }
	if (ref $reference eq "HASH") { $values = values %$reference }
	foreach my $value ($values) {
		if (ref($value) eq ("ARRAY" or "HASH")) {
			normalize($value, $divisor);
		}
		else {
			$value /= $divisor;
		}
	}
}

sub newNodePairs {
	my ($nodes) = @_;
	print "making node pairs";
	my $nodePairs = [];
	my $nodeIds = [sort {$a <=> $b } keys %$nodes];

	foreach my $firstID (@$nodeIds){
		foreach my $secondID (@$nodeIds){
			next
			  if not checkPair($nodes, $firstID, $secondID)
			; # this removes about 1/2 pairs: 797353 vs 1630728 pairs

			my $nodePair = {

				# These values don't change. int() makes this slightly smaller
				firstID => $firstID,
				secondID => $secondID,
				MRCA =>
				  $nodes->findMRCA($nodes->{$firstID}, $nodes->{$secondID}),

				# These values must be updated every generation
				converge => 0,

				#				convergeSites => {},
				pseudoConverge => 0,

				#				pseudoConvergeSites => {},
				diverge => 0,
				total => 0,
				distance => 0,
			};

			push(@$nodePairs, $nodePair);
		}
	}
	print "Made ", $#$nodePairs+1, " node pairs";
	return $nodePairs;
}

sub checkPair {
	my ($nodes, $firstNode, $secondNode) = @_;
	return 0 if $firstNode == $secondNode;
	return 0 if $firstNode > $secondNode;

	# False if nodes are sisters, parents are the same
	# Don't use numeric compare because parent of root = ""
	return 0
	  if $nodes->{$firstNode}->{parent} eq $nodes->{$secondNode}->{parent};

	# False if one node is an ancestor of the other
	die "checkPair requires ancestors"
	  if not defined $nodes->{$firstNode}->{ancestors};
	die "checkPair requires ancestors"
	  if not defined $nodes->{$secondNode}->{ancestors};

	return 0 if $firstNode ~~ @{$nodes->{$secondNode}->{ancestors}};
	return 0 if $secondNode ~~ @{$nodes->{$firstNode}->{ancestors}};

	# True else
	return 1;
}

sub mean {
	my (@array) = @_;
	return sum(@array) / (
		$#array + 1
	  );
}

# I might rather use 'map' or 'do' than a for loop.
sub sum {
	my (@array) = @_;
	my $sum = 0;
	foreach my $element (@array) { $sum += $element }
	return $sum;
}

sub nodePairDistances {
	my ($nodePairs, $nodes) = @_;
	foreach my $nodePair (@$nodePairs) {
		$nodePair->{distance} =
		  $nodes->distanceBetween($nodes->{$nodePair->{firstID}},
			$nodes->{$nodePair->{secondID}});
	}
}

##   subroutines 	##

## Stephen's subroutines ##
#Subroutines
# Notes
# I want to separate input from output of a function. When the funtion only
# takes arguments and does not return anything, inputs and outputs are
# ambiguous.

# Constructor for array of generations
sub newGenerations {
	my ($nGens) = @_;
	my $gens = [];
	for (my $g = 0; $g < $nGens; $g++){
		$gens->[$g] = {};
		$gens->[$g]->{id} = $g;
	}
	return $gens;
}

# Attaches substitutions to gen
sub findAllSubs {
	my ($nodes, $gen) = @_;
	my ($child, $childId, $subs, $parent, $states1, $states2);
	hashifempty($gen, "nodes");
	foreach $child (values %$nodes) {
		$childId = $child->{id};
		hashifempty($gen->{nodes}, $childId);
		$subs = arrayifempty($gen->{nodes}->{$childId}, "subs");
		next if ($child->{isRoot});
		$states1 = $child->{states};
		$states2 = $nodes->{$child->{parent}}->{states};
		for (my $s=0; $s <= $#$states1; $s++){
			if ($states1->[$s] ne $states2->[$s]) {
				push(@$subs, $s);
			}
		}
	}
}

## Print subroutines
sub printGenerations {
	my ($fh, $generations, $format) = @_;
	print "printing generations";
	foreach my $generation (@{$generations}) {
		print $fh "Generation $generation->{'id'} exists";
		print $fh "//";
	}
}

# Print all sequences in all generations
sub printSequences {
	my ($fh, $generations, $format) = @_;
	foreach my $generation (@{$generations}) {
		print $fh "\*Generation $generation->{id}";
		foreach my $node (@{$generation->{nodes}}){
			next if not defined $node;
			print $fh ">$node->{name}\n$node->{state}";
		}
		print $fh "//";
	}
}

sub printStates {
	my ($fh, $generations, $format) = @_;
	foreach my $generation (@{$generations}) {
		print $fh "\*Generation $generation->{'id'}";
		my $states = $generation->{"states"};
		foreach my $state (keys %{$states}){
			print $fh ">$state";
			print $fh @{$states->{$state}};
		}
		print $fh "//";
	}
}

# Get branch info (id, etc), save to branches
sub parseTree {
	my ($line) = @_;

	my $nodes = NodeHash->new();

	my $list = [split(/[(),;]+/, $line)];
	shift @$list;

	foreach (@$list) {
		$_ =~ /(\w*) (\d*)..(\d*):(\d*\.?\d*)/;
		my $name = $1;
		my $parent = $2;
		my $id = $3;
		my $distance = $4;
		my $isTip = $1 ? 1 : 0;
		$nodes->addNode(Node->new($id, $name, $isTip, $parent, $distance, 0));
	}

	foreach my $node (keys %$nodes){
		if (ref($nodes->{$node}) eq "HASH") {
			$nodes->addNode(Node->new($node, "Node_" . $node, 0, '', 0, 1));
		}
	}
	return $nodes;
}

## End Stephen's subroutines ##

sub splitifempty{
	my ($inhash,$inkey,$kmer)=@_;
	if (!$inhash->{$inkey}) {
		$inhash->{$inkey} = [];
		$inhash->{kmer}=$kmer;
		for(my $i=0; $i<length $kmer;$i++){
			$inhash->{$inkey}->[$i] = substr($kmer,$i,1);
		}
	}
}

sub makekeys{
	my ($target,$key,$source) = @_;
	arrayifempty($target,$key);
	@{$target->{$key}} = sort keys %$source;
}

sub makekeys2{
	my ($target,$key,$source) = @_;
	arrayifempty($target,$key);
	@{$target->{$key}} = sort {$a <=> $b} keys %$source;
}

sub arrayifempty {
	my ($inhash, $inkey) = @_;
	if (!$inhash->{$inkey}) { $inhash->{$inkey} = []; }
	return $inhash->{$inkey};
}

sub hashifempty {
	my ($inhash, $inkey) = @_;
	if (!$inhash->{$inkey}) { $inhash->{$inkey} = {}; }
	return $inhash->{$inkey};
}

# Read the sequences nodes
# This could be robust against different end of line characters (windows, unix, osx)

# Really all I need this function is update the states of the nodes. I should
# combine reading the sequences with translating them to states.
sub readNextFastaPlex {
	my ($fh) = @_;
	my $name;
	my $line;
	my $seqs = {};
	while($line = <$fh>)	{
		chomp $line;
		if ($line =~ /\/\// )  { return ($seqs, $fh); }
		elsif ($line =~ />/ )  { $name = getName($line)}
		else {$seqs->{$name} .= uc($line);}
	}
	die "end of seq file reached; burnin might be too large";
}

# Extract the name from a fasta sequence
sub getName()	{
	my ($name, $format) = @_;
	$name =~ s/^>//;
	my @words;
	if (!$format)  { (@words) = split(/ /, $name);}
	elsif ($format eq "lava")  { (@words) = split(/\s+/, $name);}
	$name = $words[0];
	$name =~ s/ //g;
	chomp($name);
	return $name;
}
