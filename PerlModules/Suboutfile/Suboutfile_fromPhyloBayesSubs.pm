# Convert a phylobayes substitution file into a PLEX suboutfile.

package Suboutfile;

use strict;
use warnings;

use Data::Dumper;

use Tree;

open my $tree_out, ">", "tree_out";

unless (caller) {
	my $phyloBayesSubs = $ARGV[0]
	  // "../AC_cCDS_10heuristic.trimal_sample.sub";
	my $suboutfile = $ARGV[2] // ($phyloBayesSubs . ".suboutfile");

	FromPhyloBayesSubs($phyloBayesSubs, $suboutfile);
}

sub ToSuboutfile {
	my ($substitutions, $filename) = @_;
	open my $fh, ">", $filename;

	print $fh "Name\tSubstitutions\n";
	while (my ($name, $subs) = each %$substitutions) {
		print join("\t", $name, @$subs) . "\n";
	}
}

sub FromPhyloBayesSubs {
	my ($phyloBayesSubs, $suboutfile) = @_;

	open my $fh_in, "<",$phyloBayesSubs or die $!;
	open my $fh_out, ">", $suboutfile;

	print $fh_out "Branch\tSubstitutions\n";
	my $substitutions = {};

	while (<$fh_in>) {
		if ($_ =~ /^point (\d+)/) {
			print "point $1\n";
			next if ($_ =~ /^point 1$/); # Don't print for first generation

			#			print Dumper $substitutions;
			while (my ($name, $subs) = each %$substitutions) {
				print $fh_out join("\t", $name, @$subs) . "\n";
			}
			print $fh_out "//\n";
			$substitutions = {};
			next;
		}
		next if $_ =~ /^rep/;
		chomp;
		next if not $_; # Skip empty lines.

		my $fields = [split];
		my $site = $fields->[0];
		my $tree = $fields->[1];
		UpdateSubstitutionsFromTree($substitutions, $site, $tree);
	}

	# Print last generation
	while (my ($name, $subs) = each %$substitutions) {
		print $fh_out join("\t", $name, @$subs) . "\n";
	}
	print $fh_out "//\n";
}

sub UpdateSubstitutionsFromTree {
	my ($substitutions, $site, $tree_string) = @_;

	print $tree_out "Original tree: $tree_string\n";
	my $tree_string_without_subs =
	  Tree::StripPhylobayesSubstitutions($tree_string);
	print $tree_out "No subs tree: $tree_string_without_subs\n";

	# The next line will generate lots of warnings about branches without
	# lengths
	my $tree = Tree::FromString($tree_string_without_subs);

	# The internals nodes don't have any names
	$tree = Tree::GenerateNamesForInternalNodes_orderIndependent($tree);

	$tree_string = Tree::ApplyInternalNamesToString($Tree::names, $tree_string);

	# Have to handle the root
	$tree_string =~ s/([A-Z]);/$1:0:$1;/;
	print $tree_out "With subs tree: $tree_string\n";
	my $subs = [];

	# (Cfol_v3_20843_31-35kDa_M:0:M,(Cf_Cfol_v3_20842_M:0.108931:M,(((
	@$subs = ($tree_string =~ /([^\(\),;]+)/g);

	for my $sub (@$subs) {

		# Notice this ignores multiple substitutions and only considers
		# the node above and the node below.
		$sub =~ /^([^:]+)_([A-Z]):.*:([A-Z])$/;

		my $name = $1;
		my $to = $2;
		my $from = $3;

		if ($from ne $to) {
			if (defined $substitutions->{$name}
				and $substitutions->{$name}->[$#{$substitutions->{$name}}] eq
				$from.$site.$to)
			{
				die "Found a duplicate!\n"
				  . "name: $name\n"
				  .$from
				  .$site
				  .$to."\n"
				  . $sub ."\n";

			}
			push @{$substitutions->{$name}}, $from.$site.$to;
		}
	}
}
