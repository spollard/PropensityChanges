
package Substitutions;

use strict;
use warnings;

use Data::Dumper;

unless (caller) {
	my $suboutfile ="name1\tA1C\tC2G\tG3T\tT4A\n";
	my $substitutions = {
		name1 => [qw(A1C C2G G3T T4A)],
		name2 => [qw(A5C C6G G7T T8A)],
	};
	$substitutions = KeepSitesBetweenRange($substitutions, 3, 6);

	print Dumper($substitutions);
}

sub FromFile {
	my ($suboutfile) = @_;

	open my $fh, $suboutfile or die $!;
	<$fh>; # Burn header line
	my $substitutions = {};
	while (<$fh>) {
		chomp;
		if ($_ =~ /\/\//) {
			warn "Give me a suboutfile with a single "
			  . "generation, not the whole suboutfile.";
			last;
		}
		my $fields = [split];
		my $name = shift @$fields;
		$substitutions->{$name} = $fields;
	}
	return $substitutions;
}

sub FromSuboutfile {
	return FromFile(@_);	
}

sub ToFile {
	my ($substitutions, $filename) = @_;
	open my $fh, ">", $filename;

	print $fh "Name\tSubstitutions\n";
	while (my ($name, $subs) = each %$substitutions) {
		print $fh join("\t", $name, @$subs) . "\n";
	}
}

sub ToSuboutfile {
	ToFile(@_);	
}

sub FromPhyloBayesFile {
	my ($filename) = @_;
	die "Substitutions::FromPhyloBayesFile does not work yet";
	open my $fh, $filename or die $!;

	<$fh>; # Burn header line
	my $substitutions = {};
	while (<$fh>) {
		if ($_ =~ /^point/) {
			warn "Give me a phylobayes substitution file with a single "
			  . "generation, not the whole file.";
			last;
		}
		next if $_ =~ /^rev/;

	}
}

sub ManyGensFromSuboutfile {
	my ($suboutfile, $burnin, $gens) = @_;

	sub initSubs {
		my ($file, $burnin) = @_;
		open my $fh, "<", $file or die $! . $file;
		my $g = 0;
		my $line;
		<$fh>; # burn header line
		while ($g != $burnin and defined($line = <$fh>)) {
			if ($line =~ /\/\//) {
				$g++;
			}

		}
		die "End of burnin not reached; check burnin" if $g != $burnin;
		return $fh;
	}
	my $fh = initSubs($suboutfile, $burnin);

	scalar <$fh>; # Burn header line
	my $substitutions = {};
	my $g = $burnin;
	while (<$fh>) {
		chomp;
		if ($_ =~ /\/\//) {
			$g++;
			if ($g >= $burnin + $gens){
				last;
			}
			else {
				next;
			}

		}
		my $fields = [split];
		my $name = shift @$fields;
		push @{$substitutions->{$name}}, $fields;
	}
	return $substitutions;
}

sub Combine {
	my ($substitutions_1, $substitutions_2) = @_;
	
	while (my ($name, $subs) = each (%$substitutions_2)) {
		die if not defined $substitutions_1->{$name};
		push @{$substitutions_1->{$name}}, @{$substitutions_2->{$name}}; 
	}
	
	return $substitutions_1;	
}

# Legacy interface
sub SplitSubstitutions {
	Split(@_);
}

sub Split {
	my ($substitutions) = @_;

	for my $subs (values %$substitutions) {
		my $subs_hash = {};
		for my $sub (@$subs) {
			$subs_hash->{substr $sub,1,-1} =
			  { from => (substr $sub, 0, 1), to => substr($sub, -1)};
		}
		$subs = $subs_hash;
	}
}

# Legacy interface
sub Join {
	my ($substitutions) = @_;

	for my $subs (values %$substitutions) {
		my $subs_array = [];
		while (my ($site, $sub) = each %$subs) {
			push @$subs_array, $sub->{from} . $site . $sub->{to};
		}
		$subs = $subs_array;
	}
}

sub JoinSubstitutions{
	Join(@_);
}

sub CountNumberOfSubsPerSite {
	my ($suboutfile) = @_;
	open my $fh, "<", $suboutfile or die $!;
	<$fh>; # Burn header line
	my $sites = [];
	while (<$fh>) {
		chomp;
		if ($_ =~ /\/\//) {
			warn "Give me a suboutfile with a single "
			  . "generation, not the whole suboutfile.";
			last;
		}
		my $fields = [split];

		# Get rid of the name
		my $name = shift @$fields;
		my $substitutions = $fields;
		for my $sub (@$substitutions) {
			$sites->[substr $sub, 1, length($sub)-2]++;
		}
	}
	return $sites;
}

sub PrintNumberOfSubsPerSite {
	my ($sites, $sitesFile) = @_;
	open my $fh_out, ">", $sitesFile or die $!;
	print $fh_out "Site\tNumber_of_substitutions\n";

	while (my ($site, $number_of_subs) = each @$sites) {
		print $fh_out $site, "\t", ($number_of_subs // 0), "\n";
	}
}

sub KeepSitesBetweenRange {

	# the ends of the range are inclusive so the ends are kept
	# Sites are one indexed

	my ($substitutions, $lower_site, $upper_site) = @_;

	for my $subs (values %$substitutions) {
		$subs = [
			grep {
				my $site = substr $_, 1, -1;
				$site >= $lower_site and $site <= $upper_site
			  } @$subs
		];
	}

	return $substitutions;
}
