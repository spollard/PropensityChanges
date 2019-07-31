# Convert suboutfile to a R dataframe format holding all substitutions

use strict; 
use warnings;
use 5.010.001;


local $\ = "\n";


my $suboutfile = "suboutfile";

mkdir "subs" if not -d "subs";

my $fh = initSubs($suboutfile, 250);

my $branches = {};
foreach (1 .. 100) {
	print "working on generation $_";
	open my $fh_out, ">subs/" . $_ or die $!;
	readWriteSubs($fh, $fh_out);
}

sub readWriteSubs {
	my ($fh, $fh_out) = @_;
	return $fh if not $fh;
	local $\ = "";
	print $fh_out "parent\tchild\tsite\tAA_parent\tAA_child\n";
	while (my $line = <$fh>) {
		if ($line =~ /\/\//) { return $fh } 
		if ($line =~ "parent") { next }
		chomp $line;
		print $fh_out writeBranch(readBranch($line))
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

sub writeAllSubs { 
	my ($branches) = @_;
	
	my $string = "";
	
	foreach my $branch (values %$branches) {
		$string .= writeBranch($branch); 	
	}
	
	
	return $string;
}

sub writeBranch { 
	my ($branch) = @_;
	my $bs = "$branch->{from}\t$branch->{to}\t";
	my $string = ""; 
	
	if (exists $branch->{subs}) {
		foreach my $site (sort { $a <=> $b } keys %{$branch->{subs}}) {
			$string .= $bs . writeSub($branch->{subs}->{$site}) . "\n";
		} 
	}
#	chop $string; # removed last "\n"
	return $string
}


sub updateSubs {
	my ($branches, $fh) = @_;
	return $fh if not $fh;
	local $\ = "";
	while (my $line = <$fh>) {
		if ($line =~ /\/\//) { return $fh } 
		chomp $line;
		my $branch = readBranch($line);
		$branches->{$branch->{to}} = $branch;
	}
}




sub readBranch { 
	my ($string) = @_;
	my $fields = [split " ", $string];
	
	(shift @$fields) =~ /(\d*)..(\d*)/;
	
	my $branch = {from => $1, to => $2};
	
	foreach my $field (@$fields) {
		my $sub = readSub($field);
		$branch->{subs}->{$sub->{site}} = $sub;
	}
	
	return $branch	
}


sub readSub { 
	my ($string) = @_;
	
	$string =~ /([ACDEFGHIKLMNPQRSTVWY])(\d*)([ACDEFGHIKLMNPQRSTVWY])/;

	return newSub($2, $1, $3);
}

sub writeSub {
	my ($sub) = @_;
	return "$sub->{site}\t$sub->{from}\t$sub->{to}"
}

sub newSub { 
	my ($site, $from, $to) = @_;
		# Check for good amino acid
	die "unknown aa found: $from" if not ($from =~ /[ACDEFGHIKLMNPQRSTVWY]/);
	die "unknown aa found: $to" if not ($to =~ /[ACDEFGHIKLMNPQRSTVWY]/);
	
	my $sub = {
		site => $site,
		from => $from,
		to => $to,
	};
	
	return $sub;
}

