# Statistics functions


use strict;
use warnings FATAL => "all";
use diagnostics;

use List::Util qw( sum );

unless (caller) {
	
}


sub average {
	my ($xs) = @_;

	return sum(@$xs) / scalar @$xs
}

sub variance {
	my ($xs) = @_;
	my $mean = average($xs);
	# could use reduce here instead
	# my $sum = reduce {$a + ($b - $mean)**2}, 0, @$xs;
	my $sum = 0;
	for my $i (@$xs) {
		$sum += ($i - $mean)**2;
	}
	return $sum/scalar(@$xs)
}

sub sd {
	my ($xs) = @_;

	return sqrt(variance($xs));
}

sub variance_unbiased {
    my ($xs) = @_;
	my $mean = average($xs);
	# could use reduce here instead
	# my $sum = reduce {$a + ($b - $mean)**2}, 0, @$xs;
	my $sum = 0;
	for my $i (@$xs) {
		$sum += ($i - $mean)**2;
	}
	return $sum/(scalar(@$xs)-1)
}

sub sd_unbiased {
	my ($xs) = @_;

	return sqrt(variance_unbiased($xs))
}

1;

