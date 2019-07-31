# Distributions

my $factorials = [1];
for (1 .. 10) {
	$factorials->[$_] = $factorials->[$_ - 1] * $_;
}

my $log_nCk = [];
for my $n (0 .. 10) {
	for my $k (0 .. 10) {
		$log_nCk->[$n]->[$k] =
		  log($factorials->[$n]) - log($factorials->[$k]) - log($factorials->[$n-$k]);
	}
}

# log_pbinom(k, n, p)
sub log_pbinom {
	my ($k, $n, $p) = @_;
	die "Probability is out of bounds: $p" if $p <= 0 or $p >= 1; 
	
	my $log_likelihood = 0;
	if ($n <= 10) {
		# use binomial equation
		#			print "using Binomial for $_\n";
		$log_likelihood +=
		  $log_nCk->[$n]->[$k] +
		  $k * log($p) +
		  ($n - $k) * log(1 - $p);
	}
	else {

		# Use normal approximation
		my $variance = $n * $p * (1 - $p);
		$log_likelihood += -($k - $n * $p)**2/(2 * $variance);
		$log_likelihood += -log(sqrt($variance));

		# I'm ignoring the 1/(sqrt(2)*pi) because all points will have that
		# contribution and it will be the same for every model and so it
		# will divide (or subtract when we're using logs) out of the
		# likelihood ratio.
	}
	return $log_likelihood;
}

1;