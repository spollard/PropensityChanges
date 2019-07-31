# Simulates a pair of propensities and s the amount of convergence
# between the two.

package Propensities;

use List::Util 'sum';

#require "samplers.pl";
use Samplers;

unless (caller) {
#	use Math::Random 'random_gamma';
	Test();
	exit;

	my $prop1 = [map {1/4} 1..4];
	my $prop2 = [@$prop1];

	Normalize($prop1);
	print("@$prop1");

	print("@{Mutate($prop1, 0.1)}");

	print("Constraint:" . Constraint($prop1) . "\n");
	print("CD:" . CD($prop1, $prop2) . "\n");
	print("Distance:" . Distance($prop1, $prop2) . "\n");
}

sub Test {
	my $output_dir = "Propensities/";
	mkdir $output_dir if not -d $output_dir;

	my $propensities = GenerateRandomPropensities();
	print "Starting propensities @$propensities\n";

	# test Mutate
#	open my $normal_p_out,">", $output_dir . "normal_p";
#	my $normal_p_propensities = [@$propensities];
#
#		for (1 .. 100000) {
#			Mutate2($normal_p_propensities, 1/100);
#			print $normal_p_out "@$normal_p_propensities\n";
#		}

	# test Mutate2
#	open my $binom_step_out, ">", $output_dir . "binom_step";
#	my $binom_step_propensities = [@$propensities];
#	for (1 .. 100000) {
#		Mutate($binom_step_propensities, 0.5);
#		print $binom_step_out "@$binom_step_propensities\n";
#	}
	
	
#	open my $random2, ">", $output_dir . "random";
#	
#	for (1 .. 100000) {
#		print $random2 join(" ", @{GenerateRandomPropensities(4)}) . "\n";
#	}
}

sub CD {
	my ($prop1, $prop2) = @_;

	die "Lengths are not equal" if not $#$prop1 == $#$prop2;

	my $prob_conv = sum(map {$prop1->[$_] * $prop2->[$_]} 0..$#$prop1);

	my $CD = $prob_conv / (1-$prob_conv);
	return $CD;
}

sub Constraint {
	my ($prop) = @_;

	# ternary operator required because cannot  log of
	# propensities with value zero
	my $constraint = exp( -sum(map {$_ ? $_ * log $_ : 0} @$prop));

  #to which a site is constrained in its amino acid composition was determined
  #	by using A , the effective size of the alphabet of amino acids possible at a
  #	 given location. This is defined as the exponential of the sequence entropy
  #	 at this location.

	# A = exp ( - sum over x (pi_x ln pi_x))

	return $constraint;
}

sub Distance {

	# Euclidian distance
	my ($prop1, $prop2) = @_;
	die "Lengths are not equal" if not $#$prop1 == $#$prop2;

	my $distance =sqrt(sum(map {($prop1->[$_] - $prop2->[$_])**2} 0..$#$prop1));

	return $distance;

}

# This DOES NOT produce dirichlet distributed propensities
sub Mutate2 {
	my ($prop, $sampling_variance) = @_;

	my $new_prop = [map {Samplers::normal_p($_, $sampling_variance)} @$prop];
	Normalize($new_prop);
	@$prop = @$new_prop;
}

# This DOES produce dirichlet distributed propensities
sub Mutate {
	my ($prop, $max_step_size) = @_;

	my $first = Samplers::random_choice([0 .. $#$prop]);
	my $second;
	do { $second = Samplers::random_choice([0 .. $#$prop])}
	  until ( $second != $first);

	($prop->[$first], $prop->[$second]) =
	  Samplers::binomial_step($prop->[$first], $prop->[$second], $max_step_size);
}

sub Normalize {
	my ($prop) = @_;

	my $sum = sum @$prop;
	for (values @$prop) {$_ = $_ / $sum }
}

sub MutateKeepConstraint {
	my ($prop, $sampling_variance, $constraint_wobble) = @_;
	die "MutateKeepConstraint does not work yet.";
	my $constraint = Constraint($prop);

	my $new_constraint; # must be declared outside the loop for scoping
	my $new_prop;
	do {

		#copy old props
		@$new_prop = @$prop;
		Mutate($new_prop, $sampling_variance);
		$new_constraint = Constraint($new_prop);
		print(  "old constraint:"
			  . $constraint
			  . "\n proposed constraint: "
			  . $new_constraint
			  . "\n");
	  } while $new_constraint < $constraint - $constraint_wobble
		  or $new_constraint > $constraint + $constraint_wobble;

	@$prop = @$new_prop;
}

sub log10 {
	my $n = shift;
	return log($n)/log(10);
}

# This does not generate Dirichlet distributed random propensities
# but it probably is random enough for most use cases 
sub GenerateRandomPropensities {
	my ($number_of_states) = @_;
	$number_of_states = $number_of_states // 4;

	my $prop = [map {rand} 1..$number_of_states];
	Propensities::Normalize($prop);
	return $prop;
}

# This does not generate Dirichlet distributed random propensities
sub GenerateRandomPropensities2 {
	my ($number_of_states) = @_;
	$number_of_states = $number_of_states // 4;
	my $prop = [];
	my $sum = 0;
	for (1 .. ($number_of_states-1)) {
		my $num = rand(1 - $sum);
		$sum += $num;
		push @$prop, $num;
	}
	push @$prop, (1-$sum);

	return $prop;
}

1;

