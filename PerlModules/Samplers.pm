# I want my samplers to be separate from my MCMCs. This is the script that
# holds the samplers. Right now it is simply a sampler that draws from a normal
# distribution.

package Samplers;

# testing
unless (caller) {
    my $output_dir = "Samplers/";

}

use Math::Random::OO::Normal;

use List::Util 'sum', 'min';

my $rng = Math::Random::OO::Normal->new();

# Pass in the old value and the sigma for the distribution as arguments.
# It returns the new value.
sub normal {
    if ( @_ > 1 ) {
        $rng->mean( $_[0] );
        $rng->stdev( abs( $_[1] ) );
    }
    elsif ( @_ == 1 ) {
        $rng->mean( $_[0] );
        $rng->stdev(1);
    }
    else {
        $rng->mean(0);
        $rng->stdev(1);
    }
    return $rng->next();
}

sub normal_p {
    if ( @_ > 1 ) {
        $rng->mean( $_[0] );
        $rng->stdev( abs( $_[1] ) );
    }
    elsif ( @_ == 1 ) {
        $rng->mean( $_[0] );
        $rng->stdev(1);
    }
    else {
        $rng->mean(0);
        $rng->stdev(1);
    }
    my $prob = 0;
    do { $prob = $rng->next()}while ($prob <= 0 or $prob >= 1);
    return $prob;
}

sub normal_positive {
    if ( @_ > 1 ) {
        $rng->mean( $_[0] );
        $rng->stdev( abs( $_[1] ) );
    }
    elsif ( @_ == 1 ) {
        $rng->mean( $_[0] );
        $rng->stdev(1);
    }
    else {
        $rng->mean(0);
        $rng->stdev(1);
    }
    my $prob = 0;
    do { $prob = $rng->next()}while ($prob <= 0);
    return $prob;
}

sub random_step {
    my $old_value = $_[0] // 0;
    my $step_size = $_[1] // 1;

    # return rand < 0.5 ? $old_value + $step_size : $old_value - $step_size

    my $new_value = $old_value;
    if (rand() < 0.5) {
        $new_value += $step_size;
    }
    else{
        $new_value -= $step_size;
    }

    return  $new_value;
}

sub random_choice {
    my ($choices) = @_;

    return $choices->[int(rand(scalar(@$choices)))];
}

sub weighted_random_choice {
    my ($weights) = @_;

    # cmf = cumulative mass function
    my $accumulation = 0;
    my $weights_cmf = [map {$accumulation += $_; $accumulation} @$weights];
    die "Weights sum to zero @$weights" if not $accumulation;
    my $r = rand(1) * $accumulation;

    while (my ($index, $cumulative_weight) = each @$weights_cmf) {
        if ($r < $cumulative_weight) {
            return $index;
        }
    }

    die "weighted_random_choice broke";
}

sub log_weighted_random_choice {
    my ($log_weights) = @_;

    # Add the minimum log likelihood to all log likelihoods
    # exponentiate
    # use weighted_random_choice
    my $min = min(@$log_weights);
    my $weights = [map {exp($_ - $min)} @$log_weights];
    return weighted_random_choice($weights);

    die "log_weighted_random_choice broke";
}

sub binomial_step {
    my ($first, $second, $max_step_size) = @_;
    my $total = $first + $second;

    # max step size is a fraction of the total, not a fixed step size
    my $max_step = $total * $max_step_size;
    my $step = rand(2 * $max_step) - $max_step;

    if (($first-$step) < 0 or ($first-$step) > 1 or  ($second+$step) < 0 or ($second+$step) > 1) {
        return $first, $second;
    }
    else {
        return $first-$step, $second+$step;
    }
}

1;
