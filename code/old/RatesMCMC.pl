# Uses all the instances of a certain amino acid and tries to
# find classes
#
#
# Stephen Pollard
# 5/11/2015

use strict;
use warnings FATAL => 'all';

use List::Util 'sum';
use List::MoreUtils 'any';

use Data::Dumper;

use Propensities;
use Samplers;

require "MCMC.pl";
require "UniqueDirectory.pl";


# Model options
my $amino_acids_considered = ["S", "I", "M", "A"];
my $subs_minimum = 1;
my $number_of_classes = 4;


my $out_directory = "../data/$number_of_classes"."classes_$subs_minimum"."subs/";
#$out_directory = UniqueDirectoryAppendNumber($out_directory);
mkdir $out_directory if not -d $out_directory;

# Files
my $instances_filename = "../data/AllInstancesExtraDataOfT";
my $mcmc_out_filename = $out_directory . "mcmc";
my $current_model_out_filename = $out_directory . "current_model";
my $proposed_model_out_filename = $out_directory . "proposed_model";
my $class_assignments_out_filename = $out_directory . "class_assignments";


# MCMC options
my $generations = 10000;

# Sampling options
# Perhaps number of propensities to sample per generation
#my $maximum_means_to_sample_per_generation = $number_of_windows/4;
my $max_step_size = 1 / 100;
my $number_of_parameters_sampled = 0;
my $propensities_sampling_frequency = 1;
my $class_assignments_sampling_frequency = 0.3;

# data = [distance, C/D]
print "Reading Data\n";
my $data = ReadData($instances_filename, $out_directory . "used.instances");

my $number_of_points = scalar @$data;
print "Finished reading data. Number of points: $number_of_points\n";

# The following code departs from the standard MCMC code because of the
# class assignments. Class assignments are Gibbs sampled, not Metropolis-
# Hastings sampled.
# Initialize output files
my $mcmc_out = InitializeMcmcOutput($mcmc_out_filename);
my $current_model_out = InitializeModelOutput($current_model_out_filename);
my $proposed_model_out = InitializeModelOutput($proposed_model_out_filename);

my $current_model = GenerateRandomModel($number_of_points, $number_of_classes);

my $current_log_likelihood =
  LogLikelihoodOfModelGivenData($current_model, $data);

print "Beginning MCMC loop\n";

my $time_between_prints = 10;
my $next_time = time + $time_between_prints;
my $last_printed_gen = 0;
for my $generation (1 .. $generations) {

	if (time > $next_time) {
		UpdateScreen($generation);
	}

	RecordModelState($current_model, $current_model_out);
	$number_of_parameters_sampled = 0;
	if (rand() < $propensities_sampling_frequency) {
		my $proposed_model = GenerateNewModelFromCurrent($current_model);
		RecordModelState($proposed_model, $proposed_model_out);

		my $proposed_log_likelihood =
		  LogLikelihoodOfModelGivenData($proposed_model,$data);

		my $is_accepted =
		  IsProposedLogLikelihoodAccepted($current_log_likelihood,
			$proposed_log_likelihood);
		RecordMcmcState($mcmc_out, $generation, $current_log_likelihood,
			$proposed_log_likelihood, $is_accepted,
			$number_of_parameters_sampled);
		if ($is_accepted) {
			$current_model = $proposed_model;
			$current_log_likelihood = $proposed_log_likelihood;
		}
	}

	if (rand() < $class_assignments_sampling_frequency) {
		UpdateClassAssignments($data, $current_model);
		$current_log_likelihood =
		  LogLikelihoodOfModelGivenData($current_model,$data);
	}

}
print "MCMC Finished\n";

sub InitializeModelOutput {
	my ($model_out_filename) = @_;
	open my $model_out, ">", $model_out_filename;
	my $model_header = [];

	for my $class_number (1 .. $number_of_classes){
		push @$model_header, map {$_ . $class_number} @$amino_acids_considered;
	}

	push @$model_header, (1 .. $number_of_points);

	$model_header = join "\t", @$model_header;

	print $model_out $model_header . "\n";

	return $model_out;
}

sub RecordModelState {
	my ($model, $fh) = @_;
	my $flat_model = [];
	for my $class (@{$model->{classes}}) {
		push @$flat_model, @$class;
	}

	push @$flat_model, @{$model->{assignments}};
	print $fh join("\t", @$flat_model) . "\n";
}

sub LogLikelihoodOfModelGivenData {
	my ($model, $data) = @_;

	my $log_likelihood = 0;
	while (my ($frequency_index, $frequencies) = each @$data) {
		$log_likelihood +=
		  LogLikelihoodOfFrequenciesGivenPropensities($frequencies,
			$model->{classes}->[$model->{assignments}->[$frequency_index]]);
	}

	return $log_likelihood;
}

sub GenerateNewModelFromCurrent {
	my ($current_model) = @_;

	my $new_model = {
		classes=>[],
		# Copying the reference to the assignment is ok here because the 
		# assignments are not changed for this part of the sampling
		assignments=>$current_model->{assignments},
	};

	for (@{$current_model->{classes}}) {
		my $new_class = [@$_];
		Propensities::Mutate($new_class, $max_step_size);
		push @{$new_model->{classes}}, $new_class;

		$number_of_parameters_sampled += scalar @$amino_acids_considered;
	}

	return $new_model;
}

sub GenerateRandomModel {
	my ($number_of_points, $number_of_classes) = @_;
	my $model = {};
	for (1 .. $number_of_classes) {
		push @{$model->{classes}},
		  Propensities::GenerateRandomPropensities(
			scalar @$amino_acids_considered);
	}

	$model->{assignments} =
	  GenerateRandomClassAssignments($number_of_points, $number_of_classes);

	return $model;

}

sub GenerateRandomClassAssignments{
	my ($number_of_points, $number_of_classes) = @_;

	my $class_assignments = [];
	for (1 .. $number_of_points) {
		push @$class_assignments,
		  Samplers::random_choice([0 .. ($number_of_classes-1)]);
	}

	return $class_assignments;
}

sub UpdateClassAssignments {
	my ($data, $model) = @_;

	$model->{assignments} = [];
	for my $frequencies (@$data) {
		my $log_likelihoods = [];
		for my $propensities (@{$model->{classes}}) {
			push @$log_likelihoods,
			  LogLikelihoodOfFrequenciesGivenPropensities($frequencies,
				$propensities);
		}
		push @{$model->{assignments}},
		  Samplers::log_weighted_random_choice($log_likelihoods);
	}
}

sub LogLikelihoodOfFrequenciesGivenPropensities {
	my ($frequencies, $propensities) = @_;

	my $log_likelihood = 0;
	for (0 .. $#$frequencies) {
		if ($propensities->[$_]){
			$log_likelihood +=$frequencies->[$_] * log($propensities->[$_]);
		}
	}
	return $log_likelihood;
}

sub ReadData {
	my ($filename, $filename_out) = @_;
	open my $fh, "<", $filename or die $!;
	open my $fh_out, ">", $filename_out or die $!;
	my $data = [];

	print $fh_out scalar <$fh>; # Burn header line
	while (<$fh>) {
		chomp;
		my $fields = [split ' ', $_];

		my $frequencies = {};
		for my $residue (split '', $fields->[4]) {
			$frequencies->{$residue}++;
		}

		# Hash slice to get the frequencies in order
		$frequencies = [@$frequencies{@$amino_acids_considered}];
		if ( any {$_} @$frequencies) {

			$frequencies = [map {$_ // 0} @$frequencies];
			if (sum(@$frequencies) >= $subs_minimum){
				push $data,$frequencies;

				print $fh_out $_."\n";

			}
		}
	}

	return $data;
}

sub UpdateScreen {
	my ($generation) = @_;
	my $rate = ($generation - $last_printed_gen )  / $time_between_prints;
	print "Generation $generation\nRate "
	  . sprintf("%.1f",$rate)
	  . " generations per second\n";
	$next_time = time + $time_between_prints;
	$last_printed_gen = $generation;
	my $time_left = ($generations - $generation) / $rate;
	my $hours = int($time_left / 3600);
	my $mins = int(($time_left % 3600) / 60);
	my $seconds = ($time_left % 3600) % 60;
	print "Time remaining: "
	  .sprintf("%i:%02i:%04.1f",$hours, $mins, $seconds) . "\n";
}
