# Takes a tree and a site and estimates the likelihood of a fixed covarion
# fitting the data.
#
# I want to turn this into a generic phylogenetic MCMC which can take arbitrary
# models. 
# 

# Stephen Pollard
# 2015-6-13

use strict;
use warnings FATAL => 'all';

use List::Util 'sum';
use List::MoreUtils 'any';

use Data::Dumper;

use Clone qw(clone);

use lib 'C:\Users\Stephen\Documents\Lab work\PerlModules';

use Tree;
use Sequences;
use Propensities;
use Samplers;


require "MCMC.pl";
require "Progress.pl";

my $site = $ARGV[0] // 1; # sites are 1 indexed

# Data Files
my $tree_filename = $ARGV[1] // "../data/test_tree3.newick";
my $sequences_filename = $ARGV[2] // "../data/test_seqs3.fasta";
my $tree_out_filename = $ARGV[3] // "../data/test/tree.newick";

# Globals for now
my $ancestors_known = 1;
my ($ancestral_node_names, $all_node_names);
my $amino_acids_considered = [split "", "AB"];
my $max_branch_length = 0.6;

FixedCovarionMCMC($site, $tree_filename, $sequences_filename, $tree_out_filename);

sub FixedCovarionMCMC {
    my ($site, $tree_filename, $sequences_filename, $tree_out_filename) = @_;

    # Model options
    #my $amino_acids_considered = [split "", "ACDEFGHIKLMNPQRSTVWY"];
    
    my $out_directory = "../data/test/";
    #$out_directory = UniqueDirectoryAppendNumber($out_directory);
    mkdir $out_directory if not -d $out_directory;


    # MCMC files
    my $mcmc_out_filename = $out_directory . "mcmc";
    my $model_out_filename = $out_directory . "model";

    # MCMC options
    my $generations = 10000;

    # Sampling options
    # Perhaps number of propensities to sample per generation
    #my $maximum_means_to_sample_per_generation = $number_of_windows/4;
    # my $max_step_size = 1 / 100;
    # my $number_of_parameters_sampled = 0;

    # data = [distance, C/D]
    print "Reading Data\n";


    my $data = ReadData($tree_filename, $sequences_filename);

    # globals for now
    $ancestral_node_names = Tree::GetAncestorNames($data);
    $all_node_names = Tree::GetAllNames($data);




    # Initialize output files
    my $mcmc_out = InitializeMcmcOutput($mcmc_out_filename);
    my $model_out = InitializeModelOutput($model_out_filename, $amino_acids_considered);

    my $model = InitializeModel($amino_acids_considered, $data);
    my $log_likelihood = LogLikelihoodOfModelGivenData($model, $data);


    print "Beginning MCMC loop\n\n";

    my $time_between_prints = 10;
    my $progress = Progress($generations, $time_between_prints);


    my $resample_probabilities = 0.3;

    for my $generation (1 .. $generations) {
        $progress->($generation);

        RecordModelState($model, $model_out);
        
        Tree::recurse_pre($data, \&GibbsSampleHiddenStates, $model);
        GibbsSampleProbabilities($model, $data);
        $log_likelihood = LogLikelihoodOfModelGivenData($model, $data);
        RecordMcmcState($mcmc_out, $generation, $log_likelihood, 1, 1, 0);
    }
    print "MCMC Finished\n";

}
sub InitializeModelOutput {
	my ($model_out_filename, $amino_acids_considered) = @_;
	open my $model_out, ">", $model_out_filename;
    
    print $model_out join("\t", qw(Fixed_sub Variable_sub Fixed_switch Variable_switch)), "\t";
    
    if ($ancestors_known){
        print $model_out join("\t", @$all_node_names),"\n";
      }
      else {
          print $model_out join("\t",
            #@$amino_acids_considered,
            # R for residue
            map({"R".$_} @$ancestral_node_names),
            # H for hidden state
            map({"H".$_} @$all_node_names)),
          "\n";
      }

	return $model_out;
}

sub RecordModelState {
	my ($model, $fh) = @_;
    
    print $fh join("\t", (@{$model->{SubProbs}}, @{$model->{SwitchProbs}})), "\t";
    
    if ($ancestors_known) {
	    print $fh join("\t", @{$model->{HiddenStates}}{@$all_node_names} ) . "\n";
    } else {
        print $fh join("\t", @{$model->{AncestralResidues}}{@$ancestral_node_names}, @{$model->{HiddenStates}}{@$all_node_names} ) . "\n";
    }
}


sub LogLikelihoodNode {
    my ($node, $model, $log_likelihood) = @_;
    if ($node->{Up}) {
        # residue sub
        my $ancestral_residue = $node->{Up}->{Sequence} // $model->{AncestralResidues}->{$node->{Up}->{Name}};
        my $node_residue = $node->{Sequence} // $model->{AncestralResidues}->{$node->{Name}};
        
        if (not $node_residue){print "Node " . $node->{Name} . " has no sequence\n";exit;}
        
        my $anc_hidden_state = $model->{HiddenStates}->{$node->{Up}->{Name}};
        my $node_hidden_state = $model->{HiddenStates}->{$node->{Name}};
        
        my $sub = $ancestral_residue ne $node_residue;
        
        my $anc_subprob = $sub ? $model->{SubProbs}->[$anc_hidden_state] : 1 - $model->{SubProbs}->[$anc_hidden_state];
        my $node_subprob = $sub ? $model->{SubProbs}->[$node_hidden_state] : 1 - $model->{SubProbs}->[$node_hidden_state];
        
        # Since we don't know if the sub happened on the top half or bottom half of the branch
        #$$log_likelihood += log(0.5 * $anc_subprob + 0.5 * $node_subprob);
        $$log_likelihood += log($node_subprob);
        
        # hidden state switch
        my $switch_prob = $anc_hidden_state != $node_hidden_state ? $model->{SwitchProbs}->[$anc_hidden_state] 
                                                                 : 1-$model->{SwitchProbs}->[$anc_hidden_state];
        $$log_likelihood += log($switch_prob); 
    }
}
    
sub LogLikelihoodOfModelGivenData {
	my ($model, $data) = @_;

	my $log_likelihood = 0;
    Tree::recurse_pre($data, \&LogLikelihoodNode, $model, \$log_likelihood);

	return $log_likelihood;
}
    

sub GibbsSampleHiddenStates {
    my ($node, $model) = @_;
    my $probs = [];
    for my $state (0,1) {
        $probs->[$state] = 0;
        $model->{HiddenStates}->{$node->{Name}} = $state;
        LogLikelihoodNode($node, $model, \$probs->[$state]);
        LogLikelihoodNode($node->{Left}, $model, \$probs->[$state]) if $node->{Left};
        LogLikelihoodNode($node->{Right}, $model, \$probs->[$state]) if $node->{Right};
        $probs->[$state] = exp($probs->[$state]);
    }
    #print("Probs: @$probs\n");
    $model->{HiddenStates}->{$node->{Name}} = Samplers::weighted_random_choice($probs);
}
    
sub GibbsSampleProbabilities {
   my ($model, $data) = @_;
   
   my $number_of_branches = scalar(@{Tree::GetAllNames($data)});
   
   my $number_of_variable_subs = 0;
   my $number_of_fixed_subs = 0;
   my $number_of_variable_switches = 0;
   my $number_of_fixed_switches = 0;
   
   Tree::recurse_pre($data, \&VariableSub, \$number_of_variable_subs, $model);
   Tree::recurse_pre($data, \&FixedSub, \$number_of_fixed_subs, $model);
   Tree::recurse_pre($data, \&VariableSwitch, \$number_of_variable_switches, $model);
   Tree::recurse_pre($data, \&FixedSwitch, \$number_of_fixed_switches, $model);
   
   my $optimal_variable_sub = $number_of_variable_subs / $number_of_branches;
   my $optimal_fixed_sub = $number_of_fixed_subs / $number_of_branches;
   my $optimal_variable_switch = $number_of_variable_switches / $number_of_branches;
   my $optimal_fixed_switch  = $number_of_fixed_switches / $number_of_branches;
   
    
    my $sub_stdev = 0.1;
    my $switch_stdev = 0.1;
   
   #$model->{SubProbs} = [Samplers::normal_p($optimal_fixed_sub,$sub_stdev), Samplers::normal_p($optimal_variable_sub, $sub_stdev)];
   $model->{SwitchProbs} = [Samplers::normal_p($optimal_fixed_switch, $switch_stdev), Samplers::normal_p($optimal_variable_switch, $switch_stdev)];
}

sub SampleProbabilities {
    my ($model) = @_;
    
    my $sub_stdev = 0.1;
    my $switch_stdev = 0.1;
    
    $model->{SubProbs} = [Samplers::normal_p($model->{SubProbs}->[0], $sub_stdev), Samplers::normal_p($model->{SubProbs}->[1], $sub_stdev)];
    $model->{SwitchProbs} = [Samplers::normal_p($model->{SwitchProbs}->[0], $switch_stdev), Samplers::normal_p($model->{SwitchProbs}->[1], $switch_stdev)];
}

sub VariableSub {
    my ($node, $count,$model) = @_;
    
    if ($node->{Up}) {
        my $ancestral_residue = $node->{Up}->{Sequence} // $model->{AncestralResidues}->{$node->{Up}->{Name}};
        my $node_residue = $node->{Sequence} // $model->{AncestralResidues}->{$node->{Name}};
        
        my $node_hidden_state = $model->{HiddenStates}->{$node->{Name}};
        
        my $sub = $ancestral_residue ne $node_residue;
        
        $$count++ if (($node_hidden_state  == 1) and $sub);
    }
}

sub FixedSub {
    my ($node, $count, $model) = @_;
    
    if ($node->{Up}) {
        my $ancestral_residue = $node->{Up}->{Sequence} // $model->{AncestralResidues}->{$node->{Up}->{Name}};
        my $node_residue = $node->{Sequence} // $model->{AncestralResidues}->{$node->{Name}};
        
        my $node_hidden_state = $model->{HiddenStates}->{$node->{Name}};
        
        my $sub = $ancestral_residue ne $node_residue;
        
        $$count++ if (($node_hidden_state  == 0) and $sub);
    }
}

sub VariableSwitch {
    my ($node, $count,$model) = @_;
    
    if ($node->{Up}) {
        my $anc_hidden_state = $model->{HiddenStates}->{$node->{Up}->{Name}};
        my $node_hidden_state = $model->{HiddenStates}->{$node->{Name}};
        
        my $switch = $anc_hidden_state != $node_hidden_state;
        $$count++ if (($node_hidden_state  == 1) and $switch);
    }
}

sub FixedSwitch {
    my ($node, $count, $model) = @_;
    
    if ($node->{Up}) {
        my $anc_hidden_state = $model->{HiddenStates}->{$node->{Up}->{Name}};
        my $node_hidden_state = $model->{HiddenStates}->{$node->{Name}};
        
        my $switch = $anc_hidden_state != $node_hidden_state;
        $$count++ if (($node_hidden_state  == 0) and $switch);
    }
}


sub InitializeModel {
	my ($amino_acids_considered, $data) = @_;
	my $model = {
		#Propensities=>Propensities::GenerateRandomPropensities(
		#	scalar @$amino_acids_considered
		#),
		HiddenStates => {},
		AncestralResidues => {},
        # The following parameters could be sampled, but I'd rather 
        # have them be input by the user. This also reduces the 
        # amount of sampling required.
        # for fixed or variable hidden state, prob of not sub or sub
        #SubProbs => [[.25, 0.05], [0.75, 0.5]],
        #SubProbs => [[1, 1], 
        #            [.1, 1]],
        # The first sub prob is prob of sub given off state. This should
        # be very low to force switching to on state for subs. 
        # The second number is sub prob of sub given on state. This will
        # adjust how long you are in the on state before and after a sub.
        SubProbs => [0.0005, 0.75],
        # probs [[ fixed->fixed, fixed->variable], [variable->fixed, variable->variable]]
        #SwitchProbs => [[.99, 0.01], [0.01, 0.99]],
        #SwitchProbs => [[1, 1], 
        #                [1, 1]],
        # The first number is how likely you are to switch from off to on. 
        # The second is from on to off. 
        # I'm leaving this at 0.5 so that it can switch freely for testing.
        SwitchProbs => [0.01, 0.01], 
	};

    if (not $ancestors_known){
        for my $name (@$ancestral_node_names) {
            # The ancestral residues will be read in eventually
            $model->{AncestralResidues}->{$name} = Samplers::random_choice($amino_acids_considered);
        }
    }
    for my $name (@$all_node_names) {
        # 0 is 'off' and 1 is 'on'
        #$model->{HiddenStates}->{$name} = Samplers::random_choice([0,1]);
        $model->{HiddenStates}->{$name} = 0;
    }
	return $model;
}


sub ReadData {
	my ($tree_filename,$sequences_filename) = @_;
	my $tree = Tree::FromFile($tree_filename);
    
	my $sequences = Sequences::FromFasta($sequences_filename);
    $sequences = Sequences::Select($sequences,[$site]);
	Tree::AttachSequences($tree, $sequences);
    
    # B1 allows us to ignore double flip probability
    Tree::recurse_pre($tree, \&Tree::ToB1, $max_branch_length);
    
    # In case the root has been segmented
    while (defined $tree->{Up}){ $tree = $tree->{Up} }
    
    Tree::ToFile($tree, $tree_out_filename);
    Tree::recurse_post($tree, \&PropagateSequences);
    
	my $data = $tree;

	return $data;
}

sub PropagateSequences {
    my ($tree) = @_;
    
    if (not defined $tree->{Sequence}) {
        $tree->{Sequence} = $tree->{Left}->{Sequence};
    }
}
