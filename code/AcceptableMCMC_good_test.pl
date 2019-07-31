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

use List::Util qw(sum shuffle);
use List::MoreUtils 'any';

use Set::Tiny;

use Data::Dumper;

use Clone qw(clone);

use lib 'C:\Users\Stephen\Documents\Lab work\PerlModules';

use Tree;
use Sequences;
use Propensities;
use Samplers;


require "MCMC.pl";
require "Progress.pl";



# Globals for now
my $ancestors_known = 1;
my ($ancestral_node_names, $all_node_names);
my $amino_acids_considered = [split "", "ABCDE"];
my $max_branch_length = 0.4;

unless(caller){
    my $site = $ARGV[0] // 2; # sites are 1 indexed

    # Data Files
    my $tree_filename = $ARGV[1] // "../data/test_tree3.newick";
    my $sequences_filename = $ARGV[2] // "../data/test_seqs4.fasta";
    my $tree_out_filename = $ARGV[3] // "../data/test/tree.newick";


    AcceptableMCMC($site, $tree_filename, $sequences_filename, $tree_out_filename);
}

sub AcceptableMCMC {
    my ($site, $tree_filename, $sequences_filename, $tree_out_filename) = @_;

    # Model options
    #my $amino_acids_considered = [split "", "ACDEFGHIKLMNPQRSTVWY"];
    
    my $out_directory = "../data/acceptable/";
    #$out_directory = UniqueDirectoryAppendNumber($out_directory);
    mkdir $out_directory if not -d $out_directory;


    # MCMC files
    my $mcmc_out_filename = $out_directory . "mcmc";
    my $model_out_filename = $out_directory . "model";

    # MCMC options
    my $generations = 2000;

    # data = [distance, C/D]
    print "Reading Data\n";


    my $data = ReadData($tree_filename, $sequences_filename,$site, $tree_out_filename);

    # globals for now
    $ancestral_node_names = Tree::GetAncestorNames($data);
    $all_node_names = Tree::GetAllNames($data);

    # Initialize output files
    my $mcmc_out = InitializeMcmcOutput($mcmc_out_filename);
    my $model_out = InitializeModelOutput($model_out_filename, $amino_acids_considered);

    my $model = InitializeModel($amino_acids_considered, $data);
    my $log_likelihood = LogLikelihoodOfModelGivenData($model, $data);


    print "Beginning MCMC loop\n\n";

    my $time_between_prints = 3;
    my $progress = Progress($generations, $time_between_prints);

    my $resample_probabilities = 0.3;

    for my $generation (1 .. $generations) {
        $progress->($generation);

        RecordModelState($model, $model_out);
        
        Tree::recurse_pre($data, \&GibbsSampleAcceptables, $model);
        GibbsSampleProbabilities($model, $data);
        $log_likelihood = LogLikelihoodOfModelGivenData($model, $data);
        RecordMcmcState($mcmc_out, $generation, $log_likelihood, 1, 1, 0);
    }
    print "MCMC Finished\n";

}


sub InitializeModelOutput {
	my ($model_out_filename, $amino_acids_considered) = @_;
	open my $model_out, ">", $model_out_filename;
    
    print $model_out join("\t", qw(Prior_A SwitchProb FewAcceptables )), "\t";
    
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
    
    print $fh join("\t", ($model->{Prior_A}, $model->{SwitchProb}, $model->{FewAcceptables})), "\t";
    
    if ($ancestors_known) {
	    print $fh join("\t", map({$model->{Acceptables}->{$_}->as_string()} @$all_node_names) ) . "\n";
    } else {
        print $fh join("\t", @{$model->{AncestralResidues}}{@$ancestral_node_names}, map({$model->{Acceptables}->{$_}->as_string()} @$all_node_names) ) . "\n";
    }
}


sub LogLikelihoodNode {
    my ($node, $model, $log_likelihood) = @_;
    
    # residue sub
    #my $ancestral_residue = $node->{Up}->{Sequence} // $model->{AncestralResidues}->{$node->{Up}->{Name}};
    my $node_residue = $node->{Sequence} // $model->{AncestralResidues}->{$node->{Name}};
    
    if (not $node_residue){die "Node " . $node->{Name} . " has no sequence\n";}
    
    my $node_acceptables = $model->{Acceptables}->{$node->{Name}};
    
    #my $sub = $ancestral_residue ne $node_residue;
    
    #my $anc_subprob = $sub ? $model->{SubProbs}->[$anc_hidden_state] : 1 - $model->{SubProbs}->[$anc_hidden_state];
    #my $node_subprob = $sub ? $model->{SubProbs}->[$node_hidden_state] : 1 - $model->{SubProbs}->[$node_hidden_state];
    
    # Since we don't know if the sub happened on the top half or bottom half of the branch
    #$$log_likelihood += log(0.5 * $anc_subprob + 0.5 * $node_subprob);
    #$$log_likelihood += log($node_subprob);
    
    # Prior of being in acceptable or not
    $$log_likelihood += log($node_acceptables->has($node_residue) ? $model->{Prior_A} : 1-$model->{Prior_A});
    
    # Prior against having a high number of acceptables
    $$log_likelihood += log(ProbOfNumberAcceptables($node_acceptables->size()));
    
    if ($node->{Up}) {
                    my $anc_acceptables = $model->{Acceptables}->{$node->{Up}->{Name}};
        # How many residues switched to/from acceptable
        my $differences = $anc_acceptables->symmetric_difference($node_acceptables)->size();
        if ($differences){
            $$log_likelihood += log($model->{SwitchProb} * $differences);
        }
    }
}

sub ProbOfNumberAcceptables {
    my ($number_of_acceptables) = @_;
    
    my $probs = [0.00001, 1, 0.1, 0.01];
    return 0.001 if $number_of_acceptables > $#$probs;
    
    return $probs->[$number_of_acceptables];
}
    
sub LogLikelihoodOfModelGivenData {
	my ($model, $data) = @_;

	my $log_likelihood = 0;
    Tree::recurse_pre($data, \&LogLikelihoodNode, $model, \$log_likelihood);

	return $log_likelihood;
}
    


sub GibbsSampleAcceptables {
    my ($node, $model) = @_;
    
    my $acceptables = $model->{Acceptables}->{$node->{Name}};
    # Need to loop over residues in a random order
    for my $residue (shuffle @$amino_acids_considered) {
    
        my $prob_N = 0;
        $acceptables->remove($residue);
        LogLikelihoodNode($node, $model, \$prob_N);
        LogLikelihoodNode($node->{Left}, $model, \$prob_N) if $node->{Left};
        LogLikelihoodNode($node->{Right}, $model, \$prob_N) if $node->{Right};
        
        my $prob_A = 0;
        $acceptables->insert($residue);
        LogLikelihoodNode($node, $model, \$prob_A);
        LogLikelihoodNode($node->{Left}, $model, \$prob_A) if $node->{Left};
        LogLikelihoodNode($node->{Right}, $model, \$prob_A) if $node->{Right};
        
        # These are log probabilities
        # Here is where you would put a prior against being acceptable
        if (rand() < exp($prob_A)/(exp($prob_N)+exp($prob_A))) {
            $acceptables->insert($residue);
        }
        # It has already been inserted so you don't have to do this
        else {
            $acceptables->remove($residue);
        }
    }
}
    
sub GibbsSampleProbabilities {
    my ($model, $data) = @_;

    my $number_of_branches = scalar(@{Tree::GetAllNames($data)});

    my $aggregate = {
        sites_with_acceptable_residue => 0,
        switches => 0,
        # Don't use this yet. Ideally I would estimate the probability
        # distribution of the number of acceptale amino acids. Until
        # then simply use fixed probabilities.
        acceptables_sizes => [0,0,0,0],
    };

    # You could do this in a single pass
    Tree::recurse_pre($data, \&Aggregator, $aggregate, $model);

    my $optimal_acceptable_prior = $aggregate->{sites_with_acceptable_residue} / $number_of_branches;
    my $optimal_switch_prob = $aggregate->{switches} / $number_of_branches / scalar(@$amino_acids_considered);

    #print("Prior A $optimal_acceptable_prior\n  Switch prob $optimal_switch_prob\n");
    #print("Acceptable sizes " . join("\t",map{$_//0}@{$aggregate->{acceptables_sizes}}) . "\n");
    my $prior_stdev = 0.01;
    my $switch_stdev = 0.01;
   
   #$model->{Prior_A} = Samplers::normal_p($optimal_acceptable_prior,$prior_stdev);
   $model->{SwitchProb} = Samplers::normal_p($optimal_switch_prob, $switch_stdev);
}


sub Aggregator {
    my ($node, $aggregate,$model) = @_;
    
    if ($node->{Up}) {
        #my $ancestral_residue = $node->{Up}->{Sequence} // $model->{AncestralResidues}->{$node->{Up}->{Name}};
        my $node_residue = $node->{Sequence} // $model->{AncestralResidues}->{$node->{Name}};
        
        my $anc_acceptables = $model->{Acceptables}->{$node->{Up}->{Name}};
        my $node_acceptables = $model->{Acceptables}->{$node->{Name}};
        
        $aggregate->{sites_with_acceptable_residue}++ if $node_acceptables->has($node_residue);
        
        $aggregate->{switches} += $anc_acceptables->symmetric_difference($node_acceptables)->size();
        
        $aggregate->{acceptables_sizes}->[$node_acceptables->size()]++;
    }
}


sub InitializeModel {
	my ($amino_acids_considered, $data) = @_;
	my $model = {
		#Propensities=>Propensities::GenerateRandomPropensities(
		#	scalar @$amino_acids_considered
		#),
		Acceptables => {},
		AncestralResidues => {},
        # The following parameters could be sampled, but I'd rather 
        # have them be input by the user. This also reduces the 
        # amount of sampling required.
        # Priors for a site to have a residue which is Not-acceptable and Acceptable
        Prior_A => 0.9,
        # prob of being acceptable at the end of the branch
        # why penalize substitutions?
        #SubProb => 0.9,
        # Prob of an amino acid switching to/from acceptable
        SwitchProb => 0.01, 
        # Prior for a residue to be acceptable
        # does nothing at the moment
        FewAcceptables => 0.01,
	};

    if (not $ancestors_known){
        for my $name (@$ancestral_node_names) {
            # The ancestral residues will be read in eventually
            $model->{AncestralResidues}->{$name} = Samplers::random_choice($amino_acids_considered);
        }
    }
    for my $name (@$all_node_names) {
        #Randomly choose 1 amino acid to be acceptable
        $model->{Acceptables}->{$name} = Set::Tiny->new(Samplers::random_choice($amino_acids_considered));
    }
	return $model;
}


sub ReadData {
	my ($tree_filename,$sequences_filename, $site, $tree_out_filename) = @_;
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
