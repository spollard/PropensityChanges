use strict;
use warnings;

use Data::Dumper;

use Set::Tiny;

our $PerlModulesDir;
BEGIN { require "PerlModulesDir.pl"; }
use lib $PerlModulesDir;

use Sequences;
use Samplers;

require "amino_acids.pl";
our $amino_acids;

unless (caller){
    #~ my $acceptabilities_filename_1 = $ARGV[0] // "../data/acceptabilities_actual.fasta";
    #~ my $acceptabilities_filename_2 = $ARGV[1] // "../data/acceptabilities_correct_tightness.fasta";
    
    #my $acceptabilities_filename_1 = $ARGV[0] // "../data/1914x100.acceptables";
    #my $acceptabilities_filename_2 = $ARGV[1] // "../data/1914x100/acceptabilities2.fasta";
    
    my $acceptabilities_filename_1 = $ARGV[0] // "../data/100x100.acceptables";
    my $acceptabilities_filename_2 = $ARGV[1] // "../data/simulated_100x100_g10k_sampled/acceptabilities_3kburnin.fasta";


    compare_acceptabilities($acceptabilities_filename_1, $acceptabilities_filename_2);
}

sub compare_acceptabilities {
    # Print a matrix of the lengths of the first acceptabilities compared to the 
    # lengths of the second acceptabilities 
    my ($acceptabilities_filename_1, $acceptabilities_filename_2) = @_;

    my $acc_1 = Sequences::FromFasta($acceptabilities_filename_1);
    my $acc_2 = Sequences::FromFasta($acceptabilities_filename_2);
    
    my $names_1 = [keys %$acc_1];
    my $names_2 = [keys %$acc_2];
    
    for my $name_1 (@$names_1) {
            if (not( $name_1 ~~ @$names_2)) {
                    warn "$name_1 not found in second set of acceptabilities";      
            }
    }
        
        # Use a hash because the comparison matrix might be sparse
    my $acceptabilities = {};
    my $max_1 = 0;
    my $max_2 = 0;
    
    my $correct = 0;
    my $missing = 0;
    my $added = 0;
    my $total = 0;
    
    my $acc_1_totals = [];
    my $acc_2_totals = [];

    for my $name_2 (@$names_2) {
        if (not( $name_2 ~~ @$names_1)) {
                warn "$name_2 not found in first set of acceptabilities";       
        }
       # must include a negative split limit
        my $accs_1 = [split /\|/, $acc_1->{$name_2}, -1];
        my $accs_2 = [split /\|/, $acc_2->{$name_2}, -1];
        
        if (scalar(@$accs_1) != scalar(@$accs_2)){
            die "Number of sites for $name_2 does not match in both files: " . 
            scalar(@$accs_1) . ":" . scalar(@$accs_2);  
        }
        for (0..$#$accs_1){
            my $ac1 = Set::Tiny::set(split //, $accs_1->[$_]);
            my $ac2 = Set::Tiny::set(split //, $accs_2->[$_]);
            
            $total += $ac1->size();
            $correct += $ac1->intersection($ac2)->size();
            $missing += $ac1->difference($ac2)->size();
            $added += $ac2->difference($ac1)->size();
            
            $acceptabilities->{length($accs_1->[$_])}->{length($accs_2->[$_])}++;
            $acc_1_totals->[length($accs_1->[$_])]++;
            $acc_2_totals->[length($accs_2->[$_])]++;
            if(length($accs_1->[$_]) > $max_1){$max_1 = length($accs_1->[$_])};
            if(length($accs_2->[$_]) > $max_2){$max_2 = length($accs_2->[$_])};
        }
    }
    
    print join("\t", "total:", $total, "correct:", $correct, "missing:", $missing, 
          "added:", $added), "\n";
    for my $i (0..$#$acc_1_totals){
        print("Correct total $i: " . ($acc_1_totals->[$i]//0) . "\n")
    }
    for my $i (0..$#$acc_2_totals){
        print("Estimated total $i: " . ($acc_2_totals->[$i]//0) . "\n")
    }
    for my $i (0..$max_1) {
        print $i, "\t";
        for my $j (0..$max_2){
            print $acceptabilities->{$i}->{$j} // 0, "\t";
        }
        print "\n";
    }
}
