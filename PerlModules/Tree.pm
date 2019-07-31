# Tree Functions
# I only use Newick trees and so I don't need Newick in the name

package Tree;

#tree = {
#       Name
#       Distance
#       Left
#       Right
#}

use strict;
use warnings FATAL => 'all';

use Math::Random qw(random_exponential);
use Array::Utils qw(:all);

unless (caller) {
        my $string = "((a:1,b:1):1,(d:1,c:1):1):1;";
    
        
        my $tree = FromString($string);
    GenerateSimilarRandomTree($tree);
        print ToString($tree);
}


# Modified from www.perlmonks.org/?node_id=717769
sub FromString {
    my ($string) = @_;

    $string =~ s/;//g;
    $string =~ s/\s+//g;
    
    my $tree = {};

    # This does not read scientific notation
#       $string =~ s/:(\d*\.?\d*)$//;

    $string =~ s/:([^:]+)$//;
    if (defined $1) {$tree->{Distance} = $1}
    else { warn "Cannot find distance for subtree $string"; }
    
    if ( $string !~ /[(),]/ ) {
        $tree->{Name} = $string;
    }
    else {
        $string =~ s/([^)]*)$//;
        if (defined $1) { $tree->{Name} = $1 }
        else { warn "Cannot find name for subtree $string" }

        my $depth = 0;

        for my $position ( 0 .. length($string) - 1 ) {
            my $c = substr $string, $position, 1;
            $depth++ if $c eq "(";
            $depth-- if $c eq ")";
            if ( $c eq "," and $depth == 1 ) {
                    $tree->{Left} = FromString(substr $string, 1, $position - 1);
                    $tree->{Left}->{Up} = $tree;
                    $tree->{Right} = FromString(substr $string, $position + 1, -1);
                    $tree->{Right}->{Up} = $tree;
                    last;
                }
            # found a segment
            # ((((...)a:1)b:1),...)c:1
            # b only has 1 child = a
            # If we make it here, we have scanned the whole 
            if ($depth == 0) {
                $tree->{Left} = FromString(substr $string, 1, $position - 1);
                                $tree->{Left}->{Up} = $tree;
            }
        }
    }
    return $tree;
}

# David's way...
sub FromString2 {
        my ($string) = @_;

        $string =~ s/;//g;
        $string =~ s/\s+//g;

        my $tree = {};
        $string =~ s/:(\d*\.?\d*)$//;
        if (defined $1) {$tree->{Distance} = $1}
        else {print $1 . "\n"; warn "Cannot find distance for subtree $string"; }

        if ( not $string =~ /[(),]/ ) {
                $tree->{Name} = $string;
        }
        else {
                $string =~ s/([^)]*)$//;
                if (defined $1) { $tree->{Name} = $1 }
                else { warn "Cannot find name for subtree $string" }

                ($tree->{Left}, $string) = FromString(substr $string, 1);

                ($tree->{Right},$string) = FromString(substr $string, 1, -1);
        }
        return $tree, $string;
}

sub FromFile {
        my ($file) = @_;
        open my $fh, $file or die $!;
        my $string = <$fh>;
    
        $string = CleanUpTreeLine($string);
        my $tree = FromString($string);
        return $tree;
}

sub ToString {
    my ($tree) = @_;

    my $string = "";

    if ($tree->{Left} and $tree->{Right}) {
        $string .= "(";
        $string .= ToString($tree->{Left});
        $string .= ",";
        $string .= ToString($tree->{Right});
        $string .= ")";
    }
    if ($tree->{Left} and not $tree->{Right}) {
        $string .= "(";
        $string .= ToString($tree->{Left});
        $string .= ")";
    }

    # Remove semicolons added by subtrees
    $string =~ s/;//g;

    $string .= $tree->{Name} // "";
    if (defined $tree->{Distance}) {
            $string .= ":";
            $string .= sprintf("%01.3f",$tree->{Distance});
    }

    # Add semicolon for final tree
    $string .= ";";

    return $string;
}

sub ToFile {
        my ($tree, $file) = @_;
        open my $fh, ">", $file or die $!;
    
    print $fh ToString($tree);
}

sub PrintTopology {
        my ($tree, $depth) = @_;
        print "\t" x $depth . $tree->{Name} . "\n";
        PrintTopology($tree->{Left}, $depth + 1) if $tree->{Left};
        PrintTopology($tree->{Right}, $depth + 1) if $tree->{Right};
}

# new method
sub Length {
        my ($tree) = @_;

        my $sum = $tree->{Distance};
        $sum += CalculateLength($tree->{Left}) if $tree->{Left};
        $sum += CalculateLength($tree->{Right}) if $tree->{Right};
        return $sum;
}

# old method
sub CalculateLength {
        return Length(@_);      
}

sub FindDeepestDescendant {
        my ($tree) = @_;

        my $longest_distance = 0;
        if ($tree->{Left} and $tree->{Right}) {
                my $l_dist = FindDeepestDescendant($tree->{Left});
                my $r_dist = FindDeepestDescendant($tree->{Right});

                if ($l_dist > $r_dist) {
                        $longest_distance = $l_dist + $tree->{Distance};
                }
                else {
                        $longest_distance = $r_dist + $tree->{Distance};
                }
        }
    elsif ($tree->{Left} and not $tree->{Right}) {
        my $l_dist = FindDeepestDescendant($tree->{Left});
        $longest_distance = $tree->{Distance} + $l_dist;
    }
        else {
                $longest_distance = $tree->{Distance};
        }

        return $longest_distance;
}

sub AttachSequences {
        my ($tree, $sequences)  = @_;
    
    recurse_pre($tree, \&AttachSequence, $sequences)
}

sub AttachSequence {
        my ($tree, $sequences)  = @_;
    if ($tree->{Name}) {
                if (defined $sequences->{$tree->{Name}}) {
                        $tree->{Sequence} = $sequences->{$tree->{Name}};
                }
                else{
                        warn "Cannot find a sequence for node with name " . $tree->{Name};
                }
        }
    }

sub AllSequencesToString {
        my ($tree) = @_;
        my $string = ">$tree->{Name}\n$tree->{Sequence}\n";

        $string .= AllSequencesToString($tree->{Left}) if $tree->{Left};
        $string .= AllSequencesToString($tree->{Right}) if $tree->{Right};

        return $string;
}

sub ToSequences {
        my ($tree) = @_;
    
    my $sequences = {};
    recurse_pre($tree, sub {my ($node, $sequences) = @_; $sequences->{$node->{Name}}=$node->{Sequence}//""}, $sequences);

        return $sequences;
}

sub AllHiddenStatesToString {
        my ($tree) = @_;
        my $string = ">$tree->{Name}\n$tree->{HiddenState}\n";

        $string .= AllHiddenStatesToString($tree->{Left}) if $tree->{Left};
        $string .= AllHiddenStatesToString($tree->{Right}) if $tree->{Right};

        return $string;
}

sub AllHiddenStatesToString2 {
        my ($tree) = @_;
        my $string = ">$tree->{Name}\n" . join("",$tree->{HiddenState}->members()) . "\n";

        $string .= AllHiddenStatesToString2($tree->{Left}) if $tree->{Left};
        $string .= AllHiddenStatesToString2($tree->{Right}) if $tree->{Right};

        return $string;
}

sub AttachSubstitutions {
        my ($tree, $substitutions)  = @_;

        if ($tree->{Name}) {
                if (defined $substitutions->{$tree->{Name}}) {
                        $tree->{Substitutions} = $substitutions->{$tree->{Name}};
                }
                else{
                        warn "Cannot find substitutions for node with name ". $tree->{Name};
                }
        }

        if ($tree->{Left}) {
                AttachSubstitutions($tree->{Left}, $substitutions);
        }
        if ($tree->{Right}) {
                AttachSubstitutions($tree->{Right}, $substitutions);
        }
}

sub CleanUpTreeLine {
        my ($line) = @_;
        $line =~ s/^tree gen\d* = //; # this is for a PLEX tree out file

        $line =~ s/\) \d+\.\./\)Node_/g; # This only affects internal branches
        $line =~ s/ \d+\.\.\d+//g; # Remove PLEX's old label from all branches

        $line =~ s/\);/\)Node_0:0;/; # Append the root label and distance
        return $line;
}

sub SequencesToSubstitutions {
        my ($tree) = @_;

        if ($tree->{Left}) {
                $tree->{Left}->{Substitutions} =
                  GetSubstitutions($tree->{Sequence}, $tree->{Left}->{Sequence});
                SequencesToSubstitutions($tree->{Left});
        }
        if ($tree->{Right}) {
                $tree->{Right}->{Substitutions} =
                  GetSubstitutions($tree->{Sequence}, $tree->{Right}->{Sequence});
                SequencesToSubstitutions($tree->{Right});
        }

}

sub GetSubstitutions {
        my ($ancestral_sequence, $descendant_sequence) = @_;

        my $subs = {};

        # Using 1 indexed counting
        foreach my $site (1 .. length($ancestral_sequence)) {
                my $from = substr $ancestral_sequence, $site-1, 1;
                my $to = substr $descendant_sequence, $site-1, 1;
                if ($from ne $to) {
                        $subs->{$site} = {from => $from, to => $to};
                }
        }
        return $subs;
}

sub GenerateNamesForInternalNodes {

        #       my ($filename, $filename_out) = @_;
        #       open my $fh, $filename or die $!;
        #       my $string = <$fh>;
        my ($string) = @_;
        my $current_node_id = 0;

        # must use the old construction because the length needs to be recalculated
        # every time.
        for(my $position = 0; $position < length($string); $position++) {
                my $c = substr $string, $position, 1;
                if($c eq ")") {
                        if (substr($string, $position+1, 1) eq ":"){

                                #insert name
                                substr($string, $position+1, 0, "Node_$current_node_id");
                        }
                        else {

                                # If there is a state there, as in the phylobayes substitution
                                # files
                                substr($string, $position+1, 0, "Node_$current_node_id" . "_");
                        }
                        $current_node_id++;
                }
        }
        return $string;
}

sub ApplyInternalNamesToString {
        my ($internal_names, $string) = @_;

        # must use the old construction because the length needs to be recalculated
        # every time.
        for(my $position = 0; $position < length($string); $position++) {
                my $c = substr $string, $position, 1;
                if($c eq ")") {
                        die "Ran out of internale names" if not @$internal_names;
                        if (substr($string, $position+1, 1) eq ":"){
                                #insert name
                                substr($string, $position+1, 0, shift @$internal_names);
                        }
                        else {

                                # If there is a state there, as in the phylobayes substitution
                                # files
                                substr($string, $position+1, 0, shift(@$internal_names) . "_");
                        }
                }
        }
        return $string;
}

our $names = [];

sub GenerateNamesForInternalNodes_orderIndependent {
        my ($tree) = @_;

        if ($tree->{Left} and $tree->{Right}) {
                GenerateNamesForInternalNodes_orderIndependent($tree->{Left});
                GenerateNamesForInternalNodes_orderIndependent($tree->{Right});
                $tree->{Name} =
                  join("_", sort ($tree->{Left}->{Name}, $tree->{Right}->{Name}));
                push @$names, $tree->{Name};
        }
    # print "Finished with " . $tree->{Name} . "\n";
        # If a leaf, do nothing because it already has a name.

}

sub StripPhylobayesSubstitutions {
        my ($string) = @_;
        $string =~ s/([A-Z]);/$1:0:$1;/;
        $string =~ s/_?[A-Z]:([^:]+)[^,)]*:[A-Z]([,);])/:$1$2/g;
        return $string;
}

sub GetAllNames {
        my ($tree) = @_;
        
        my $names = [$tree->{Name}];

        push @$names, @{GetAllNames($tree->{Left})} if ($tree->{Left});
        push @$names, @{GetAllNames($tree->{Right})} if ($tree->{Right});
    
        return $names;
}

sub GetAncestorNames {
        my ($tree) = @_;
        
        my $names = [];
        if ($tree->{Left} and $tree->{Right}) {
                push @$names, $tree->{Name};
                push @$names, @{GetAncestorNames($tree->{Left})};
                push @$names, @{GetAncestorNames($tree->{Right})};
        }
        return $names;
}

sub GenerateSimilarRandomTree {
    my ($tree) = @_;
    my $names = GetAllNames($tree);
    my $ancs = GetAncestorNames($tree);
    
    # Average branch length
    my $len = Length($tree) / (2 * (scalar @$names)-1);
    my $lens = [random_exponential(2 * (scalar @$names)-1, $len)];
    
    $names = [array_minus(@$names, @$ancs)]; # Only keep leaf names
    $names = [grep {$_} @$names]; # get rid of empty or nill names
    
    
    my $nodes = [map {{Name=>$_, Distance=>shift @$lens}} @$names];
    
    
    while (@$nodes > 1) {
        my$r1 = int(rand @$nodes);
        my$r2 = int(rand @$nodes);
        next if $r1 == $r2;
        # Must splice the larger of the two random nodes first
        if ($r1 < $r2) {($r1, $r2) = ($r2, $r1)}
        print "$r1 $r2\n";
        my $n1 = splice @$nodes, $r1, 1;
        my $n2 = splice @$nodes, $r2, 1;
        
        # Randomly insert new node
        my $r3 = int(rand @$nodes);
        splice @$nodes, $r3, 0, {Distance=>shift @$lens, Left=>$n1, Right=>$n2};
    }
    
    return $nodes->[0]
}

sub ToSubstitutions {
        my ($tree) = @_;
        
        my $subs_l = $tree->{Left} ? ToSubstitutions($tree->{Left}) : {};
        my $subs_r =$tree->{Right} ?  ToSubstitutions($tree->{Right}) : {};
        
        my $subs = {%$subs_l , %$subs_r};
        $subs->{$tree->{Name}} = [map { $tree->{Substitutions}->{$_}->{from} . $_ . $tree->{Substitutions}->{$_}->{to}} keys %{$tree->{Substitutions}}];

        return $subs;
}

sub FindLongBranches {
    my ($tree, $threshold) = @_;
    
    my $long_branches = [];
    push @$long_branches, $tree->{Name} if $tree->{Distance} > $threshold;
    
    if ($tree->{Left}){
        push @$long_branches, @{FindLongBranches($tree->{Left}, $threshold)};
        }
    if ($tree->{Right}){
        push @$long_branches, @{FindLongBranches($tree->{Right}, $threshold)};
        }
    return $long_branches
}

sub ToB1 {
    my ($tree, $threshold) = @_;
    
    if ($tree->{Distance} > $threshold) {
        
        # Add new nodes above the current node
        my $number_of_new_nodes = int($tree->{Distance} / $threshold);
        my $new_length = $tree->{Distance} / ($number_of_new_nodes + 1);
        $tree->{Distance} = $new_length;
        
        my $Up = $tree->{Up};
        my $old_descendant = $tree;
        
        for my $i (1 .. $number_of_new_nodes) {
            # Here I could copy the $tree so that all the attached 
            # information (sequences, models, etc) would propagate up
            my $new_node = {"Name" => $tree->{Name} . "_" . $i,
                            "Distance" => $new_length,
                            "Left" => $old_descendant,
                           };
            
            if ($old_descendant->{Sequence}){
                if($Up and $Up->{Sequence} and $i > ($number_of_new_nodes/2)){
                    $new_node->{Sequence} = $Up->{Sequence};
                }
                else {
                    $new_node->{Sequence} = $old_descendant->{Sequence};
                }
            }
            $old_descendant->{Up} = $new_node;
            $old_descendant = $new_node;
        }
        # old_descendant is now the node furthest from the $tree node
        
        # Have to tell the Up that it has a new node attached
        if ($Up){
            $old_descendant->{Up} = $Up;
            if ($tree == $Up->{Left}){
                $Up->{Left} = $old_descendant;
            } else {
                $Up->{Right} = $old_descendant;
            }
        }
    }
}

sub recurse {
        my ($tree, $pre_order_sub, $in_order_sub, $post_order_sub, @args) = @_;
        
        $pre_order_sub->($tree, @args);
        
        recurse($tree->{Left}, $pre_order_sub, $in_order_sub, $post_order_sub, @args) if $tree->{Left};
        $in_order_sub->($tree, @args);
        recurse($tree->{Right}, $pre_order_sub, $in_order_sub, $post_order_sub, @args) if $tree->{Right};
        
        $post_order_sub->($tree, @args);                
}

sub recurse_pre {
        my ($tree, $pre_order_sub, @args) = @_;
        recurse($tree, $pre_order_sub, sub{}, sub{}, @args);    
}

sub recurse_in {
        my ($tree, $in_order_sub, @args) = @_;
        recurse($tree, sub{}, $in_order_sub,  sub{}, @args);    
}

sub recurse_post {
        my ($tree, $post_order_sub, @args) = @_;
        recurse($tree, sub{}, sub{}, $post_order_sub, @args);   
}

1;
