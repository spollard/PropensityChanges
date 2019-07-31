# Progress indicator.
#
# Stephen Pollard
# 2016-11-18

use strict;
use warnings FATAL => 'all';
use diagnostics;

use Time::HiRes 'sleep', 'time';


unless (caller) {
print(join(":",time_to_hours_mins_seconds(86973)));exit;
    my $n = 20;
    my $progress = Progress($n,3);
    for (1..$n) {
        my $time = 1 + (rand(10) + 1)/100.0;
        #print "Sleeping for $time.\n";
        sleep $time;
        $progress->($_);
    }
}


sub Progress {
    my ($total_iterations, $update_interval) = @_;
    
    $update_interval //= 10; #10 seconds default interval

    my $start_time = time();
    my $previous_time = $start_time;
    my $previous_iteration = 0;
    print "Iteration 1\n";
    
    return sub {
        my ($at) = @_;
        if (((time() - $previous_time) > $update_interval) or $at == $total_iterations) {
            my $time_taken = (time() - $start_time) | 1;
            my $rate = $at / $time_taken;
            
            print "Iteration $at - "
            . int($at / $total_iterations * 100) . "% - "
            . sprintf("%.2f",$rate) . " iterations/s\n";
              
            my ($hours_taken, $mins_taken, $seconds_taken) = time_to_hours_mins_seconds($time_taken);
            print "Time taken (H:MM:SS): "    
              .sprintf("%i:%02i:%02i",$hours_taken, $mins_taken, $seconds_taken) . "\n";
              
            my $time_left = ($total_iterations - $at) / $rate;
            my ($hours_left, $mins_left, $seconds_left) = time_to_hours_mins_seconds($time_left);
            
              
            my ($hours_tot, $mins_tot, $seconds_tot) = time_to_hours_mins_seconds($time_left + $time_taken);
            print "Estimated total time (H:MM:SS): "  
              .sprintf("%i:%02i:%02i",$hours_tot, $mins_tot, $seconds_tot) . "\n";

            print "Estimated time remaining (H:MM:SS): "  
              .sprintf("%i:%02i:%02i",$hours_left, $mins_left, $seconds_left) . "\n\n";
              
            
            $previous_time = time();
            $previous_iteration = $at;
        }
    }
}

sub time_to_hours_mins_seconds {
    my ($time) = @_;
    my $hours = int($time / 3600);
    my $mins = int(($time % 3600) / 60);
    my $seconds = ($time % 3600) % 60;
    return $hours, $mins, $seconds;
}

1;
