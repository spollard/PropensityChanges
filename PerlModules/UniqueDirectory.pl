# Generates a unique directory

use strict;
use warnings;

use DateTime;

unless (caller) {
	my $base = $ARGV[0] // "";
	print UniqueDirectory($base);	
}


sub UniqueDirectory {
	my ($base) = @_;
	my $dir;

	my $i = 0;
	do {
		my $dt = DateTime->now(time_zone => 'America/Denver');
		$dir = $base . $dt->ymd("-") . 'T' . $dt->hms("_") . "_$i/";
		$i++;
	}
	while (-d $dir);
	
	return $dir;
}


sub UniqueDirectoryAppendNumber {
	my ($base) = @_;
	my $dir;

	my $i = 0;
	do {
		$dir = $base . "_$i/";
		$i++;
	}
	while (-d $dir);
    
	return $dir;
}


1;