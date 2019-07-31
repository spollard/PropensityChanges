# README.txt generator

use strict;
use warnings;

unless (caller) {
	GenerateReadme("Readme.template",
		{"DATE" => "2016-10-04", "AUTHOR" => "Stephen T. Pollard"}, "Readme_example.txt");
}

sub GenerateReadme {
	my ($template_filename, $options, $readme_filename) = @_;

	local $/ = undef; # Turn on slurp mode
	open my $template, "<", $template_filename or die $!;
	my $string = <$template>;
	while (my ($option, $value) = each %$options) {
		if ($string !~ s/$option/$value/g) {
			warn "Could not find a place for option $option";
		}
	}
	open my $readme, ">", $readme_filename or die $!;
	print $readme $string;
}

1;
