#!/usr/bin/env perl 

use strict; use warnings;
use jlagardeJsonschema;
#natively "use App::jsonschema;"; this was modified in order to load a customized version of this module able to check json string format
		
&usage if @ARGV < 2;
my $schema_file = shift;

my $app = App::jsonschema->new(	schema_file => $schema_file );
$app->validate(@ARGV);

sub usage {
	print STDERR "Usage: $0 schema.json file1.json [file2.json ...]\n";
	exit 1;
}


# ABSTRACT: Validate JSON files using JSON Schema
# PODNAME: jsonschema.pl

__END__

=pod

=encoding utf-8

=head1 NAME

jsonschema.pl - Validate JSON files using JSON Schema

=head1 VERSION

version 0.03

=head1 SYNOPSIS

jsonschema.pl schema.json file1.json [file2.json ...]

=head1 SEE ALSO

L<App::jsonschema>, L<JSON>, L<JSON::Schema>

=head1 AUTHOR

Andre Santos <andrefs@cpan.org>

=head1 COPYRIGHT AND LICENSE

This software is copyright (c) 2013 by Andre Santos.

This is free software; you can redistribute it and/or modify it under
the same terms as the Perl 5 programming language system itself.

=cut
