
use JSON;
#use Try::Tiny;
sub processJsonToHash{
	open JSON, "$_[0]" or die "$_[0] : $!\n";
	my $whole_json_file='';
#	{
#		local $/;
#		$whole_json_file=<JSON>;
#	}
	while (my $line = <JSON>){
		$line=~s/(\r)|(\n)|(\t)//g; # must process line by line as opposed to slurp whole file, since perl can't process > 1GB long string ("Error: substitution loop" blahblah)
	    $whole_json_file .= $line;
	}
	close JSON;
	#$whole_json_file=~s/(\r)|(\n)|(\t)//g;
	my $tree=undef;
	eval {
		$tree = decode_json($whole_json_file);
		1;
	} or do  {
#	unless (defined $tree){
		print STDERR "processJsonToHash: Caught error: $@ \nprocessJsonToHash: Returning empty hash ('{}')\n";
#		$tree=undef;
	};

	return $tree;
}

1;
