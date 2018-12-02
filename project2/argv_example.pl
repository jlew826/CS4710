## subroutine takes multiple scalar input
print("## subroutine takes multiple scalar input ##\n");
sub printArg{
	my ($a, $b, $c) = @_;

	print "Inside sub: $a, $b, $c\n\n";
	return($a);

}

$from_sub = printArg(100, 200, 300);
print "$from_sub\n=======\n";
