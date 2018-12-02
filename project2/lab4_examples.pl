### lab4 subroutines

## subroutine takes single scalar input
print("## subroutine takes single scalar input ##\n");
sub addACGT {
    my ($dna) = @_;
    $dna .= 'ACGT';
    return $dna;
}

print addACGT('AAA'), "\n";
print $dna,"\n=======\n";


## subroutine takes multiple scalar input
print("## subroutine takes multiple scalar input ##\n");
sub printArg{
	my ($a, $b, $c) = @_;

	print "Inside sub: $a, $b, $c\n\n";
	return($a);

}

$from_sub = printArg(100, 200, 300);
print "$from_sub\n=======\n";

## subroutine takes single array
print("## subroutine takes single array ##\n");
sub addone{
	my @a = @_;
	print "Inside sub: $a[0], $a[1]\n\n";
	return($a[0]+1);

}

$from_sub = addone((10, 20));
print "$from_sub\n=======\n";

## subroutine takes multiple arrays
print("## subroutine takes multiple arrays ##\n");

sub add {
	my @a, @b = @_;
	print(@a,"\n");
	print(@b,"\n");
   print "Inside sub: @_\n=====\n";  # 2 3 7 8 5
}
my @first  = (2, 3);
my @second = (7, 8, 5);
add(@first, @second);

print("## subroutine takes multiple arrays by reference ##\n");

sub add_ref {

	my ($one_ref, $two_ref) = @_;
    my @one = @{ $one_ref };       # dereferencing and copying each array
    my @two = @{ $two_ref };

    print "@one\n";    # 2 3
    print "@two\n";    # 7 8 5
    print "${$one_ref}[0]\n======\n"

}
my @first  = (2, 3);
my @second = (7, 8, 5);
add_ref(\@first, \@second);

print("## mixed situations ##\n");

sub mixed {

	my ($one_ref, $two_ref, $scalar) = @_;
    my @one = @{ $one_ref };       # dereferencing and copying each array
    my @two = @{ $two_ref };

    print "@one\n";    # 2 3
    print "@two\n";    # 7 8 5
    print "$scalar\n======\n";

}
my @first  = (2, 3);
my @second = (7, 8, 5);
mixed(\@first, \@second, "hello world");


#practice making two arrays, and passing them as references, and find the minimum.

sub findmin {
  my($one_ref, $two_ref, $scalar) = @_;
    my @one = @{ $one_ref};
    my @two = @{ $two_ref};
    
    print "@one\n";
    print "@two\n";
    print "$scalar\n======\n";
}
my @first = (2, 3, 4, 5, 6, 7, 8, 9);
my @second = (1, 3, 4, 5, 6, 7, 8, 9);
findmin(\@first, \@second, "hello world");
