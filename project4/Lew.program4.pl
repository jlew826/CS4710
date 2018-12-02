#!/usr/bin/perl
use Math::Round;
use Data::Dumper;

#storing command line argument as variable
$filename = $ARGV[0];

#opening .sif file
if ($filename =~ /.*\.sif/){
  open (FILE, $filename)
    or die "Couldn't open \"$filename\"\n";

    $node1;
    @node1 = ();
    $node2;
    @node2 = ();

    while (<FILE>){
      chomp $_;
      @columns = split /\s/, $_; #splits line of string by any whitespace
      $node1[@node1] = $columns[0]; #node1 is in the first position of the array
      $node2[@node2] = $columns[2]; #node2 is in the second position of the array
    }

    #close the file;
    close $filename;
}
else {
  print "Not a .sif file\n";
  exit;
}

#creates hash of key = unique nodes, value = array of interacting nodes
sub CreateAdjacencyMatrix {
    my ($node1_ref, $node2_ref) = @_;
  my @node1 = @{ $node1_ref };
  my @node2 = @{ $node2_ref };

  #finding total unique nodes
  @totalnodes = @node1;
  @totalnodes = (@totalnodes, @node2);
  %totalnodes = map { $_ => 1 } @totalnodes;
  @totaluniqnodes = keys %totalnodes;

  #finding number of unique interactions
  $numinteractions = scalar @node1;

  #making a hash %totalppi in which every unique node is a key, and values are the interacting proteins
  %totalppi = ();
  foreach my $node (@totaluniqnodes){
    for (my $i = 0; $i < $numinteractions; $i++){
      if (($node eq $node1[$i]) && ($node ne $node2[$i])) {
        push @{$totalppi{$node}}, $node2[$i];
      }
      elsif (($node eq $node2[$i]) && ($node ne $node1[$i])) {
        push @{$totalppi{$node}}, $node1[$i];
      }
    }
  }
  return %totalppi;
  return @totaluniqnodes;
}
CreateAdjacencyMatrix(\@node1, \@node2);


#Computes and prints the degree distribution of the network
sub  degreedist {
    my ($hash_ref, $array_ref) = @_;
  my %hash = %{ $hash_ref };
  my @array = @{ $array_ref };
  my $arraylen = scalar @array;
  %degreek = ();

  #creating new hash of key = degrees(k), values = nodes that have that degree (protein name)
  foreach my $node (@array){
    push @{$degreek{scalar @{$hash{$node}}}}, $node;
  }

  #prints out the nuerically sorted key(node1), and size of value array containing all of the interacting proteins(node2)
  print "\nDegree Distribution\n";
  print "k\tpk\n";
  print("-------------------\n");
  foreach my $key (sort {$a<=>$b} keys %degreek){
	   print "$key\t";
     for (my $i = 0; $i < scalar @{$degreek{$key}}; $i++){
       print "#";
     }
     print " (".scalar @{$degreek{$key}}.")";
     print "\n";
  }
  print "\n\n";
  return %degreek;
}
print "\n**********  For $filename  **********\n";
degreedist(\%totalppi, \@totaluniqnodes);


#Computes and prints the average clustering coefficient
sub AVG_C {
    my ($hash_ref, $array_ref, $deghash_ref) = @_;
  my %hash = %{ $hash_ref };
  my @uniqarray = @{ $array_ref };
  my %deghash = %{ $deghash_ref };
  my %nuhash = ();

  #calculating the number of neighbors
  foreach my $firstnode (@uniqarray){
    my $ku1 = scalar @{$hash{$firstnode}};
    my @sourcenodearray = @{$hash{$firstnode}};
    my $nu = 0;
    if ($ku1 > 1){ #for nodes that have a degree greater than 1
      foreach my $secondnode (@sourcenodearray){
        my @possibleneighbors = @{$hash{$secondnode}};
        my $ku2 = scalar @possibleneighbors;
        for (my $i = 0; $i < $ku1; $i++){
          for (my $j = 0; $j < $ku2; $j++){
            if ($possibleneighbors[$j] eq @sourcenodearray[$i]){
              $nu++
            }
          }
        }
      }
      $nu = $nu/2; #removes duplicates (i.e a->b, b->a)
      $nuhash{$firstnode} = $nu;
    }
    else { #has 0 neighbors, since they only have one interacting node (k=1)
      $nuhash{$firstnode} = $nu;
    }
  }

  #finding network average clustering coefficient
  my $totalC = 0;
  foreach my $firstnode (@uniqarray){
    my $nu = $nuhash{$firstnode};
    my $ku = scalar @{$hash{$firstnode}};
    if($ku > 1){
      my $Cu = ((2 * $nu)/($ku * ($ku - 1)));
      $totalC += $Cu;
    }
    else { #sets clustering coeffiecient to 0 for nodes with no neighbors
      $Cu = 0;
    }
  }
  print "The network AVG_C = ".$totalC/(scalar @uniqarray)."\n\n\n";

  #finding average clustering coefficient for each k value in the network
  foreach $k (sort {$a <=> $b} keys %deghash){
    my $totalkC = 0;
    foreach $protein (@{$deghash{$k}}){
      my $nu = $nuhash{$protein};
      if($k > 1){ #ignoring nodes that have k = 1
        my $Cu = ((2 * $nu)/($k * ($k - 1)));
        $totalkC += $Cu;
      }
      else {
        $Cu = 0;
        $totalkC += $Cu;
      }
    }
    print "The AVG_C (k=$k) is ".$totalkC/(scalar @{$deghash{$k}})."\n";
  }
  print "\n\n";
}
AVG_C(\%totalppi, \@totaluniqnodes, \%degreek);

#creating a random graph using the same number of nodes as kshv.sif (n = 59)
sub randomGraph {
    my ($numprob_ref, $numnodes_ref) = @_;
  my $numprob = 10 * shift;
  my $numnodes = shift;
  my %randnetwork = ();
  my $randcounter = 0;
  @randnode1 = ();
  @randnode2 = ();

  #set the seed (optional: if you want randomly generated numbers everytime, take this out)
  srand(1000);

  #creating random node interactions (using input probability) without duplicate or self interactions
  #format: [node 1, node 2] [1||0]
  for(my $i = 0; $i < $numnodes; $i++){
    for (my $j = $i; $j < $numnodes; $j++){ #avoids creating duplicate interactions
      if ($i ne $j){ #avoids creating self interactions
        if ((1 + int rand(10)) <= $numprob){
          $randnetwork{($i+1).",".($j+1)} = 1;
        }
        else{
          $randnetwork{($i+1).",".($j+1)} = 0;
        }
      }
    }
   }

  #organizing interacting nodes (1 is interacting, 0 is not interacting) from random graph in the format to be used in CreateAdjacencyMatrix subroutine
  foreach my $nodepair (keys %randnetwork){
    if ($randnetwork{$nodepair} == 1){
      my @nodes = (split /,/, $nodepair); #spliting interacting nodes into sepearate arrays
      $randnode1[$randcounter] = $nodes[0];
      $randnode2[$randcounter] = $nodes[1];
      $randcounter++;
    }
  }
  return @randnode1;
  return @randnode2;
}

#for probability = 0.5
print "\n******  For random graph p=0.5  ******\n";
randomGraph(0.5, scalar @totaluniqnodes);
CreateAdjacencyMatrix (\@randnode1, \@randnode2);
degreedist(\%totalppi, \@totaluniqnodes);
AVG_C(\%totalppi, \@totaluniqnodes, \%degreek);

#for probability = 0.1
print "\n******  For random graph p=0.1  ******\n";
randomGraph(0.1, scalar @totaluniqnodes);
CreateAdjacencyMatrix (\@randnode1, \@randnode2);
degreedist(\%totalppi, \@totaluniqnodes);
AVG_C(\%totalppi, \@totaluniqnodes, \%degreek);

#for probability = 0.8
print "\n******  For random graph p=0.8  ******\n";
randomGraph(0.8, scalar @totaluniqnodes);
CreateAdjacencyMatrix (\@randnode1, \@randnode2);
degreedist(\%totalppi, \@totaluniqnodes);
AVG_C(\%totalppi, \@totaluniqnodes, \%degreek);
