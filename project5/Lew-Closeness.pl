#!/usr/bin/perl
#use Data::Dumper;
#use Data::Dumper qw(Dumper);

$filename = $ARGV[0];

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

#The subroutine Create-Adjacency-Matrix takes as input argument the .sif representation
#of a graph and generates the adjacency matrix representation of the graph. It uses hashes
#(as defined in perl) to map the protein ids into indexes of the adjacency matrix.
sub createAdjacencyMatrix {
  my ($node1_ref, $node2_ref) = @_;
my @node1 = @{ $node1_ref };
my @node2 = @{ $node2_ref };

#finding total unique nodes
@totalnodes = (@node1, @node2);
%totalnodes = map { $_ => 1 } @totalnodes;
@totaluniqnodes = sort keys %totalnodes; #creating array of unique nodes

# hash of uniq nodes and their position
%totaluniqnodes;
my $nodecounter = 0;
foreach my $uniqnode(@totaluniqnodes){
  $totaluniqnodes{$uniqnode} = $nodecounter;
  $nodecounter++;
}

#finding number of edges
$numinteractions = scalar @node1;
my $counter = 0;

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

#intializing matrix
@matrix;
push @matrix, [("x") x (scalar @totaluniqnodes)] for 1 .. (scalar @totaluniqnodes);

#creating matrix
for (my $row = 0; $row < (scalar @totaluniqnodes); $row++){
  for (my $col = 0; $col < (scalar @totaluniqnodes); $col++){
    for (my $node = 0; $node < $numinteractions; $node++){
      if ($totaluniqnodes[$row] eq $totaluniqnodes[$col]) { #if row and col is the same protein
        $matrix [$col][$row] = "-";
      }
      if (($totaluniqnodes[$row] eq $node1[$node]) && ($totaluniqnodes[$col] eq $node2[$node])) { #if row and col have an interaction
        $matrix [$col][$row] = 1;
      }
      if (($totaluniqnodes[$row] eq $node1[$node]) && ($totaluniqnodes[$col] eq $node2[$node])) { #if row and col have an interaction
        $matrix [$col][$row] = 1;
      }
      if (($totaluniqnodes[$row] eq $node2[$node]) && ($totaluniqnodes[$col] eq $node1[$node])) { #if row and col have an interaction
        $matrix [$col][$row] = 1;
      }
      if (($totaluniqnodes[$row] eq $node2[$node]) && ($totaluniqnodes[$col] eq $node1[$node])) { #if row and col have an interaction
        $matrix [$col][$row] = 1;
      }
      if ($matrix [$col][$row] eq x){ #if row and col do not have an interaction
          $matrix [$col][$row] = 0;
      }
    }
  }
}

return @matrix;
return %totalppi;
return @totaluniqnodes;
}
createAdjacencyMatrix(\@node1, \@node2);

#Computes the degree of all nodes, then prints the top 10 nodes of highest degrees.
sub  Degree {
    my ($hash_ref, $array_ref) = @_;
  my %hash = %{ $hash_ref };
  my @array = @{ $array_ref };
  my $arraylen = scalar @array;
  %degreek;
  my $counter = 0;

  #creating new hash of key = degrees(k), values = nodes that have that degree (protein name)
  foreach my $node (@array){
    $degreek{$node} = scalar @{$hash{$node}};
  }

  #prints out the top 10 nodes (values) with highest degrees (key)
  print "Top 10 Nodes with Highest Degree\n";
  print "----------------------------------\n";
  foreach my $key (sort {$degreek{$b} <=> $degreek{$a}} keys %degreek){
    if ($counter == 10){
      last;
    }
    print "Node: $key --- Degree: $degreek{$key}\n";
    $counter++;
  }
  print "\n\n";
  return %degreek;
}
Degree(\%totalppi, \@totaluniqnodes);


#The subroutine AllPairsShortestPaths has as input argument the adjacency matrix
#representation of a graph. The subroutine returns a two-dimensional matrix, D, with D[u,v]
#giving the length of the shortest path between nodes u and v.
sub AllPairsShortestPaths {
    my ($matrixarray_ref, $totaluniqnodes_ref) = @_;
  my @matrix = @{ $matrixarray_ref };
  my @totaluniqnodes = @{ $totaluniqnodes_ref };
  my $counter = 0;

  #intializing new path matrix
  @path;
  push @path, [("inf") x (scalar @totaluniqnodes)] for 1 .. (scalar @totaluniqnodes);

  #initializing new predecessor matrix (for subroutine GetShortestPath)
  @predecessor;
  push @predecessor, [("inf") x (scalar @totaluniqnodes)] for 1 .. (scalar @totaluniqnodes);

  for (my $k = 0; $k < scalar @totaluniqnodes; $k++){
    for (my $i = 0; $i < scalar @totaluniqnodes; $i++){
      for (my $j = 0; $j < scalar @totaluniqnodes; $j++){
        if ($matrix[$i][$j] eq "-"){
          $path[$i][$j] = 0;
        }
        if ($matrix[$i][$j] eq 1){
          $path[$i][$j] = 1;
        }
        if ($path[$i][$j] > $path[$i][$k] + $path[$k][$j]) {
          $path[$i][$j] = ($path[$i][$k] + $path[$k][$j]);
          $predecessor[$i][$j] = "$totaluniqnodes[$k]";
        }
      }
    }
  }
  return @path;
  return @predecessor;
}
AllPairsShortestPaths(\@matrix, \@totaluniqnodes);

#subroutine NodeCentrality computes the closeness centrality of all nodes and returns the top 10
sub NodeCentrality {
    my ($path_ref, $totaluniqnodesarr_ref) = @_;
  my @path = @{ $path_ref };
  my @totaluniqnodes = @{ $totaluniqnodesarr_ref };

  %cc;

  for (my $i = 0; $i < scalar @totaluniqnodes; $i++){
    my $sumdist = 0;
    for (my $j = 0; $j < scalar @totaluniqnodes; $j++){
      if ($path[$i][$j] ne "-"){
        $sumdist += $path[$i][$j];
      }
    }
    $cc{$totaluniqnodes[$i]} = 1/$sumdist;
  }

  #finding nodes with highest and lowest CC (to use for subroutine GetShortestPath)
  @sortedCCnodes = sort { $cc{$b} <=> $cc{$a}} keys %cc;
  $highestCCnode = $sortedCCnodes[0];
  $lowestCCnode = $sortedCCnodes[-1];

  #finding top 10 nodes with highest CC
  print "Top 10 Nodes with Highest Closeness Centrality (CC)\n";
  print "----------------------------------------------------\n";
  foreach my $key (sort {$cc{$b} <=> $cc{$a}} keys %cc){
    if ($counter == 10){
      last;
    }
    print "Node: $key --- CC: $cc{$key}\n";
    $counter++;
  }
  print "\n\n";
  return %cc;
  return $highestCCnode;
  return $lowestCCnode;
}
NodeCentrality(\@path, \@totaluniqnodes);

#The subroutine Get-Shortest-Path takes in input a pair of nodes and the adjacency matrix
#representation of a graph and returns a shortest path between the two nodes. Call this
#subroutine from the main using as input parameters the two nodes with highest and lowest
#value of node-centrality.
sub GetShortestPath {
    my ($node1_ref, $node2_ref, $path_ref, $hash_ref, $path2_ref) = @_;
  my $node1 = $node1_ref;
  my $node2 = $node2_ref;
  my @predecessor = @{ $path_ref };
  my %nodehash = %{ $hash_ref };
  my @path = @{ $path2_ref };

  $node1index = $nodehash{$node1};
  $node2index = $nodehash{$node2};

  my @pathway; #new array to hold the pathway

  #finding the shortest path from node1 to node 2
  until ($path[$node1index][$node2index] == 1){
    $copynode2 = $node2;
    unshift @pathway, "$predecessor[$node1index][$node2index] -> ";
    $copynode2 = $predecessor[$node1index][$node2index];
    $node2index = $nodehash{$copynode2};
  }

  print "Shortest path from $node1 to $node2:\n$node1 -> ";
  print "@pathway";
  print "$node2\n\n";
}
GetShortestPath ($highestCCnode, $lowestCCnode, \@predecessor, \%totaluniqnodes, \@path);
