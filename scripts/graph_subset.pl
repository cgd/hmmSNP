#!/usr/bin/perl -w


# Programmer:  Glen Beane
# Purpose:     perl script to trim hmmSNP graphviz output to a subset of SNPs
# Overview:    reads in graphviz .dot file and writes new .dot file with only 
#              the specified SNPs
# Usage:       graph_subset.pl -h <num_haplotypes> -s <start_SNP> -e <end_SNP> \
#                              -i inFile.dot -o outFile.dot
# Inputs:      -h - number of haplotypes used in hmm
#              -s - starting SNP #
#              -e - end SNP #
#              -i - path to .dot file output from hmmSNP
#              -o - path to output file (will be created)
#              
# Assumptions: Assumes that there will be whitespace between node number 
#              and node attributes, and on each side of "->" in edge list
#
#

#  Copyright (c) 2010 The Jackson Laboratory
# 
#  This software was developed by Gary Churchill's Lab at The Jackson
#  Laboratory (see http://research.jax.org/faculty/churchill).
# 
#  This is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
# 
#  This software is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with this software. If not, see <http://www.gnu.org/licenses/>.

use Getopt::Std;

%options=();
getopts("h:s:e:i:o:", \%options);

if (defined $options{h} && defined $options{s} && defined $options{e}
    && defined $options{i} && defined $options{o}) {
   print "num haplotypes:  $options{h}\n";
   print "starting SNP#:   $options{s}\n";
   print "ending SNP#:     $options{e}\n";
   print "Input File:      $options{i}\n";
   print "Output File:     $options{o}\n";
}
else {
   print "you must specify -h -s -e -i and -o\n";
}


$start_node = 2 + $options{h} * ($options{s} - 1);
$end_node =  1 + $options{h} * ($options{e} );




print "\n\ncopying nodes $start_node through $end_node\n";

open(IN, "<$options{i}") or die("Unable to open input file");
open(OUT, ">$options{o}") or die("Unable to open output file");

while ($line = <IN>) {

   @split_line = split(/\s+/, $line);
   $node = $split_line[0];
   if ($node && $node =~ /^\d+$/) {
      if ($split_line[1] =~ /^->$/) {
         $node2 = $split_line[2];
         if ($node >= $start_node && $node <= $end_node
             && $node2 >= $start_node && $node2 <= $end_node) {
            print OUT $line; 
         }
      }
      elsif ($node >= $start_node && $node <= $end_node) {
         print OUT $line;
      }
   }
   else {
      print OUT $line;
   }
}

close OUT;
close IN;

exit 0;
