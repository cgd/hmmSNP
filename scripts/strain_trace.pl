#!/usr/bin/perl 

# Programmer:  Glen Beane
# Purpose:     perl script to add strain path to graphviz file
# Overview:    reads in graphviz .dot file and hmm path file and writes 
#               new .dot file with the specified strain(s) path(s) highlighted
# Usage:       strain_trace.pl -h <num_haplotypes> -g <graph.dot> \ 
#                              -p <path.csv> -o <out.dot> [-c] <strain nums>
# Inputs:      -h - number of haplotypes used in hmm
#              -g - input graphviz file
#              -p - path file from hmmSNP
#              -o - path to output file (will be created)
#              -n (no argument) - keep original node coloring
#              -e (no argument) - keep original edge coloring
#              <strain nums> space separated list of strain numbers to trace
#                  strain numbers start at 1
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


@colors = ("cyan", "pink2", "khaki", "sienna", "skyblue", "magenta", "orange", "blue4", "red", "green", "mediumpurple", "purple4", "steelblue", "black", "cyan4", "darkgreen");

#@colors = ("cyan4", "blue4", "red", "yellow", "orange", "mediumpurple3", "skyblue3", "pink");

%options=();
getopts("g:p:o:h:ne", \%options);

if (defined $options{g} && defined $options{p} && defined $options{o} && 
    defined $options{h}) {
   print "graphviz file: :   $options{g}\n";
   print "Path File:         $options{p}\n";
   print "Output File:       $options{o}\n";
   print "#Haplotypes:       $options{h}\n";
}
else {
   print "you must specify -h -g -p and -o\n";
   exit(1);
}

if (defined $options{n}) {
   $keep_node_colors = 1;
}
else {
   $keep_node_colors = 0;
}

if (defined $options{e}) {
   $keep_edge_colors = 1;
}
else {
   $keep_edge_colors = 0;
}

open(GRAPH, "<$options{g}") or die("Unable to open graph file");
open(OUT, ">$options{o}") or die("Unable to open output file");

$haplotypes  = $options{h};

if (@ARGV > @colors) {
   print "current version of script only supports up to $colors paths\n";
   exit(1);
}


#first copy the nodes and edges from the original graph viz file
# (everything but the closing curly brace)
while ($line = <GRAPH> ) {

   if ($line !~ /\s*}\s*/) {

      
     
      if (!$keep_node_colors) {
         $line =~ s/style=filled//;
         $line =~ s/color=grey[0-9]+//;
      }
      if (!$keep_edge_colors) {
         $line =~ s/"setlinewidth\([0-9]+\.[0-9]+\)"/invis/;
      }
     

     
      print OUT $line;

      if ($line =~ /\s*{\s*/) {
         print OUT "node [fontsize=18];node [style=\"setlinewidth(2.0)\"]\n";
      }


   }
   else
   {
      last;
   }
}

print OUT "\n";

#now we need to put the edges for each strain

for ($strain_num = 0; $strain_num < @ARGV; $strain_num++) {
   print "Doing Strain $ARGV[$strain_num]\n";
   open(PATH, "<$options{p}") or die("Unable to open path file");
   @path = ();

   $strain = $ARGV[$strain_num] - 1;

   $header = <PATH>;

   @split_header = split(",", $header);

   print "$split_header[$strain + 6]\n";

   while ($line = <PATH>)
   {

      @split_line = split(",", $line);
      $hap = $split_line[$strain + 6]; 
      $hap =~ s/\s*//;
      if ($hap =~ /^[\d]$/)
      {
         push(@path, $hap);
      }
   }
   
   close(PATH);


   $style = "[style=\"setlinewidth(4.0)\" color=$colors[$strain_num]]";

   $prev = $path[0] + 2;
   for ($i = 1; $i < @path; $i++) {
      $current =  2 + $path[$i] + ($i  * $haplotypes);
      print OUT "$prev -> $current $style;\n";
      $prev = $current;
   }

   print OUT "\n";

}

print OUT "}";
print "done\n";

close(OUT);
close(GRAPH);
