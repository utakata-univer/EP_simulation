#!/usr/local/bin/perl

# File: find_beta.pl
# Created: Fri Jun  9 1995
# Author: J. E. Hetrick <hetrick@corelli.physics.Arizona.EDU>
#
# Description: utility for picking out measurements at a particular
#              beta from multi-input set data.
#
# Usage: find_beta.pl datafile betavalue

$DATA_ON = 0;
open(DATA, $ARGV[0]);
$beta = $ARGV[1];

while(<DATA>) {
	if(($DATA_ON == 1) && (/plaq. beta = /)) { 
		$DATA_ON = 0;
		next;
	}
	if(($DATA_ON == 0) && (/plaq. beta = $beta/)) { $DATA_ON = 1; }
	if($DATA_ON == 1) {
		print if(/GMES/);
	}
}

