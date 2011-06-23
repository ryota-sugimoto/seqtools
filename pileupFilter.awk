#!/usr/bin/env awk -f

BEGIN {
	FS = "\t"
}
( $6 >= 20 ) && ( $8 >= 20 ) { 
	print
}
