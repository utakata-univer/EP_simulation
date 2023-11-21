#!/bin/csh

alias ave "awk '/GMES/{print .5\*(\$5+\$6)}' \!* | ~/preon/../doug/bin/averblock 10 -t"