#!/bin/bash

# Needs to be run for each chromosome

# Variables
CHR=Chr1;
IN1=/path/to/pop_${CHR}_phased.xpehh.out;

# Software
NORM=/path/to/norm

${NORM} --xpehh --files ${IN1} --bp-win --winsize 10000 > ${CHR}_norm_10Kb_win.log 2>&1
