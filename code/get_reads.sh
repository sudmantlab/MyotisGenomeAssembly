#!/bin/bash

ls data/PacBio-HiFi/M_*/data2/pb/*/*_*/*.subreads.bam | sed "s/data\/PacBio-HiFi\///; s/data2\/pb\///; s/.subreads.bam//; s/\//\t/g"
