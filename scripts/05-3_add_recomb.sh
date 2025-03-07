#!/usr/bin/env sh

# LocusZoom directory on Eurac server
app=/usr/local/stow/locuszoom-1.4

# Downloaded lifted over merged map file
map=~/projects/gwas/pairways_LD/recomb-hg38

echo python2 ${app}/bin/dbmeister.py   \
	--db ${app}/data/database/locuszoom_hg38.db   \
	--recomb_rate ${map}/genetic_map_GRCh38_merged.tab.headered

