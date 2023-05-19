#!/usr/bin/env python3

import sys

with open("data/input/uniprot-myoLuc-transtable.tsv") as tt: 
    tdict = {l[0]: l[1] for l in [line.strip().split("\t") for line in tt]}

with open("output/myoLuc3/AvA-pcScore0.1_pcIdent0.8_pcQuerySpan0.5_reverse-myoLuc2-RecSearchReady/output/myoLuc3_RecBlastOutput.bed") as inbed, open("output/myoLuc3/AvA-pcScore0.1_pcIdent0.8_pcQuerySpan0.5_reverse-myoLuc2-RecSearchReady/output/myoLuc3_RecBlastOutput.newname.bed", "w") as outbed:
    for line in inbed:
        line = line.strip().split("\t")
        id_full = line[3].split("_")
        if id_full[1] in tdict.keys():
            line[3] = tdict[id_full[1]] + "_" + id_full[-1]
        else:
            error = "Item " + line[3] + " is missing in translation table!"
            print(error, file=sys.stderr)
        line = "\t".join(line) + "\n"
        outbed.write(line)

