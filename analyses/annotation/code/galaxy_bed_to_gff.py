#!/usr/bin/env python

# this is a heavily changed version of the script from Galaxy to deal with TOGA.

from __future__ import print_function

import sys

assert sys.version_info[:2] >= (2, 6)


def __main__():
    input_name = sys.argv[1]
    output_name = sys.argv[2]
    skipped_lines = 0
    first_skipped_line = 0
    i = 0
    with open(output_name, "w") as out, open(input_name) as fh_in:
        for i, line in enumerate(fh_in):
            complete_bed = False
            line = line.rstrip("\r\n")
            if line and not line.startswith("#") and not line.startswith("track") and not line.startswith("browser"):
                try:
                    elems = line.split("\t")
                    if len(elems) == 12:
                        complete_bed = True
                    chrom = elems[0]
                    if complete_bed:
                        feature = "mRNA"
                    else:
                        try:
                            feature = elems[3]
                        except Exception:
                            feature = "feature%d" % (i + 1)
                    start = int(elems[1]) + 1
                    end = int(elems[2])
                    try:
                        score = elems[4]
                    except Exception:
                        score = "0"
                    try:
                        strand = elems[5]
                    except Exception:
                        strand = "+"
                    try:
                        group = elems[3]
                    except Exception:
                        group = "group%d" % (i + 1)
                    if complete_bed:
                        out.write(
                            "%s\tTOGA\t%s\t%d\t%d\t%s\t%s\t.\tID=%s;Name=%s;\n"
                            % (chrom, "gene", start, end, score, strand, group + ".gene", group + '.TOGA')
                        )
                        out.write(
                            "%s\tTOGA\t%s\t%d\t%d\t%s\t%s\t.\tID=%s;Parent=%s;\n"
                            % (chrom, feature, start, end, score, strand, group, group + ".gene")
                        )
                    if complete_bed:
                        # We have all the info necessary to annotate exons for genes and mRNAs
                        block_count = int(elems[9])
                        block_sizes = elems[10].split(",")
                        block_starts = elems[11].split(",")
                        for j in range(block_count):
                            exon_start = int(start) + int(block_starts[j])
                            exon_end = exon_start + int(block_sizes[j]) - 1
                            out.write(
                                "%s\tTOGA\texon\t%d\t%d\t%s\t%s\t.\tID=%s_exon_%s;Parent=%s;\n"
                                % (chrom, exon_start, exon_end, score, strand, group, j, group)
                            )
                            out.write(
                                "%s\tTOGA\tCDS\t%d\t%d\t%s\t%s\t.\tID=%s_exon_%s_cds;Parent=%s;\n"
                                % (chrom, exon_start, exon_end, score, strand, group, j, group)
                            )
                except Exception:
                    skipped_lines += 1
                    if not first_skipped_line:
                        first_skipped_line = i + 1
            else:
                skipped_lines += 1
                if not first_skipped_line:
                    first_skipped_line = i + 1
    info_msg = "%i lines converted to GFF version 3.  " % (i + 1 - skipped_lines)
    if skipped_lines > 0:
        info_msg += "Skipped %d blank/comment/invalid lines starting with line #%d." % (
            skipped_lines,
            first_skipped_line,
        )
    print(info_msg)


if __name__ == "__main__":
    __main__()


