#!/usr/bin/env python3
"""
Fetch selected chromosomes from a remote FASTA using its .fai index, with HTTP
range requests. No samtools required.

A .fai line is:  NAME  LENGTH  OFFSET  LINEBASES  LINEWIDTH
so the byte span of a contig is computable exactly, and the server is asked for
only that span. The output keeps the source line wrapping, which means a correct
.fai for the new file can be written directly rather than re-derived.

Usage
  python3 scripts/02_fetch_ref.py --chroms chr1 chr17 --out ref/GRCh38.panel_chroms.fa
"""

import argparse
import os
import sys
import urllib.request

REF_URL = ("https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/"
           "GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa")


def fetch(url, start=None, end=None, retries=4):
    """GET url, optionally a byte range [start, end] inclusive."""
    for attempt in range(retries):
        req = urllib.request.Request(url)
        if start is not None:
            req.add_header("Range", "bytes=%d-%d" % (start, end))
        try:
            with urllib.request.urlopen(req, timeout=120) as r:
                if start is not None and r.status != 206:
                    raise RuntimeError(
                        "server ignored the Range header (status %d). "
                        "Range requests are required for this script." % r.status)
                return r.read()
        except Exception as e:
            if attempt == retries - 1:
                raise
            sys.stderr.write("  retry %d after %s\n" % (attempt + 1, e))


def contig_bytes(length, linebases, linewidth):
    """Byte span of a contig's sequence block, newlines included."""
    nl = linewidth - linebases
    full, rem = divmod(length, linebases)
    return full * linewidth + (rem + nl if rem else 0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--url", default=REF_URL)
    ap.add_argument("--chroms", nargs="+", required=True)
    ap.add_argument("--out", default="ref/GRCh38.panel_chroms.fa")
    ap.add_argument("--cache", default="ref/chroms")
    args = ap.parse_args()

    os.makedirs(args.cache, exist_ok=True)
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)

    sys.stderr.write("reading remote .fai\n")
    fai_txt = fetch(args.url + ".fai").decode()
    index = {}
    for line in fai_txt.splitlines():
        f = line.split("\t")
        if len(f) >= 5:
            index[f[0]] = (int(f[1]), int(f[2]), int(f[3]), int(f[4]))

    missing = [c for c in args.chroms if c not in index]
    if missing:
        sys.exit("not in the reference index: %s" % missing)

    # fetch each chromosome once, cached on disk so a failed run resumes
    for c in args.chroms:
        path = os.path.join(args.cache, c + ".seq")
        length, offset, lb, lw = index[c]
        want = contig_bytes(length, lb, lw)
        if os.path.exists(path) and os.path.getsize(path) == want:
            sys.stderr.write("  %-6s cached\n" % c)
            continue
        sys.stderr.write("  %-6s fetching %.1f Mb ...\n" % (c, want / 1e6))
        data = fetch(args.url, offset, offset + want - 1)
        if len(data) != want:
            sys.exit("%s: expected %d bytes, got %d" % (c, want, len(data)))
        bases = len(data) - data.count(b"\n")
        if bases != length:
            sys.exit("%s: expected %d bases, got %d" % (c, length, bases))
        with open(path, "wb") as fh:
            fh.write(data)

    # concatenate and write our own .fai as we go
    sys.stderr.write("writing %s\n" % args.out)
    fai_lines = []
    with open(args.out, "wb") as out:
        pos = 0
        for c in args.chroms:
            header = (">%s\n" % c).encode()
            out.write(header)
            pos += len(header)
            length, _, lb, lw = index[c]
            fai_lines.append("%s\t%d\t%d\t%d\t%d" % (c, length, pos, lb, lw))
            with open(os.path.join(args.cache, c + ".seq"), "rb") as fh:
                while True:
                    chunk = fh.read(1 << 22)
                    if not chunk:
                        break
                    out.write(chunk)
                    pos += len(chunk)

    with open(args.out + ".fai", "w") as fh:
        fh.write("\n".join(fai_lines) + "\n")

    total = sum(index[c][0] for c in args.chroms)
    sys.stderr.write("done: %d contigs, %.2f Gb, %s\n"
                     % (len(args.chroms), total / 1e9, args.out))
    sys.stderr.write("wrote %s.fai\n" % args.out)


if __name__ == "__main__":
    main()
