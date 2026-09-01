#!/usr/bin/env python3
"""Generate small GRCh38 BAMs with WGS-like Illumina evidence for teaching/IGV.

The BAMs contain only the requested loci.  They are not whole-genome BAMs, but
their read length, depth, qualities and pairing emulate a standard 2x150 WGS.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import random
import tempfile
import time
import urllib.parse
import urllib.request
from dataclasses import dataclass
from pathlib import Path

import pysam


SEED = 20260830
READ_LEN = 150
CONTIGS = {
    "chr2": ("NC_000002.12", 242193529),
    "chr3": ("NC_000003.12", 198295559),
    "chr6": ("NC_000006.12", 170805979),
    "chr8": ("NC_000008.11", 145138636),
    "chr11": ("NC_000011.10", 135086622),
    "chr13": ("NC_000013.11", 114364328),
}


@dataclass(frozen=True)
class Region:
    chrom: str
    start0: int
    end0: int


REGIONS = {
    # RefSeq gene bounds plus ~1 kb flanks (0-based, half-open).
    "brca2": Region("chr13", 32314076, 32401268),
    "cyp21a2": Region("chr6", 32037414, 32042644),
    "cyp21a1p": Region("chr6", 32004600, 32009450),
    "mnv": Region("chr8", 93753843, 93819121),       # TMEM67
    "cyp1b1": Region("chr2", 38066508, 38077151),
    "btd": Region("chr3", 15600360, 15723516),
    "del214": Region("chr11", 66256086, 66268860),  # KLC2 + deletion flank
    "dup": Region("chr6", 128882137, 129517566),    # LAMA2
}

# Synthetic teaching SNVs (0-based positions). KLC2 variants inside the
# heterozygous deletion are placed on the intact HP2 allele; flanking variants
# are placed on HP1. LAMA2 variants are placed on duplicated HP1, yielding 2/3
# ALT copies in a total copy number of three.
KLC2_HEMIZYGOUS_SNVS = [66257121, 66257191, 66257261]
KLC2_FLANKING_HET_SNVS = [
    66256336, 66256666, 66256906, 66257640, 66257720, 66258050,
]
DUPLICATION_BAF_SNVS = [
    129060122, 129075320, 129090516, 129105710,
    129120906, 129136102, 129151298, 129166494,
    129181690, 129196886, 129212082, 129227278,
    129242474, 129257670, 129272866, 129288062,
]


def alternate_base(ref: str) -> str:
    return {"A": "C", "C": "A", "G": "T", "T": "G"}[ref]


def snv_map(genome: "Genome", chrom: str, positions: list[int]) -> dict[int, str]:
    result = {}
    for pos0 in positions:
        ref = genome.fetch(chrom, pos0, pos0 + 1)
        if ref not in "ACGT":
            raise ValueError(f"Cannot create SNV at {chrom}:{pos0 + 1}; REF={ref}")
        result[pos0] = alternate_base(ref)
    return result


def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


class Genome:
    def __init__(self, cache_dir: Path, fasta: str | None):
        self.cache_dir = cache_dir
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.fasta = pysam.FastaFile(fasta) if fasta else None
        self.loaded: dict[str, tuple[Region, str]] = {}

    def preload(self) -> None:
        for key, region in REGIONS.items():
            self.loaded[key] = (region, self._load_region(key, region))

    def _load_region(self, key: str, region: Region) -> str:
        if self.fasta:
            seq = self.fasta.fetch(region.chrom, region.start0, region.end0).upper()
        else:
            cache = self.cache_dir / f"{key}.{region.chrom}.{region.start0 + 1}-{region.end0}.fa.gz"
            if cache.exists():
                with gzip.open(cache, "rt") as handle:
                    seq = "".join(x.strip() for x in handle if not x.startswith(">"))
            else:
                accession = CONTIGS[region.chrom][0]
                params = urllib.parse.urlencode({
                    "db": "nuccore", "id": accession,
                    "seq_start": region.start0 + 1, "seq_stop": region.end0,
                    "rettype": "fasta", "retmode": "text",
                })
                url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + params
                last_error = None
                for attempt in range(5):
                    try:
                        with urllib.request.urlopen(url, timeout=120) as response:
                            text = response.read().decode()
                        seq = "".join(x.strip() for x in text.splitlines() if not x.startswith(">"))
                        if len(seq) != region.end0 - region.start0:
                            raise RuntimeError(f"NCBI returned {len(seq)} bp, expected {region.end0-region.start0}")
                        with gzip.open(cache, "wt") as handle:
                            handle.write(f">{region.chrom}:{region.start0 + 1}-{region.end0}\n{seq}\n")
                        break
                    except Exception as exc:  # retry transient NCBI errors
                        last_error = exc
                        time.sleep(2 ** attempt)
                else:
                    raise RuntimeError(f"Could not fetch {key} from NCBI: {last_error}")
        if len(seq) != region.end0 - region.start0:
            raise ValueError(f"Wrong sequence length for {key}")
        return seq

    def fetch(self, chrom: str, start0: int, end0: int) -> str:
        if self.fasta:
            return self.fasta.fetch(chrom, start0, end0).upper()
        for region, seq in self.loaded.values():
            if region.chrom == chrom and region.start0 <= start0 and end0 <= region.end0:
                return seq[start0 - region.start0:end0 - region.start0]
        raise KeyError(f"No cached sequence covers {chrom}:{start0 + 1}-{end0}")


def header(sample: str) -> dict:
    return {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": c, "LN": n, "AS": "GRCh38"} for c, (_, n) in CONTIGS.items()],
        "RG": [{"ID": sample, "SM": sample, "PL": "ILLUMINA", "LB": "WGS", "PU": "SYNTHETIC"}],
        "PG": [{"ID": "synthetic_variant_bams", "PN": "generate_bams.py", "VN": "1.0"}],
        "CO": ["Synthetic teaching data; continuous full-gene WGS-like evidence, not a whole-genome BAM."],
    }


def mutate_substitutions(seq: str, start0: int, subs: dict[int, str]) -> str:
    chars = list(seq)
    for pos0, alt in subs.items():
        idx = pos0 - start0
        if 0 <= idx < len(chars):
            chars[idx] = alt
    return "".join(chars)


def make_read(
    genome: Genome, chrom: str, start0: int, qname: str, sample: str,
    reverse: bool = False, read2: bool = False, mate_start0: int | None = None,
    template_length: int = 0, mapq: int = 60, length: int = READ_LEN,
    subs: dict[int, str] | None = None, hp: int | None = None,
    cigar: list[tuple[int, int]] | None = None, aligned_seq: str | None = None,
    variant_bq: dict[int, int] | None = None, extra_tags: list[tuple[str, object]] | None = None,
    proper: bool = True, nm: int | None = None,
) -> pysam.AlignedSegment:
    r = pysam.AlignedSegment()
    r.query_name = qname
    flag = 0x1 | (0x2 if proper else 0) | (0x80 if read2 else 0x40)
    if reverse:
        flag |= 0x10
    if mate_start0 is not None and not reverse:
        flag |= 0x20
    r.flag = flag
    r.reference_id = list(CONTIGS).index(chrom)
    r.reference_start = start0
    r.mapping_quality = mapq
    r.cigartuples = cigar or [(0, length)]
    if aligned_seq is None:
        aligned_seq = genome.fetch(chrom, start0, start0 + length)
        aligned_seq = mutate_substitutions(aligned_seq, start0, subs or {})
    # BAM/SAM stores SEQ in the orientation shown against the reference; flag
    # 0x10 records the original molecule orientation.
    query_seq = aligned_seq
    r.query_sequence = query_seq
    # Typical NovaSeq-like profile: high central qualities with a mild decline
    # at both read ends. Values vary per base/read but remain mostly Q30-Q40.
    quals = []
    for i in range(len(query_seq)):
        edge = min(i, len(query_seq) - 1 - i)
        penalty = max(0, 5 - edge) + (2 if i > len(query_seq) - 12 else 0)
        quals.append(max(22, min(40, int(round(random.gauss(37 - penalty, 1.4))))))
    for pos0, q in (variant_bq or {}).items():
        idx = pos0 - start0
        if 0 <= idx < len(aligned_seq):
            quals[idx] = q
    r.query_qualities = quals
    if mate_start0 is not None:
        r.next_reference_id = r.reference_id
        r.next_reference_start = mate_start0
    r.template_length = template_length
    if nm is None:
        ref_consumed = sum(n for op, n in r.cigartuples if op in (0, 2, 3, 7, 8))
        nm = sum(start0 <= p < start0 + ref_consumed for p in (subs or {}))
    tags: list[tuple[str, object]] = [("RG", sample), ("NM", nm)]
    if hp is not None:
        tags.extend([("HP", hp), ("PS", 1)])
    if extra_tags:
        tags.extend(extra_tags)
    r.set_tags(tags)
    return r


def add_standard_pairs(
    out: pysam.AlignmentFile, genome: Genome, sample: str, chrom: str,
    center0: int, pairs: int, alt_fraction: float, subs: dict[int, str],
    fragment: int = 350, mapq: int = 60, low_bq: dict[int, int] | None = None,
) -> None:
    for i in range(pairs):
        # Keep the requested locus inside R1 while preserving a normal insert.
        start0 = center0 - READ_LEN // 2 + (i % 11) - 5
        mate0 = start0 + fragment - READ_LEN
        alt = i < round(pairs * alt_fraction)
        hp = 1 if alt else 2
        chosen = subs if alt else {}
        qname = f"{sample}_HP{hp}_{i:04d}"
        out.write(make_read(genome, chrom, start0, qname, sample, False, False, mate0,
                            fragment, mapq, subs=chosen, hp=hp, variant_bq=low_bq))
        out.write(make_read(genome, chrom, mate0, qname, sample, True, True, start0,
                            -fragment, mapq, subs=chosen, hp=hp, variant_bq=low_bq))


def best_hamming(query: str, target: str) -> tuple[int, int]:
    best_i, best_d = 0, len(query) + 1
    for i in range(0, len(target) - len(query) + 1):
        d = sum(a != b for a, b in zip(query, target[i:i + len(query)]))
        if d < best_d:
            best_i, best_d = i, d
    return best_i, best_d


def write_snv_good(path: Path, genome: Genome) -> None:
    sample, chrom, pos0 = "01_SNV_good_BRCA2", "chr13", 32319079
    assert genome.fetch(chrom, pos0, pos0 + 1) == "T"
    with pysam.AlignmentFile(path, "wb", header=header(sample)) as out:
        add_standard_pairs(out, genome, sample, chrom, pos0, 32, 0.5, {pos0: "G"})


def write_snv_homology(path: Path, genome: Genome) -> None:
    sample, chrom, pos0 = "02_SNV_low_quality_CYP21A2", "chr6", 32038609
    assert genome.fetch(chrom, pos0, pos0 + 1) == "A"
    pseudo_region, pseudo_seq = genome.loaded["cyp21a1p"]
    with pysam.AlignmentFile(path, "wb", header=header(sample)) as out:
        for i in range(28):
            start0 = pos0 - 70 + (i % 9) - 4
            mate0 = start0 + 210
            alt = i < 9  # ~0.32 apparent VAF
            hp = 1 if alt else 2
            subs = {pos0: "T"} if alt else {}
            qname = f"{sample}_AMBIG_HP{hp}_{i:04d}"
            for reverse, read2, read_start, mate_start, tlen in [
                (False, False, start0, mate0, 360), (True, True, mate0, start0, -360)
            ]:
                aligned = mutate_substitutions(genome.fetch(chrom, read_start, read_start + READ_LEN), read_start, subs)
                # Add one deterministic sequencing mismatch away from the variant.
                if i % 4 == 0:
                    j = 20 + (i % 30)
                    if read_start + j != pos0:
                        bases = "ACGT"
                        aligned = aligned[:j] + bases[(bases.index(aligned[j]) + 1) % 4] + aligned[j + 1:]
                primary = make_read(genome, chrom, read_start, qname, sample, reverse, read2,
                                    mate_start, tlen, 8 if alt else 12, aligned_seq=aligned,
                                    hp=hp, variant_bq={pos0: 9 if alt else 18},
                                    extra_tags=[("XA", "CYP21A1P_homolog")])
                primary.set_tag("NM", 1 if i % 4 == 0 else 0)
                out.write(primary)
                aligned_forward = aligned
                idx, nm = best_hamming(aligned_forward, pseudo_seq)
                secondary = make_read(genome, chrom, pseudo_region.start0 + idx, qname, sample,
                                      reverse, read2, None, 0, 0, aligned_seq=aligned,
                                      hp=hp, proper=False, extra_tags=[("NH", 2)])
                secondary.flag |= 0x100
                secondary.set_tag("NM", nm)
                out.write(secondary)


def deletion_read(genome: Genome, chrom: str, start0: int, read_len: int,
                  del_start0: int, del_end0: int, subs: dict[int, str] | None = None) -> tuple[str, list[tuple[int, int]]]:
    left = del_start0 - start0
    if not (1 <= left < read_len):
        raise ValueError("Deletion does not fall inside read")
    right = read_len - left
    seq = genome.fetch(chrom, start0, del_start0) + genome.fetch(chrom, del_end0, del_end0 + right)
    seq = mutate_substitutions(seq, start0, subs or {})
    return seq, [(0, left), (2, del_end0 - del_start0), (0, right)]


def write_mnv(path: Path, genome: Genome) -> None:
    sample, chrom = "03_complex_MNV_hom_cis_chr8", "chr8"
    del_start0, del_end0, snv0 = 93797350, 93797351, 93797351
    assert genome.fetch(chrom, 93797349, 93797352) == "TGC"
    with pysam.AlignmentFile(path, "wb", header=header(sample)) as out:
        for i in range(32):
            start0 = 93797300 + (i % 13) - 6
            seq, cigar = deletion_read(genome, chrom, start0, READ_LEN, del_start0, del_end0)
            # The SNV is the first reference base after the deletion.
            idx = del_start0 - start0
            seq = seq[:idx] + "T" + seq[idx + 1:]
            qname = f"{sample}_HOM_CIS_{i:04d}"
            out.write(make_read(genome, chrom, start0, qname, sample, False, False, None, 0,
                                60, aligned_seq=seq, cigar=cigar, hp=1,
                                extra_tags=[("MC", "components:93797350_TG>T;93797352_C>T")], nm=2))


def write_small_del(path: Path, genome: Genome) -> None:
    sample, chrom = "04_small_del_hom_CYP1B1", "chr2"
    d0, d1 = 38071277, 38071290
    assert genome.fetch(chrom, d0, d1) == "TCTGCCTGCACTC"
    lengths = [100, 125, 150, 150, 175, 250]
    with pysam.AlignmentFile(path, "wb", header=header(sample)) as out:
        for i in range(36):
            length = lengths[i % len(lengths)]
            start0 = d0 - (length // 2) + (i % 7) - 3
            seq, cigar = deletion_read(genome, chrom, start0, length, d0, d1)
            qname = f"{sample}_DEL_{length}bp_{i:04d}"
            out.write(make_read(genome, chrom, start0, qname, sample, i % 2 == 1,
                                i % 2 == 1, None, 0, 55, length=length,
                                aligned_seq=seq, cigar=cigar, hp=1, nm=13))
        # Deliberately clipped breakpoint evidence of several read sizes.
        for i, length in enumerate([90, 110, 130, 150, 175, 200]):
            left = length // 2
            start0 = d0 - left
            aligned = genome.fetch(chrom, start0, d0) + genome.fetch(chrom, d1, d1 + length - left)
            qname = f"{sample}_SOFTCLIP_{length}bp_{i:02d}"
            out.write(make_read(genome, chrom, start0, qname, sample, False, False, None, 0,
                                20, length=length, aligned_seq=aligned,
                                cigar=[(0, left), (4, length - left)], hp=1,
                                extra_tags=[("SC", "intentional_breakpoint_softclip")]))


def write_btd(path: Path, genome: Genome, cis: bool) -> None:
    sample = "06_BTD_cis_pair_view" if cis else "05_BTD_trans_pair_view"
    chrom = "chr3"
    left0 = 15644856 if cis else 15644916
    right0 = 15645185
    expected_left = "T" if cis else "C"
    assert genome.fetch(chrom, left0, left0 + 1) == expected_left
    assert genome.fetch(chrom, right0, right0 + 1) == "G"
    # Put the left locus in R1 and the right locus in R2 of the same fragment.
    fragment = (right0 - left0) + 145
    with pysam.AlignmentFile(path, "wb", header=header(sample)) as out:
        for i in range(30):
            hp = 1 if i < 15 else 2
            if cis:
                subs = {left0: "A", right0: "C"} if hp == 1 else {}
            else:
                subs = {left0: "T"} if hp == 1 else {right0: "C"}
            start0 = left0 - 70 + (i % 9) - 4
            mate0 = start0 + fragment - READ_LEN
            qname = f"{sample}_HP{hp}_{i:04d}"
            out.write(make_read(genome, chrom, start0, qname, sample, False, False,
                                mate0, fragment, 60, subs=subs, hp=hp))
            out.write(make_read(genome, chrom, mate0, qname, sample, True, True,
                                start0, -fragment, 60, subs=subs, hp=hp))


def write_del214(path: Path, genome: Genome) -> None:
    sample, chrom = "07_deletion_214bp_het", "chr11"
    d0, d1 = 66257086, 66257300
    with pysam.AlignmentFile(path, "wb", header=header(sample)) as out:
        # Reference allele reads.
        add_standard_pairs(out, genome, sample, chrom, (d0 + d1) // 2, 18, 0.0, {}, fragment=420)
        # Alternate junction-spanning reads.
        for i in range(18):
            start0 = d0 - 75 + (i % 11) - 5
            seq, cigar = deletion_read(genome, chrom, start0, READ_LEN, d0, d1)
            qname = f"{sample}_HP1_DEL_{i:04d}"
            out.write(make_read(genome, chrom, start0, qname, sample, i % 2 == 1,
                                i % 2 == 1, None, 0, 60, aligned_seq=seq,
                                cigar=cigar, hp=1, proper=False, nm=214))


def add_coverage_pairs(out: pysam.AlignmentFile, genome: Genome, sample: str,
                       chrom: str, start0: int, end0: int, depth: int, hp: int) -> None:
    # Approximate depth for paired 2x150 reads with a 350-bp fragment.
    step = max(1, int((2 * READ_LEN) / depth))
    i = 0
    for s0 in range(start0, end0 - 350, step):
        mate0 = s0 + 350 - READ_LEN
        qname = f"{sample}_COV_HP{hp}_{i:07d}"
        out.write(make_read(genome, chrom, s0, qname, sample, False, False, mate0, 350, 60, hp=hp))
        out.write(make_read(genome, chrom, mate0, qname, sample, True, True, s0, -350, 60, hp=hp))
        i += 1


def supplementary_segment(genome: Genome, sample: str, qname: str, chrom: str,
                          start0: int, aligned_seq: str, reverse: bool, sa: str) -> pysam.AlignedSegment:
    r = make_read(genome, chrom, start0, qname, sample, reverse, False, None, 0,
                  20, length=len(aligned_seq), aligned_seq=aligned_seq,
                  cigar=[(0, len(aligned_seq))], hp=1, proper=False,
                  extra_tags=[("SA", sa)])
    r.flag |= 0x800
    return r


def write_dup(path: Path, genome: Genome) -> None:
    sample, chrom = "08_tandem_dup_251kb_het", "chr6"
    d0, d1 = 129049897, 129300892  # assumed exact breakpoints from ClinVar inner bounds
    with pysam.AlignmentFile(path, "wb", header=header(sample)) as out:
        # 30x diploid baseline across the entire region, plus 15x inside = CN3.
        add_coverage_pairs(out, genome, sample, chrom, REGIONS["dup"].start0,
                           REGIONS["dup"].end0, 30, 2)
        add_coverage_pairs(out, genome, sample, chrom, d0, d1, 15, 1)
        # Tandem junction: right end followed by left start, represented as split reads.
        for i in range(16):
            flank = 75
            left_part = genome.fetch(chrom, d1 - flank, d1)
            right_part = genome.fetch(chrom, d0, d0 + flank)
            qname = f"{sample}_JUNCTION_HP1_{i:03d}"
            primary = make_read(genome, chrom, d1 - flank, qname, sample, False, False,
                                None, 0, 20, aligned_seq=left_part + right_part,
                                cigar=[(0, flank), (4, flank)], hp=1, proper=False,
                                extra_tags=[("SA", f"{chrom},{d0 + 1},+,75S75M,20,0;")])
            out.write(primary)
            supp = supplementary_segment(genome, sample, qname, chrom, d0,
                                         left_part + right_part, False,
                                         f"{chrom},{d1-flank+1},+,75M75S,20,0;")
            supp.cigartuples = [(4, flank), (0, flank)]
            out.write(supp)


# ---------------------------------------------------------------------------
# Gene-window WGS simulator (v2). These definitions intentionally replace the
# earlier locus-only writers above while retaining their low-level BAM helpers.
# ---------------------------------------------------------------------------

def basic_read_model(genome: Genome, chrom: str, subs_by_hp: dict[int, dict[int, str]]):
    def model(hp: int, start0: int):
        seq = genome.fetch(chrom, start0, start0 + READ_LEN)
        subs = subs_by_hp.get(hp, {})
        seq = mutate_substitutions(seq, start0, subs)
        nm = sum(start0 <= p < start0 + READ_LEN for p in subs)
        return seq, [(0, READ_LEN)], nm
    return model


def indel_read_model(
    genome: Genome, chrom: str, del_start0: int, del_end0: int,
    deleted_hps: set[int], subs_by_hp: dict[int, dict[int, str]],
):
    def model(hp: int, start0: int):
        subs = subs_by_hp.get(hp, {})
        if hp not in deleted_hps:
            seq = mutate_substitutions(genome.fetch(chrom, start0, start0 + READ_LEN), start0, subs)
            nm = sum(start0 <= p < start0 + READ_LEN for p in subs)
            return seq, [(0, READ_LEN)], nm
        if del_start0 <= start0 < del_end0:
            return None
        left = del_start0 - start0
        if 1 <= left < READ_LEN:
            right = READ_LEN - left
            seq = genome.fetch(chrom, start0, del_start0) + genome.fetch(
                chrom, del_end0, del_end0 + right
            )
            chars = list(seq)
            sub_nm = 0
            for pos0, alt in subs.items():
                if start0 <= pos0 < del_start0:
                    idx = pos0 - start0
                elif del_end0 <= pos0 < del_end0 + right:
                    idx = left + pos0 - del_end0
                else:
                    continue
                chars[idx] = alt
                sub_nm += 1
            cigar = [(0, left), (2, del_end0 - del_start0), (0, right)]
            return "".join(chars), cigar, (del_end0 - del_start0) + sub_nm
        seq = mutate_substitutions(genome.fetch(chrom, start0, start0 + READ_LEN), start0, subs)
        nm = sum(start0 <= p < start0 + READ_LEN for p in subs)
        return seq, [(0, READ_LEN)], nm
    return model


def add_wgs_gene_pairs(
    out: pysam.AlignmentFile, genome: Genome, sample: str, region: Region,
    model, depth: int = 30, hps: tuple[int, ...] = (1, 2),
    mapq_model=None, name_prefix: str = "WGS",
) -> int:
    """Tile independent PE150 fragments with jitter to approximate uniform WGS.

    A start every 300/depth bases yields the requested physical read depth.
    Insert sizes follow a truncated N(350,45), and every fragment has a unique
    start/insert combination so IGV does not show artificial read islands.
    """
    step = max(1, round((2 * READ_LEN) / depth))
    written = 0
    pair_index = 0
    for nominal in range(region.start0, region.end0 - 250, step):
        start0 = max(region.start0, nominal + random.randint(-4, 4))
        fragment = max(250, min(500, int(round(random.gauss(350, 45)))))
        if start0 + fragment > region.end0:
            continue
        mate0 = start0 + fragment - READ_LEN
        hp = hps[pair_index % len(hps)]
        r1_model = model(hp, start0)
        r2_model = model(hp, mate0)
        pair_index += 1
        if r1_model is None or r2_model is None:
            continue
        mapq1 = mapq_model(start0, pair_index, hp) if mapq_model else 60
        mapq2 = mapq_model(mate0, pair_index, hp) if mapq_model else 60
        qname = f"{sample}_{name_prefix}_HP{hp}_{pair_index:07d}"
        seq1, cigar1, nm1 = r1_model
        seq2, cigar2, nm2 = r2_model
        out.write(make_read(genome, region.chrom, start0, qname, sample, False, False,
                            mate0, fragment, mapq1, aligned_seq=seq1, cigar=cigar1,
                            hp=hp, nm=nm1))
        out.write(make_read(genome, region.chrom, mate0, qname, sample, True, True,
                            start0, -fragment, mapq2, aligned_seq=seq2, cigar=cigar2,
                            hp=hp, nm=nm2))
        written += 2
    return written


def add_forced_phase_pairs(
    out: pysam.AlignmentFile, genome: Genome, sample: str, chrom: str,
    left0: int, right0: int, subs_by_hp: dict[int, dict[int, str]], pairs: int = 12,
) -> None:
    model = basic_read_model(genome, chrom, subs_by_hp)
    fragment = (right0 - left0) + 145
    for i in range(pairs):
        hp = 1 if i < pairs // 2 else 2
        start0 = left0 - 70 + (i % 5) - 2
        mate0 = start0 + fragment - READ_LEN
        qname = f"{sample}_PHASE_HP{hp}_{i:04d}"
        seq1, cigar1, nm1 = model(hp, start0)
        seq2, cigar2, nm2 = model(hp, mate0)
        out.write(make_read(genome, chrom, start0, qname, sample, False, False,
                            mate0, fragment, 60, aligned_seq=seq1, cigar=cigar1,
                            hp=hp, nm=nm1, extra_tags=[("PV", "forced_pair_view")]))
        out.write(make_read(genome, chrom, mate0, qname, sample, True, True,
                            start0, -fragment, 60, aligned_seq=seq2, cigar=cigar2,
                            hp=hp, nm=nm2, extra_tags=[("PV", "forced_pair_view")]))


def add_softclip_pairs(out: pysam.AlignmentFile, genome: Genome, sample: str,
                       chrom: str, d0: int, d1: int) -> None:
    for i in range(8):
        left = 55 + (i % 5) * 10
        start0 = d0 - left
        mate0 = d1 + 170 + i * 3
        fragment = mate0 + READ_LEN - start0
        seq1 = genome.fetch(chrom, start0, d0) + genome.fetch(chrom, d1, d1 + READ_LEN - left)
        seq2 = genome.fetch(chrom, mate0, mate0 + READ_LEN)
        qname = f"{sample}_SOFTCLIP_PAIR_{i:03d}"
        out.write(make_read(genome, chrom, start0, qname, sample, False, False,
                            mate0, fragment, 20, aligned_seq=seq1,
                            cigar=[(0, left), (4, READ_LEN - left)], hp=1,
                            proper=False, nm=0,
                            extra_tags=[("SC", "intentional_breakpoint_softclip")]))
        out.write(make_read(genome, chrom, mate0, qname, sample, True, True,
                            start0, -fragment, 50, aligned_seq=seq2,
                            cigar=[(0, READ_LEN)], hp=1, proper=False, nm=0))


def add_cyp21_ambiguous_pairs(out: pysam.AlignmentFile, genome: Genome, sample: str,
                              pos0: int) -> None:
    pseudo_region, pseudo_seq = genome.loaded["cyp21a1p"]
    for i in range(10):
        hp = 1 if i < 5 else 2
        start0 = pos0 - 70 + (i % 5) - 2
        fragment = max(300, min(430, int(round(random.gauss(350, 35)))))
        mate0 = start0 + fragment - READ_LEN
        subs = {pos0: "T"} if hp == 1 else {}
        qname = f"{sample}_MULTIMAP_HP{hp}_{i:03d}"
        for reverse, read2, read_start, mate_start, tlen in [
            (False, False, start0, mate0, fragment),
            (True, True, mate0, start0, -fragment),
        ]:
            seq = mutate_substitutions(genome.fetch("chr6", read_start, read_start + READ_LEN),
                                       read_start, subs)
            primary = make_read(genome, "chr6", read_start, qname, sample, reverse,
                                read2, mate_start, tlen, 5, aligned_seq=seq, hp=hp,
                                nm=int(read_start <= pos0 < read_start + READ_LEN and hp == 1),
                                variant_bq={pos0: 12}, extra_tags=[("NH", 2)])
            out.write(primary)
            idx, mismatches = best_hamming(seq, pseudo_seq)
            secondary = make_read(genome, "chr6", pseudo_region.start0 + idx, qname,
                                  sample, reverse, read2, None, 0, 0,
                                  aligned_seq=seq, hp=hp, proper=False,
                                  nm=mismatches, extra_tags=[("NH", 2)])
            secondary.flag |= 0x100
            out.write(secondary)


def add_dup_junction_pairs(out: pysam.AlignmentFile, genome: Genome, sample: str,
                           chrom: str, d0: int, d1: int) -> None:
    for i in range(16):
        flank = 75
        qname = f"{sample}_TANDEM_JUNCTION_HP1_{i:03d}"
        junction = genome.fetch(chrom, d1 - flank, d1) + genome.fetch(chrom, d0, d0 + flank)
        mate0 = d0 + 220 + i
        mate_seq = genome.fetch(chrom, mate0, mate0 + READ_LEN)
        primary = make_read(genome, chrom, d1 - flank, qname, sample, False, False,
                            mate0, 0, 20, aligned_seq=junction,
                            cigar=[(0, flank), (4, flank)], hp=1, proper=False, nm=0,
                            extra_tags=[("SA", f"{chrom},{d0 + 1},+,75S75M,20,0;")])
        out.write(primary)
        supplementary = make_read(genome, chrom, d0, qname, sample, False, False,
                                  mate0, 0, 20, aligned_seq=junction,
                                  cigar=[(4, flank), (0, flank)], hp=1,
                                  proper=False, nm=0,
                                  extra_tags=[("SA", f"{chrom},{d1-flank+1},+,75M75S,20,0;")])
        supplementary.flag |= 0x800
        out.write(supplementary)
        out.write(make_read(genome, chrom, mate0, qname, sample, True, True,
                            d1 - flank, 0, 50, aligned_seq=mate_seq,
                            cigar=[(0, READ_LEN)], hp=1, proper=False, nm=0))


def circular_dup_read(genome: Genome, chrom: str, d0: int, d1: int,
                      rel_start: int, subs: dict[int, str]):
    """Render one read from a tandem copy treated as a circular interval.

    Reads crossing d1->d0 become complementary primary/supplementary segments.
    This keeps third-copy coverage stationary up to both breakpoints instead of
    creating start/end ramps or stacks at one exact coordinate.
    """
    length = d1 - d0
    rel_start %= length
    genomic_start = d0 + rel_start
    if rel_start + READ_LEN <= length:
        seq = genome.fetch(chrom, genomic_start, genomic_start + READ_LEN)
        seq = mutate_substitutions(seq, genomic_start, subs)
        nm = sum(genomic_start <= p < genomic_start + READ_LEN for p in subs)
        return seq, genomic_start, [(0, READ_LEN)], None, nm

    left = length - rel_start
    right = READ_LEN - left
    seq = genome.fetch(chrom, genomic_start, d1) + genome.fetch(chrom, d0, d0 + right)
    chars = list(seq)
    nm = 0
    for pos0, alt in subs.items():
        if genomic_start <= pos0 < d1:
            idx = pos0 - genomic_start
        elif d0 <= pos0 < d0 + right:
            idx = left + pos0 - d0
        else:
            continue
        chars[idx] = alt
        nm += 1
    primary_cigar = [(0, left), (4, right)]
    supplementary = (d0, [(4, left), (0, right)])
    return "".join(chars), genomic_start, primary_cigar, supplementary, nm


def add_circular_dup_copy_pairs(
    out: pysam.AlignmentFile, genome: Genome, sample: str, chrom: str,
    d0: int, d1: int, subs: dict[int, str], depth: int = 15,
) -> None:
    length = d1 - d0
    step = max(1, round((2 * READ_LEN) / depth))
    for i, nominal in enumerate(range(0, length, step), start=1):
        rel1 = (nominal + random.randint(-4, 4)) % length
        fragment = max(250, min(500, int(round(random.gauss(350, 45)))))
        rel2 = (rel1 + fragment - READ_LEN) % length
        r1 = circular_dup_read(genome, chrom, d0, d1, rel1, subs)
        r2 = circular_dup_read(genome, chrom, d0, d1, rel2, subs)
        seq1, start1, cigar1, supp1, nm1 = r1
        seq2, start2, cigar2, supp2, nm2 = r2
        wraps = supp1 is not None or supp2 is not None or rel2 < rel1
        proper = not wraps
        tlen = fragment if proper else 0
        qname = f"{sample}_DUP_COPY_HP1_{i:07d}"
        sa1 = [] if supp1 is None else [("SA", f"{chrom},{d0 + 1},+,{cigar1[0][1]}S{cigar1[1][1]}M,20,{nm1};")]
        sa2 = [] if supp2 is None else [("SA", f"{chrom},{d0 + 1},-,{cigar2[0][1]}S{cigar2[1][1]}M,20,{nm2};")]
        out.write(make_read(genome, chrom, start1, qname, sample, False, False,
                            start2, tlen, 60 if proper else 20,
                            aligned_seq=seq1, cigar=cigar1, hp=1,
                            proper=proper, nm=nm1, extra_tags=sa1))
        out.write(make_read(genome, chrom, start2, qname, sample, True, True,
                            start1, -tlen, 60 if proper else 20,
                            aligned_seq=seq2, cigar=cigar2, hp=1,
                            proper=proper, nm=nm2, extra_tags=sa2))
        for read2, seq, supplementary, nm, mate_start in [
            (False, seq1, supp1, nm1, start2),
            (True, seq2, supp2, nm2, start1),
        ]:
            if supplementary is None:
                continue
            supp_start, supp_cigar = supplementary
            supp = make_read(genome, chrom, supp_start, qname, sample, read2,
                             read2, mate_start, 0, 20, aligned_seq=seq,
                             cigar=supp_cigar, hp=1, proper=False, nm=nm)
            supp.flag |= 0x800
            out.write(supp)


def write_snv_good(path: Path, genome: Genome) -> None:
    sample, pos0 = "01_SNV_good_BRCA2", 32319079
    assert genome.fetch("chr13", pos0, pos0 + 1) == "T"
    with pysam.AlignmentFile(path, "wb", header=header(sample)) as out:
        model = basic_read_model(genome, "chr13", {1: {pos0: "G"}, 2: {}})
        add_wgs_gene_pairs(out, genome, sample, REGIONS["brca2"], model)


def write_snv_homology(path: Path, genome: Genome) -> None:
    sample, pos0 = "02_SNV_low_quality_CYP21A2", 32038609
    assert genome.fetch("chr6", pos0, pos0 + 1) == "A"
    def low_mapq(start0: int, i: int, hp: int) -> int:
        if start0 < 32041644 and start0 + READ_LEN > 32038414:
            return (0, 5, 10, 20)[i % 4]
        return 40 if i % 3 else 60
    with pysam.AlignmentFile(path, "wb", header=header(sample)) as out:
        model = basic_read_model(genome, "chr6", {1: {pos0: "T"}, 2: {}})
        add_wgs_gene_pairs(out, genome, sample, REGIONS["cyp21a2"], model,
                           mapq_model=low_mapq)
        add_cyp21_ambiguous_pairs(out, genome, sample, pos0)


def write_mnv(path: Path, genome: Genome) -> None:
    sample, chrom = "03_complex_MNV_hom_cis_chr8", "chr8"
    d0, d1, snv0 = 93797350, 93797351, 93797351
    assert genome.fetch(chrom, 93797349, 93797352) == "TGC"
    with pysam.AlignmentFile(path, "wb", header=header(sample)) as out:
        model = indel_read_model(genome, chrom, d0, d1, {1, 2},
                                 {1: {snv0: "T"}, 2: {snv0: "T"}})
        add_wgs_gene_pairs(out, genome, sample, REGIONS["mnv"], model)


def write_small_del(path: Path, genome: Genome) -> None:
    sample, chrom = "04_small_del_hom_CYP1B1", "chr2"
    d0, d1 = 38071277, 38071290
    assert genome.fetch(chrom, d0, d1) == "TCTGCCTGCACTC"
    with pysam.AlignmentFile(path, "wb", header=header(sample)) as out:
        model = indel_read_model(genome, chrom, d0, d1, {1, 2}, {1: {}, 2: {}})
        add_wgs_gene_pairs(out, genome, sample, REGIONS["cyp1b1"], model)
        add_softclip_pairs(out, genome, sample, chrom, d0, d1)


def write_btd(path: Path, genome: Genome, cis: bool) -> None:
    sample = "06_BTD_cis_pair_view" if cis else "05_BTD_trans_pair_view"
    left0 = 15644856 if cis else 15644916
    right0 = 15645185
    if cis:
        subs = {1: {left0: "A", right0: "C"}, 2: {}}
    else:
        subs = {1: {left0: "T"}, 2: {right0: "C"}}
    with pysam.AlignmentFile(path, "wb", header=header(sample)) as out:
        model = basic_read_model(genome, "chr3", subs)
        add_wgs_gene_pairs(out, genome, sample, REGIONS["btd"], model)
        add_forced_phase_pairs(out, genome, sample, "chr3", left0, right0, subs)


def write_del214(path: Path, genome: Genome) -> None:
    sample, chrom = "07_deletion_214bp_het", "chr11"
    d0, d1 = 66257086, 66257300
    flank_snvs = snv_map(genome, chrom, KLC2_FLANKING_HET_SNVS)
    inside_snvs = snv_map(genome, chrom, KLC2_HEMIZYGOUS_SNVS)
    with pysam.AlignmentFile(path, "wb", header=header(sample)) as out:
        # HP1 carries the deletion and flanking heterozygous SNVs. HP2 is intact
        # and carries the SNVs inside the deletion, so those are observed at 100%.
        model = indel_read_model(genome, chrom, d0, d1, {1},
                                 {1: flank_snvs, 2: inside_snvs})
        add_wgs_gene_pairs(out, genome, sample, REGIONS["del214"], model)


def write_dup(path: Path, genome: Genome) -> None:
    sample, chrom = "08_tandem_dup_251kb_het", "chr6"
    d0, d1 = 129049897, 129300892
    dup_snvs = snv_map(genome, chrom, DUPLICATION_BAF_SNVS)
    with pysam.AlignmentFile(path, "wb", header=header(sample)) as out:
        # HP1 is duplicated and carries the ALT alleles. Baseline HP1 + duplicated
        # HP1 contribute two ALT copies, while HP2 contributes one REF copy.
        model = basic_read_model(genome, chrom, {1: dup_snvs, 2: {}})
        add_wgs_gene_pairs(out, genome, sample, REGIONS["dup"], model, depth=30)
        add_circular_dup_copy_pairs(out, genome, sample, chrom, d0, d1,
                                    dup_snvs, depth=15)


def write_control_reference(path: Path, genome: Genome) -> None:
    """Write a variant-free control over every teaching gene window.

    The control deliberately uses one sample/read group and the same baseline
    PE150 WGS model as the scenario BAMs. Both haplotypes are identical to
    GRCh38, so the KLC2 interval remains diploid and the LAMA2 interval remains
    copy-neutral (CN=2) at approximately 30x.
    """
    sample = "00_control_reference"
    control_windows = (
        ("BRCA2", REGIONS["brca2"]),
        ("CYP21A2", REGIONS["cyp21a2"]),
        ("TMEM67", REGIONS["mnv"]),
        ("CYP1B1", REGIONS["cyp1b1"]),
        ("BTD", REGIONS["btd"]),
        ("KLC2", REGIONS["del214"]),
        ("LAMA2", REGIONS["dup"]),
    )
    with pysam.AlignmentFile(path, "wb", header=header(sample)) as out:
        for gene, region in control_windows:
            model = basic_read_model(genome, region.chrom, {1: {}, 2: {}})
            add_wgs_gene_pairs(out, genome, sample, region, model, depth=30,
                               name_prefix=f"CONTROL_{gene}")


def pileup_counts(bam: pysam.AlignmentFile, chrom: str, pos0: int,
                  ref: str, alt: str) -> tuple[int, int, int]:
    ref_count = alt_count = 0
    for col in bam.pileup(chrom, pos0, pos0 + 1, truncate=True,
                          min_base_quality=20, stepper="nofilter"):
        for pr in col.pileups:
            read = pr.alignment
            if (pr.is_del or pr.is_refskip or read.is_secondary
                    or read.is_supplementary or read.is_duplicate):
                continue
            base = read.query_sequence[pr.query_position]
            ref_count += base == ref
            alt_count += base == alt
    return ref_count, alt_count, ref_count + alt_count


def write_duplication_baf(outdir: Path, genome: Genome) -> None:
    bam_path = outdir / "08_tandem_dup_251kb_het.bam"
    snvs = snv_map(genome, "chr6", DUPLICATION_BAF_SNVS)
    with pysam.AlignmentFile(bam_path) as bam, \
            (outdir / "08_LAMA2_duplication_BAF.tsv").open("w") as tsv, \
            (outdir / "08_LAMA2_duplication_BAF.bedgraph").open("w") as bedgraph:
        tsv.write("CHROM\tPOS\tREF\tALT\tREF_COUNT\tALT_COUNT\tDP\tBAF\tEXPECTED_BAF\n")
        bedgraph.write("track type=bedGraph name=\"LAMA2 duplication BAF\" description=\"ALT allele fraction; expected 0.667\"\n")
        for pos0, alt in sorted(snvs.items()):
            ref = genome.fetch("chr6", pos0, pos0 + 1)
            ref_count, alt_count, dp = pileup_counts(bam, "chr6", pos0, ref, alt)
            baf = alt_count / dp if dp else float("nan")
            tsv.write(f"chr6\t{pos0 + 1}\t{ref}\t{alt}\t{ref_count}\t{alt_count}\t{dp}\t{baf:.4f}\t0.6667\n")
            bedgraph.write(f"chr6\t{pos0}\t{pos0 + 1}\t{baf:.4f}\n")


def write_vcf_and_manifest(outdir: Path, genome: Genome) -> None:
    vcf = outdir / "expected_variants.grch38.vcf"
    records = []
    def snv(chrom: str, pos1: int, ref: str, alt: str, sample: str, gt: str, info: str):
        assert genome.fetch(chrom, pos1 - 1, pos1) == ref
        records.append((chrom, pos1, ".", ref, alt, ".", "PASS", info, "GT:PS", f"{gt}:1", sample))
    snv("chr13", 32319080, "T", "G", "01_SNV_good_BRCA2", "0|1", "CLINVAR=266991")
    snv("chr6", 32038610, "A", "T", "02_SNV_low_quality_CYP21A2", "0|1", "CLINVAR=12183;LOW_QUALITY=HOMOLOGY")
    # Component records for the homozygous cis complex allele.
    records.append(("chr8", 93797350, ".", "TG", "T", ".", "PASS", "SOURCE=gnomAD_r4;COMPONENT=1", "GT:PS", "1|1:1", "03_complex_MNV_hom_cis_chr8"))
    snv("chr8", 93797352, "C", "T", "03_complex_MNV_hom_cis_chr8", "1|1", "SOURCE=gnomAD_r4;COMPONENT=2")
    anchor = genome.fetch("chr2", 38071276, 38071290)
    records.append(("chr2", 38071277, ".", anchor, anchor[0], ".", "PASS", "CLINVAR=282564", "GT:PS", "1|1:1", "04_small_del_hom_CYP1B1"))
    snv("chr3", 15644917, "C", "T", "05_BTD_trans_pair_view", "1|0", "CLINVAR=2230099")
    snv("chr3", 15645186, "G", "C", "05_BTD_trans_pair_view", "0|1", "CLINVAR=1900")
    snv("chr3", 15644857, "T", "A", "06_BTD_cis_pair_view", "1|0", "CLINVAR=2203317")
    snv("chr3", 15645186, "G", "C", "06_BTD_cis_pair_view", "1|0", "CLINVAR=1900")
    records.append(("chr11", 66257086, ".", genome.fetch("chr11", 66257085, 66257086), "<DEL>", ".", "PASS", "END=66257300;SVTYPE=DEL;SVLEN=-214;CLINVAR=1684657", "GT:PS", "0|1:1", "07_deletion_214bp_het"))
    for pos0, alt in snv_map(genome, "chr11", KLC2_FLANKING_HET_SNVS).items():
        ref = genome.fetch("chr11", pos0, pos0 + 1)
        snv("chr11", pos0 + 1, ref, alt, "07_deletion_214bp_het", "1|0",
            "SYNTHETIC_TYPE=FLANKING_HET;EXPECTED_BAF=0.5")
    for pos0, alt in snv_map(genome, "chr11", KLC2_HEMIZYGOUS_SNVS).items():
        ref = genome.fetch("chr11", pos0, pos0 + 1)
        snv("chr11", pos0 + 1, ref, alt, "07_deletion_214bp_het", ".|1",
            "SYNTHETIC_TYPE=HEMIZYGOUS_INSIDE_DEL;EXPECTED_BAF=1.0")
    records.append(("chr6", 129049898, ".", genome.fetch("chr6", 129049897, 129049898), "<DUP:TANDEM>", ".", "PASS", "END=129300892;SVTYPE=DUP;SVLEN=250995;CLINVAR=543888;BREAKPOINTS=ASSUMED_INNER_BOUNDS", "GT:PS", "0|1:1", "08_tandem_dup_251kb_het"))
    for pos0, alt in snv_map(genome, "chr6", DUPLICATION_BAF_SNVS).items():
        ref = genome.fetch("chr6", pos0, pos0 + 1)
        snv("chr6", pos0 + 1, ref, alt, "08_tandem_dup_251kb_het", "1|0",
            "SYNTHETIC_TYPE=DUPLICATED_HAPLOTYPE;CN=3;EXPECTED_BAF=0.6667")
    with vcf.open("w") as h:
        h.write("##fileformat=VCFv4.3\n##reference=GRCh38\n")
        for chrom, (_, length) in CONTIGS.items():
            h.write(f"##contig=<ID={chrom},length={length}>\n")
        h.write("##INFO=<ID=END,Number=1,Type=Integer,Description=End coordinate>\n")
        h.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=SV type>\n")
        h.write("##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=SV length>\n")
        h.write("##INFO=<ID=SCENARIO,Number=1,Type=String,Description=Corresponding BAM scenario>\n")
        h.write("##INFO=<ID=CLINVAR,Number=1,Type=String,Description=ClinVar Variation ID>\n")
        h.write("##INFO=<ID=LOW_QUALITY,Number=1,Type=String,Description=Reason for low-quality evidence>\n")
        h.write("##INFO=<ID=SOURCE,Number=1,Type=String,Description=Variant source>\n")
        h.write("##INFO=<ID=COMPONENT,Number=1,Type=Integer,Description=Complex-allele component number>\n")
        h.write("##INFO=<ID=BREAKPOINTS,Number=1,Type=String,Description=Breakpoint modeling note>\n")
        h.write("##INFO=<ID=SYNTHETIC_TYPE,Number=1,Type=String,Description=Teaching SNV category>\n")
        h.write("##INFO=<ID=EXPECTED_BAF,Number=1,Type=Float,Description=Expected ALT allele fraction>\n")
        h.write("##INFO=<ID=CN,Number=1,Type=Integer,Description=Expected total copy number>\n")
        h.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>\n")
        h.write("##FORMAT=<ID=PS,Number=1,Type=Integer,Description=Phase set>\n")
        h.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYNTHETIC\n")
        for rec in records:
            row = list(rec[:10])
            row[7] += f";SCENARIO={rec[10]}"
            h.write("\t".join(map(str, row)) + "\n")
    manifest_rows = [
        ["00_control_reference.grch38.bam", "controle negativo GRCh38", "todas as janelas do BED", "0|0 em todos os sítios", "GRCh38", "um RG/SM; PE150; ~30x; CN=2; sem variantes, indels, soft-clips ou reads suplementares"],
        ["00_all_scenarios.grch38.bam", "todos os oito cenários", "múltiplos loci", "preserva genótipos/fase individuais", "arquivos 01-08", "BAM combinado com oito RG/SM; recomendado para IGV"],
        ["01_SNV_good_BRCA2.bam", "SNV boa qualidade", "chr13:32319080 T>G", "0|1", "ClinVar 266991", "BRCA2 completo; PE150; ~30x; MAPQ60; BQ variável"],
        ["02_SNV_low_quality_CYP21A2.bam", "SNV em região homóloga", "chr6:32038610 A>T", "0|1", "ClinVar 12183", "CYP21A2 completo; MAPQ0-20 no gene; alinhamentos secundários CYP21A1P"],
        ["03_complex_MNV_hom_cis_chr8.bam", "alelo complexo/MNV", "chr8:93797350 TG>T + chr8:93797352 C>T", "1|1 / cis", "gnomAD r4", "TMEM67 completo; ambos os componentes nos dois haplótipos"],
        ["04_small_del_hom_CYP1B1.bam", "deleção pequena", "chr2:38071278-38071290 del13", "1|1", "ClinVar 282564", "CYP1B1 completo; PE150; CIGAR 13D e soft-clips intencionais"],
        ["05_BTD_trans_pair_view.bam", "variantes em trans", "chr3:15644917 C>T / chr3:15645186 G>C", "1|0 / 0|1", "ClinVar 2230099 / 1900", "BTD completo; insert N(350,45); pares de fase adicionais e tags HP"],
        ["06_BTD_cis_pair_view.bam", "variantes em cis", "chr3:15644857 T>A / chr3:15645186 G>C", "1|0 / 1|0", "ClinVar 2203317 / 1900", "BTD completo; insert N(350,45); pares de fase adicionais e tags HP"],
        ["07_deletion_214bp_het.bam", "deleção heterozigótica", "chr11:66257087-66257300 del214", "0|1", "ClinVar 1684657", "KLC2 completo; SNVs flanqueadores ~50% e SNVs internos hemizigóticos ~100%"],
        ["08_tandem_dup_251kb_het.bam", "duplicação tandem", "chr6:129049898-129300892 dup", "0|1 assumido", "ClinVar 543888", "LAMA2 completo; terceira cópia circularizada; SNVs no haplótipo duplicado ~66%; TSV/bedGraph BAF"],
    ]
    with (outdir / "manifest.tsv").open("w", newline="") as h:
        w = csv.writer(h, delimiter="\t")
        w.writerow(["bam", "scenario", "grch38_variant", "genotype_phase", "source", "synthetic_evidence"])
        w.writerows(manifest_rows)
    bed_rows = [
        ("chr13", 32314076, 32401268, "BRCA2_window"),
        ("chr6", 32037414, 32042644, "CYP21A2_window"),
        ("chr8", 93753843, 93819121, "TMEM67_window"),
        ("chr2", 38066508, 38077151, "CYP1B1_window"),
        ("chr3", 15600360, 15723516, "BTD_window"),
        ("chr11", 66256086, 66268860, "KLC2_and_deletion_window"),
        ("chr6", 128882137, 129517566, "LAMA2_window"),
    ]
    with (outdir / "gene_windows.grch38.bed").open("w") as h:
        for row in bed_rows:
            h.write("\t".join(map(str, row)) + "\n")


def finalize_bam(unsorted: Path, final: Path) -> None:
    pysam.sort("-o", str(final), str(unsorted))
    unsorted.unlink()
    pysam.index(str(final))
    # Fast structural validation.
    pysam.quickcheck(str(final))


def merge_scenarios(individual_bams: list[Path], merged: Path, tmpdir: Path) -> None:
    merged_unsorted = tmpdir / "all_scenarios.merged.bam"
    pysam.merge("-f", str(merged_unsorted), *map(str, individual_bams))
    finalize_bam(merged_unsorted, merged)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="output", help="Output directory")
    ap.add_argument("--reference", help="Indexed hg38 FASTA (optional; otherwise fetches NCBI slices)")
    args = ap.parse_args()
    random.seed(SEED)
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    genome = Genome(outdir / "reference_cache", args.reference)
    genome.preload()
    jobs = [
        ("01_SNV_good_BRCA2.bam", write_snv_good),
        ("02_SNV_low_quality_CYP21A2.bam", write_snv_homology),
        ("03_complex_MNV_hom_cis_chr8.bam", write_mnv),
        ("04_small_del_hom_CYP1B1.bam", write_small_del),
        ("05_BTD_trans_pair_view.bam", lambda p, g: write_btd(p, g, False)),
        ("06_BTD_cis_pair_view.bam", lambda p, g: write_btd(p, g, True)),
        ("07_deletion_214bp_het.bam", write_del214),
        ("08_tandem_dup_251kb_het.bam", write_dup),
    ]
    with tempfile.TemporaryDirectory(prefix="synthetic_variant_bams_") as tmp:
        individual_bams = []
        for filename, writer in jobs:
            final = outdir / filename
            unsorted = Path(tmp) / f"{filename}.unsorted"
            writer(unsorted, genome)
            finalize_bam(unsorted, final)
            individual_bams.append(final)
            print(final.name)
        merged = outdir / "00_all_scenarios.grch38.bam"
        merge_scenarios(individual_bams, merged, Path(tmp))
        print(merged.name)
        # Generate the negative control after the scenarios so adding it does
        # not alter their seeded random geometry or previously validated BAFs.
        control = outdir / "00_control_reference.grch38.bam"
        control_unsorted = Path(tmp) / "00_control_reference.grch38.bam.unsorted"
        write_control_reference(control_unsorted, genome)
        finalize_bam(control_unsorted, control)
        print(control.name)
    write_duplication_baf(outdir, genome)
    write_vcf_and_manifest(outdir, genome)
    print("Done")


if __name__ == "__main__":
    main()
