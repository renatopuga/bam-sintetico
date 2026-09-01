#!/usr/bin/env python3
"""Validate BAM integrity, WGS geometry, coverage and expected variant evidence."""

from collections import Counter, defaultdict
from pathlib import Path
import argparse
import csv
import statistics
import pysam


REGIONS = {
    "01_SNV_good_BRCA2.bam": ("chr13", 32314076, 32401268, []),
    "02_SNV_low_quality_CYP21A2.bam": ("chr6", 32037414, 32042644, []),
    "03_complex_MNV_hom_cis_chr8.bam": ("chr8", 93753843, 93819121, [(93797350, 93797351)]),
    "04_small_del_hom_CYP1B1.bam": ("chr2", 38066508, 38077151, [(38071277, 38071290)]),
    "05_BTD_trans_pair_view.bam": ("chr3", 15600360, 15723516, []),
    "06_BTD_cis_pair_view.bam": ("chr3", 15600360, 15723516, []),
    "07_deletion_214bp_het.bam": ("chr11", 66256086, 66268860, []),
    "08_tandem_dup_251kb_het.bam": ("chr6", 128882137, 129517566, []),
}

CONTROL_WINDOWS = [
    ("chr13", 32314076, 32401268),
    ("chr6", 32037414, 32042644),
    ("chr8", 93753843, 93819121),
    ("chr2", 38066508, 38077151),
    ("chr3", 15600360, 15723516),
    ("chr11", 66256086, 66268860),
    ("chr6", 128882137, 129517566),
]

CONTROL_REFERENCE_SITES = [
    ("chr13", 32319080, "T"),
    ("chr6", 32038610, "A"),
    ("chr8", 93797352, "C"),
    ("chr2", 38071278, "T"),
    ("chr3", 15644917, "C"),
    ("chr3", 15645186, "G"),
    ("chr3", 15644857, "T"),
    ("chr11", 66257122, None),
    ("chr6", 129060123, None),
]


def allele_by_hp(path: Path, chrom: str, pos1: int):
    total, hp = Counter(), defaultdict(Counter)
    with pysam.AlignmentFile(path) as bam:
        for col in bam.pileup(chrom, pos1 - 1, pos1, truncate=True,
                              min_base_quality=0, stepper="nofilter"):
            for pr in col.pileups:
                read = pr.alignment
                if (pr.is_del or pr.is_refskip or read.is_secondary
                        or read.is_supplementary):
                    continue
                base = read.query_sequence[pr.query_position]
                total[base] += 1
                hp[read.get_tag("HP") if read.has_tag("HP") else 0][base] += 1
    return total, hp


def depths(bam: pysam.AlignmentFile, chrom: str, start0: int, end0: int):
    cov = bam.count_coverage(chrom, start0, end0, quality_threshold=0,
                             read_callback="all")
    return [sum(x) for x in zip(*cov)]


def mean_depth(bam: pysam.AlignmentFile, chrom: str, start0: int, end0: int) -> float:
    return statistics.mean(depths(bam, chrom, start0, end0))


def assert_haplotype(hp, hap: int, expected: str):
    assert hp[hap], f"No observations for HP{hap}"
    assert set(hp[hap]) == {expected}, (hap, hp[hap], expected)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("outdir", nargs="?", default="output")
    out = Path(ap.parse_args().outdir)

    bam_paths = sorted(out.glob("*.bam"))
    assert len(bam_paths) == 10, f"Expected 10 BAMs, found {len(bam_paths)}"
    for path in bam_paths:
        pysam.quickcheck(str(path))
        assert Path(str(path) + ".bai").exists(), f"Missing index: {path}"

    for filename, (chrom, start0, end0, allowed_zero_intervals) in REGIONS.items():
        path = out / filename
        inserts, primary = [], []
        with pysam.AlignmentFile(path) as bam:
            for read in bam.fetch(until_eof=True):
                if read.is_secondary or read.is_supplementary:
                    continue
                primary.append(read)
                assert read.is_paired and read.query_length == 150
                if read.is_read1 and read.template_length:
                    inserts.append(abs(read.template_length))
        proper_fraction = sum(r.is_proper_pair for r in primary) / len(primary)
        assert proper_fraction > 0.98, (filename, proper_fraction)
        assert 335 < statistics.mean(inserts) < 365, (filename, statistics.mean(inserts))
        assert 30 < statistics.pstdev(inserts) < 60, (filename, statistics.pstdev(inserts))

        inner_start, inner_end = start0 + 500, end0 - 500
        with pysam.AlignmentFile(path) as bam:
            cov = depths(bam, chrom, inner_start, inner_end)
        for i, value in enumerate(cov, start=inner_start):
            allowed_zero = any(a <= i < b for a, b in allowed_zero_intervals)
            assert value > 0 or allowed_zero, (filename, chrom, i + 1)

    individual = sorted(p for p in bam_paths if not p.name.startswith("00_"))
    merged_path = out / "00_all_scenarios.grch38.bam"
    with pysam.AlignmentFile(merged_path) as merged:
        merged_count = sum(1 for _ in merged.fetch(until_eof=True))
        rg = merged.header.to_dict().get("RG", [])
        assert len(rg) == 8 and len({x["SM"] for x in rg}) == 8
    individual_count = 0
    for path in individual:
        with pysam.AlignmentFile(path) as bam:
            individual_count += sum(1 for _ in bam.fetch(until_eof=True))
    assert merged_count == individual_count

    # Negative control: same WGS geometry and full set of windows, but one
    # copy-neutral reference sample with no alternate bases or SV signatures.
    control_path = out / "00_control_reference.grch38.bam"
    with pysam.AlignmentFile(control_path) as control:
        rg = control.header.to_dict().get("RG", [])
        assert len(rg) == 1 and rg[0]["SM"] == "00_control_reference"
        primary = [r for r in control.fetch(until_eof=True)
                   if not r.is_secondary and not r.is_supplementary]
        assert primary and all(r.is_paired and r.is_proper_pair for r in primary)
        assert all(r.query_length == 150 and r.cigarstring == "150M" for r in primary)
    for chrom, start0, end0 in CONTROL_WINDOWS:
        with pysam.AlignmentFile(control_path) as control:
            cov = depths(control, chrom, start0 + 500, end0 - 500)
        assert min(cov) > 0, ("control uncovered base", chrom)
        assert 27 <= statistics.mean(cov) <= 33, (chrom, statistics.mean(cov))
    for chrom, pos1, expected_ref in CONTROL_REFERENCE_SITES:
        total, hp = allele_by_hp(control_path, chrom, pos1)
        assert total, ("control missing locus", chrom, pos1)
        if expected_ref is not None:
            assert set(total) == {expected_ref}, (chrom, pos1, total)
        assert set(hp) == {1, 2}, (chrom, pos1, hp)

    total, hp = allele_by_hp(out / "01_SNV_good_BRCA2.bam", "chr13", 32319080)
    assert_haplotype(hp, 1, "G"); assert_haplotype(hp, 2, "T")
    assert 0.35 < total["G"] / sum(total.values()) < 0.65

    total, hp = allele_by_hp(out / "02_SNV_low_quality_CYP21A2.bam", "chr6", 32038610)
    assert_haplotype(hp, 1, "T"); assert_haplotype(hp, 2, "A")
    with pysam.AlignmentFile(out / "02_SNV_low_quality_CYP21A2.bam") as bam:
        local = list(bam.fetch("chr6", 32038550, 32038670))
        assert max(r.mapping_quality for r in local if not r.is_secondary) <= 20
    with pysam.AlignmentFile(out / "02_SNV_low_quality_CYP21A2.bam") as bam:
        assert any(r.is_secondary for r in bam.fetch(until_eof=True))

    total, hp = allele_by_hp(out / "03_complex_MNV_hom_cis_chr8.bam", "chr8", 93797352)
    assert set(total) == {"T"} and set(hp) == {1, 2}

    _, trans_left = allele_by_hp(out / "05_BTD_trans_pair_view.bam", "chr3", 15644917)
    _, trans_right = allele_by_hp(out / "05_BTD_trans_pair_view.bam", "chr3", 15645186)
    assert_haplotype(trans_left, 1, "T"); assert_haplotype(trans_right, 1, "G")
    assert_haplotype(trans_left, 2, "C"); assert_haplotype(trans_right, 2, "C")

    _, cis_left = allele_by_hp(out / "06_BTD_cis_pair_view.bam", "chr3", 15644857)
    _, cis_right = allele_by_hp(out / "06_BTD_cis_pair_view.bam", "chr3", 15645186)
    assert_haplotype(cis_left, 1, "A"); assert_haplotype(cis_right, 1, "C")
    assert_haplotype(cis_left, 2, "T"); assert_haplotype(cis_right, 2, "G")

    with pysam.AlignmentFile(out / "03_complex_MNV_hom_cis_chr8.bam") as bam:
        crossing = [r for r in bam.fetch("chr8", 93797340, 93797360) if "1D" in r.cigarstring]
        assert crossing and all(r.get_tag("NM") == 2 for r in crossing)
    with pysam.AlignmentFile(out / "04_small_del_hom_CYP1B1.bam") as bam:
        cigars = [r.cigarstring for r in bam.fetch("chr2", 38071200, 38071350)]
        assert any("13D" in c for c in cigars) and any("S" in c for c in cigars)
    with pysam.AlignmentFile(out / "07_deletion_214bp_het.bam") as bam:
        cigars = [r.cigarstring for r in bam.fetch("chr11", 66257000, 66257400)]
        assert any("214D" in c for c in cigars) and any(c == "150M" for c in cigars)
    with pysam.AlignmentFile(out / "08_tandem_dup_251kb_het.bam") as bam:
        outside = mean_depth(bam, "chr6", 128900000, 128902000)
        inside = mean_depth(bam, "chr6", 129060000, 129062000)
        assert 1.4 < inside / outside < 1.7, (outside, inside)
        assert any(r.is_supplementary for r in bam.fetch("chr6", 129049800, 129050100))

    # KLC2: flanking SNVs remain heterozygous; SNVs hidden by the deleted HP1
    # allele are observed only on intact HP2 and therefore appear hemizygous.
    klc_flank_vafs, klc_inside_vafs = [], []
    with pysam.VariantFile(out / "expected_variants.grch38.vcf") as vcf:
        for rec in vcf:
            kind = rec.info.get("SYNTHETIC_TYPE")
            if rec.chrom != "chr11" or not kind:
                continue
            total, _ = allele_by_hp(out / "07_deletion_214bp_het.bam", rec.chrom, rec.pos)
            vaf = total[rec.alts[0]] / sum(total.values())
            if kind == "FLANKING_HET":
                klc_flank_vafs.append(vaf)
            elif kind == "HEMIZYGOUS_INSIDE_DEL":
                klc_inside_vafs.append(vaf)
    assert len(klc_flank_vafs) == 6 and all(0.35 <= x <= 0.65 for x in klc_flank_vafs)
    assert len(klc_inside_vafs) == 3 and all(x >= 0.95 for x in klc_inside_vafs)

    # The circular third-copy model must not create breakpoint coverage spikes.
    d0, d1 = 129049897, 129300892
    with pysam.AlignmentFile(out / "08_tandem_dup_251kb_het.bam") as bam:
        plateau = mean_depth(bam, "chr6", d0 + 5000, d0 + 7000)
        boundary_means = [
            mean_depth(bam, "chr6", d0, d0 + 75),
            mean_depth(bam, "chr6", d0 + 75, d0 + 150),
            mean_depth(bam, "chr6", d1 - 150, d1 - 75),
            mean_depth(bam, "chr6", d1 - 75, d1),
        ]
    assert all(abs(x - plateau) < 6 for x in boundary_means), (plateau, boundary_means)

    baf_path = out / "08_LAMA2_duplication_BAF.tsv"
    with baf_path.open() as handle:
        baf_rows = list(csv.DictReader(handle, delimiter="\t"))
    bafs = [float(r["BAF"]) for r in baf_rows]
    assert len(bafs) == 16
    assert 0.63 <= statistics.mean(bafs) <= 0.70, statistics.mean(bafs)
    assert all(0.55 <= x <= 0.78 for x in bafs)
    assert (out / "08_LAMA2_duplication_BAF.bedgraph").exists()

    with pysam.VariantFile(out / "expected_variants.grch38.vcf") as vcf:
        assert sum(1 for _ in vcf) == 36
    assert (out / "gene_windows.grch38.bed").exists()
    print("OK: variant-free control plus PE150 FR scenarios, continuous gene-window coverage, variants, phase and SV evidence validated")


if __name__ == "__main__":
    main()
