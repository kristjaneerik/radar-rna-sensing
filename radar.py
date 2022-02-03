import re
import tqdm
import fire
from Bio.Seq import Seq
from Bio import SeqIO
from collections import Counter, defaultdict
import numpy as np
import pandas as pd
import gzip
import hashlib
import tempfile
import subprocess
import multiprocessing
import functools
pd.set_option("display.max_rows", 1000)
pd.set_option("display.max_columns", 100)
pd.set_option("display.max_colwidth", 1000)
pd.set_option("display.expand_frame_repr", False)


GENOME = "/home/eerik/Downloads/genomics/hg38.fa.gz"
FA_GENOME ="/home/eerik/Downloads/genomics/for_blast/hg38.fa"  # masked
MOUSE_FA_GENOME = "/home/eerik/Downloads/genomics/for_blast/GRCm39.fa"

MART_EXPORT = "/home/eerik/Downloads/genomics/UTR/mart_export.txt.gz"
MART_EXPORT_CDS = "/home/eerik/Downloads/genomics/human_cDNA_mart_export.txt.gz"

NONPOLAR = {"F", "L", "I", "M", "V", "P", "A", "W", "G"}
POLAR = {"S", "T", "Y", "Q", "N", "C"}
BASIC = {"H", "K", "R"}
ACIDIC = {"D", "E"}
AA_TO_TYPE = {"*": "*"}
for a in NONPOLAR: AA_TO_TYPE[a] = "N"
for a in POLAR: AA_TO_TYPE[a] = "P"
for a in BASIC: AA_TO_TYPE[a] = "B"
for a in ACIDIC: AA_TO_TYPE[a] = "A"


def count_residues(tln):
    return (
        len([x for x in tln if x in NONPOLAR]),
        len([x for x in tln if x in POLAR]),
        len([x for x in tln if x in BASIC]),
        len([x for x in tln if x in ACIDIC]),
    )


def eval_for_homopolymers(seq):
    nt = seq[0]
    nt_count = 0
    longest_homopolymer = 0
    total_homopolymer = 0
    for this_nt in seq:
        if this_nt == nt:
            nt_count += 1
        else:
            longest_homopolymer = max(longest_homopolymer, nt_count)
            if nt_count >= 4:
                total_homopolymer += nt_count
            nt = this_nt
            nt_count = 1
    longest_homopolymer = max(longest_homopolymer, nt_count)
    return longest_homopolymer, total_homopolymer


def find_arRNAs(
    seq, variant="CCA", version="stop", left=42, right=45, verbose=True,
    max_longest_homopolymer=5, max_total_homopolymer=10, min_gc=0.35, max_gc=0.8,
    max_candidates=None, max_candidates_rand=None, allow_extra_stops=0,
):
    """
    v1: Just search for CCA sequences. Exclude any candidates that have an in-frame stop in the
    -42 ... +45 bases around the CCA (total sensor length 42 (left) + 3 + 45 (right) = 90).
    This also works for TCA sequences: one would use TAA as the corresponding stop codon.
    left and right could be reduced to find more potential candidates,
    since e.g., a 72-base sensor (left = 33, right = 36) would also work,
    just less effectively, but would have a smaller chance of in-frame stop codons.
    CTA is a perfect match for UAG, and still shows editing activity, so is another variant.

    Note: for UGA, you would use either CCA (with C:A mismatch), or TCA (no mismatch).
    """
    assert variant in {
        "CCA", "TCA", "CTA", "CAA", "GCA",  # stop codon based
        "NNN",  # for split sensors -- just look for stretch without stop
    }
    assert left % 3 == 0
    assert right % 3 == 0
    seq = seq.upper()
    candidates = []
    starts = []
    ends = []
    if variant != "NNN":
        for m in re.finditer(variant, seq):
            # a candidate is a subsequence of seq
            starts.append(m.start() - left)
            ends.append(m.end() + right)
    else:
        starts = list(range(len(seq)))
        ends = [s + left + 3 + right for s in starts]
        if max_candidates_rand and len(starts) > max_candidates_rand:
            prng = np.random.RandomState(42)
            ixs = prng.randint(0, len(starts), max_candidates_rand)
            starts = [starts[i] for i in ixs]
            ends = [ends[i] for i in ixs]
    for candidate_start, candidate_end in zip(starts, ends):
        # a candidate is a subsequence of seq
        if candidate_start < 0 or candidate_end > len(seq):
            # candidate out of bounds
            continue
        candidate = seq[candidate_start:candidate_end]
        rc = Seq(candidate).reverse_complement()
        tln = rc.translate()  # the amino acids that are part of the sensor
        stops = [i for i, a in enumerate(tln) if a == "*"]
        if stops:
            if verbose >= 3:
                print(f"{len(stops)}x* at positions {stops} in rc: {tln}")
            # candidate has in-frame stop
            if variant == "CTA" or variant == "TCA":
                # for the CTA and TCA variants we actually have the stop codon there already,
                # so we just can't have more than *one* stop
                if len(stops) > 1 + allow_extra_stops:
                    continue
            elif len(stops) > allow_extra_stops:
                continue
        sensor = list(rc)  # this is the nucleotide sequence of the sensor as a list
        if variant != "NNN":
            sensor[right + 1] = "a"  # this is the A that gets edited
            if version == "stop":
                sensor[right + 0] = "T"
                sensor[right + 2] = "G"
            # in the NNN case we are not looking to edit anything
        sensor = "".join(sensor)  # now we are back to a string
        nt_counts = Counter(candidate)
        gc = round((nt_counts["G"] + nt_counts["C"]) / len(candidate), 2)
        if gc < min_gc:
            if verbose >= 3:
                print(f"low GC: {gc}")
            continue
        elif gc > max_gc:
            if verbose >= 3:
                print(f"high GC: {gc}")
            continue
        longest_homopolymer, total_homopolymer = eval_for_homopolymers(sensor)
        if longest_homopolymer > max_longest_homopolymer:
            if verbose >= 3:
                print(f"long homopolymer: {longest_homopolymer}")
            continue
        if total_homopolymer > max_total_homopolymer:
            if verbose >= 3:
                print(f"too much total homopolymer: {total_homopolymer}")
            continue
        sensor_tln = str(Seq(sensor).translate())
        nonp, p, b, a = count_residues(sensor_tln)
        candidates.append({
            "start": candidate_start,
            "end": candidate_end,
            "candidate": candidate,
            "rc": str(rc),
            "sensor": sensor,
            "sensor_tln": sensor_tln,
            "translation": str(tln),
            "GC": gc,
            "longest_homopolymer": longest_homopolymer,
            "total_homopolymer": total_homopolymer,
            "nonpolar": nonp, "polar": p, "basic": b, "acidic": a,
            "variant": variant, "left": left, "right": right,
        })
        if max_candidates and len(candidates) >= max_candidates:
            break
    if max_candidates_rand and len(candidates) > max_candidates_rand:
        prng = np.random.RandomState(42)
        candidates = [candidates[i] for i in prng.randint(0, len(candidates), max_candidates_rand)]
    if verbose:
        print(f"{len(candidates)} candidate(s)")
        for candidate in candidates:
            print(f"{candidate['start']}-{candidate['end']} (0-index, right-exclusive)")
            print("\tmRNA\t" + candidate["candidate"].replace(variant, variant.lower()))
            print(f"\tsensor\t{candidate['sensor']}")
            print(f"\tcoding\t{candidate['sensor_tln']}")
            print(f"\ttype\t{''.join([AA_TO_TYPE[a] for a in candidate['sensor_tln']])}")
            nonp, p, b, a = count_residues(candidate["sensor_tln"])
            print(f"\t\t{nonp} Nonpolar, {p} Polar, {b} Basic, {a} Acidic")
            # NOTE: 1 nonpolar added from the * -> W edit
            # with the 2A tags, the nature of the AAs is not really a concern
            print(f"\tGC: {candidate['GC']:.1%}")
            print("")
    return candidates


def reverse_complement(s):
    return str(Seq(s).reverse_complement())


def design_oligos(candidate, left, right, overhang=4, prefix="o"):
    """Annotates the candidate dictionary with four oligos: 1 & 3 are annealed, 2 & 4 are annealed.
    The 2 & 4 pair contains the stop codon.
    To make the trigger, you need oligos 5 & 6 annealed: put them together with the 1 & 3 anneal.
    1 & 3 + 5 & 6 also work for making the positive control.
    """
    bad_overhangs = {"GATC", "CCGG"}  # used in cloning
    # the three below comes from e.g. "CCA" being in the middle of "left" and "right"
    o1 = "gatcc" + candidate["sensor"][:left - overhang + 3]
    o2 = candidate["sensor"][left - overhang + 3:] + "a"
    if o2[:overhang] in bad_overhangs or reverse_complement(o2[:overhang]) in bad_overhangs:
        # make the overhang longer if it happens to be the thing used in cloning; this is still
        # not great, it would probably be better to shift it over or something like that.
        return design_oligos(candidate, left, right, overhang + 1, prefix)
    o3 = reverse_complement(candidate["sensor"][:left + 3]) + "g"
    o4 = "ccggt" + reverse_complement(candidate["sensor"][left + 3:])
    o5 = list(o2)
    o5[5] = "g"
    o5 = "".join(o5)
    o6 = list(o4)
    o6[-2] = "c"
    o6 = "".join(o6)
    candidate["oligos"] = {
        prefix + "1": o1, prefix + "2": o2, prefix + "3": o3, prefix + "4": o4,
        prefix + "5": o5, prefix + "6": o6,
    }
    return candidate


def get_by_gene_from_mart_export(
    mart_export=MART_EXPORT,
    utr_data=None,
):
    gene_header_ix = 4
    is_mouse = False
    if "musmusculus" in mart_export:
        gene_header_ix = 5
        is_mouse = True
    if utr_data is None:
        with gzip.open(mart_export, "rt") as handle:
            utr_data = SeqIO.to_dict(SeqIO.parse(handle, "fasta"), lambda rec: rec.description)
    by_gene = {}
    for k in utr_data:
        splits = k.split("|")
        try:
            gene = splits[gene_header_ix]
        except IndexError:  # happens with the mouse export? maybe with the human one too, silently?
            continue
        if gene == "":  # no gene name -> ignored
            continue
        if is_mouse:
            desc = splits[gene_header_ix - 1].lower()
        else:
            desc = splits[gene_header_ix + 1].lower()
        if (  # many show up on the CDS dataset, some in the UTR one as well
            desc.startswith("novel")
            or desc.startswith("uncharacterized")
            or "pseudogene" in desc
            or "antisense" in desc
            or "readthrough" in desc
            or "overlapping transcript" in desc
            or "predicted gene" in desc
            or "riken cdna" in desc
            or "cdna riken" in desc
            or "expressed sequence" in desc
            or "cdna sequence" in desc
            or "opposite strand" in desc
            or desc.startswith("est ")
            or desc.startswith("dna segment")
        ):
            continue
        if gene not in by_gene:
            by_gene[gene] = {}  # so that we can keep track of no-transcript genes
        if utr_data[k].seq == "Sequenceunavailable":
            continue
        by_gene[gene][k] = utr_data[k]
    return by_gene


def design_based_on_mart_export(
    gene="GAPDH",
    mart_export=MART_EXPORT,
    min_transcripts=0,
    variant="CCA", left=42, right=45, min_gc=0.35, max_gc=0.8,
    utr_data=None, by_gene=None, longest_homopolymer=5, total_homopolymer=10,
    max_candidates=None, max_candidates_rand=None, allow_extra_stops=0,
    fa_genome=FA_GENOME, threads=12,
    verbose=True,
):
    arguments = locals()
    del arguments["by_gene"]
    by_gene = by_gene or get_by_gene_from_mart_export(mart_export, utr_data)

    candidates = []
    utrs = by_gene[gene].values()
    if verbose:
        print(f"{len(utrs)} RNA variants recorded:")
    for utr in utrs:
        info = utr.description.split("|")
        if verbose:
            print(f"{len(utr.seq)} bases, {info[3]}")
        utr.annotations = info  # normally a dict
    seen_seqs = set()
    for rec in utrs:
        seq = str(rec.seq)
        if seq in seen_seqs:  # exactly this sequence seen before
            continue
        seen_seqs.add(seq)
        this_candidates = find_arRNAs(
            seq, variant=variant, left=left, right=right, min_gc=min_gc, max_gc=max_gc,
            max_longest_homopolymer=longest_homopolymer, max_total_homopolymer=total_homopolymer,
            max_candidates=max_candidates, max_candidates_rand=max_candidates_rand,
            allow_extra_stops=allow_extra_stops, verbose=max(verbose - 1, 0),
        )
        for cand in this_candidates:
            cand["source_description"] = rec.description
            # cand["source_seq"] = seq
        candidates.extend(this_candidates)

    unique_arRNAs = set()
    unique_candidates = []
    for cand in candidates:
        if cand["sensor"] in unique_arRNAs:
            continue
        unique_arRNAs.add(cand["sensor"])
        unique_candidates.append(cand)

    if verbose:
        print(len(candidates), "candidates,", len(unique_candidates), "unique")
        print("#" * 50)
        for cand in unique_candidates:
            good_for = [rec for rec in utrs if cand["candidate"] in rec.seq]
            if min_transcripts and len(good_for) < min_transcripts:
                continue
            blast_results = filter_blast(
                blast_candidate(cand["sensor"], fa_genome), left=left, right=right,
            )
            cand["n_blast"] = blast_results.shape[0]
            cand = design_oligos(cand, left, right)
            print(cand)
            print(
                len(good_for), "/", len(utrs), ":",
                ", ".join([f"{rec.annotations[3]} ({len(rec.seq)} bases)" for rec in good_for]),
                f"{blast_results.shape[0]} blast hits",
            )
            if blast_results.shape[0]:
                for ix, row in blast_results.drop(columns=["seq"]).iterrows():
                    print(blast_results.loc[[ix]])
                    print_blast_alignment(row.query_seq, row.subject_seq)
                    print("")
            else:
                print("no blast results")
            print("#" * 50)
        print("arguments =", arguments)
        print("-" * 50)

    return unique_candidates


def blast_candidates(candidates, fa_genome=FA_GENOME, threads=12):
    """Input: list of strings to blast hg38 against"""
    with tempfile.TemporaryDirectory() as tempdir:
        seq2seq = {}
        with open(tempdir + "/query.seq", "w") as fp:
            for i, candidate in enumerate(candidates, start=1):
                fp.write(f"> Seq_{i}\n{candidate}\n")
                seq2seq[f"Seq_{i}"] = candidate
        subprocess.run([
            "blastn", "-query", tempdir + "/query.seq",
            "-db", fa_genome, "-task", "blastn",
            "-out", tempdir + "/out.tsv", "-num_threads", f"{threads}",
            "-evalue", "1",  # unlikely to be a problem if evalue higher
            "-outfmt", (
                "6 qaccver saccver pident length nident mismatch gapopen qstart qend "
                "sstart send evalue bitscore qseq sseq qframe sframe sstrand"),
        ])
        results = pd.read_csv(tempdir + "/out.tsv", sep="\t", names=[
            "query_acc", "subject_acc", "percent_identity", "alignment_length",
            "matches", "mismatches", "gap_opens", "qstart", "qend", "sstart", "send",
            "evalue", "bit_score", "query_seq", "subject_seq", "query_frame",
            "subject_frame", "subject_strand",
        ])
    results["seq"] = results.query_acc.map(seq2seq.get)
    results = results.drop(columns=["query_acc"])
    # deal with multiple matches on standard & non-standard chromosome (e.g. _alt).
    # We want to know about multiple mappings on the same subject_acc,
    # but we don't care about multiple mappings on the additional chromosomes (e.g., alt)
    # if there is an equivalent mapping on the main chromosome already.
    # this is pretty slow.
    if results.shape[0]:
        results["chrom"] = results.subject_acc.str.split("_", n=1, expand=True)[0]
    else:
        results["chrom"] = results.subject_acc
    results["is_nonstd_chrom"] = results.chrom != results.subject_acc
    ixs = []
    for x, sdf in results.groupby(
        ["seq", "alignment_length", "matches", "mismatches",
         "gap_opens", "qstart", "qend", "chrom"],
        sort=False,
    ):
        if sdf.is_nonstd_chrom.sum() == 0:  # only standard chromosome matches
            ixs.extend(sdf.index)
            continue
        std_chrom = sdf[~sdf.is_nonstd_chrom]
        if std_chrom.shape[0] == 0:  # no matches on standard chromosome
            nonstd_chrom = sdf[sdf.is_nonstd_chrom]
            if nonstd_chrom.shape[0]:  # add first match from a non-standard chrom
                ixs.append(nonstd_chrom.index[0])
        else:
            ixs.extend(std_chrom.index)  # add all matches on standard chromosome
    ixs = np.sort(ixs)
    results = results.loc[ixs]
    return results


def blast_candidates_in_batches(
    candidates, fa_genome=FA_GENOME, batch_size=10000, threads=12, verbose=True,
):
    all_results = []
    for i in tqdm.trange(
        0, len(candidates), batch_size, desc=f"blasting {len(candidates)} sequences",
        disable=not verbose,
    ):
        this_candidates = candidates[i:i + batch_size]
        all_results.append(blast_candidates(this_candidates, fa_genome, threads=threads))
    results = pd.concat(all_results, ignore_index=True)
    return results


def blast_candidate(candidate, fa_genome=FA_GENOME, threads=12):
    return blast_candidates([candidate], fa_genome, threads=threads)


def filter_blast(blast_results, left, right, min_match=20):
    """Filter out some blast results. Particularly, if the qstart...qend range does not really
    match around the CCA, then it's not relevant.
    So if qstart > left + 3 (matching starts after CCA), ignore it,
    or if qend < left (matching ends before CCA), also ignore it.
    NOTE: we are not looking at whether this sequence may be transcribed.
    NOTE: existence of pseudogenes means some good candidates may be filtered out."""
    blast_results = blast_results[~(blast_results.qstart > left + 3)]
    blast_results = blast_results[~(blast_results.qend < left)]
    blast_results = blast_results[blast_results.alignment_length >= min_match]
    return blast_results.sort_values(["bit_score", "matches"], ascending=False)


def print_blast_alignment(query, subject):
    print(query)
    comparison = []
    comparison = ["*" if q != s else " " for q, s in zip(query, subject)]
    print("".join(comparison))
    print(subject)


def look_at_all_mart_genes(
    mart_export=MART_EXPORT,
    variants=("CCA", "CTA"), left=42, right=45, min_match=30,
    max_candidates=None, max_candidates_rand=None,
    utr_data=None, by_gene=None, fa_genome=FA_GENOME, threads=12, verbose=True,
):
    """Goes through all entries in a BioMart export, and calls:
    - design_based_on_mart_export
    It then annotates all results via blast of the sensor (not trigger) sequences using:
    - annotate_candidates_by_gene_with_blast.
        - blast_candidates_in_batches
            - blast_candidates
    """
    by_gene = by_gene or get_by_gene_from_mart_export(mart_export, utr_data)

    candidates_by_gene = {}  # gene -> list of candidates (dictionary)

    for gene in tqdm.tqdm(
        by_gene, total=len(by_gene), desc="finding sensor candidates", disable=not verbose,
    ):
        # about 1.5 minutes for just CCA; ~2:15 for CCA + CTA
        candidates_by_gene[gene] = []
        for variant in variants:
            candidates = design_based_on_mart_export(
                gene=gene, by_gene=by_gene, verbose=False, variant=variant,
                left=left, right=right,
                max_candidates=max_candidates, max_candidates_rand=max_candidates_rand,
            )
            candidates_by_gene[gene].extend(candidates)

    # make sure sensors are unique
    blast_results, candidates_by_gene = annotate_candidates_by_gene_with_blast(
        candidates_by_gene, fa_genome, 1, min_match, threads=threads, verbose=verbose,
    )

    return by_gene, candidates_by_gene, blast_results


def annotate_candidates_by_gene_with_blast(
    candidates_by_gene, fa_genome=FA_GENOME, max_hits=1, min_match=30, suffix="", threads=12,
    verbose=True,
):
    candidate_seqs = set()
    for gene, candidates in candidates_by_gene.items():
        candidate_seqs.update(cand["sensor"] for cand in candidates)
    candidate_seqs = list(candidate_seqs)
    blast_results = blast_candidates_in_batches(
        candidate_seqs, fa_genome, threads=threads, verbose=verbose,
    )
    gb = blast_results.groupby("seq")
    for gene, candidates in tqdm.tqdm(
        candidates_by_gene.items(), total=len(candidates_by_gene),
        desc="annotating with blast", disable=not verbose,
    ):  # ~2 mins
        for candidate in candidates:
            try:
                sdf = gb.get_group(candidate["sensor"])
                candidate["blast_results" + suffix] = sdf
                candidate["bad_blast" + suffix] = bool(
                    (sdf.matches >= min_match).sum() > max_hits
                )
            except KeyError:
                # 0 hits: this can happen if the sequence has been masked as repetitive in hg38
                candidate["blast_results" + suffix] = None
                candidate["bad_blast" + suffix] = bool(max_hits >= 1)
                # ^ max_hits = 0 for blasting human sensors against mouse genome,
                # so no hits is actually a good blast
    return blast_results, candidates_by_gene  # NOTE: mutates candidates_by_gene



def print_transcripts(gene, mart_export=MART_EXPORT):
    """Helper function to print transcripts from a Biomart export on a single line"""
    by_gene = get_by_gene_from_mart_export(mart_export)
    for name, transcript in by_gene.get(gene, {}).items():
        print(name)
        print(transcript.seq)
        print("#" * 50)


if __name__ == "__main__":
    fire.Fire()
