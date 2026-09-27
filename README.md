# ARGenus

**ARG detection with honest genus / species / replicon-context reporting from metagenomes**

[![Crates.io](https://img.shields.io/crates/v/argenus.svg)](https://crates.io/crates/argenus)
[![License](https://img.shields.io/crates/l/argenus.svg)](LICENSE)

ARGenus detects antibiotic resistance genes (ARGs) in metagenomic reads and, for each
ARG, reports the bacterial **source** it was found in — using the DNA that *flanks*
the gene on the assembled contig. Unlike tools that only detect ARGs, ARGenus links
each ARG to a genus (and, when the flanking is specific enough, a species) and tells
you how much to trust that link.

## What's new in 0.4.1

- **ARGs the assembly lost are now reported (on by default).** When a gene is broken
  across contigs it fails the coverage floor and used to vanish from `results.tsv`, even
  though the read→ARG alignment had already seen it. Those genes now appear as
  **align-path rows** (`Detection_Path` = `align`). They carry no flanking, so they are
  never attributed — `Genus`/`Species`/`Context` stay `Unknown` and `Limited_By` is
  `no_context`; the read evidence is in `Top_Matches`. This means a 0.4.1 run can report
  more rows than 0.4.0 on the same data; `--align-path off` restores the old behaviour.
- **`Limited_By` says why a row is `Unknown`.** Previously every unresolved row
  reported `none` — the same value a confident single-genus call reports. Unknown rows
  now name their cause (`gene_not_in_reference`, `context_unmatched`,
  `flank_too_short`, …), so you can tell whether more sequencing would help or only a
  larger reference would.
- **Per-variant SNP/INDEL counts.** `N_SNP`, `N_INDEL`, `INDEL_bp` and `Variants` are
  parsed from the blastn traceback and reported for every hit. ARGenus does not use
  them to accept or reject a call; they are there so you can judge point-mutation
  resistance genes yourself.

## What's new in 0.4.0

- **Self-contained binary.** The GTDB genus distance/lineage tables and the conformal
  calibration are embedded in the binary, so only the flanking DB is downloaded
  separately and the kernel constants always match the table scale they were calibrated
  to. A `genus_dist.tsv` / `genus_lineage.tsv` / `conformal.tsv` in the db dir still
  overrides the embedded copy.
- **Kernel-posterior genus/family classification.** Genus/family scoring became a
  posterior over lineages rather than a count-weighted identity score, so a sparse genus
  borrows evidence from close relatives instead of losing to a deeply-sampled one.
  Detection itself is unchanged.

## What's new in 0.3.1

- **`--db-dir` picks one flanking DB unambiguously.** If a database folder holds more
  than one `.fdb` (e.g. both the 1 kbp and 5 kbp resolutions), ARGenus now stops and
  asks you to choose one with `-f/--flanking-db`, instead of silently loading the
  alphabetically-first file.

## What's new in 0.3.0

ARGs frequently sit on **plasmids and other mobile elements** that move between
genera, so forcing a single "source genus" is often wrong. 0.3.0 replaces the single
genus call with an **honest 4-axis report** so an ambiguous locus looks ambiguous
instead of confidently wrong:

- **Genus** — a single genus, or `multi-genus(N):A/B/C/…` when several genera share
  the flanking within a small identity margin (a promiscuous gene), or `Unknown`.
- **Species** — a stricter, separately-thresholded call (`--species-identity`);
  `multi-species(N):…` / `Unknown` when the flanking isn't near-identical.
- **Context** — replicon of the matched flanking: `plasmid` / `chromosome` /
  `ambiguous` / `NA`, derived from PLSDB provenance of the flanking references.
  `chromosome` + single-genus is trustworthy; `plasmid` + multi-genus is not.
- **Specificity** — how gene-specific (breadth) the flanking evidence is.

Also new:
- **Per-locus reassembly** (`--reassemble`, opt-in): core/flank read split + SPAdes
  to recover a classifiable flanking for stalled loci.
- **Per-locus exports** (`--emit-*`): write gene / flanking sequences (assembled or
  as reads) for the resolved / no-flank-match / gene-not-in-DB classes.
- **Pluggable read filter** (`--mapper strobealign|minimap2|bwa-mem2`).
- **Contig-only mode** (`--classify-contigs`) to classify a pre-assembled FASTA.
- **Bounded-memory extension** — extension consensus is accumulated as per-position
  base counts, and each contig end is capped (`--max-extension`), so runtime and RAM
  no longer blow up on high-coverage/repetitive loci.

## Features

- **Direct ARG→genus/species linkage** through flanking sequence analysis
- **Honest reporting**: multi-genus / multi-species / replicon context, not a forced call
- **Targeted assembly**: read filtering → MEGAHIT → k-mer extension (optional SPAdes reassembly)
- **SNP verification**: confirms resistance-conferring mutations for point-mutation ARGs
- **Compact database**: custom FDB format (zstd-compressed, on-demand gene blocks)
- **Scales with depth**: bounded extension memory; multi-threaded

## Dependencies

Nothing is vendored. ARGenus needs a Rust toolchain to build, a few external aligners on
`PATH` at run time, and two databases.

### Build

Rust **1.88** or newer — the dependency tree already requires it.

### External tools

Resolved at startup by searching `PATH`. If a required one is missing ARGenus aborts with
`<tool> not found in PATH` before doing any work.

| Tool | Needed for | Role | Path override |
|---|---|---|---|
| **`blastn`** + **`makeblastdb`** (BLAST+) | **every run**, `--classify-contigs` included | contig→ARG detection (`dc-megablast`) | `--blastn-path` (blastn only) |
| **`minimap2`** | **every run** | flanking alignment for genus/species classification; also the read filter under `--mapper minimap2` | — |
| `megahit` | the read pipeline (not `--classify-contigs`) | assembly of the filtered reads | — |
| `strobealign` | `--mapper strobealign` (**default**) | read filtering (Step 1) | `--strobealign-path` |
| `bwa-mem2` | `--mapper bwa-mem2` | read filtering (Step 1) | `--bwa-mem2-path` |
| `paftools.sh` | optional | SAM→PAF for strobealign / bwa-mem2; a built-in converter is used when absent | `--paftools-path` |
| `spades.py` | `--reassemble` | per-locus reassembly | `--spades-path` |
| `blastdbcmd` | `--build-db --mode long` only | pulling 5,000 bp flanks out of `nt_prok` | `--blastdbcmd-path` |

So a standard read run needs **blastn, makeblastdb, minimap2, megahit and strobealign**;
`--classify-contigs` needs only **blastn, makeblastdb and minimap2**.

`--build-db` looks up none of the above — it returns before the run-time tools are resolved.
Building the 5,000 bp flanking DB (`--mode long`) instead takes `blastn` and `blastdbcmd` from
`--blastn-path` / `--blastdbcmd-path`, which are mandatory there and are **not** searched on
`PATH`.

### Data

An ARG reference FASTA and a flanking `.fdb`; neither ships with the crate. See
[Databases](#databases) for the pre-built download or how to build your own.

## Installation

### From crates.io

```bash
cargo install argenus
```

### From source

```bash
git clone https://github.com/necoli1822/ARGenus.git
cd ARGenus
cargo build --release
```

## Databases

ARGenus needs an **ARG reference** (a FASTA) and a **flanking database**
(`.fdb`). These are **not** shipped with the crate (too large) — download the
pre-built set or build your own.

### Download the pre-built database (recommended)

The **1,000 bp** PanRes database is published as a GitHub Release asset. Download,
extract, and point ARGenus at the folder with `-d`:

```bash
# 297 MB download, ~312 MB extracted into ./db/
curl -L -o argenus-db-1kbp-v0.3.0.tar.gz \
  https://github.com/necoli1822/ARGenus/releases/download/db-1kbp-v0.3.0/argenus-db-1kbp-v0.3.0.tar.gz
tar xzf argenus-db-1kbp-v0.3.0.tar.gz       # -> db/

# -d auto-discovers everything inside the folder
argenus -d db/ -1 R1.fq.gz -2 R2.fq.gz -o results/
```

ARGenus ships two flanking-DB resolutions:

- **1,000 bp** (genus-level, high coverage) — the GitHub Release bundle above.
- **5,000 bp** (species-level, higher resolution) — much larger (~9 GB), archived on
  **Zenodo**: <https://doi.org/10.5281/zenodo.21321983> (CC BY-NC 4.0).

**Picking a resolution.** With `-d`, ARGenus loads the single `.fdb` it finds in the
folder. If the folder holds more than one — e.g. both resolutions — it stops and asks
you to choose one **explicitly with `-f`**:

```bash
# 1,000 bp (genus-level)
argenus -f flanking_1kbp.fdb -a db/PanRes_genes_v1.0.0.fa -1 R1.fq.gz -2 R2.fq.gz -o results/
# 5,000 bp (species-level)
argenus -f flanking_5kbp.fdb -a db/PanRes_genes_v1.0.0.fa -1 R1.fq.gz -2 R2.fq.gz -o results/
```

You can also build either database yourself (see below).

> **Database license:** the pre-built database is derived from
> [PanRes v1.0.0](https://doi.org/10.5281/zenodo.8055116) and is licensed
> **[CC BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/) (non-commercial)** —
> separate from the MIT-licensed ARGenus code. See [NOTICE](NOTICE) and
> [License & citation](#license--citation).

The bundle contains:

| File | Size | Role |
|------|------|------|
| `PanRes_genes_v1.0.0.fa` | 13 MB | ARG reference (FASTA) |
| `flanking_1kbp.fdb` | 293 MB | Flanking database (**required**) |
| `plasmid_contigs.txt` | 372 KB | Context axis (plasmid/chromosome) |
| `contig_species.tsv` | 6.9 MB | Species axis |

`-d/--db-dir` auto-discovers the ARG reference, the flanking `.fdb`, and the optional
`plasmid_contigs.txt` / `contig_species.tsv` inside the folder. Explicit
`-a/-f/--plasmid-contigs/--species-map` still override what's found there.

> **ARG reference format.** The reference ships as a **FASTA**, which every aligner can
> use: ARGenus indexes it for minimap2 on the fly (the correct `-x asm20` / `-x sr`
> preset per step) and lets strobealign / bwa-mem2 build their own index. You may pass a
> prebuilt minimap2 `.mmi` to `-a` for a small speed-up, but it must be built with
> ARGenus's presets or minimap2 will override them — so the FASTA is the recommended,
> foolproof choice.

### Build your own

```bash
# Build the ARG reference (e.g. from PanRes)
argenus -b arg -x panres -o databases/

# Build a 1,000 bp flanking DB (GenBank/PLSDB)
argenus -b flank --mode short -a databases/PanRes_genes_v1.0.0.fa -o databases/ -e you@email.com

# Compress an existing flanking TSV to FDB (external-sort build)
argenus -b fdb -a flanking.tsv -o databases/flanking_1kbp.fdb
```

The two side files (`plasmid_contigs.txt`, `contig_species.tsv`) are auto-loaded from
beside the `.fdb` if present, and enable the Context and Species axes respectively. You
can also pass them explicitly with `--plasmid-contigs` / `--species-map`.

## Usage

```bash
# Single sample (with a database folder — see Databases below)
argenus -d db/ -1 R1.fq.gz -2 R2.fq.gz -o results/

# ...or point at each database file explicitly
argenus -1 R1.fq.gz -2 R2.fq.gz \
    -a databases/PanRes_genes_v1.0.0.fa \
    -f databases/flanking_1kbp.fdb \
    -o results/

# Batch: a directory (auto-detect *_R[12].fastq.gz) or an ID-list file
argenus -l fastq_dir/ -d db/ -o results/
```

Results are written to `results/results.tsv`.

### Common options

| Option | Default | Description |
|--------|---------|-------------|
| `-1, -2 <FILE>` | — | Paired FASTQ(.gz); comma-separated for multiple samples |
| `-l, --samples <PATH>` | — | Batch: ID-list file or directory of FASTQs |
| `-a, --arg-db <FILE>` | — | ARG reference (`.mmi` or FASTA) |
| `-f, --flanking-db <FILE>` | — | Flanking database (`.fdb`) |
| `-d, --db-dir <DIR>` | — | Auto-discover ARG ref + flanking `.fdb` (+ side files) in one folder |
| `-o, --outdir <DIR>` | `.` | Output directory |
| `-t, --threads <N>` | auto | Total threads |
| `--mapper <TOOL>` | `strobealign` | Read filter: `strobealign` / `minimap2` / `bwa-mem2` |
| `--ref-fasta <FILE>` | derived | FASTA reference for strobealign/bwa-mem2 |
| `-i, --arg-identity <F>` | `0.80` | Min identity for ARG detection |
| `-c, --arg-coverage <F>` | `0.70` | Min coverage for ARG detection |
| `--align-path <on\|off>` | `on` | Also report ARGs seen only in the read alignment (never attributed) |
| `--align-min-breadth <F>` | `0.80` | Min read-covered reference breadth for an align-path row |
| `--align-min-reads <N>` | `10` | Min aligned reads for an align-path row |
| `-n, --max-flanking <BP>` | `1000` | Flanking length used for classification |
| `-u, --keep-temp` | off | Keep per-sample intermediates |
| `-v, --verbose` | off | Progress to stderr |

### Honest-reporting options

| Option | Default | Description |
|--------|---------|-------------|
| `--genus-identity <F>` | `0.90` | Min flanking identity to separate genera |
| `--species-identity <F>` | `0.96` | Min flanking identity to call species (0 disables) |
| `--context-plasmid-frac <F>` | `0.5` | Plasmid-fraction ≥ this → Context `plasmid` |
| `--context-chromosome-frac <F>` | `0.1` | Plasmid-fraction ≤ this → Context `chromosome` |
| `--plasmid-contigs <FILE>` | auto | Plasmid accessions for the Context axis |
| `--species-map <FILE>` | auto | `contig<TAB>species` for the Species axis |

### Assembly / reassembly / exports

| Option | Default | Description |
|--------|---------|-------------|
| `--max-extension <BP>` | `0` (=2×max-flanking) | Cap bp added to each contig end (stops runaway extension; no effect on classification) |
| `--reassemble` | off | Per-locus core/flank reassembly (SPAdes) for stalled loci |
| `--spades-path <FILE>` | `spades.py` | SPAdes for `--reassemble` |
| `--reassemble-jobs <N>` | `4` | Concurrent SPAdes jobs |
| `--classify-contigs <FASTA>` | — | Classify a pre-assembled contig FASTA (skip read pipeline) |
| `--emit-class/-part/-state <LIST>` | — | Per-locus FASTA exports (opt-in) |

Run `argenus --help` for the complete list.

## Output format (`results.tsv`)

Tab-delimited, one row per ARG locus (29 columns):

| Column | Description |
|--------|-------------|
| Sample | Sample identifier |
| Contig_ID | Contig identifier |
| ARG_Name | ARG gene name |
| ARG_Class | Antimicrobial class |
| **Genus** | Source genus, or `multi-genus(N):…`, or `Unknown` |
| **Species** | Source species, or `multi-species(N):…`, or `Unknown` |
| Confidence | Mean flanking identity of the call |
| Specificity | Gene-specificity (breadth) of the flanking evidence |
| **Context** | `plasmid` / `chromosome` / `ambiguous` / `NA` |
| ARG_Identity | ARG sequence identity (`0.0` on align-path rows — there is no contig alignment to score) |
| ARG_Coverage | ARG sequence coverage; on align-path rows this is read-covered reference breadth |
| Contig_Len | Assembled contig length |
| ARG_Start / ARG_End | ARG position on the contig |
| Upstream_Len / Downstream_Len | Flanking length recovered on each side |
| Extension_Method | `strict` / `flexible` / `reassemble` |
| Top_Matches | Top genus candidates with scores |
| **Credible_Set** | Genera whose posterior mass reaches the conformal threshold |
| **Resolution_Rank** | Rank the evidence supports — `genus`, `family`, `order`, … |
| **Resolution_Taxon** | Taxon name at `Resolution_Rank` |
| Support | Posterior mass of the credible set |
| Resolution_Distance | Mash radius of the credible set |
| **Limited_By** | Why resolution stopped — see below |
| Detection_Path | `assembly`, or `align` for a read-only row that carries no flanking |
| **N_SNP / N_INDEL / INDEL_bp** | Substitutions and gaps against the reference allele |
| **Variants** | The substitutions themselves, in reference coordinates |

### Align-path evidence (`Top_Matches`)

An align-path row has no flanking to classify, so `Top_Matches` carries its read evidence
instead of genus candidates:

```text
read_only:0.94;alleles:3;support:18/22;margin:0.28
```

| Field | Meaning |
|---|---|
| `read_only` | Reference breadth covered by reads |
| `alleles` | Sibling alleles in this cluster that cleared the thresholds |
| `support` | Winner's perfect-match reads / total reads on it |
| `margin` | Winner's lead over the runner-up, **on perfect-match reads only** |

`margin` measures only the separation *between sibling alleles*. It is not a confidence that
the gene is present, nor that the winning allele is right in absolute terms:

- a single surviving candidate scores `margin:1.00` by definition — there is nothing to
  compare against — however weak its own support;
- `margin` is undefined when the winner has zero perfect-match reads, and is also reported as
  `1.00` in that case if it is the only candidate.

`support` is what tells you a call is thin. A row reading `alleles:1;support:0/22;margin:1.00`
means the gene is covered by reads but **no** read matches this reference exactly, so the
allele label is a nearest neighbour, not an identification. Always read `margin` together with
`alleles` and `support`.

### `Limited_By`

On rows that called a genus:

| Value | Meaning |
|---|---|
| `none` | A single genus — nothing limited it |
| `flank_truncated` | The contig edge cut the flank short; deeper data may help |
| `flank_shared` | Full flank recovered, but the context is genuinely shared across genera |

> `flank_truncated` and `flank_shared` were called `query` and `biology` up to 0.4.0.

On rows where `Genus` is `Unknown`, the cause is named, because the remedies differ:

| Value | Meaning |
|---|---|
| `no_context` | Align-path row — no contig, so context was never attempted |
| `flank_too_short` | Contig gave < 50 bp of flank on both sides |
| `gene_not_in_reference` | The flanking DB holds nothing for this gene — only a larger reference helps |
| `context_unmatched` | The DB holds the gene, but no reference flank matched this sample — possibly a host absent from the reference |
| `alignment_failed` | The aligner errored on this locus |
| `no_flanking_db` | Run without a flanking database |

**Reading it:** `chromosome` + single Genus + single Species = trustworthy.
`multi-genus(N)` and/or `plasmid` = a promiscuous / mobile gene — the genus is not a
reliable single source. An `Unknown` genus is not a failure to report: check
`Limited_By` to see whether more data would help.

## License & citation

ARGenus uses **two different licenses** — one for the code, one for the database:

| Component | License | Notes |
|-----------|---------|-------|
| ARGenus source code | [MIT](LICENSE) | Free for any use, including commercial |
| Pre-built database (`argenus-db-*.tar.gz`) | [CC BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/) | Derived from PanRes — **non-commercial only** |

The pre-built database is derived from the **PanRes v1.0.0** resistance-gene collection
([Zenodo, DOI 10.5281/zenodo.8055116](https://doi.org/10.5281/zenodo.8055116), CC BY-NC 4.0), which aggregates
ResFinder, ResFinderFG, CARD, MEGARes, NCBI AMRFinderPlus, ARG-ANNOT, and BacMet.
See [NOTICE](NOTICE) for full attribution and the modifications ARGenus makes.

If you use ARGenus, please cite this tool (citation to be added upon publication)
**and** the PanRes / ARGprofiler database:

```
Martiny H-M, et al. ARGprofiler—a pipeline for large-scale analysis of antimicrobial
resistance genes and their flanking regions in metagenomic datasets.
Bioinformatics 40(3):btae086 (2024). https://doi.org/10.1093/bioinformatics/btae086
```

## Contact

Created and maintained by **Sunju Kim** ([ORCID 0000-0002-2384-2425](https://orcid.org/0000-0002-2384-2425)).
Questions and bug reports: open an issue at https://github.com/necoli1822/ARGenus
