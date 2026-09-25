# Changelog

All notable changes to ARGenus will be documented in this file.

## [0.4.1] - 2026-09-21

### Added

- **Align path: ARGs the assembly lost are now reported (on by default).** A gene that
  MEGAHIT breaks across contigs fails the coverage floor and disappeared from
  `results.tsv` entirely, even though the read→ARG alignment that selected the reads for
  assembly had already seen it (observed for `aph(3'')-Ib` on the Abramova spike-in:
  0.32–0.54 contig coverage against a 0.70 floor). Those genes are now emitted as
  **align-path rows**, gated on reference breadth and read count.

  These rows carry no flanking, so they are **never attributed**: `Genus`, `Species` and
  `Context` stay `Unknown`, `Limited_By` is `no_context`, and the new `Detection_Path`
  column marks them `align` (assembly rows are `assembly`). They must not be counted as
  genus attributions — they say only "the gene is present, its context could not be
  reconstructed". `Top_Matches` carries the evidence:
  `read_only:<breadth>;alleles:<n>;support:<perfect>/<total>;margin:<lead>`. Read
  `margin` together with `alleles` and `support` — a single surviving allele scores
  `margin:1.00` by definition, however thin its own support.

  Dedup against the assembly path is by **gene**, not reference id: PanRes holds many
  near-identical alleles of one gene (327 TEM, 49 NDM entries), so id-level comparison
  would re-report the same gene under a sibling allele.

  New flags: `--align-path on|off` (default `on`), `--align-min-breadth` (default `0.80`),
  `--align-min-reads` (default `10`). `--align-path off` restores 0.4.0 behaviour exactly.

- **Per-variant SNP/INDEL counts from the blastn traceback.** Detection has used
  blastn dc-megablast since 0.4.0, so each hit already carries a BTOP string; it is
  now parsed instead of discarded. `parse_btop` walks the traceback in reference
  coordinates (handling reverse hits where `send < sstart`) and separates
  substitutions from gaps, counting a gap-open only on the first base of a gap run.
  `results.tsv` gains `N_SNP`, `N_INDEL`, `INDEL_bp` and `Variants`.

  On a minus-strand hit blastn prints the subject reverse-complemented, so the traceback
  letters are the complement of the reference's own bases while the coordinates stay on
  the plus strand. `Variants` reports both on the reference strand, so a token can be
  compared against a point-mutation catalogue directly whichever way the contig ran.
  Align-path rows have no traceback and report `NA`/`-`.

  This exposes the evidence for point-mutation resistance genes without asserting a
  verdict. `snp.rs` never fires on PanRes (no parseable mutation names) and CARD
  separates homologs by curated per-gene bit-score cutoffs rather than mutation
  lists, so the counts are reported for the user to judge.

### Fixed

- **`Limited_By` now says why a row is `Unknown`.** Every `Unknown` row on the
  assembly path reported `none`, which is also what a confident single-genus call
  reports — a row with no answer was indistinguishable from the most resolved row in
  the file. Unknown rows now name their cause: `flank_too_short` (contig gave < 50 bp
  of flank on both sides), `gene_not_in_reference` (the flanking DB holds nothing for
  this gene), `context_unmatched` (the DB holds the gene but no reference flank
  matched this sample), `alignment_failed`, and `no_flanking_db`. The align path uses
  `no_context`; answered rows report `none`, `flank_truncated` or `flank_shared` (the
  latter two renamed from `query`/`biology` — see Compatibility).

  The causes have different remedies — deeper sequencing fixes `flank_too_short`,
  only a larger reference fixes `gene_not_in_reference`, and `context_unmatched` may
  indicate a host absent from the reference.

### Compatibility

- `results.tsv` gains **five** columns at the end — `Detection_Path`, `N_SNP`, `N_INDEL`,
  `INDEL_bp`, `Variants` — taking the row from 24 to 29 fields. Readers that index by
  column name are unaffected; readers that assume a fixed column count need updating.
- **The align path adds rows by default.** A run against the same data can now report more
  ARGs than 0.4.0 did. Every added row has `Detection_Path` = `align`, `Genus` = `Unknown`
  and `Limited_By` = `no_context`, so filtering on any of those three reproduces the old
  row set exactly; `--align-path off` does the same at the source. Note that align rows
  report `ARG_Identity` `0.0` (there is no contig-vs-reference alignment to score) and put
  reference breadth in `ARG_Coverage`, so an identity filter silently drops them.
- **`Limited_By` on answered rows was renamed.** `query` → `flank_truncated` and
  `biology` → `flank_shared`. The old names said how the cause was categorised rather than
  what the reader would look for. Code matching on `"query"`/`"biology"` must be updated —
  it will not error, it will just stop matching.
- `Limited_By` values on `Unknown` rows change from `none` to the specific causes
  listed above. Code that tested `limited_by == "none"` to mean "resolved" should
  test the `Genus` column instead.

## [0.4.0] - 2026-07-17

### Added

- **Self-contained binary: embedded GTDB phylogenetic-context tables.** The GTDB
  genus distance/lineage tables and the conformal calibration are now embedded in
  the binary (zstd-compressed, `src/embedded/`), so a built binary is self-contained
  — only the flanking DB (FDB) is downloaded separately. This eliminates "old table
  + new constant" version drift: the kernel constants (λ=0.3, coherence radius 0.5,
  absent-patristic 3.0) are calibrated to the embedded GTDB patristic scale and ship
  together. A `genus_dist.tsv`/`genus_lineage.tsv`/`conformal.tsv` in the db dir still
  overrides the embedded copy when present.

### Changed

- **Kernel-posterior genus/family classification.** Reworked the phylogenetic-context
  scoring toward a posterior over lineages. Validated with no regression and net
  improvement: Zymo genus-wrong 5.8→4.3%, GTDB (7,077 genomes) family 91.1→92.2% /
  genus-wrong 14.4→12.6%, RAPID (775 MAGs) detection 98.7% byte-identical with host
  Bracken corroboration 90.3→93.4%. Detection specificity unchanged (ResFinder golden:
  recall 99.9%, ARG-negative specificity 99.5%).

## [0.3.1] - 2026-07-12

### Changed

- **`--db-dir` flanking-DB discovery**: when a database folder contains more than one
  `*.fdb` file (e.g. both `flanking_1kbp.fdb` and `flanking_5kbp.fdb`), ARGenus now
  errors and asks you to choose one with `-f/--flanking-db`, instead of silently
  picking the alphabetically-first file.

## [0.3.0] - 2026-07-12

### Added

- **Honest 4-axis classification report**: replaces the single forced "source genus" call
  - **Genus**: single genus, `multi-genus(N):A/B/C/…` for promiscuous genes, or `Unknown`
  - **Species**: stricter, separately-thresholded call (`--species-identity`)
  - **Context**: replicon of the matched flanking (`plasmid` / `chromosome` / `ambiguous` / `NA`), from PLSDB provenance
  - **Specificity**: how gene-specific the flanking evidence is
- **Per-locus reassembly** (`--reassemble`, opt-in): core/flank read split + SPAdes to recover classifiable flanking for stalled loci (new `reassemble` module)
- **Per-locus exports** (`--emit-*`): write gene / flanking sequences for resolved / no-flank-match / gene-not-in-DB classes
- **Pluggable read filter** (`--mapper strobealign|minimap2|bwa-mem2`)
- **Contig-only mode** (`--classify-contigs`) to classify a pre-assembled FASTA

### Changed

- **Bounded-memory extension**: extension consensus accumulated as per-position base counts, each contig end capped (`--max-extension`), so runtime/RAM no longer blow up on high-coverage/repetitive loci

## [0.2.3] - 2026-07-05

### Fixed

- **Dependencies**: Updated the transitive `lz4_flex` lockfile entry from the yanked 0.11.5 to 0.11.6 (no source changes).

## [0.2.2] - 2026-07-05

### Fixed

- **Package metadata**: Corrected placeholder author to `Sunju Kim <n.e.coli.1822@gmail.com>` and LICENSE copyright holder

## [0.2.1] - 2026-02-14

### Added

- **Contig_ID column**: New `Contig_ID` column in output TSV after `Sample` column
  - Links each ARG detection to its source contig (e.g., contig_1, contig_2)
  - Enables tracing ARG variants back to assembled contigs

### Changed

- **Package naming**: Standardized to `ARGenus` (capital ARG) across all configurations
- **Documentation**: Updated README with Contig_ID in Output Format table

## [0.2.0] - 2026-02-13

### Added

- **Dual database mode**: New `--mode short|long` option for flanking database building
  - `short` mode (1,000 bp): High coverage (97.6%) from GenBank + PLSDB
  - `long` mode (5,000 bp): High resolution (92.8%) from NCBI nt_prok
- **New source file**: `flanking_db_ntprok.rs` for 5,000 bp database construction using BLASTN
- **Streaming FDB builder**: Memory-efficient processing for large datasets
  - `--sorted` flag for pre-sorted input (streaming mode)
  - External merge sort support
  - Works with 8-16 GB RAM for 190+ GB datasets
- **Auto-download taxdump**: Automatic download of NCBI taxonomy files if not present

### Changed

- **FDB format v2**: Enhanced binary format with improved compression
  - ~22x compression ratio (194 GB → 8.7 GB for 5,000 bp database)
  - O(1) random access via gene name index
- **Improved classifier**: Enhanced genus classification with confidence metrics
- **Updated dependencies**: All dependencies updated to latest stable versions

### Database Statistics

| Database | Records | Genes | Coverage | Genus Resolution |
|----------|---------|-------|----------|------------------|
| 1,000 bp | 1,069,848 | 11,835 | 97.6% | 83.9% |
| 5,000 bp | 23,184,244 | 11,092 | 91.5% | 92.8% |

### Performance

- Database query: sub-millisecond per ARG match
- FDB building: ~700 MB peak memory (streaming mode)
- Processing: 5-10 minutes per sample (16 threads)

## [0.1.5] - 2026-02-06

### Initial Release

- ARG detection using minimap2 alignment
- Genus classification via flanking sequence analysis
- SNP verification for point mutation ARGs
- Targeted assembly workflow with MEGAHIT
- Compressed FDB format for flanking database
