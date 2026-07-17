# Changelog

All notable changes to ARGenus will be documented in this file.

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
