# ARGenus

**Taxonomic inference of antibiotic resistance genes (ARGs) from metagenomic data using flanking sequence analysis**

[![Crates.io](https://img.shields.io/crates/v/argenus.svg)](https://crates.io/crates/argenus)
[![License](https://img.shields.io/crates/l/argenus.svg)](LICENSE)

ARGenus is a bioinformatics tool that simultaneously detects antibiotic resistance genes (ARGs) and identifies their source bacterial genera from metagenomic sequencing data. Unlike existing tools that only detect ARGs, ARGenus provides direct ARG-to-genus linkage through flanking sequence analysis.

## Features

- **Direct ARG-genus linkage**: Identifies the bacterial source of each detected ARG
- **Targeted assembly**: Efficient processing through read filtering and localized assembly
- **SNP verification**: Filters false positives by confirming resistance-conferring mutations
- **High compression**: 450-fold compressed flanking database (27GB → 60MB)
- **Fast processing**: 5-10 minutes per sample with 16 threads

## Installation

### From crates.io

```bash
cargo install argenus
```

### From source

```bash
git clone https://github.com/necoli1822/ARGenus.git
cd argenus
cargo build --release
```

### Dependencies

ARGenus requires the following tools in your PATH:
- [minimap2](https://github.com/lh3/minimap2) - sequence alignment
- [MEGAHIT](https://github.com/voutcn/megahit) - metagenomic assembly

**Note**: ARGenus has been tested on Linux and macOS. Windows is not currently supported due to minimap2 and MEGAHIT dependencies.

## Database Setup

ARGenus requires two databases: an ARG reference database and a flanking sequence database.

```bash
# Step 1: Build ARG database
argenus -b arg -o databases/

# Step 2: Build flanking database (requires NCBI GenBank data, several hours)
argenus -b flank -a databases/AMR_NCBI.mmi -o databases/ -e your@email.com
```

## Usage

### Basic usage

```bash
argenus \
    -1 sample_R1.fastq.gz \
    -2 sample_R2.fastq.gz \
    -a databases/AMR_NCBI.mmi \
    -f databases/flanking.fdb \
    -o results/
```

### Options

```
argenus [OPTIONS]

Required:
    -1, --r1 <FILE>             Forward reads (FASTQ/FASTQ.gz)
    -2, --r2 <FILE>             Reverse reads (FASTQ/FASTQ.gz)
    -a, --arg-db <FILE>         ARG database (.mmi or .fasta)
    -f, --flanking-db <FILE>    Flanking sequence database (.fdb)
    -o, --outdir <DIR>          Output directory [default: .]

Optional:
    -l, --samples <PATH>        Batch mode: directory or sample ID list
    -t, --threads <N>           Total threads [default: auto]
    -s, --threads-per-sample <N> Threads per sample [default: 8]
    -i, --arg-identity <F>      Minimum ARG identity [default: 0.80]
    -c, --arg-coverage <F>      Minimum ARG coverage [default: 0.70]
    -g, --min-contig-len <BP>   Minimum contig length [default: 200]
    -u, --keep-temp             Keep intermediate files
    --all-hits                  Include WildType/NotCovered in output
```

## Output Format

ARGenus produces a tab-delimited file with the following columns:

| Column | Description |
|--------|-------------|
| Sample | Sample identifier |
| ARG_Name | ARG gene name |
| ARG_Class | Antimicrobial drug class |
| Genus | Assigned source genus |
| Confidence | Classification confidence score |
| Specificity | Database prevalence of gene in assigned genus (%) |
| ARG_Identity | ARG sequence identity (%) |
| ARG_Coverage | ARG sequence coverage (%) |
| Contig_Len | Length of assembled contig |
| Upstream_Len | Upstream flanking sequence length |
| Downstream_Len | Downstream flanking sequence length |
| Extension_Method | Contig extension method (strict/flexible) |
| SNP_Status | SNP verification result |
| Top_Matches | Top 5 genus matches with scores |

## Workflow

1. **Read Filtering**: Align reads against ARG database using minimap2
2. **De Novo Assembly**: Assemble filtered reads with MEGAHIT
3. **Contig Extension**: Extend contigs using k-mer overlap analysis
4. **ARG Detection**: Identify ARGs in extended contigs
5. **Genus Classification**: Classify genus using flanking sequence homology
6. **SNP Verification**: Confirm resistance-conferring mutations

## Performance

| Metric | Value |
|--------|-------|
| Processing time | 5-10 min/sample (16 threads) |
| Genus classification rate | ~73% |
| False positive reduction | ~72% (via SNP filtering) |
| Database size | 60 MB (compressed) |

## Comparison with Other Tools

| Feature | ARGenus | RGI | KMA | AMRFinderPlus |
|---------|---------|-----|-----|---------------|
| ARG detection | ✓ | ✓ | ✓ | ✓ |
| Genus classification | ✓ | ✗ | ✗ | ✗ |
| SNP verification | ✓ | ✓ | ✗ | ✓ |
| Targeted assembly | ✓ | ✗ | ✗ | ✗ |

## Citation

If you use ARGenus in your research, please cite:

```
[Citation information to be added upon publication]
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

<table>
  <tbody>
    <tr>
      <td align="center">
        <a href="https://github.com/guinamwee">
          <img src="https://github.com/guinamwee.png" width="100px;" alt=""/><br />
          <sub><b>GNWEE</b></sub>
        </a>
      </td>
      <td align="center">
        <a href="https://github.com/kyuwon-shim-ARL">
          <img src="https://github.com/kyuwon-shim-ARL.png" width="100px;" alt=""/><br />
          <sub><b>kyuwon-shim-ARL</b></sub>
        </a>
      </td>
    </tr>
  </tbody>
</table>

## Contact

[Sunju Kim](n.e.coli.1822@gmail.com)

For questions and feedback, please open an issue on GitHub.
