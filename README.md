# rsomics-bam

`rsomics-bam` is a single command-line product for SAM, BAM, and CRAM
inspection, filtering, conversion, validation, and pileup workflows.

## Install

```sh
cargo install rsomics-bam
```

## Commands

| Command | Purpose |
|---|---|
| `depth` | Compute one per-input depth column at each position |
| `flags` | Convert numeric and symbolic SAM flags |
| `flagstat` | Count records by SAM flag category |
| `head` | Print alignment headers and leading records as SAM |
| `mpileup` | Generate per-position text pileup |
| `quickcheck` | Validate headers and format-specific end markers |
| `samples` | List samples and other read-group metadata |
| `view` | Filter records and convert to SAM or BAM |

```sh
rsomics-bam view -b -@ 4 -q 20 -o selected.bam input.bam
rsomics-bam view -c -F UNMAP,SECONDARY input.cram -T reference.fa
rsomics-bam depth -a -b targets.bed sample.bam
rsomics-bam mpileup -f reference.fa -Q 20 input.bam
```

The commands accept SAM, BAM, and CRAM where their help declares those input
formats. `view` writes SAM or BAM; CRAM output is intentionally unavailable.
Indexed region queries require a usable BAI, CSI, or CRAI. A non-zero thread
request for CRAM decoding fails explicitly because the current decoder does
not provide ordered parallel decoding.

Stable behavior is tested against samtools 1.24 across SAM, BAM, and CRAM.
Named alignment and pileup outputs are committed only after successful
processing. See [PERFORMANCE.md](PERFORMANCE.md) for the representative BAM
throughput, memory, and decoded-output gate.

License: MIT OR Apache-2.0.
