# rsomics-bam

`rsomics-bam` is a single command-line product for SAM, BAM, and CRAM
inspection, filtering, conversion, validation, collation, sorting, and pileup
workflows.

## Install

```sh
cargo install rsomics-bam
```

## Commands

| Command | Purpose |
|---|---|
| `collate` | Group all alignments for each read name with bounded memory |
| `depth` | Compute one per-input depth column at each position |
| `flags` | Convert numeric and symbolic SAM flags |
| `flagstat` | Count records by SAM flag category |
| `head` | Print alignment headers and leading records as SAM |
| `index` | Build BAI, CSI, or CRAI random-access indexes |
| `merge` | Merge ordered alignment files into BAM |
| `mpileup` | Generate per-position text pileup |
| `quickcheck` | Validate headers and format-specific end markers |
| `samples` | List samples and other read-group metadata |
| `sort` | Sort alignments with bounded memory and external runs |
| `view` | Filter records and convert to SAM or BAM |

```sh
rsomics-bam view -b -@ 4 -q 20 -o selected.bam input.bam
rsomics-bam view -c -F UNMAP,SECONDARY input.cram -T reference.fa
rsomics-bam collate -m 128M -o grouped.bam input.bam
rsomics-bam depth -a -b targets.bed sample.bam
rsomics-bam index -c -m 14 sample.bam
rsomics-bam merge lane1.bam lane2.cram --reference reference.fa -o merged.bam
rsomics-bam mpileup -f reference.fa -Q 20 input.bam
rsomics-bam sort -m 768M -o sorted.bam input.bam
```

The commands accept SAM, BAM, and CRAM where their help declares those input
formats. `view` writes SAM or BAM; CRAM output is intentionally unavailable.
Indexed region queries require a usable BAI, CSI, or CRAI. A non-zero thread
request for CRAM decoding fails explicitly because the current decoder does
not provide ordered parallel decoding. `index` uses up to four additional
workers when `-@` is omitted; pass `-@ 0` for one-thread indexing. Named index
outputs are committed only after the complete index has been built and parsed.
`collate` accepts SAM, BAM, and CRAM input and writes BAM with contiguous QNAME
groups. Its total record-memory budget and external merge fan-in are bounded;
the order between groups is intentionally unspecified. Fast-mode filtering and
early-pair buffering are not yet exposed.
`merge` accepts named, already ordered SAM, BAM, and CRAM inputs and writes BAM.
It reconciles reference, read-group, and program records, validates both the
declared and observed order, and keeps its per-input read-ahead bounded. It
supports coordinate, natural query-name, bytewise query-name, and
template-coordinate order. Named outputs replace their destination only after
the complete BGZF stream passes validation.
`sort` accepts SAM, BAM, and CRAM input and writes BAM. It supports coordinate,
natural query-name, bytewise query-name, and template-coordinate order. Its
memory option is a total record budget, external merges use bounded fan-in,
and named outputs replace their destination only after a complete BGZF stream
passes validation. It uses up to four additional workers when `-@` is omitted;
pass `-@ 0` for one-thread sorting.

Stable behavior is tested against samtools 1.24 across SAM, BAM, and CRAM.
Named alignment and pileup outputs are committed only after successful
processing. See [PERFORMANCE.md](PERFORMANCE.md) for the representative BAM
throughput, memory, and decoded-output gate.

License: MIT OR Apache-2.0.
