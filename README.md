# rsomics-bam

`rsomics-bam` is a single command-line product for SAM, BAM, and CRAM
inspection, filtering, conversion, validation, collation, mate repair, sorting,
compressed file editing, read-group editing, duplicate marking, FASTQ import,
and pileup workflows.

## Install

```sh
cargo install rsomics-bam
```

## Commands

| Command | Purpose |
|---|---|
| `addreplacerg` | Add or replace header read groups and record RG tags |
| `cat` | Concatenate BAM files without reencoding alignment blocks |
| `collate` | Group all alignments for each read name with bounded memory |
| `depth` | Compute one per-input depth column at each position |
| `flags` | Convert numeric and symbolic SAM flags |
| `flagstat` | Count records by SAM flag category |
| `fasta` | Convert name-grouped alignments to one FASTA stream |
| `fastq` | Convert name-grouped alignments to one FASTQ stream |
| `fixmate` | Repair mate fields in name-grouped alignments |
| `head` | Print alignment headers and leading records as SAM |
| `index` | Build BAI, CSI, or CRAI random-access indexes |
| `import` | Convert FASTQ reads to unmapped SAM or BAM |
| `merge` | Merge ordered alignment files into BAM |
| `markdup` | Mark or remove duplicate alignments in coordinate order |
| `mpileup` | Generate per-position text pileup |
| `quickcheck` | Validate headers and format-specific end markers |
| `reheader` | Replace a BAM header without reencoding alignment blocks |
| `samples` | List samples and other read-group metadata |
| `sort` | Sort alignments with bounded memory and external runs |
| `view` | Filter records and convert to SAM or BAM |

```sh
rsomics-bam view -b -@ 4 -q 20 -o selected.bam input.bam
rsomics-bam view -c -F UNMAP,SECONDARY input.cram -T reference.fa
rsomics-bam addreplacerg -r ID:lane1 -r SM:sample input.bam -o tagged.bam
rsomics-bam cat lane1.bam lane2.bam -o combined.bam
rsomics-bam collate -m 128M -o grouped.bam input.bam
rsomics-bam depth -a -b targets.bed sample.bam
rsomics-bam fastq -o reads.fq.bgzf name-sorted.bam
rsomics-bam fixmate -m grouped.bam -o fixed.bam
rsomics-bam index -c -m 14 sample.bam
rsomics-bam import read1.fastq read2.fastq -@ 4 -o unmapped.bam
rsomics-bam merge lane1.bam lane2.cram --reference reference.fa -o merged.bam
rsomics-bam markdup -r fixed-and-sorted.bam -o deduplicated.bam
rsomics-bam mpileup -f reference.fa -Q 20 input.bam
rsomics-bam reheader replacement.sam input.bam -o reheadered.bam
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
`fixmate` preserves name-grouped record order while repairing mate coordinates,
flags, TLEN, MC, and MQ. `-m` adds the mate-score tags consumed by `markdup`.
Sanitizer selection, template-cigar and base-modification repair, and alternate
output formats are not yet exposed.
`fasta` and `fastq` consume adjacent QNAME groups, prefer the first alignment
with qualities in each READ1, READ2, or other category, restore original read
orientation, and write one unified stream. Named `.gz`, `.bgz`, and `.bgzf`
outputs are BGZF. Split read-end files, singleton routing, copied auxiliary
tags, UMI/CASAVA headers, and soft-clip removal are not exposed.
`import` accepts one FASTQ, two paired FASTQs, or the corresponding
`-0`/`-s`/`-1`/`-2` forms. A single input derives PAIRED, READ1, and READ2 from
`/1` and `/2` suffixes. Standard output defaults to SAM; `.bam` selects BAM for
named output, and `-O` overrides format selection. Named outputs are
transactional. Read groups, input-order tags, plain/gzip/BGZF input, and BAM
compression workers are supported. Index FASTQs, CASAVA and UMI parsing,
FASTQ-comment auxiliary tags, and CRAM output are not yet exposed.
`addreplacerg` accepts SAM, BAM, and reference-backed CRAM input and writes
SAM or BAM. A new read group can be assembled with repeatable `-r` fields, an
existing ID can be selected with `-R`, and omitting both selects the first
header read group. Overwrite mode replaces every record tag and makes a new
group the only header read group; orphan mode preserves existing groups and
tags. Same-ID header replacement requires `-w`. Named output is transactional.
CRAM output, automatic indexing, and HTSlib format-option strings are not yet
exposed.
`markdup` consumes coordinate-sorted output from `fixmate -m`, marks or removes
duplicates, and supports template and sequence decision modes. It preserves
SAM, BAM, and CRAM input semantics while producing transactional BAM output.
Optical duplicate classification, barcode partitioning, duplicate chains,
non-primary propagation, and alternate output formats are not yet exposed.
`cat` accepts named BAM inputs, repeated files of filenames, and an optional
SAM, BAM, or CRAM header source. It preserves compressed alignment blocks,
record order, and compatible read groups. `reheader` accepts a replacement
header from SAM, BAM, or CRAM and preserves the input BAM records. Both write
one complete BGZF stream, reject input/output aliases, and replace named
outputs only after validation. CRAM concatenation, CRAM output, in-place
editing, and shell header transformation are not exposed.
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
Read-group editing is field-matched across new, existing, implicit, overwrite,
orphan, SAM, BAM, CRAM, and non-string-tag cases.
FASTQ import is field-matched for positional and explicit input modes, read
groups, order tags, compressed input, standard input, SAM, and BAM output.
FASTA/FASTQ extraction also has bytewise stdin and historical-fixture checks.
Named alignment and pileup outputs are committed only after successful
processing. See [PERFORMANCE.md](PERFORMANCE.md) for the representative BAM
throughput, memory, and decoded-output gate.

License: MIT OR Apache-2.0.
