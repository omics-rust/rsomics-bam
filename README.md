# rsomics-bam

`rsomics-bam` is a single command-line product for SAM, BAM, and CRAM
inspection, filtering, conversion, validation, collation, mate repair, sorting,
compressed file editing, read-group editing, duplicate marking, FASTQ import,
padded-reference projection, alignment reset, amplicon sequencing, pileup, and
content-checksum and alignment-to-interval workflows.

## Install

```sh
cargo install rsomics-bam
```

## Commands

| Command | Purpose |
|---|---|
| `addreplacerg` | Add or replace header read groups and record RG tags |
| `ampliconclip` | Clip primer sequence from coordinate-ordered BAM alignments |
| `ampliconstats` | Report amplicon coverage and template statistics |
| `bedcov` | Append alignment coverage totals to BED regions |
| `calmd` | Recalculate alignment MD and NM tags from a reference |
| `cat` | Concatenate BAM files without reencoding alignment blocks |
| `checksum` | Compute content checksums or merge prior reports |
| `collate` | Group all alignments for each read name with bounded memory |
| `coverage` | Summarize reads, breadth, depth, and quality by reference |
| `cram-size` | Report CRAM storage by content ID, codec, and data series |
| `depad` | Project padded-reference alignments into unpadded coordinates |
| `depth` | Compute one per-input depth column at each position |
| `flags` | Convert numeric and symbolic SAM flags |
| `flagstat` | Count records by SAM flag category |
| `fasta` | Convert name-grouped alignments to one FASTA stream |
| `fastq` | Convert name-grouped alignments to one FASTQ stream |
| `fixmate` | Repair mate fields in name-grouped alignments |
| `head` | Print alignment headers and leading records as SAM |
| `index` | Build BAI, CSI, or CRAI random-access indexes |
| `import` | Convert FASTQ reads to unmapped SAM or BAM |
| `idxstats` | Report mapped and unmapped counts from an index or sorted scan |
| `merge` | Merge ordered alignment files into BAM |
| `markdup` | Mark or remove duplicate alignments in coordinate order |
| `mpileup` | Generate per-position text pileup |
| `quickcheck` | Validate headers and format-specific end markers |
| `reheader` | Replace a BAM header without reencoding alignment blocks |
| `reset` | Restore primary alignments to unaligned reads |
| `samples` | List samples and other read-group metadata |
| `sort` | Sort alignments with bounded memory and external runs |
| `stats` | Produce comprehensive alignment statistics |
| `to-bed` | Convert alignments to BED6, BED12, or BEDPE |
| `view` | Filter records and convert to SAM or BAM |

```sh
rsomics-bam view -b -@ 4 -q 20 -o selected.bam input.bam
rsomics-bam view -c -F UNMAP,SECONDARY input.cram -T reference.fa
rsomics-bam addreplacerg -r ID:lane1 -r SM:sample input.bam -o tagged.bam
rsomics-bam ampliconclip -b primers.bed input.bam -o clipped.bam
rsomics-bam ampliconstats primers.bed clipped.bam -o amplicons.txt
rsomics-bam bedcov targets.bed sample.bam
rsomics-bam calmd -b -o recalculated.bam input.bam reference.fa
rsomics-bam cat lane1.bam lane2.bam -o combined.bam
rsomics-bam checksum -a -o sample.chk input.bam
rsomics-bam collate -m 128M -o grouped.bam input.bam
rsomics-bam coverage -q 20 -Q 20 sample.bam
rsomics-bam cram-size -e sample.cram
rsomics-bam depad -T padded.fa -@ 4 -o unpadded.bam padded.bam
rsomics-bam depth -a -b targets.bed sample.bam
rsomics-bam fastq -o reads.fq.bgzf name-sorted.bam
rsomics-bam fixmate -m grouped.bam -o fixed.bam
rsomics-bam index -c -m 14 sample.bam
rsomics-bam import read1.fastq read2.fastq -@ 4 -o unmapped.bam
rsomics-bam idxstats sample.bam
rsomics-bam merge lane1.bam lane2.cram --reference reference.fa -o merged.bam
rsomics-bam markdup -r fixed-and-sorted.bam -o deduplicated.bam
rsomics-bam mpileup -f reference.fa -Q 20 input.bam
rsomics-bam reheader replacement.sam input.bam -o reheadered.bam
rsomics-bam reset --keep-tag RG,BC aligned.bam -o unmapped.bam
rsomics-bam sort -m 768M -o sorted.bam input.bam
rsomics-bam stats -r reference.fa -o sample.bamstat sample.bam
rsomics-bam to-bed --split-d -o blocks.bed alignments.bam
```

The commands accept SAM, BAM, and CRAM where their help declares those input
formats. `view` writes SAM or BAM; CRAM output is intentionally unavailable.
`cram-size` accepts CRAM 2.1, 3.0, and 3.1 from a file or standard input and
reports physical block sizes without decoding alignment records. Default,
verbose, and encoding-map output match samtools 1.24. Named output is
transactional; `--json` writes the compatibility report to a named file and
returns the typed report through the shared envelope.
Indexed region queries require a usable BAI, CSI, or CRAI. A non-zero thread
request for CRAM decoding fails explicitly because the current decoder does
not provide ordered parallel decoding. `index` uses up to four additional
workers when `-@` is omitted; pass `-@ 0` for one-thread indexing. Named index
outputs are committed only after the complete index has been built and parsed.
`checksum` accepts one or more SAM, BAM, no-reference CRAM, FASTA, or FASTQ
inputs, as well as standard input. It supports flag and tag selection, strand
normalization, order-aware checksums, position, CIGAR, mate and QC columns,
sanitization, tabular output,
and biobambam2-compatible reports. Merge mode validates native and
bamseqchksum versions, schemas, totals, and repeated rows before combining
them. Named reports are transactional; `--json` requires one and returns typed
groups through the shared envelope. Additional workers apply to named BAM
input. Reference-backed CRAM and HTSlib input-format option strings are not
exposed.
`collate` accepts SAM, BAM, and CRAM input and writes BAM with contiguous QNAME
groups. Its total record-memory budget and external merge fan-in are bounded;
the order between groups is intentionally unspecified. Fast-mode filtering and
early-pair buffering are not yet exposed.
`ampliconclip` consumes coordinate-ordered BAM and a three- or six-column
primer BED, then writes clipped BAM plus optional statistics, rejects, and
per-primer counts. `ampliconstats` consumes coordinate-ordered clipped BAM and
the six-column BED to emit samtools-compatible per-file and combined reports.
Both use the shared CLI and JSON envelope, reject aliased outputs, and commit
named outputs only after successful processing. SAM, CRAM, and plot generation
are not exposed for this workflow.
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
`coverage` emits the complete nine-column reference summary for SAM, BAM, and
CRAM, including filtered read counts, covered bases, mean depth, and mean base
and mapping qualities. It accepts multiple inputs, input lists, standard input,
and indexed regions. `bedcov` preserves BED row order while appending one
coverage total per input, with optional depth-threshold and read-count columns;
it accepts gzip-compressed BED and explicit BAI, CSI, or CRAI paths. `idxstats`
uses index metadata when available and otherwise validates a coordinate-sorted
scan. Named outputs are transactional, and machine summaries use the shared
`rsomics-help` JSON envelope. Histogram-only coverage modes are not exposed.
`calmd` recalculates MD and NM for mapped records from an indexed reference.
It accepts SAM, BAM, CRAM, and standard input, writes SAM or BAM, preserves
mapped records without query sequence with an explicit warning, and can rewrite
reference-matching query bases to `=`. Existing correct tags retain their
position; corrected tags follow samtools replacement order. Named outputs are
transactional. BAQ, mapping-quality adjustment, and CRAM output are not exposed.
`depad` removes padded-reference columns from alignment coordinates, CIGARs,
mate positions, and reference lengths. It accepts SAM, BAM, no-reference CRAM,
and standard input and writes SAM or BAM. A padded FASTA supplied with `-T` is
indexed in memory without creating a sidecar; without it, embedded reference
records drive projection and header lengths remain padded. Named output is
transactional. Reference-backed CRAM input and CRAM output are not exposed.
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
`stats` emits the complete samtools 1.24 alignment-statistics report for SAM,
BAM, CRAM, and standard input. It supports references, target intervals,
indexed regions and custom indexes, flag and read-group filters, mate-overlap
removal, tagged split reports, reference statistics, and threaded BAM or CRAM
decompression. Coverage and cycle state are streamed or sparse, split report
cardinality is bounded, and named main and split outputs commit as one group.
`--json` requires a named compatibility report and returns the typed report
through the shared `rsomics-help` envelope.
`reset` accepts SAM, BAM, CRAM, and standard input and writes SAM, BAM, or
CRAM. It drops secondary and supplementary alignments, restores reverse reads
to read orientation, clears alignment and mate coordinates, removes the
samtools default alignment tags, and supports explicit remove, keep, read-group,
program-chain, and duplicate-flag policies. Named output is transactional;
CRAM is fully decoded before commit. HTSlib input and output format-option
strings are not exposed.
`to-bed` accepts SAM, BAM, CRAM, ordinary gzip-compressed SAM, and standard
input. It emits BED6 by default, split BED6 at CIGAR `N` or `D` boundaries,
blocked BED12, or BEDPE from adjacent name-grouped pairs. Scores come from
MAPQ, `NM`, or another integer tag; CIGAR text can be appended to unsplit
BED6. Named outputs are transactional, reference-backed CRAM uses `-T`, and
machine summaries stay separate through the shared JSON envelope.

Stable behavior is tested against samtools 1.24 across SAM, BAM, and CRAM.
Coverage summaries, BED coverage totals, custom indexes, filtering, depth
limits, regions, and CRAM quality semantics are byte-matched to samtools 1.24.
Read-group editing is field-matched across new, existing, implicit, overwrite,
orphan, SAM, BAM, CRAM, and non-string-tag cases.
FASTQ import is field-matched for positional and explicit input modes, read
groups, order tags, compressed input, standard input, SAM, and BAM output.
FASTA/FASTQ extraction also has bytewise stdin and historical-fixture checks.
The full `stats` section order and stable body are byte-matched across the
upstream CIGAR, target, overlap, barcode, split, reference-statistics, indexed
region, and CRAM fixtures.
Named alignment and pileup outputs are committed only after successful
processing. See [PERFORMANCE.md](PERFORMANCE.md) for the representative BAM
throughput, memory, and decoded-output gate.

License: MIT OR Apache-2.0.
