# Performance

The release gate compares `rsomics-bam view` with `samtools view` while both
read and write BAM with four additional I/O threads. Every timing round is
accepted only after samtools verifies both files and their decoded headers and
records match.

## Release benchmark

The 2026-08-02 benchmark used `rsomics-bam` revision
`2e3781c5eaa54ec565289445da1c9a6793827d15` and samtools/HTSlib 1.24 on an
Ubuntu 22.04 machine with two Intel Xeon Gold 6238R processors. Both commands
were restricted to CPUs 40–44. One warm-up preceded five alternating timed
rounds.

```text
rsomics-bam view -b -@ 4 --no-pg -o ours.bam input.bam
samtools view -b -@ 4 --no-PG -o samtools.bam input.bam
```

| Tool | Median wall time | Mean wall time | Median CPU time | Median peak RSS |
|---|---:|---:|---:|---:|
| `rsomics-bam view` | 2.39 s | 2.40 s | 7.12 s | 4,480 KiB |
| `samtools view` | 4.14 s | 4.36 s | 9.02 s | 10,752 KiB |

| Round | `rsomics-bam` | samtools |
|---:|---:|---:|
| 1 | 2.53 s | 3.81 s |
| 2 | 2.31 s | 3.85 s |
| 3 | 2.29 s | 4.53 s |
| 4 | 2.39 s | 5.47 s |
| 5 | 2.50 s | 4.14 s |

On this fixture, `rsomics-bam view` was 1.73 times as fast by median wall time,
used 21.1% less median CPU time, and used 58.3% less median peak RSS. These
ratios describe this workload and machine rather than all alignment data.

The coordinate-sorted fixture contained 3,000,000 records and was 188,400,612
bytes. Its SHA-256 was
`fcd794dd0f865d5c80ede9bdb6e8048c7994ef480c1a32c504824926c815420d`.
The decoded header SHA-256 was
`55a9ddc7b5c667a2d61cd96b4b93e5283376300fbcab82acf7a92e0335cbff21`,
and the decoded record-stream SHA-256 was
`63330b1405bb2544a038c44567c92681d94fa98c55bb313ec7f71cb9d248e885`
for the warm-up and all five rounds. BAM byte streams differed because the
encoders used different BGZF block layouts.

The release-mode `rsomics-bam` binary was built with rustc 1.95.0 and had
SHA-256
`47bd36dedf73ee5793ac119eeac169ed896f2a997f44a8dca5569ea4ddeb69af`.
The samtools 1.24 binary had SHA-256
`0ea9344e09afa7dcde414adc3e5dae2a139a49bbe519dc0d05c3bac034de85bb`.

## Depth benchmark

The 2026-08-03 depth gate used revision
`ebf7f9606db853d1a191f4110519981fbd2885d2` and samtools/HTSlib 1.24 on an
Apple M2 Mac mini with four performance cores, four efficiency cores, and 8
GiB of memory. Twenty timed rounds alternated command order after one warm-up.
Both commands used their default single-threaded behavior and formatted the
full depth stream to `/dev/null`. A separate untimed pass wrote named files and
compared the complete streams byte for byte.

```text
rsomics-bam depth input.bam > /dev/null
samtools depth input.bam > /dev/null
```

| Tool | Mean wall time | Mean user time | Mean system time | Mean peak RSS |
|---|---:|---:|---:|---:|
| `rsomics-bam depth` | 0.4400 s | 0.3735 s | 0.0200 s | 3,978,854 bytes |
| `samtools depth` | 0.4900 s | 0.4020 s | 0.0315 s | 6,572,442 bytes |

The paired mean wall-time difference was -0.0500 seconds with a 0.0801-second
sample standard deviation. `rsomics-bam` won 16 of 20 pairs; the paired
t-statistic was -2.790. On this fixture it reduced mean wall time by 10.20%,
user time by 7.09%, system time by 36.51%, and peak RSS by 39.46%.

The coordinate-sorted 5,000,000-base fixture was 35 MiB at approximately 30x
coverage. Its SHA-256 was
`33b6780ec3758a8ccde746935366dec441e89aaafb5b0253a19cfa1af350282c`.
The complete correctness pass matched samtools byte for byte; the output
SHA-256 was
`f9bbc936dab1d5e7ef17c834a505187edefcfb854935be7021614fd51aaf2a69`.
The complete 20-pair timing ledger has SHA-256
`ca2cad4d5d55aa0e492f1916c1d7c3ccf646308c63323b6615b99fa5afaec193`;
its generated summary has SHA-256
`dccf26d1c6b115e7a86811ef0a7114ba998ade2a94a76f8c99936a5cca25ea45`.
The measured `rsomics-bam` and samtools binaries had SHA-256 values
`00fbb30027eaed8817bc6e2d47e4cab69523e1e58dd6f5dca67640b9164054c4`
and `c265b440b09c4b21d1f25a65963cf907b0d9f9d18caa9382c31104158f89d027`,
respectively.
These measurements apply to the default BAM depth path on this fixture and do
not establish the same advantage for quality filtering, region queries,
multiple inputs, SAM, or CRAM.

## Index benchmark

The 2026-08-03 index gate used revision
`dce21e7341cfc0a39ac66f9d82a027bda2c4cbc2` and samtools/HTSlib 1.24 on an
Apple M2 Mac mini with 8 GiB of memory. Twenty timed pairs alternated command
order after one warm-up. `rsomics-bam index` selected its default four
additional workers; samtools used its default one-thread indexing path.

```text
rsomics-bam index -o ours.bai input.bam
samtools index -o samtools.bai input.bam
```

| Tool | Mean wall time | Mean user time | Mean system time | Mean peak RSS |
|---|---:|---:|---:|---:|
| `rsomics-bam index` | 0.4280 s | 0.5810 s | 0.1780 s | 5,968,691 bytes |
| `samtools index` | 0.7675 s | 0.4570 s | 0.0580 s | 6,887,834 bytes |

The default `rsomics-bam` path was 1.79 times as fast by mean wall time and
used 13.34% less mean peak RSS. It won 19 of 20 timed pairs; the paired mean
wall-time difference was -0.3395 seconds, its sample standard deviation was
0.3517 seconds, and the paired t-statistic was -4.317. Automatic parallelism
increased mean CPU time by 47.38%; pass `-@ 0` when one-thread behavior is
preferred over the default latency target.

The coordinate-sorted WGBS fixture contained 4,000,000 records and was
77,438,045 bytes. Its SHA-256 was
`fe4f1977a9eb9352faafec62f5ab44e77f93757fd5557917d83b4558bc5530d6`.
Every round produced byte-identical BAI files, and samtools `idxstats -X`
accepted both files with identical output. Their SHA-256 was
`9c904c043df9e2252bcb527a571ac46d8947882e6a3e4c53abc0fe6e01c0bb7f`.

The timing ledger, generated summary, environment record, and JSON index
summary had SHA-256 values
`c8d618a06fac44013f944f03007352b0e66565ba050b9e127504b3fd99758d4c`,
`8dbbc698c710d1c1f9bd30db045143b48065bb5247b77903efed09d8bc300782`,
`47ac9fad280c29bfd7d1927954806c6ce7228dc9ba3e86738b68f843570aa618`,
and `38204d348bcee017fd714ef77760f62b950de290a8be66905e17966e6455f677`,
respectively. The measured `rsomics-bam` and samtools binaries had SHA-256
values `e96947024854724583af0e688d936b23b0918c275d9dea3a8a8789b5820dab37`
and `c265b440b09c4b21d1f25a65963cf907b0d9f9d18caa9382c31104158f89d027`.
These results establish the default BAM/BAI gate on this fixture; they do not
claim the same advantage for explicit equal thread counts, CSI, BGZF SAM, or
CRAM.

## Sort benchmark

The 2026-08-03 sort gate used revision
`8433aea711d59d9977f9c2734f5396c09bd6de32` and samtools/HTSlib 1.24 on an
Apple M2 Mac mini with 8 GiB of memory. Twenty timed pairs alternated command
order after one warm-up. Both commands received a 768 MiB memory setting.
`rsomics-bam sort` selected its default four additional workers; samtools used
its default one-thread path.

```text
rsomics-bam sort --no-PG -m 768M -o ours.bam input.bam
samtools sort --no-PG -m 768M -o samtools.bam input.bam
```

| Tool | Mean wall time | Mean user time | Mean system time | Mean peak RSS |
|---|---:|---:|---:|---:|
| `rsomics-bam sort` | 6.7300 s | 9.0265 s | 2.8475 s | 848,793,600 bytes |
| `samtools sort` | 11.5080 s | 7.1435 s | 1.8650 s | 960,046,694 bytes |

The default `rsomics-bam` path was 1.71 times as fast by mean wall time and
used 11.59% less mean peak RSS. It won all 20 timed pairs; the paired mean
wall-time difference was -4.7780 seconds, its sample standard deviation was
1.0959 seconds, and the paired t-statistic was -19.499. Automatic parallelism
increased mean CPU time by 31.81%; pass `-@ 0` when one-thread behavior is
preferred over the default latency target.

The query-name-sorted WGBS fixture contained 4,000,000 records and was
77,438,055 bytes. Its SHA-256 was
`057e5ff8c46f5870d7c925d28f429a3bb61745a2448c0f0c948e110d131e452e`.
It was derived from the index benchmark fixture with samtools 1.24 natural
query-name sorting and no added program record. The rsomics warm-up produced
two temporary runs and one merge pass. Every warm-up and timed pair had an
identical complete header and order-sensitive `samtools checksum -a -O`
report; the final checksum artifact SHA-256 was
`44eeec9a781436072463cc707feb982e727278ce5b33125742c86483947550a8`.

The timing ledger, generated summary, environment record, and JSON sort
summary had SHA-256 values
`077fcc9b4183162834c837ea38b136f92a4c690084760dc7d35d6e20927f4abb`,
`922f3fb081d964ef6613d4c080c40b981bd0c0a964756ed3ddb22aa96877945a`,
`fa7926c350d86604f45365d78a0d07be93cb299fc0c286042737d833f5d9f761`,
and `046cc7ab2eddeff0fc98e4cd0d090a82bf66d067653c91504b35ba48e27ba39b`.
The measured `rsomics-bam` and samtools binaries had SHA-256 values
`3ce9d6a56d0bc7511416a5e3f2ac0f7147ae5c81eeca179249c5ad243b3934af`
and `c265b440b09c4b21d1f25a65963cf907b0d9f9d18caa9382c31104158f89d027`.
These results establish the default coordinate-sort latency gate on this
fixture; they do not claim the same advantage for explicit equal worker
counts, other sort orders, SAM or CRAM input, or different compression and
memory settings.

## Reproduction

`benchmarks/view-vs-samtools.sh` records the machine, tool versions, binary and
input checksums, GNU time results, and decoded-output checksums. It fails on a
malformed output or any semantic disagreement:

```sh
RSOMICS_COMMIT=2e3781c benchmarks/view-vs-samtools.sh \
  target/release/rsomics-bam \
  /path/to/samtools \
  /path/to/input.bam \
  /path/to/results \
  5 4 40-44
```

`benchmarks/depth-vs-samtools-macos.sh` performs the corresponding bytewise
depth comparison, records macOS resource usage, and alternates command order:

```sh
RSOMICS_COMMIT=ebf7f96 benchmarks/depth-vs-samtools-macos.sh \
  target/release/rsomics-bam \
  /path/to/samtools \
  /path/to/input.bam \
  /path/to/results \
  20
```

`benchmarks/index-vs-samtools-macos.sh` records the default index comparison,
alternates command order, and rejects any BAI or `idxstats` disagreement:

```sh
RSOMICS_COMMIT=dce21e7 benchmarks/index-vs-samtools-macos.sh \
  target/release/rsomics-bam \
  /path/to/samtools \
  /path/to/input.bam \
  /path/to/results \
  20
```

`benchmarks/sort-vs-samtools-macos.sh` records the default coordinate-sort
comparison, alternates command order, and rejects any header or full-record
checksum disagreement:

```sh
RSOMICS_COMMIT=8433aea benchmarks/sort-vs-samtools-macos.sh \
  target/release/rsomics-bam \
  /path/to/samtools \
  /path/to/queryname-sorted.bam \
  /path/to/results \
  20 default
```

Pass `equal-workers` as the final argument to compare four additional workers
for both tools while dividing samtools' per-thread memory setting to keep the
requested sort-record budget near 768 MiB.
