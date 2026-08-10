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

## Merge benchmark

The 2026-08-03 merge gate used feature revision
`83b73a0c727436d9e69924f365af66d689c4f3aa` and samtools/HTSlib 1.24 on an
Apple M2 Mac mini with 8 GiB of memory. It merged the same natural
query-name-sorted 4,000,000-record BAM twice into an 8,000,000-record BAM.
Twelve default-mode pairs alternated command order after one warm-up.

```text
rsomics-bam merge --no-PG -c -p -n input.bam input.bam -o ours.bam
samtools merge --no-PG -c -p -n -@ 0 -f -o samtools.bam input.bam input.bam
```

| Tool | Mean wall time | Mean user time | Mean system time | Mean peak RSS |
|---|---:|---:|---:|---:|
| `rsomics-bam merge` | 3.0158 s | 10.6917 s | 0.8658 s | 20,959,232 bytes |
| `samtools merge` | 9.3550 s | 8.0958 s | 0.5117 s | 13,873,152 bytes |

The automatic `rsomics-bam` path was 3.10 times as fast by mean wall time.
Its mean peak RSS was 1.51 times samtools' and its mean CPU time was 1.34
times samtools'. The paired samtools-minus-rsomics wall-time difference was
6.3392 seconds, with sample standard deviation 0.2091 seconds and paired
t-statistic 105.00.

An additional eight-pair gate gave both tools four additional workers. Mean
wall time was 3.3613 seconds for `rsomics-bam` and 3.3288 seconds for samtools;
the paired difference was -0.0325 seconds with a -0.157 t-statistic. Mean peak
RSS was 20,170,752 versus 17,508,352 bytes. This does not establish an
equal-worker throughput advantage.

Every warm-up and timed pair passed samtools quickcheck and produced identical
complete headers and order-sensitive `samtools checksum -a -O -T` reports.
The header and record-checksum artifact SHA-256 values were
`3337edeca3b276da02efcc727b3892ed201a54664a2491ea51971d9802e7e198` and
`4661c6ad4c3488583659a7e19cf18589f5a8ab23f2614177589679b0f24d8491`.
The measured rsomics and samtools binaries had SHA-256 values
`7cfd2a587d8fd2ce6e0ba3e18d50aa492c8513c7c94724640ec110d855ed9381` and
`c265b440b09c4b21d1f25a65963cf907b0d9f9d18caa9382c31104158f89d027`.
The default environment, timing, and summary files had SHA-256 values
`68d31dde9866c82faf15a24d425d3d838a298a129ded8f9f30bb8e97007ffa85`,
`35af7a6edff77ebb5539b1ed42cb01db522d0006d8c3802eaf897b4a18c2b89e`, and
`9f89349184dee6d0728ef8371251a84287f6e5fd423457ff502206123711fe3c`.
The equal-worker files had SHA-256 values
`1971065ebfb2c830c19a67cf9cd6db0b49737f15c2585b893b991d548c308495`,
`59023213de23264404a3f2f6e6b36d2d8a60d957127ef16ef8d14da1dc40b540`, and
`1b658de5b94e6d1ea1f8cde7c1a79b8ce691a92893dd0a61e69b62f8abef20ef`.
These results establish the natural query-name BAM merge gate only; they do
not claim the same performance for coordinate or template order, disjoint
inputs, SAM or CRAM, or materially different input counts.

## Collate benchmark

The 2026-08-03 collate gate used feature revision
`24095b8650c26cde05cbc6684c076bc28adc1ea5` and samtools/HTSlib 1.24 on an
Apple M2 Mac mini with 8 GiB of memory. The coordinate-sorted 4,000,000-record
BAM was collated with a 128 MiB rsomics record budget. Twelve default-mode
pairs alternated command order after one warm-up.

```text
rsomics-bam collate --no-PG -m 128M input.bam -o ours.bam
samtools collate --no-PG -@ 0 -o samtools.bam input.bam
```

| Tool | Mean wall time | Mean CPU time | Mean peak RSS |
|---|---:|---:|---:|
| `rsomics-bam collate` | 8.3400 s | 15.5483 s | 155,277,995 bytes |
| `samtools collate` | 13.5742 s | 7.0508 s | 40,013,824 bytes |

`rsomics-bam` won all 12 pairs and was 1.63 times as fast by mean wall time.
The paired samtools-minus-rsomics difference was 5.2342 seconds, with a
1.8032-second sample standard deviation and paired t-statistic 10.055. The
throughput advantage cost 2.21 times samtools' CPU time and 3.88 times its
peak RSS. The rsomics warm-up created nine temporary runs and one merge pass.

An eight-pair equal-worker gate passed four additional workers to both tools.
Mean wall time was 9.2463 seconds for `rsomics-bam` and 11.6563 seconds for
samtools, a 1.26 times advantage; rsomics won all eight pairs. The paired
difference was 2.4100 seconds with a 1.0738-second standard deviation and
paired t-statistic 6.348. Mean RSS was 154,777,600 versus 62,615,552 bytes,
and rsomics used 1.82 times samtools' CPU time.

Every warm-up and timed pair passed samtools quickcheck, had byte-identical
complete headers, and produced the same order-independent multiset fingerprint
over all four million complete SAM records after canonical auxiliary-field
ordering. The header and fingerprint artifact SHA-256 values were
`4f5236583648f1c96db66603dd2efdd062dc73ea95425e602121795d1705dfc8` and
`a37f38f50018235bd8b85810f0721c73059d660084e879d8bff40ba4d824dd13`.
The measured rsomics and samtools binary SHA-256 values were
`2ddaffe39cfcf9f5cfea7ef191de52286a830fdd89b81e00cd97fa0338988527` and
`c265b440b09c4b21d1f25a65963cf907b0d9f9d18caa9382c31104158f89d027`.

The input was 77,438,045 bytes with SHA-256
`fe4f1977a9eb9352faafec62f5ab44e77f93757fd5557917d83b4558bc5530d6`.
The default environment, timing, summary, and JSON files had SHA-256 values
`c2facd8c65a9019cea02122fd34253a15d04beb9b8222d305c0fd52b56e5f035`,
`ea565b4be15d5a6948b67c872e4f4ecaf99d42a0963946f91cdbfc7ba5b0f56d`,
`0bdcba7be536ab7185cc66cfbdf2c4865c97154b7953704f2604fc7507a697e1`,
and `120fb5a36ebd0a8cebf9d91de08f1b04456d1a7eefe27795e5de950d12a501a5`.
The equal-worker files had SHA-256 values
`191e3ee70c3edf1f92893953d264e16a3ef7ff8bda5cbe0404bd43878a118cdf`,
`36014358fe15b436eea5be11572853b25b41faac5ae478ef08d8e80f675b6fb7`,
`f6817ca9304bcb5943b0703d4bbd6a5651980b119b4a434c03bada44c8e5470c`,
and `5952e9281ffe46cd97be73ca0131d3a0f0e75082500bd83f1e479a1a1f79f3ea`.
These results establish a BAM wall-time advantage for standard collation on
this machine and fixture. They establish no CPU or memory advantage and do not
cover fast mode, SAM or CRAM input, or materially different name distributions.

## Fixmate benchmark

The 2026-08-10 fixmate gate used feature revision
`a8a684ba57c613d78addfee466cb7eca16937d71` and samtools/HTSlib 1.24 on an
Apple M2 Mac mini with 8 GiB of memory. The input contained 4,000,000 records
in 2,000,000 consecutive paired templates. Both tools calculated mate scores,
disabled program records, and ran with samtools' sanitizer disabled. Twelve
timed pairs alternated command order after one warm-up.

```text
rsomics-bam fixmate --no-pg -m input.name.bam -o ours.bam
samtools fixmate -z off --no-PG -m -@ 0 input.name.bam samtools.bam
```

| Tool | Mean wall time | Mean CPU time | Mean peak RSS |
|---|---:|---:|---:|
| `rsomics-bam fixmate` | 2.2708 s | 8.2300 s | 7,002,795 bytes |
| `samtools fixmate` | 7.0083 s | 6.4442 s | 7,289,515 bytes |

The default `rsomics-bam` path selected four additional compression workers,
won all 12 pairs, and was 3.09 times as fast by mean wall time. The paired
rsomics-minus-samtools difference was -4.7375 seconds, with sample standard
deviation 0.2796 seconds and paired t-statistic -58.69. Its mean peak RSS was
3.93% lower, while its mean CPU time was 27.71% higher.

An additional 12-pair gate gave both tools four additional workers. Mean wall
time was 2.3300 seconds for `rsomics-bam` and 2.5050 seconds for samtools;
`rsomics-bam` won 7 of 12 pairs. The paired difference was -0.1750 seconds
with a -1.40 t-statistic, so this gate does not establish a stable
equal-worker throughput advantage. Mean peak RSS was 7,028,736 versus
13,100,373 bytes, a 46.35% reduction.

The 82,830,601-byte fixture had SHA-256
`8582710cca9390ba70c3c219f41dc6b0dd8cd66fc5303f3a26d46f4a9ac9b9b5`.
Every warm-up and timed pair passed samtools quickcheck and produced identical
complete headers and order-sensitive record streams after normalizing only
auxiliary-tag order. The canonical record stream contained 4,000,000 records
and had SHA-256
`9653a33166b90b5121dcc56281b96e3f574d5b6d3e0e9cbed6e3673b9072dcd5`.

The measured rsomics and samtools binaries had SHA-256 values
`0f4e0015c7c39438f319765c19bcda2802bda03fa0aee7c440bd396b188d5cd3`
and `c265b440b09c4b21d1f25a65963cf907b0d9f9d18caa9382c31104158f89d027`.
The default environment, timing, summary, and JSON files had SHA-256 values
`321f45332f6083b29a098e13dd91e7d2a78cd7220e60540d33dfa95d86d6cab0`,
`c9ee0c944273963be00ea0e4d8794abf566c3ffc2fcb74f14028092b2394c9f2`,
`c71b179cdefa83c47bff1aa8d0fe792da7f99e55d8edce3cdfc53ee555c3ccc7`,
and `3e42140f1d21d3b2efc84b8d507d9cfe94242528b91df8b61957c3a69f6e21c7`.
The equal-worker files had SHA-256 values
`d7643e54d50bba937ae8fd72142d94d65cd35b7e48030424072139bb04367f3e`,
`0aad39d24c078d930e50cd98cc53486cafdf79fc71501d92968119ddfef1ce71`,
`285792ba013ec666ea716938ac74449fe2230707a45310e8b4e8f1e2f483161e`,
and `9ea00249d072debfbb03542feb46e28877de5516c78430f6e12424eb6fc3a800`.
These results cover name-grouped BAM input, BAM output, mate-score calculation,
and the measured thread settings only. They do not establish the same
performance for SAM, CRAM, standard input, other option combinations, or
different template-size distributions.

## Markdup benchmark

The 2026-08-10 markdup gate used feature revision
`5c7dc5603dabdeed212ef270474270f6abdbb47d` and samtools/HTSlib 1.24 on an
Apple M2 Mac mini with 8 GiB of memory and macOS 26.6.1. The 92,673,552-byte
coordinate-sorted BAM contained 4,000,000 records and had SHA-256
`bc2257da48b4c06da643edafbec1a383e946b7d1a0c0dd09dc21edc48dc2ef2d`.
It was generated through the release product path with `fixmate -m` followed
by coordinate sorting. Twelve timed pairs alternated command order after one
warm-up.

```text
rsomics-bam markdup --no-pg input.coordinate.bam -o ours.bam
samtools markdup --no-PG -@ 0 input.coordinate.bam samtools.bam
```

| Tool | Mean wall time | Mean CPU time | Mean peak RSS |
|---|---:|---:|---:|
| `rsomics-bam markdup` | 2.3475 s | 8.3583 s | 6,980,949 bytes |
| `samtools markdup` | 7.6617 s | 7.2467 s | 7,587,157 bytes |

The default product path selected four additional compression workers, won all
12 pairs, and was 3.26 times as fast by mean wall time. The paired
rsomics-minus-samtools difference was -5.3142 seconds, with sample standard
deviation 0.1308 seconds and paired t-statistic -140.70. Mean peak RSS was
7.99% lower, while mean CPU time was 15.34% higher.

An additional 12-pair gate gave both tools four additional workers. Mean wall
time was 2.5317 seconds for `rsomics-bam` and 2.7792 seconds for samtools;
`rsomics-bam` won 10 of 12 pairs. The paired difference was -0.2475 seconds
with sample standard deviation 0.3793 seconds and paired t-statistic -2.26.
Mean CPU time was 8.3092 versus 9.0942 seconds, an 8.63% reduction, and mean
peak RSS was 6,990,507 versus 13,660,160 bytes, a 48.82% reduction.

Every warm-up and timed pair passed samtools quickcheck and produced an
identical complete `samtools view -h --no-PG` stream. The exact SAM stream
fingerprint was
`6279ec79c152d1b2f6092b31021a32f8a62935615a0e2f3668c42e9a17011c99`.
The measured product binary had SHA-256
`beeb3ddfa7c789a0841af9594b941d6e20c2a08392abd7a6502c32a85f0ff60b`.

The default environment, timing, summary, and JSON files had SHA-256 values
`bedfb17f61486845a2e494116afd300ad3fe5cb8fe13b4d6b3c0bce08289ae62`,
`e54c8d8e80db0bcbb84cfd2ae3ae1d81c8de9e39b1ac1858417ef7364b9aded9`,
`42c7a2b418fc47c8b7ad3d9bc46e0172faf1888ca370d4ec9c651c40defda4b1`,
and `74df4b1de37bca65af449cc164f4128f29240817f05d74f0bcc17466c9f6af73`.
The equal-worker files had SHA-256 values
`3a261d01dc313f761df8472939e97dcb27d22ed469cbcf6f27ce7f4a129d0b82`,
`9fac2e4acef6586ff21b33881112498d572478ddb44630d98a13840c49ac2562`,
`ac64e18fc1616a848c995c7e61f09cafe5316bc8346d09f504f332a6f7de4be1`,
and `b902cbc0455780fac10ab191f8a8491cdce93a7581e8b3205098b873be939684`.
These results cover default-mode BAM input and BAM output on the measured
fixture and thread settings. They do not establish the same performance for
SAM, CRAM, standard input, removal or clearing, sequence mode, or different
duplicate distributions.

## FASTA/FASTQ extraction benchmark

The 2026-08-10 extraction gate used implementation revision
`d6cbf10707061c705eeabd85f7bfd40b30d9bf80` and samtools/HTSlib 1.24 on an
Apple M2 Mac mini with 8 GiB of memory and macOS 26.6.1. One warm-up preceded
10 timed pairs with alternating command order. Both tools used their default
single-threaded paths and formatted the complete output stream to `/dev/null`.

```text
rsomics-bam fasta input.queryname.bam > /dev/null
samtools fasta input.queryname.bam > /dev/null
rsomics-bam fastq input.queryname.bam > /dev/null
samtools fastq input.queryname.bam > /dev/null
```

| Format and tool | Mean wall time | Wall-time standard deviation | Mean peak RSS |
|---|---:|---:|---:|
| `rsomics-bam fasta` | 1.294 s | 0.032 s | 5,559,091 bytes |
| `samtools fasta` | 1.705 s | 0.091 s | 6,709,248 bytes |
| `rsomics-bam fastq` | 1.773 s | 0.071 s | 5,559,091 bytes |
| `samtools fastq` | 2.608 s | 0.071 s | 6,696,141 bytes |

FASTA reduced mean wall time by 24.11% and mean peak RSS by 17.14%; FASTQ
reduced them by 32.02% and 16.98%. `rsomics-bam` won all 10 pairs for each
format. The paired rsomics-minus-samtools wall-time differences were -0.411
seconds for FASTA and -0.835 seconds for FASTQ, with sample standard
deviations of 0.099 and 0.079 seconds.

The natural query-name-sorted BAM contained 4,000,000 records, was 92,673,124
bytes, and had SHA-256
`b949852de15f08a5e13d8c6d908b6d5801ef9f254eca58699aa353883cf88326`.
The complete FASTA and FASTQ streams matched samtools byte for byte, with
SHA-256 values
`7e13841514dd9137e08b0d9994afa5b4baafd0583bf7228740949bfcd6de80e3`
and
`a961649b1a0ee9beec27439fae61441e1d48694aa1486c7c74c1c64238ce4988`.
The measured product and samtools binaries had SHA-256 values
`8333eb7d4e1da96504af58182e19e7688841f2fff70ef6da18972a1eee1f9e3e`
and
`c265b440b09c4b21d1f25a65963cf907b0d9f9d18caa9382c31104158f89d027`.

The environment, complete timing ledger, and generated summary had SHA-256
values
`6c8d2c7d6d58bd7785d47f383e3521e725e6b1387d075163ccf14d82a6993b39`,
`abb92e745073024a2ac809e14e93ae5eaaffa3120122864cd5cb85abaf58ed45`,
and
`bd9833a1ce53dc15691df740bb11fedcc23318bd81d3cdf1cdc72412020e07d5`.
These performance claims cover uncompressed unified output from BAM on this
fixture. SAM, CRAM, stdin, filters, OQ replacement, and BGZF output have
compatibility coverage but no throughput claim here.

## Compressed file-operation benchmark

The 2026-08-10 file-operation gate used the final implementation tree from
revisions `1178474e4a4577f61882fed249f6f1909491ca40` and
`62023e111562d6e4bd47c26719b741499d0e1815` with samtools/HTSlib 1.24 on an
Apple M2 Mac mini with 8 GiB of memory and macOS 26.6.1. One warm-up preceded
12 alternating timed pairs. Every pair passed samtools quickcheck and produced
the same complete `samtools view -h --no-PG` stream.

The cat fixture contained four BAM shards totalling 90,038,862 bytes and
4,000,027 records. Their SHA-256 values were
`666433497d0653b907d4eb0da46f49bcc2f996eda2b3c497bf4bc2ed9f1b44e9`,
`2490918d1cc917b434dd5a69669c8daffc6a3dfe434f97ae591d638da2f9d670`,
`d3ea8e9f1f1363efbc5bd8cfc9f193209c41e521e5a1359295a3ce465b7225fd`,
and `f7bdf22d5f3f2964b297a789ddef4d8f90689307a588ee7dfea3033e0a8f8cd2`.

| Cat tool | Mean wall time | Mean system time | Mean peak RSS |
|---|---:|---:|---:|
| `rsomics-bam cat` | 1.0667 s | 0.0517 s | 5,503,659 bytes |
| `samtools cat` | 0.9225 s | 0.0525 s | 7,019,179 bytes |

Cat traded 15.63% mean wall time for a 21.59% mean peak-RSS reduction while
structurally validating every BGZF frame. It won 5 of 12 wall-time pairs; the
paired difference was 0.1442 seconds with a 0.4008-second sample standard
deviation and a 1.25 t-statistic. A 4 MiB output-buffer candidate was rejected
after a separate 24-pair gate: it averaged 1.0358 versus 0.9796 seconds, won
8 pairs, and used 9,644,032 versus 7,010,304 bytes mean peak RSS. The selected
path retains the bounded 2 MiB policy and makes a resource-use claim, not a
cat throughput claim.

The reheader fixture was a 92,673,552-byte BAM with 4,000,000 records and
SHA-256
`bc2257da48b4c06da643edafbec1a383e946b7d1a0c0dd09dc21edc48dc2ef2d`.
The replacement SAM header had SHA-256
`6b73526eb169145d79be3a71bb5e5cd190626771468d1dbc13fb25d94089df50`.

| Reheader tool | Mean wall time | Mean system time | Mean peak RSS |
|---|---:|---:|---:|
| `rsomics-bam reheader` | 0.6258 s | 0.0317 s | 5,496,832 bytes |
| `samtools reheader` | 0.7225 s | 0.0300 s | 6,980,949 bytes |

Reheader reduced mean wall time by 13.38% and mean peak RSS by 21.26%, winning
9 of 12 pairs. The paired difference was -0.0967 seconds with a 0.2940-second
sample standard deviation and a -1.14 t-statistic. The decoded cat and
reheader streams had SHA-256 values
`6c34e6be99f94979dce788b23849eec7be161179f13fec37fdce84a61c208389`
and `8bc5ca00000bfa575068de363b70bf3224cbebb3919519b58e2e01f410a19a15`.

The cat environment, timing, summary, and JSON files had SHA-256 values
`9cde1ffb7b0c87727b36c48dea8aa32cfa75a93216a51c6d56553d70cc10f502`,
`c8a28d57bc82806c146b30b4ad15cb1ea0417cde1de7d22bcc10e0c0d16b2c86`,
`c3491dc97b2c642f36ffd9b6dbd6c580232ee3a4cf93be05fefb19bb9f3f1a16`,
and `275357e146048340a597054de546b686f320f9c971b25d0bd13037dd3392f47b`.
The corresponding reheader values were
`d193ebff24fe892cdcbfaaa8714b0090dee8f5969efbd26bb85c193c3b13edab`,
`4ebc2c45a82fc824c2d9babc3c04b689a2baa11dc6956a62030682b8623c86b4`,
`85fd17d9980a3efd9cf707e295646c119269e0a66e3ea1c7213dcbc1db377b2a`,
and `bad792f1c535ca76ccf9d2655113c8a00941d835b48b9dd267282e08fcc316fd`.
These claims cover named BAM output and the measured fixtures only.

Release 0.12.0 was published from revision
`dfbc321d9dbec515b704634a104fed1680238e33` after exact-head CI
`31368507277` passed on native Linux and macOS for x86_64 and aarch64.
Publication workflow `31368957921` produced an unyanked 144,963-byte registry
archive with SHA-256
`df10fadae75e377e4a3c40244ad5bfd19d47010a874cc8063b81641fc8d1182b`.
A fresh registry install reported 0.12.0 and exposed both commands through the
shared help tree. Its two-shard cat smoke stream matched samtools 1.24 at
`0c9f5514885e469f4720858c25bef106a529091844c33c37ccc449ba45feb675`;
the four-million-record reheader smoke matched at
`8bc5ca00000bfa575068de363b70bf3224cbebb3919519b58e2e01f410a19a15`.
docs.rs serves the corresponding library documentation.

Release 0.13.0 was published from revision
`4829bbb3be06fddc0b13a2ede2cf72279044976e` after exact-head CI
`31376119277` passed on native Linux and macOS for x86_64 and aarch64. The
Linux x86_64 gate built samtools 1.24 and exercised the FASTA/FASTQ oracle for
SAM, BAM, CRAM, standard input, filtering, naming, and quality modes.
Publication workflow `31376581939` produced an unyanked 151,519-byte registry
archive with SHA-256
`97cc23593d5b92a7f3c49c19ab9b9c014e466ecc37e12e392b07ccaac27cf056`;
its VCS metadata identifies the exact release revision. A fresh registry
install reports 0.13.0 and exposes both commands through the shared help tree.
Installed FASTA and decompressed BGZF FASTQ smoke streams match samtools 1.24
at `a69efdbf4ebf740457c7df6e52112d1a56b63c388ad493c2a0f9ffbc0f8e61f8`
and `ca6ae968349466db34aa481149c0fc005689a3595cc3a3f8627139316754d733`.

## Import benchmark

The 2026-08-10 import gate used feature revision
`1df18368dd7c6e969af768bbf9c770f786ad0d1e` and samtools/HTSlib 1.24 on an
Apple M2 Mac mini with 8 GiB of memory. The paired WGSIM fixture contained
500,000 reads in each mate file. Single-input measurements used the first mate
file; paired measurements interleaved both files into one unmapped BAM. Both
tools used four additional compression workers and omitted program records.
Twelve timed pairs alternated command order after one warm-up.

```text
rsomics-bam import reads-1.fq --no-PG -@ 4 -o ours.bam
samtools import reads-1.fq --no-PG -@ 4 -o samtools.bam

rsomics-bam import -1 reads-1.fq -2 reads-2.fq --no-PG -@ 4 -o ours.bam
samtools import -1 reads-1.fq -2 reads-2.fq --no-PG -@ 4 -o samtools.bam
```

| Input mode | Tool | Mean wall time | Mean user time | Mean peak RSS |
|---|---|---:|---:|---:|
| Single | `rsomics-bam import` | 0.3150 s | 1.0292 s | 6,393,856 bytes |
| Single | `samtools import` | 0.5200 s | 1.0017 s | 11,057,835 bytes |
| Paired | `rsomics-bam import` | 0.6133 s | 1.9917 s | 6,703,787 bytes |
| Paired | `samtools import` | 0.8825 s | 1.9783 s | 11,182,080 bytes |

The rsomics single-input path reduced mean wall time by 39.42% and mean peak
RSS by 42.18%, winning all 12 pairs. Its paired mean difference was -0.2050
seconds with a 0.1354-second sample standard deviation and a -5.24 paired
t-statistic. The paired-input path reduced mean wall time by 30.50% and mean
peak RSS by 40.05%, winning all 12 pairs. Its paired mean difference was
-0.2692 seconds with a 0.1147-second standard deviation and a -8.13 paired
t-statistic. Mean user time was 2.75% higher for single input and 0.67% higher
for paired input.

The 172,708,525-byte FASTQ inputs had SHA-256 values
`664f264a7c06dba94c70c97d1a2d0a0c5ebb4fad1edc4f3d7c44bea5db651efa`
and
`13e5c4001ac58069c997761ba0ba22c813c8df39203691e587581e46dd02d4ff`.
Every warm-up and timed output decoded successfully. Stable headers matched at
`0917ebfe5bcf5d83582cb55c57ce146443dfee80a113b3885b4e454d336202c1`;
complete order-sensitive record streams matched at
`4457c6f0df0f719f7d5620aaa0ff0310532b562622691c2ca4cab5aded329262`
for single input and
`e52d26e8a33a0ac6f5f52ae66e2cf7927a8777a74159ef7885650f1cc3312351`
for paired input. The rsomics BAMs were 0.09% and 0.06% larger than the
corresponding samtools files.

The measured rsomics and samtools binaries had SHA-256 values
`b3e81cc1945cba86999d37839e44c57c16f100f04dfeb1caece7d03ddb1bfe25`
and `c265b440b09c4b21d1f25a65963cf907b0d9f9d18caa9382c31104158f89d027`.
The environment, timing, summary, and paired-statistic artifacts had SHA-256
values `a73a70cb4ea82ce6fb2e0a6bf5085f747410f0f0034b8d2bf3979f5c3c1fd585`,
`e4ba48c2f7ce14a1819f3fab02d759606663f0603ca6a6df25b068279b662b20`,
`e36c65c74f70cb6f4120744cd758c90e4c24b4051a2672f887f1f4ef7dd2061c`,
and `8e7c71c6499eacbfd64afa7cee86a8e04e07c783a87dd88df2b682bb2610bc38`.
These results establish the plain FASTQ to BAM hot path with four additional
workers on this fixture. They do not claim the same advantage for compressed
input, SAM output, auxiliary-tag extraction, other thread counts, or materially
different read lengths and entropy.

Release 0.14.0 was published from revision
`d54924462ad65d4a0781545e6511f7bc3a8becb8` after exact-head CI
`31383580026` passed on native Linux and macOS for x86_64 and aarch64. Its
Linux x86_64 gate rebuilt samtools 1.24 and ran the complete product oracle,
including the import matrix. Publication workflow `31384051315` produced an
unyanked 162,589-byte registry archive with SHA-256
`d092eb6d53b301d1e9be0d9e17671502f66d216d1e5b3eb63e4f311da442dcef`;
its VCS metadata identifies the exact release revision.

A fresh registry install reports 0.14.0 and exposes `import` through the shared
help tree. Installed single-end SAM, gzip standard input, and paired BAM output
match samtools 1.24 after removing program and comment header records. The
stable single-end SAM, paired header, and complete paired record stream have
SHA-256 values `b3a74cf2ac8815237013ca55b9d2c3c466d5f651fb1abd92c749d02a714d7e37`,
`ea48c78110bc71ed96dd0346c56b395e341e562efd046edf9393bc12823267df`,
and `2a07f119149d2b36ca6415b0735a35bb0cb1ff9fbbf34506d3933b85e8f70f64`.
The installed binary also rejects non-IUPAC FASTQ input with a nonzero exit.

## Read-group editing benchmark

The 2026-08-10 `addreplacerg` gate used implementation revision
`033a7fa6c274` and benchmark revision `07bd58c7fd37` with samtools/HTSlib 1.24
on an Apple M2 Mac mini with 8 GiB of memory. The 99,545,915-byte BAM contained
4,000,260 records: 2,000,130 with `RG:Z:old` and 2,000,130 without an `RG`
field. Both tools used four additional workers and omitted program records.
Twelve timed pairs alternated command order after one warm-up for each mode.

```text
rsomics-bam addreplacerg -r 'ID:new\tSM:after' -m overwrite_all --no-PG -@ 4 -o ours.bam input.bam
samtools addreplacerg -r 'ID:new\tSM:after' -m overwrite_all --no-PG -@ 4 -o samtools.bam input.bam

rsomics-bam addreplacerg -r 'ID:new\tSM:after' -m orphan_only --no-PG -@ 4 -o ours.bam input.bam
samtools addreplacerg -r 'ID:new\tSM:after' -m orphan_only --no-PG -@ 4 -o samtools.bam input.bam
```

| Mode | Tool | Mean wall time | Mean user time | Mean peak RSS |
|---|---|---:|---:|---:|
| Overwrite all | `rsomics-bam addreplacerg` | 1.8125 s | 6.9592 s | 7,012,352 bytes |
| Overwrite all | `samtools addreplacerg` | 2.4600 s | 7.4167 s | 12,626,603 bytes |
| Orphan only | `rsomics-bam addreplacerg` | 1.7692 s | 6.9275 s | 6,985,045 bytes |
| Orphan only | `samtools addreplacerg` | 2.4825 s | 7.3633 s | 12,749,483 bytes |

The overwrite path reduced mean wall time by 26.32% and mean peak RSS by
44.46%, winning all 12 pairs. Its paired mean difference was -0.6475 seconds
with a 0.1567-second sample standard deviation and a -14.31 paired
t-statistic. The orphan path reduced mean wall time by 28.73% and mean peak
RSS by 45.21%, also winning all 12 pairs. Its paired mean difference was
-0.7133 seconds with a 0.2596-second standard deviation and a -9.52 paired
t-statistic.

Every warm-up and timed output decoded completely. Within-header `@SQ` and
`@RG` order and all field values matched. Cross-type placement of `@CO` was
normalized because noodles serializes header record types canonically while
samtools appends a new `@RG` after existing comments. Complete order-sensitive
record streams matched at
`8d5dfae4ee1e5db4cbb9be43ea9a6c73ccaff5fccbc5c3648e2c4c2d74bffb56`
for overwrite mode and
`0bab95e6775fec54b24fd81ca9eac4f058419bddc20c6243148a2561704e1008`
for orphan mode. The normalized header hashes were
`ad51c33f58087415558e256885d678787f69faa7bec5e8864fda5b058e3cd34c`
and `e36570e85d1f6ec14a052517237c0a9b410f7e25dbf11cb10f22b74481772c13`.

The input SHA-256 was
`9f82e1faae07d53bf916689828146c6923714d08a29078df726d00284363b1b3`.
The measured rsomics and samtools binaries had SHA-256 values
`f61522baf7b7cffdb2527dcbcd348308c046decc3a56141e41494797d86e9bde`
and `c265b440b09c4b21d1f25a65963cf907b0d9f9d18caa9382c31104158f89d027`.
The environment, timing, summary, and paired-statistic artifacts had SHA-256
values `7fc0b00caadfdb45d92c6a5ede740925701320bd8769f56656b587ea84f3b7d9`,
`82ddfc4d2b5ff7e7684c7efa4166cd693db86c74152aeb7cc3a3b7299212e96c`,
`03fc68b1f42b49677de717d5d312d85f04aad1b2211fa30bbdf0f89191511a6b`,
and `7de316099404a55150ea709e580a4eb3bedfa4a42b6da26ae1694a047f8cf381`.
These claims cover the measured mixed-tag BAM-to-BAM path with four additional
workers. They do not claim the same advantage for SAM or CRAM input, SAM
output, other thread counts, or materially different auxiliary-field layouts.

Release 0.15.0 was published from revision
`fe2beb388a7565ce064ed430af8b9476b821ced9` after exact-head CI
`31387911685` passed on native Linux and macOS for x86_64 and aarch64.
Publication workflow `31388331846` produced an unyanked 171,070-byte registry
archive with SHA-256
`79dec6d6cf7deff0a27443539974bec188fba213c7d0e9485059a94ddef61527`;
its VCS metadata identifies the exact release revision.

A fresh registry install reports 0.15.0 and exposes `addreplacerg` through the
shared help tree. Its binary SHA-256 is
`d79af698eda372b22a42f5828ae3c6e5ab8ef118f8c8ab4fb053b4104647ec37`.
Installed overwrite, orphan-only, and implicit-first-read-group paths matched
samtools 1.24 for complete record streams and normalized headers. Their record
stream SHA-256 values were
`2ad1b7c463c40d7f09ae8a4176bfecb0ff79be58c181249f690c7e11d48ee103`,
`f31c577c0f73294717a035faa1a1356c1f11fbd88ccd21f04f82099519bad31d`,
and `011ceaae1ba5c9104774ad64d6429f2e4267a8bc1c6f9bb0eb23d97f31e8e125`.
Both binaries rejected a conflicting header ID without `-w`, and the installed
binary emitted the shared JSON envelope with the expected three-record
summary.

## Coverage summary benchmark

The 2026-08-10 gate used revision
`981d4e91c623684d1aa123ca3f909ee8110541cb` and samtools/HTSlib 1.24 on an
Apple M2 Mac mini with 8 GiB of memory. Commands used their default
single-threaded behavior. Timed pairs alternated command order, and every
operation first produced a complete byte-identical output file.

The coordinate-sorted BAM contained 4,000,000 records and was 92,673,552
bytes. Its SHA-256 was
`bc2257da48b4c06da643edafbec1a383e946b7d1a0c0dd09dc21edc48dc2ef2d`.
The explicit BAI SHA-256 was
`d207836008a4cb7f75384ec2e357d0eecc75eb1d7876509a190de2082c958385`.

### Coverage

Ten full-reference pairs compared the complete nine-column table.

| Tool | Mean wall time | Mean user time | Mean system time | Mean peak RSS |
|---|---:|---:|---:|---:|
| `rsomics-bam coverage` | 3.155 s | 2.956 s | 0.074 s | 4,207,411 bytes |
| `samtools coverage` | 3.469 s | 3.282 s | 0.080 s | 6,565,069 bytes |

`rsomics-bam` won 9 of 10 pairs, was 1.10 times as fast by mean wall time,
and used 35.9% less mean peak RSS. The paired rsomics-minus-samtools
difference was -0.3140 seconds with a 0.1255-second sample standard deviation
and a -7.915 t-statistic. The identical output SHA-256 was
`eacdc4e84b06f55ffe38aa2986e910b1ecd8c8c525d06243e6e3be250b68d3f5`.

### BED coverage

The dense input contained 50,001 adjacent BED rows covering the reference;
its SHA-256 was
`0a285ac162bbb3e2077cdc0ba090a77cd32f97ebd4a0e51b8a087318fd765e3c`.
Five pairs exercised the ordered sweep path. The sparse input contained 10
rows and had SHA-256
`59037abb519691a2ab28435871c1d646b3d9e50b734790553a2328246160e12d`.
Thirty pairs batched 100 command invocations per timing sample and report the
per-invocation mean.

| Workload and tool | Mean wall time | Mean user time | Mean system time | Mean peak RSS |
|---|---:|---:|---:|---:|
| dense `rsomics-bam bedcov` | 0.792 s | 0.700 s | 0.026 s | 15,253,504 bytes |
| dense `samtools bedcov` | 8.238 s | 7.650 s | 0.342 s | 7,159,808 bytes |
| sparse `rsomics-bam bedcov` | 4.750 ms | 2.393 ms | 1.570 ms | 5,268,548 bytes |
| sparse `samtools bedcov` | 6.923 ms | 3.937 ms | 2.043 ms | 7,022,182 bytes |

The dense sweep won all five pairs and was 10.40 times as fast while using
2.13 times the peak RSS. The paired difference was -7.446 seconds with a
0.1126-second sample standard deviation. The sparse indexed path won all 30
pairs, was 1.46 times as fast, and used 25.0% less peak RSS. Its paired
difference was -2.173 ms with a 0.139-ms sample standard deviation. Dense and
sparse output SHA-256 values were
`2c1ddd98ae0504e40d7fc74a9e2c440516f2a3a8d5bb6fdbb10bce390f8a9954`
and
`b6720a48828b19dc4e44502ce33c5125ed1e9f8176dd80829c415d938691bee8`.

### Index statistics

Thirty pairs batched 100 explicit-index invocations per timing sample.

| Tool | Mean wall time | Mean user time | Mean system time | Mean peak RSS |
|---|---:|---:|---:|---:|
| `rsomics-bam idxstats` | 3.773 ms | 1.573 ms | 1.443 ms | 4,992,751 bytes |
| `samtools idxstats` | 5.440 ms | 2.700 ms | 1.917 ms | 6,646,443 bytes |

`rsomics-bam` won all 30 pairs, was 1.44 times as fast, and used 24.9% less
peak RSS. The paired difference was -1.667 ms with a 0.137-ms sample standard
deviation. The identical output SHA-256 was
`3a4e5fd538257de057d5f4a8264a8cc559bc5ef91cb93811f83e5414b9e79f03`.

The measured rsomics and samtools binaries had SHA-256 values
`746be956ec1df5e3c3bbf508d76d8bc1f327b7836593ad0eff7ac2130993cbe0`
and `c265b440b09c4b21d1f25a65963cf907b0d9f9d18caa9382c31104158f89d027`.
The coverage, dense BED, sparse BED, and index-statistics timing ledgers had
SHA-256 values
`5dda03b3b910eec489be629db5f96a1e30d0a3995bfb30a26879813f3dd2c0da`,
`f8541465617dc60f7b3fa965070dae6d128c406671222e5d8b75900ebb31d65e`,
`000871b3e68cf0a635ae9c0475057185d02e47dce4e0f2e229bb68db2925f6c4`,
and `dff4375e8000743761e77733f61fd1ee07df76017551b8a072e8979d0dbaf963`.

These results establish the default single-BAM summary paths on this fixture.
SAM, CRAM, multiple-input, filtered, region, threshold, and read-count behavior
is covered by compatibility tests but carries no throughput claim here. The
dense sweep trades a small bounded memory increase for its throughput gain;
other BED layouts and unusually deep data can select the indexed or pileup
path instead.

Release 0.16.0 was published from revision
`be3cafe21867c7773f40395d3669502909c2e12b` after exact-head CI
`31398246573` passed on native Linux and macOS for x86_64 and aarch64. The
Linux x86_64 job also passed the complete samtools 1.24 compatibility oracle.
Publication workflow `31398905778` produced an unyanked 190,502-byte registry
archive with SHA-256
`47f7bf82915054ac2a1fc1b66dbed35c77b47940fdb3f6c680ce789478be3345`;
its VCS metadata identifies the exact release revision.

A fresh registry install reports 0.16.0 and exposes `bedcov`, `coverage`, and
`idxstats` through the shared help tree. Its binary SHA-256 is
`9bd2f0ac01d7f67a16e4e7af4b37e3328038f75bf69c61da7195322638f2a4d7`.
Installed smoke tests produced BED coverage, complete reference coverage,
indexed-statistics fallback, and the shared JSON envelope. Their output
SHA-256 values were
`7ad62cb366232be9a14e4ec74b4661536672065a7e81f6a1321ded4d9cf633a9`,
`480fe5dea573ecee21aee4b149343309bf222736edfc2fd613e8b4430b908708`,
`b7ce5363e6b971a51f35597167dbc649c73432561ed55e8566ca379b5492aa9c`,
and `bc776a5f60348f7e963166445a8a1d464e3a8d052ed5defbd02d87ae20084402`.

## MD/NM recalculation benchmark

The 2026-08-10 gate used revision
`5e8e28129b3d06884f41d91443df1697babafe41` and samtools/HTSlib 1.24 on an
Apple M2 Mac mini with 8 GiB of memory. Both tools read and wrote BAM with four
additional I/O workers. One complete warm-up per tool preceded 20 timed pairs;
command order alternated between pairs.

```text
rsomics-bam calmd --no-pg -b -@ 4 input.bam reference.fa > rsomics.bam
samtools calmd --no-PG -b -@ 4 input.bam reference.fa > samtools.bam
```

| Tool | Mean wall time | Median wall time | Mean CPU time | Mean peak RSS |
|---|---:|---:|---:|---:|
| `rsomics-bam calmd` | 0.589 s | 0.570 s | 2.838 s | 19,285,606 bytes |
| `samtools calmd` | 0.922 s | 0.880 s | 2.913 s | 18,333,696 bytes |

`rsomics-bam` won all 20 pairs and was 1.57 times as fast by mean wall time.
The paired samtools-minus-rsomics difference was 0.3330 seconds with a
0.1430-second sample standard deviation and a 10.413 t-statistic. It used 2.6%
less mean CPU time and 5.2% more mean peak RSS.

The coordinate-sorted fixture contained 1,000,000 records, covered a
5,000,000-base reference at approximately 30x, and was 36,459,282 bytes. The
BAM and reference SHA-256 values were
`33b6780ec3758a8ccde746935366dec441e89aaafb5b0253a19cfa1af350282c`
and `13bd65f4568d0a30bc0ee218db62223cc26d9593f2b116530aa5e0b78b5f34dc`.
Both outputs passed `samtools quickcheck`; their complete decoded headers and
one million records had the identical SHA-256
`d1e0cfd0c1f1c1c88482e7140efc505ef323b0027ef1fac89be4c0b49d978eb9`.
BAM byte streams differed because their BGZF layouts differed.

The measured rsomics binary was built with rustc 1.97.1 and had SHA-256
`311632c80e80d575bd6c26a80ab8562f616a24e31e376568f2aecba2d438a605`;
the samtools binary SHA-256 was
`c265b440b09c4b21d1f25a65963cf907b0d9f9d18caa9382c31104158f89d027`.
The timing ledger, summary, and environment record had SHA-256 values
`96907c49cb6d40ce7148c5f3bbbe840abd8f7f7b4a5945d099e2f3cff7d54242`,
`6aa3ce9e7921213d0c2d6dfd342a6c694a866c662d49ad3c41d8c29568051164`,
and `61f4518bc44e9f7240c3530d6a8c00dca1927e25cc250ca6387aa3ecde57dffc`.
These results establish the default compressed BAM hot path on this fixture.
SAM, CRAM, uncompressed BAM, other worker counts, and materially different
reference or auxiliary-tag distributions carry no throughput claim here.

Release 0.17.0 was published from revision
`0debc103993f992aaf078291f23f3414b52acb3c` after exact-head CI
`31407557237` passed on native Linux and macOS for x86_64 and aarch64. The
Linux x86_64 job included the complete samtools 1.24 oracle. Publish workflow
`31408461408` produced the unyanked 198,528-byte registry archive with SHA-256
`6fd2ef2ad1c0072b3912d606b4bf52a2ee7d841a74a8af96383f53843eb6efc2`
and exact VCS metadata. A fresh registry install reported 0.17.0 and exposed
`calmd` through the shared help tree. Its named-BAM and shared-JSON smoke over
the one-million-record fixture completed all records, passed `samtools
quickcheck`, and reproduced decoded-output SHA-256
`d1e0cfd0c1f1c1c88482e7140efc505ef323b0027ef1fac89be4c0b49d978eb9`.

## Padded-reference projection benchmark

The 2026-08-11 gate used feature revision
`e1b8f89eed742984472d7230d0d65b17b4ade5c5` and samtools/HTSlib 1.24 on an
Apple M2 Mac mini with 8 GiB of memory. Both tools read and wrote compressed
BAM with four additional I/O workers and the same padded FASTA. One complete
warm-up per tool preceded 20 timed pairs; command order alternated between
pairs.

```text
rsomics-bam depad --no-pg -@ 4 -T padded.fa -o rsomics.bam padded.bam
samtools depad --no-PG --threads 4 -T padded.fa -o samtools.bam padded.bam
```

| Tool | Mean wall time | Median wall time | Mean CPU time | Mean peak RSS |
|---|---:|---:|---:|---:|
| `rsomics-bam depad` | 0.3680 s | 0.3500 s | 0.6175 s | 60,659,302 bytes |
| `samtools depad` | 0.6115 s | 0.6000 s | 0.5660 s | 37,438,259 bytes |

`rsomics-bam` won all 20 pairs and was 1.66 times as fast by mean wall time.
The paired samtools-minus-rsomics difference was 0.2435 seconds with a
0.0469-second sample standard deviation and a 23.202 t-statistic. The faster
path used 9.1% more mean CPU time and 62.0% more mean peak RSS.

The fixture contained 1,000,000 records against a 5,000,000-column padded
reference. The 4,034,944-byte BAM and 5,000,007-byte FASTA had SHA-256 values
`97d7b42c390da771119510c018d319aac782cdcde509a4d32db7607e7c8ba31f`
and `05895c75b94b2e8229bfdfa15cc16f084aa9c07d30807c9b22e286036a8e20a7`.
Both outputs passed `samtools quickcheck`; their complete decoded headers and
one million records had the identical SHA-256
`b56d7863308db97b0b081782d1bc39a8805c8c1086b00c6ff72dee68e46de904`.
BAM byte streams differed because their BGZF layouts differed.

The measured rsomics binary was built with rustc 1.91.0 and had SHA-256
`4d019db7157211f06b4d93e5e54960c93e00bccbc5ae77f2e1729b6a01616400`;
the samtools binary SHA-256 was
`c265b440b09c4b21d1f25a65963cf907b0d9f9d18caa9382c31104158f89d027`.
The timing ledger, summary, environment record, and decoded-output checksum
had SHA-256 values
`c386254a38be9d681b08635ae43f2bc15a53ebaa8b48e70a81ee8047eb6e39b2`,
`d8628dfd1e038c8c6b30417acfc38ed83d414fdb5cedb8e640506703467b6c50`,
`b3ae52071871766ec36547f8a88926ce6926354ae968a4eb7b681b01d60ff2dd`,
and `d09d1c3a484c19f98b77ade3072b5de5670cdcedb3c9e39ec6a65da22518d9c8`.
These results establish the compressed BAM path with a supplied padded FASTA
on this fixture. SAM, CRAM, embedded-reference projection, other worker
counts, and materially different gap distributions carry no throughput claim
here.

Release 0.18.0 was published from revision
`5304f278bfaaa8ca6c7d20fbcf3fb2900662884a` after exact-head CI
`31414206433` passed on native Linux and macOS for x86_64 and aarch64. The
Linux x86_64 job included package verification and the complete samtools 1.24
oracle. Publication workflow `31415017446` produced the unyanked 210,403-byte
registry archive with SHA-256
`e2a5f63c3cd11cdd8c8666883029879272467ccf1a8fa0efcd20a65675bee4f9`
and exact VCS metadata. A fresh registry install reports 0.18.0 and exposes
`depad` through the shared help tree. Its binary SHA-256 is
`5e8a960d9b79a3a5593f7960de813cc270d2ed0998a6af8653b5dd8c28c5f9fb`.
The named-BAM and shared-JSON smoke processed all one million records without
creating a FASTA sidecar, passed `samtools quickcheck`, and reproduced the
decoded-output SHA-256 above.

## CRAM storage diagnostic benchmark

The `cram-size` release gate used implementation revision `74c7a7fa8f06`,
samtools and HTSlib 1.24, and a 14,354,392-byte no-reference CRAM with SHA-256
`8b37d7ef3e2ac30236bb5b5c4bba27335b1ec2b71356e376db25e7864195d5c0`.
It contains 100 containers, 100 slices, 1,000,000 sequences, and 150,000,000
bases. The default reports were byte-identical with SHA-256
`e430528d73de3086be9032b811243138247d5ecca45f8f2f113b8e42a7570903`.
The encoding-map reports were also byte-identical.

The machine was an Apple M2 Mac mini with 8 GiB RAM and macOS 26.6.1. One
warm-up preceded 20 alternating paired rounds. Each timed process executed ten
complete commands serially and the CPU and wall values below are per command.

| Tool | Mean wall | Median wall | Mean user | Mean system | Mean peak RSS |
|---|---:|---:|---:|---:|---:|
| `rsomics-bam 74c7a7f cram-size` | 13.00 ms | 13.00 ms | 8.90 ms | 3.00 ms | 5,512,397 bytes |
| `samtools 1.24 cram-size` | 9.05 ms | 9.00 ms | 4.05 ms | 3.10 ms | 7,611,187 bytes |

The Rust implementation used 27.58% less mean peak RSS and was 43.65% slower
by mean wall time. This is a strict resource-use advantage, not a throughput
claim. The executable SHA-256 values were
`01907571489b4859d8d00160a0861d80e261b53cd32be57b7238e077d3095ee1`
for rsomics and
`c265b440b09c4b21d1f25a65963cf907b0d9f9d18caa9382c31104158f89d027`
for samtools.

The retained artifacts are under
`/Volumes/Zane's HDD/rsomics-fixtures/bam/cram-size-1m-complex-20260811/exact-74c7a7f`.
The environment, timing ledger, summary, and output-digest file have SHA-256
values `136a8d3deca43fc072748b537f60df981eeb4d6bab0d10f58f5c2fb466f2b8cf`,
`591cbbf18d290b0417bb8a461b53edb665728cc018f2b66153a714c177a2098e`,
`9728232c1eefd6a96bc7dbb7430c616e17e00e71e1b9829c2ff81e29fa099494`,
and `481fb6fe8f02e462a701a9f022fb1c35c239172a546b3301ae44e33bf250287d`.

## Reproduction

`benchmarks/cram-size-vs-samtools-macos.sh` compares complete reports
byte-for-byte, batches short runs for timer resolution, alternates command
order, and records macOS resource usage:

```sh
RSOMICS_COMMIT=74c7a7f benchmarks/cram-size-vs-samtools-macos.sh \
  target/release/rsomics-bam /path/to/samtools input.cram \
  /path/to/results 20 10 default
```

`benchmarks/depad-vs-samtools-macos.sh` compares padded-reference projection,
alternates command order, and rejects malformed output or any complete decoded
stream difference:

```sh
RSOMICS_COMMIT=e1b8f89 benchmarks/depad-vs-samtools-macos.sh \
  target/release/rsomics-bam /path/to/samtools input.bam padded.fa \
  /path/to/results 20 4
```

`benchmarks/calmd-vs-samtools-macos.sh` compares default compressed BAM
recalculation, alternates command order, and rejects malformed output or any
complete decoded-stream difference:

```sh
RSOMICS_COMMIT=5e8e281 benchmarks/calmd-vs-samtools-macos.sh \
  target/release/rsomics-bam /path/to/samtools input.bam reference.fa \
  /path/to/results 20 4
```

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

`benchmarks/merge-vs-samtools-macos.sh` records the corresponding ordered
merge comparison and rejects any complete-header or full-record checksum
disagreement:

```sh
RSOMICS_COMMIT=83b73a0 benchmarks/merge-vs-samtools-macos.sh \
  target/release/rsomics-bam \
  /path/to/samtools \
  /path/to/queryname-sorted.bam \
  /path/to/results \
  12 default natural
```

Use `equal-workers` as the mode to pass four additional workers to both tools.

`benchmarks/stats-vs-samtools-macos.sh` compares `coverage`, `idxstats`, or
`bedcov`, rejects any bytewise output difference, records checksums and macOS
resource usage, and alternates command order. Fast operations use a batch so
each timing sample exceeds the timer resolution:

```sh
RSOMICS_COMMIT=981d4e9 benchmarks/stats-vs-samtools-macos.sh \
  target/release/rsomics-bam /path/to/samtools coverage input.bam \
  /path/to/results/coverage 10

BENCH_BATCH=100 RSOMICS_COMMIT=981d4e9 \
  benchmarks/stats-vs-samtools-macos.sh \
  target/release/rsomics-bam /path/to/samtools idxstats input.bam \
  /path/to/results/idxstats 30 '' input.bam.bai

BENCH_BATCH=100 RSOMICS_COMMIT=981d4e9 \
  benchmarks/stats-vs-samtools-macos.sh \
  target/release/rsomics-bam /path/to/samtools bedcov input.bam \
  /path/to/results/bedcov 30 targets.bed input.bam.bai
```
The final order argument is one of `coordinate`, `natural`,
`lexicographical`, or `template-coordinate`.

`benchmarks/import-vs-samtools-macos.sh` compares single and paired FASTQ
import, alternates command order, and rejects any stable-header or complete
record-stream disagreement:

```sh
RSOMICS_COMMIT=1df1836 benchmarks/import-vs-samtools-macos.sh \
  target/release/rsomics-bam \
  /path/to/samtools \
  /path/to/reads-1.fq \
  /path/to/reads-2.fq \
  /path/to/results \
  12 4
```

`benchmarks/addreplacerg-vs-samtools-macos.sh` compares overwrite and orphan
read-group editing, alternates command order, and rejects any normalized-header
or complete record-stream disagreement:

```sh
RSOMICS_COMMIT=07bd58c benchmarks/addreplacerg-vs-samtools-macos.sh \
  target/release/rsomics-bam \
  /path/to/samtools \
  /path/to/mixed-read-group.bam \
  /path/to/results \
  12 4 'ID:new\tSM:after'
```

`benchmarks/file-operations-vs-samtools-macos.sh` compares either BAM
concatenation or BAM reheadering, alternates command order, and rejects any
complete decoded-stream disagreement:

```sh
RSOMICS_COMMIT=1178474 benchmarks/file-operations-vs-samtools-macos.sh \
  target/release/rsomics-bam /path/to/samtools cat /path/to/results 12 \
  shard-1.bam shard-2.bam shard-3.bam shard-4.bam

RSOMICS_COMMIT=62023e1 benchmarks/file-operations-vs-samtools-macos.sh \
  target/release/rsomics-bam /path/to/samtools reheader /path/to/results 12 \
  replacement.sam input.bam
```

`benchmarks/fastx-vs-samtools-macos.sh` compares unified FASTA and FASTQ
extraction, alternates command order, and rejects any bytewise stream
disagreement:

```sh
RSOMICS_COMMIT=d6cbf10 benchmarks/fastx-vs-samtools-macos.sh \
  target/release/rsomics-bam /path/to/samtools \
  /path/to/queryname-sorted.bam /path/to/results 10
```

`benchmarks/collate-vs-samtools-macos.sh` compares standard BAM collation,
alternates command order, and rejects any complete-header or order-independent
full-record multiset disagreement:

```sh
RSOMICS_COMMIT=24095b8 benchmarks/collate-vs-samtools-macos.sh \
  target/release/rsomics-bam \
  /path/to/samtools \
  /path/to/input.bam \
  /path/to/results \
  12 default
```

Use `equal-workers` as the mode to pass four additional workers to both tools.

`benchmarks/markdup-vs-samtools-macos.sh` compares default duplicate marking,
alternates command order, and rejects any complete-header, field, auxiliary
tag, or record-order disagreement:

```sh
RSOMICS_COMMIT=5c7dc56 benchmarks/markdup-vs-samtools-macos.sh \
  target/release/rsomics-bam \
  /path/to/samtools \
  /path/to/coordinate-sorted-fixmate.bam \
  /path/to/results \
  12 default
```

Use `equal-workers` as the mode to pass four additional workers to both tools.

`benchmarks/fixmate-vs-samtools-macos.sh` compares mate repair with mate-score
calculation, alternates command order, and rejects any complete-header or
order-sensitive full-record disagreement after normalizing auxiliary-tag
order:

```sh
RSOMICS_COMMIT=a8a684b benchmarks/fixmate-vs-samtools-macos.sh \
  target/release/rsomics-bam \
  /path/to/samtools \
  /path/to/queryname-sorted-pairs.bam \
  /path/to/results \
  12 default
```

Use `equal-workers` as the mode to pass four additional workers to both tools.
