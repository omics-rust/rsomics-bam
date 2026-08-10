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
The final order argument is one of `coordinate`, `natural`,
`lexicographical`, or `template-coordinate`.

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
