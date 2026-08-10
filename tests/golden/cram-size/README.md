# Version fixtures

The `version-*` CRAM files exercise the physical layouts written by samtools
1.24. They contain the three records in `../records.sam` and omit `@PG` lines.
The 2.1 fixture is reference-backed and uses the deterministic
`urn:rsomics:fixture:9889878875bfc855a532253c415dceb6:xxxxxxx` locator; pass
`-T ../reference.fa` when decoding it. The 3.0 and 3.1 fixtures store bases with
`no_ref=1`.

Their SHA-256 digests are:

```text
d708070cfdc253836c02ab5bf257ee105d7791cfba58c8b8ca9344250c585d5f  version-2.1.cram
5292e847141d31f662b2de55352ce8f01212dd23e599219fd59cfb2e9272f976  version-3.0.cram
0ecc49d4de424de61107b1c8f4d91bd4fe963cd4108b94118fdfd467339cb6df  version-3.1.cram
```

The adjacent text reports are direct `samtools cram-size` and
`samtools cram-size -e` output.
