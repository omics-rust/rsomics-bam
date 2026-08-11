#!/usr/bin/env python3

import sys


sequences = (
    "ACGTACGTAA",
    "ACGTACGTAA",
    "TCGTACGTAA",
    "TCGTACGTAA",
    "ACGTTCGTAA",
    "ACGTTCGTAA",
    "TCGTTCGTAA",
    "TCGTTCGTAA",
    "ACGTACGTAA",
    "TCGTTCGTAA",
    "ACGTACGTAA",
    "TCGTTCGTAA",
)


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit(f"usage: {sys.argv[0]} GROUPS")
    try:
        groups = int(sys.argv[1])
    except ValueError as error:
        raise SystemExit("GROUPS must be an integer") from error
    if not 1 <= groups <= 10_000_000:
        raise SystemExit("GROUPS must be between 1 and 10000000")
    write = sys.stdout.write
    write(f"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:{groups * 200 + 100}\n")
    for group in range(groups):
        position = group * 200 + 1
        for read, sequence in enumerate(sequences, 1):
            write(
                f"g{group:07d}r{read:02d}\t0\tchr1\t{position}\t60\t10M\t*\t0\t0\t"
                f"{sequence}\tIIIIIIIIII\n"
            )


if __name__ == "__main__":
    main()
