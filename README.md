# Nucleic Acid Sequence tools

[![Crate](https://img.shields.io/crates/v/na_seq.svg)](https://crates.io/crates/na_seq)
[![Docs](https://docs.rs/na_seq/badge.svg)](https://docs.rs/na_seq)

This small library contains types and functions used for performing operations on DNA sequences. Its most fundamental type is the `Nucleotide` enum, representing a single DNA nucleotide. This library is general, and intended to be used by any program or library that uses DNA sequences.

It includes functions to convert between `&[Nucleotide]` to string and vice-versa, and convert to and from u8 integer representations. It includes functions to serialize and deserialize in a compact binary format, with 2 bits per nucleotide.

It also includes Amino Acid types.

See [the docs](https://docs.rs/na_seq) for details on data structures and functions available.


## Utility functionality
- Sequence and nucleotide complements
- Sequence and nucleotide weight
- GC content

We may add alignment, sequence searches, and other utility features in the future.

Also includes restriction enzyme functionality basics.



This library is used by the [PlasCAD](https://github.com/David-OConnor/plascad) plasmid editor.