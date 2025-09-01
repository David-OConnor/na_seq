# Nucleic Acid Sequence tools

[![Crate](https://img.shields.io/crates/v/na_seq.svg)](https://crates.io/crates/na_seq)
[![Docs](https://docs.rs/na_seq/badge.svg)](https://docs.rs/na_seq)
[![PyPI](https://img.shields.io/pypi/v/na-seq.svg)](https://pypi.org/project/na-seq)

This Rust and Python library contains types and functions used for performing operations on DNA and amino acid sequences. 
Its most fundamental types are the `Nucleotide` and `AminoAcid` enums, representing a single DNA nucleotide, 
and single amino acid respectively. This library is general, and intended to be used by any program or library 
that uses DNA sequences. It also includes an `Element` enum, with parameters associated with each element.

It includes functions to convert between `&[Nucleotide]` to string and vice-versa, and convert to and from u8 
representations of the UTF-8 characters. It includes functions to serialize and deserialize in a compact binary
format, with 2 bits per nucleotide.

It includes forcefield-parameter amino acid variants, as used by Amber as the `AminoAcidProtenationVariant` enum.

Basic types impl `Display` and `ToStr`, with variants in some cases, e.g. to display an AA as a single letter,
or 3-letter string.

See [the docs](https://docs.rs/na_seq) for details on data structures and functions available.

Also includes restriction enzyme and ligation basics.

Minmal python example. See the rust docs for full details:

```python
from na_seq import *

aa = AminoAcid.from_codons([Nucleotide.A, Nucleotide.T, Nucleotide.G])
aa = AminoAcid.from_str("HIS")

aa.to_str(AaIdent.OneLetter)
aa.to_str(AaIdent.ThreeLetters)

Element.C.atomic_weight()
```


## Utility functionality
- Sequence and nucleotide complements
- Sequence and nucleotide weight
- GC content, and misc general properties of elements and amino acids. ( e.g. hydrophobicity, weight, CPK colors)
- A small restriction enzyme library


We may add Sequence searches, and other utility features in the future.

This library is used by the [PlasCAD](https://github.com/David-OConnor/plascad) plasmid editor and [Daedalus](https://github.com/David-OConnor/daedalus) 
molecule viewer.

Note: The Python version currently does not include restriction-enzyme and ligation functionality.