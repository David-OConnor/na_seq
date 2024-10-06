# Seq

This small library contains types and functions used for performing operations on DNA sequences. Its most fundamental type is the `Nucleotide` enum, representing a single DNA nucleotide. This library is general, and intended to be used by any program or library that uses DNA sequences.

It includes functions to convert between `&[Nucleotide]` to string and vice-versa, and convert to and from u8 integer representations.

This library is used by the [PlasCAD](https://github.com/David-OConnor/plascad) plasmid editor.

