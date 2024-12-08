use std::{io, io::ErrorKind};

use bincode::{Decode, Encode};
use num_enum::TryFromPrimitive;

pub use crate::amino_acids::{AaIdent, AminoAcid, CodingResult};
use crate::Nucleotide::*;

pub mod amino_acids;
pub mod ligation;
pub mod re_lib;
pub mod restriction_enzyme;

// Index 0: 5' end.
pub type Seq = Vec<Nucleotide>;

pub struct IndexError {}

/// A DNA nucleotide. The u8 repr is for use with a compact binary format.
/// This is the same nucleotide mapping as [.2bit format](http://genome.ucsc.edu/FAQ/FAQformat.html#format7).
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug, Encode, Decode, TryFromPrimitive)]
#[repr(u8)]
pub enum Nucleotide {
    T = 0b00,
    C = 0b01,
    A = 0b10,
    G = 0b11,
}

impl Nucleotide {
    /// For interop with FASTA, GenBank, and SnapGene formats.
    pub fn from_u8(val: u8) -> io::Result<Self> {
        match val {
            b'A' | b'a' => Ok(A),
            b'T' | b't' => Ok(T),
            b'G' | b'g' => Ok(G),
            b'C' | b'c' => Ok(C),
            _ => Err(io::Error::new(ErrorKind::InvalidData, "Invalid nucleotide")),
        }
    }

    /// Returns `b'A'` etc. For interop with FASTA, GenBank, and SnapGene formats.
    pub fn to_u8_upper(&self) -> u8 {
        match self {
            A => b'A',
            T => b'T',
            G => b'G',
            C => b'C',
        }
    }

    /// Returns `b'a'` etc. For interop with FASTA, GenBank, and SnapGene formats.
    pub fn to_u8_lower(&self) -> u8 {
        match self {
            A => b'a',
            T => b't',
            G => b'g',
            C => b'c',
        }
    }

    pub fn to_str_upper(&self) -> String {
        match self {
            A => "A".to_owned(),
            T => "T".to_owned(),
            C => "C".to_owned(),
            G => "G".to_owned(),
        }
    }

    pub fn to_str_lower(&self) -> String {
        match self {
            A => "a".to_owned(),
            T => "t".to_owned(),
            C => "c".to_owned(),
            G => "g".to_owned(),
        }
    }

    pub fn complement(self) -> Self {
        match self {
            A => T,
            T => A,
            G => C,
            C => G,
        }
    }

    /// Molecular weight, in Daltons, in a DNA strand.
    /// [Weight source: NorthWestern](http://biotools.nubic.northwestern.edu/OligoCalc.html)
    pub fn weight(&self) -> f32 {
        match self {
            A => 313.21,
            T => 304.2,
            G => 329.21,
            C => 289.18,
        }
    }

    /// Optical density of a 1mL solution, in a cuvette with 1cm pathlength.
    /// Result is in nm.
    /// http://biotools.nubic.northwestern.edu/OligoCalc.html
    pub fn a_max(&self) -> f32 {
        match self {
            A => 259.,
            T => 267.,
            G => 253.,
            C => 271.,
        }
    }

    /// Optical density of a 1mL solution, in a cuvette with 1cm pathlength.
    /// Result is in 1/(Moles x cm)
    /// http://biotools.nubic.northwestern.edu/OligoCalc.html
    pub fn molar_density(&self) -> f32 {
        match self {
            A => 15_200.,
            T => 8_400.,
            G => 12_010.,
            C => 7_050.,
        }
    }
}

/// Reverse direction, and swap C for G, A for T.
pub fn seq_complement(seq: &[Nucleotide]) -> Seq {
    let mut result = seq.to_vec();
    result.reverse();

    for nt in &mut result {
        *nt = nt.complement();
    }

    result
}

/// Create a nucleotide sequence from a string. (Case insensitive)
pub fn seq_from_str(str: &str) -> Seq {
    let mut result = Vec::new();

    for char in str.to_lowercase().chars() {
        match char {
            'a' => result.push(A),
            't' => result.push(T),
            'c' => result.push(C),
            'g' => result.push(G),
            _ => (),
        };
    }

    result
}

/// Create an amino-acid sequence from a string of single-letter identifiers. (Case insensitive)
pub fn seq_aa_from_str(str: &str) -> Vec<AminoAcid> {
    let mut result = Vec::new();

    for char in str.chars() {
        let letter = char.to_string(); // Convert `char` to `String`
        if let Ok(aa) = AminoAcid::from_str(&letter) {
            result.push(aa);
        }
    }

    result
}

/// Convert a nucleotide sequence to string.
pub fn seq_to_str_lower(seq: &[Nucleotide]) -> String {
    let mut result = String::new();

    for nt in seq {
        result.push_str(&nt.to_str_lower());
    }

    result
}

/// Convert a nucleotide sequence to string.
pub fn seq_to_str_upper(seq: &[Nucleotide]) -> String {
    let mut result = String::new();

    for nt in seq {
        result.push_str(&nt.to_str_upper());
    }

    result
}

/// Convert an amino acid sequence to string of single-letter idents.
pub fn seq_aa_to_str(seq: &[AminoAcid]) -> String {
    let mut result = String::new();

    for aa in seq {
        result.push_str(&aa.to_str(AaIdent::OneLetter));
    }

    result
}

/// Convert a sequence to bytes associated with UTF-8 letters. For compatibility with external libraries.
pub fn seq_to_u8_upper(seq: &[Nucleotide]) -> Vec<u8> {
    seq.iter().map(|nt| nt.to_u8_upper()).collect()
}

/// Convert a sequence of amino acids to bytes associated with UTF-8 letters. For compatibility with external libraries.
pub fn seq_to_u8_lower(seq: &[Nucleotide]) -> Vec<u8> {
    seq.iter().map(|nt| nt.to_u8_lower()).collect()
}

/// Convert a sequence of amino acids to bytes associated with UTF-8 letters. For compatibility with external libraries.
pub fn seq_aa_to_u8_upper(seq: &[AminoAcid]) -> Vec<u8> {
    seq.iter().map(|aa| aa.to_u8_upper()).collect()
}

/// Convert a string to bytes associated with UTF-8 letters. For compatibility with external libraries.
pub fn seq_aa_to_u8_lower(seq: &[AminoAcid]) -> Vec<u8> {
    seq.iter().map(|aa| aa.to_u8_lower()).collect()
}

/// Sequence weight, in Daltons. Assumes single-stranded.
pub fn seq_weight(seq: &[Nucleotide]) -> f32 {
    let mut result = 0.;

    for nt in seq {
        result += nt.weight();
    }

    result -= 61.96;

    result
}

/// Calculate portion of a sequence that is either the G or C nucleotide, on a scale of 0 to 1.
pub fn calc_gc(seq: &[Nucleotide]) -> f32 {
    let num_gc = seq.iter().filter(|&&nt| nt == C || nt == G).count();
    num_gc as f32 / seq.len() as f32
}

/// A compact binary serialization of our sequence. Useful for file storage.
/// The first four bytes is sequence length, big endian; we need this, since one of our nucleotides necessarily serializes
/// to 0b00.
///
/// MSB. Nucleotides are right-to-left in a given byte. Example: A byte containing
/// nucleotides TCAG is `0b1110_0100`.
pub fn serialize_seq_bin(seq: &[Nucleotide]) -> Vec<u8> {
    let mut result = Vec::new();
    result.extend(&(seq.len() as u32).to_be_bytes());

    for i in 0..seq.len() / 4 + 1 {
        let mut val = 0;
        for j in 0..4 {
            let ind = i * 4 + j;
            if ind + 1 > seq.len() {
                break;
            }
            let nt = seq[ind];
            val |= (nt as u8) << (j * 2);
        }
        result.push(val);
    }
    result
}

/// A compact binary deserialization of our sequence. Useful for file storage.
/// The first four bytes is sequence length, big endian; we need this, since one of our nucleotides necessarily serializes
/// to 0b00.
/// todo: Is this MSB or LSB?
pub fn deser_seq_bin(data: &[u8]) -> io::Result<Seq> {
    let mut result = Vec::new();

    if data.len() < 4 {
        return Err(io::Error::new(
            ErrorKind::InvalidData,
            "Bin nucleotide sequence is too short.",
        ));
    }

    let seq_len = u32::from_be_bytes(data[0..4].try_into().unwrap()) as usize;

    for byte in &data[4..] {
        for i in 0..4 {
            // This trimming removes extra 00-serialized nucleotides.
            if result.len() >= seq_len {
                break;
            }

            let bits = (byte >> (2 * i)) & 0b11;
            result.push(Nucleotide::try_from(bits).map_err(|_| {
                io::Error::new(
                    ErrorKind::InvalidData,
                    format!("Invalid NT serialization: {}, {}", byte, bits),
                )
            })?);
        }
    }

    Ok(result)
}

#[derive(Clone, Copy, PartialEq, Encode, Decode)]
pub enum SeqTopology {
    Linear,
    Circular,
}

impl Default for SeqTopology {
    fn default() -> Self {
        Self::Circular
    }
}

/// Insert a segment of one sequence into another. For example, for cloning.
/// Note that `insert_loc` uses 1-based indexing.
pub fn insert_into_seq(
    seq_vector: &mut Seq,
    insert: &[Nucleotide],
    insert_loc: usize,
) -> Result<(), IndexError> {
    if insert_loc == 0 || insert_loc > seq_vector.len() {
        eprintln!("Error: Insert location out of bounds: {insert_loc}");
        return Err(IndexError {});
    }

    let insert_i = insert_loc - 1; // 1-based indexing.
    seq_vector.splice(insert_i..insert_i, insert.iter().cloned());

    Ok(())
}
