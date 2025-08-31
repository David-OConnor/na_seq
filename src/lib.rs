use std::{io, io::ErrorKind};

use bincode::{Decode, Encode};

use crate::Nucleotide::*;
pub use crate::{
    amino_acids::{
        AaCategory, AaIdent, AminoAcid, AminoAcidGeneral, AminoAcidProtenationVariant, CodingResult,
    },
    element::{AtomTypeInRes, Element},
    nucleotide::{Nucleotide, NucleotideGeneral},
    restriction_enzyme::RestrictionEnzyme,
};

pub mod amino_acids;
pub mod element;
pub mod ligation;
pub mod nucleotide;
pub mod re_lib;
pub mod restriction_enzyme;

// Index 0: 5' end.
pub type Seq = Vec<Nucleotide>;

pub struct IndexError {}

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
        if let Ok(aa) = letter.parse::<AminoAcid>() {
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
                    format!("Invalid NT serialization: {byte}, {bits}"),
                )
            })?);
        }
    }

    Ok(result)
}

#[derive(Clone, Copy, PartialEq, Debug, Encode, Decode)]
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
