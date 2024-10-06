use std::{io, io::ErrorKind};

use bincode::{Decode, Encode};
use num_enum::TryFromPrimitive;

use crate::Nucleotide::*;

// Index 0: 5' end.
pub type Seq = Vec<Nucleotide>;

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
    pub fn from_u8_letter(val_u8: u8) -> io::Result<Self> {
        match val_u8 {
            b'A' | b'a' => Ok(A),
            b'T' | b't' => Ok(T),
            b'G' | b'g' => Ok(G),
            b'C' | b'c' => Ok(C),
            _ => Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Invalid nucleotide",
            )),
        }
    }

    /// For interop with FASTA, GenBank, and SnapGene formats.
    pub fn to_u8_letter(&self) -> u8 {
        match self {
            A => b'A',
            T => b'T',
            G => b'G',
            C => b'C',
        }
    }

    pub fn as_str(&self) -> &str {
        match self {
            A => "a",
            T => "t",
            C => "c",
            G => "g",
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

/// Convert a nucleotide sequence to string.
pub fn seq_to_str(seq: &[Nucleotide]) -> String {
    let mut result = String::new();

    for nt in seq {
        result.push_str(nt.as_str());
    }

    result
}

/// Convert a string to bytes associated with ASCII letters. For compatibility with external libraries.
pub fn seq_to_letter_bytes(seq: &[Nucleotide]) -> Vec<u8> {
    seq.iter().map(|nt| nt.to_u8_letter()).collect()
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
/// The first byte is sequence length; we need this, since one of our nucleotides necessarily serializes
/// to 0b00.
/// todo: Is this MSB or LSB?
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
/// The first byte is sequence length; we need this, since one of our nucleotides necessarily serializes
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
