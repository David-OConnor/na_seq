//! This module contains types and functions for working with nucleotides.

use std::{io, io::ErrorKind};

use bincode::{Decode, Encode};
use num_enum::TryFromPrimitive;
use Nucleotide::*;

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
    /// E.g. For interop with FASTA, GenBank, and SnapGene formats.
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

/// This includes both normal nucleotides, and "either" combinations of nucleotides.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum NucleotideGeneral {
    A,
    T,
    C,
    G,
    /// Any
    N,
    /// A or T
    W,
    /// C or G
    S,
    /// Pyrimidines: C or T
    Y,
    /// Purines: A or G
    R,
    /// A or C
    M,
    /// G or T
    K,
}

impl NucleotideGeneral {
    /// Which nucleotides this symbol matches with.
    fn nt_matches(&self) -> Vec<Nucleotide> {
        match self {
            Self::A => vec![A],
            Self::T => vec![T],
            Self::C => vec![C],
            Self::G => vec![G],
            Self::N => vec![A, C, T, G],
            Self::W => vec![A, T],
            Self::S => vec![C, G],
            Self::Y => vec![C, T],
            Self::R => vec![A, G],
            Self::M => vec![A, C],
            Self::K => vec![T, T],
        }
    }

    pub fn matches(&self, nt: Nucleotide) -> bool {
        self.nt_matches().contains(&nt)
    }

    pub fn from_u8(val: u8) -> io::Result<Self> {
        match val {
            b'A' | b'a' => Ok(Self::A),
            b'T' | b't' => Ok(Self::T),
            b'G' | b'g' => Ok(Self::G),
            b'C' | b'c' => Ok(Self::C),
            b'N' | b'n' => Ok(Self::N),
            b'W' | b'w' => Ok(Self::W),
            b'S' | b's' => Ok(Self::S),
            b'Y' | b'y' => Ok(Self::Y),
            b'R' | b'r' => Ok(Self::T),
            b'M' | b'm' => Ok(Self::M),
            b'K' | b'k' => Ok(Self::K),
            _ => Err(io::Error::new(ErrorKind::InvalidData, "Invalid nucleotide")),
        }
    }

    pub fn to_u8_lower(&self) -> u8 {
        match self {
            Self::A => b'a',
            Self::T => b't',
            Self::C => b'c',
            Self::G => b'g',
            Self::N => b'n',
            Self::W => b'w',
            Self::S => b's',
            Self::Y => b'y',
            Self::R => b'r',
            Self::M => b'm',
            Self::K => b'k',
        }
        .to_owned()
    }

    pub fn to_u8_upper(&self) -> u8 {
        match self {
            Self::A => b'A',
            Self::T => b'T',
            Self::C => b'C',
            Self::G => b'G',
            Self::N => b'N',
            Self::W => b'W',
            Self::S => b'S',
            Self::Y => b'Y',
            Self::R => b'R',
            Self::M => b'M',
            Self::K => b'K',
        }
        .to_owned()
    }

    pub fn to_str_lower(&self) -> String {
        match self {
            Self::A => "a",
            Self::T => "t",
            Self::C => "c",
            Self::G => "g",
            Self::N => "n",
            Self::W => "w",
            Self::S => "s",
            Self::Y => "y",
            Self::R => "r",
            Self::M => "m",
            Self::K => "k",
        }
        .to_owned()
    }

    pub fn to_str_upper(&self) -> String {
        match self {
            Self::A => "A",
            Self::T => "T",
            Self::C => "C",
            Self::G => "G",
            Self::N => "N",
            Self::W => "W",
            Self::S => "S",
            Self::Y => "Y",
            Self::R => "R",
            Self::M => "M",
            Self::K => "K",
        }
        .to_owned()
    }
}
