//! This module contains types and functions for working with nucleotides.

use std::{fmt, io, io::ErrorKind, str::FromStr};

use Nucleotide::*;
use bincode::{Decode, Encode};
use num_enum::TryFromPrimitive;

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

impl FromStr for Nucleotide {
    type Err = io::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::from_u8_letter(s.as_bytes()[0])
    }
}

impl fmt::Display for Nucleotide {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.to_str_upper())
    }
}

impl Nucleotide {
    /// E.g. For interop with FASTA, GenBank, and SnapGene formats.
    pub fn from_u8_letter(val: u8) -> io::Result<Self> {
        Ok(match val {
            b'A' | b'a' => A,
            b'T' | b't' => T,
            b'G' | b'g' => G,
            b'C' | b'c' => C,
            _ => {
                return Err(io::Error::new(
                    ErrorKind::InvalidData,
                    "Invalid nucleotide letter",
                ));
            }
        })
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
/// The u8 repr is for use with a binary format.
#[derive(Clone, Copy, Debug, PartialEq, Eq, TryFromPrimitive)]
#[repr(u8)]
pub enum NucleotideGeneral {
    T = 0,
    C = 1,
    A = 2,
    G = 3,
    /// Any
    N = 4,
    /// A or T
    W = 5,
    /// C or G
    S = 6,
    /// Pyrimidines: C or T
    Y = 7,
    /// Purines: A or G
    R = 8,
    /// A or C
    M = 9,
    /// G or T
    K = 10,
}

impl FromStr for NucleotideGeneral {
    type Err = io::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::from_u8_letter(s.as_bytes()[0])
    }
}

impl fmt::Display for NucleotideGeneral {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.to_str_upper())
    }
}

impl NucleotideGeneral {
    pub fn from_u8_letter(val: u8) -> io::Result<Self> {
        Ok(match val {
            b'T' | b't' => Self::T,
            b'C' | b'c' => Self::C,
            b'A' | b'a' => Self::A,
            b'G' | b'g' => Self::G,
            b'N' | b'n' => Self::N,
            b'W' | b'w' => Self::W,
            b'S' | b's' => Self::S,
            b'Y' | b'y' => Self::Y,
            b'R' | b'r' => Self::R,
            b'M' | b'm' => Self::M,
            b'K' | b'k' => Self::K,
            _ => {
                return Err(io::Error::new(
                    ErrorKind::InvalidData,
                    "Invalid nucleotide letter",
                ));
            }
        })
    }

    /// Which nucleotides this symbol matches with.
    fn nt_matches(&self) -> Vec<Nucleotide> {
        match self {
            Self::T => vec![T],
            Self::C => vec![C],
            Self::A => vec![A],
            Self::G => vec![G],
            Self::N => vec![A, C, T, G],
            Self::W => vec![A, T],
            Self::S => vec![C, G],
            Self::Y => vec![C, T],
            Self::R => vec![A, G],
            Self::M => vec![A, C],
            Self::K => vec![G, T],
        }
    }

    pub fn matches(&self, nt: Nucleotide) -> bool {
        self.nt_matches().contains(&nt)
    }

    pub fn to_u8_lower(&self) -> u8 {
        match self {
            Self::T => b't',
            Self::C => b'c',
            Self::A => b'a',
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
            Self::T => b'T',
            Self::C => b'C',
            Self::A => b'A',
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
            Self::T => "t",
            Self::C => "c",
            Self::A => "a",
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
