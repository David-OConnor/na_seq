use std::{fmt, io, str::FromStr};

use bincode::{Decode, Encode};

use crate::{Nucleotide, Nucleotide::*};

#[derive(Clone, Copy, PartialEq, Debug, Encode, Decode)]
pub enum AaIdent {
    OneLetter,
    ThreeLetters,
}

#[derive(Clone, Copy, PartialEq, Debug)]
pub enum CodingResult {
    AminoAcid(AminoAcid),
    StopCodon,
}

#[derive(Clone, Copy, PartialEq)]
pub enum AaCategory {
    Hydrophobic,
    // Neutral, // todo: Rem? Currently a placeholder for
    // Hydrophilic,
    Acidic,
    Basic,
    Polar,
}

#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug, Encode, Decode)]
pub enum AminoAcid {
    Arg,
    His,
    Lys,
    Asp,
    Glu,
    Ser,
    Thr,
    Asn,
    Gln,
    Cys,
    Sec,
    Gly,
    Pro,
    Ala,
    Val,
    Ile,
    Leu,
    Met,
    Phe,
    Tyr,
    Trp,
}

impl AminoAcid {
    pub fn to_str(&self, ident: AaIdent) -> String {
        match ident {
            AaIdent::OneLetter => match self {
                Self::Arg => "R",
                Self::His => "H",
                Self::Lys => "K",
                Self::Asp => "D",
                Self::Glu => "E",
                Self::Ser => "S",
                Self::Thr => "T",
                Self::Asn => "N",
                Self::Gln => "Q",
                Self::Cys => "C",
                Self::Sec => "U",
                Self::Gly => "G",
                Self::Pro => "P",
                Self::Ala => "A",
                Self::Val => "V",
                Self::Ile => "I",
                Self::Leu => "L",
                Self::Met => "M",
                Self::Phe => "F",
                Self::Tyr => "Y",
                Self::Trp => "W",
            },
            AaIdent::ThreeLetters => match self {
                Self::Arg => "Arg",
                Self::His => "His",
                Self::Lys => "Lys",
                Self::Asp => "Asp",
                Self::Glu => "Glu",
                Self::Ser => "Ser",
                Self::Thr => "Thr",
                Self::Asn => "Asn",
                Self::Gln => "Gln",
                Self::Cys => "Cys",
                Self::Sec => "Sec",
                Self::Gly => "Gly",
                Self::Pro => "Pro",
                Self::Ala => "Ala",
                Self::Val => "Val",
                Self::Ile => "Ile",
                Self::Leu => "Leu",
                Self::Met => "Met",
                Self::Phe => "Phe",
                Self::Tyr => "Tyr",
                Self::Trp => "Trp",
            },
        }
        .to_owned()
    }

    /// Convert to a byte for the associated single-letter ident.
    pub fn to_u8_upper(&self) -> u8 {
        match self {
            Self::Arg => b'R',
            Self::His => b'H',
            Self::Lys => b'K',
            Self::Asp => b'D',
            Self::Glu => b'E',
            Self::Ser => b'S',
            Self::Thr => b'T',
            Self::Asn => b'N',
            Self::Gln => b'Q',
            Self::Cys => b'C',
            Self::Sec => b'U',
            Self::Gly => b'G',
            Self::Pro => b'P',
            Self::Ala => b'A',
            Self::Val => b'V',
            Self::Ile => b'I',
            Self::Leu => b'L',
            Self::Met => b'M',
            Self::Phe => b'F',
            Self::Tyr => b'Y',
            Self::Trp => b'W',
        }
    }

    /// Convert to a byte for the associated single-letter ident.
    pub fn to_u8_lower(&self) -> u8 {
        match self {
            Self::Arg => b'r',
            Self::His => b'h',
            Self::Lys => b'k',
            Self::Asp => b'd',
            Self::Glu => b'e',
            Self::Ser => b's',
            Self::Thr => b't',
            Self::Asn => b'n',
            Self::Gln => b'q',
            Self::Cys => b'c',
            Self::Sec => b'u',
            Self::Gly => b'g',
            Self::Pro => b'p',
            Self::Ala => b'a',
            Self::Val => b'v',
            Self::Ile => b'i',
            Self::Leu => b'l',
            Self::Met => b'm',
            Self::Phe => b'f',
            Self::Tyr => b'y',
            Self::Trp => b'w',
        }
    }

    /// Used to make displaying a centered letter in a sequence easier; 3 characters.
    pub fn to_str_offset(&self) -> String {
        format!(" {} ", self.to_str(AaIdent::OneLetter))
    }

    /// Return the molecular weight, in Da.
    /// Source: https://www.promega.com/resources/tools/amino-acid-chart-amino-acid-structure/
    /// todo: This table is not very precise; consider updating with a better source.
    pub fn weight(&self) -> f32 {
        match self {
            Self::Arg => 174.,
            Self::His => 155.,
            Self::Lys => 146.,
            Self::Asp => 133.,
            Self::Glu => 147.,
            Self::Ser => 105.,
            Self::Thr => 119.,
            Self::Asn => 132.,
            Self::Gln => 146.,
            Self::Cys => 121.,
            Self::Sec => 168.06,
            Self::Gly => 75.,
            Self::Pro => 115.,
            Self::Ala => 89.,
            Self::Val => 117.,
            Self::Ile => 131.,
            Self::Leu => 131.,
            Self::Met => 149.,
            Self::Phe => 165.,
            Self::Tyr => 181.,
            Self::Trp => 204.,
        }
    }

    /// Used for determining protein hydropathy. High (eg positive) values intdicate hydrophilic
    /// AAs. (Seems to not be completely true from some example checks? Some traditionally hydrophilic
    /// proteins like Proline (-1.6) and Glycine (-4) are on the list, but the very negative values
    /// are not associated with traditionally hydrophillic AAs.
    /// [Kyte, Doolittle](https://web.expasy.org/protscale/pscale/Hydropath.Doolittle.html)
    pub fn hydropathicity(&self) -> f32 {
        match self {
            Self::Arg => -4.5,
            Self::His => -3.2,
            Self::Lys => -3.9,
            Self::Asp => -3.5,
            Self::Glu => -3.5,
            Self::Ser => -0.8,
            Self::Thr => -0.7,
            Self::Asn => -3.5,
            Self::Gln => -3.5,
            Self::Cys => 2.5,
            Self::Sec => 0., // todo?
            Self::Gly => -0.4,
            Self::Pro => -1.6,
            Self::Ala => 1.8,
            Self::Val => 4.2,
            Self::Ile => 4.5,
            Self::Leu => 3.8,
            Self::Met => 1.9,
            Self::Phe => 2.8,
            Self::Tyr => -1.3,
            Self::Trp => -0.9,
        }
    }

    /// https://en.wikipedia.org/wiki/DNA_and_RNA_codon_tables#/media/File:Aminoacids_table.svg
    /// If a codon has less than 3 nucleotides, it means the third can be any; this may have both conciseness,
    /// and performance advantages.
    pub fn codons(&self) -> Vec<Vec<Nucleotide>> {
        match self {
            // todo: Should we do wildcards etc, to speed up matching? Ie Arg is just [C, G].
            Self::Arg => vec![vec![C, G]],
            Self::Gln => vec![vec![C, A, G], vec![C, A, A]],
            Self::His => vec![vec![C, A, C], vec![C, A, T]],
            Self::Pro => vec![vec![C, C]],
            Self::Leu => vec![vec![C, T]],
            Self::Met => vec![vec![A, T, G]],
            _ => Vec::new(),
        }
    }

    // todo: Move to CodingResult?
    pub fn from_codons(codons: [Nucleotide; 3]) -> CodingResult {
        // Handle cases that are defined entirely by the first two codons.
        match codons[0..2] {
            [C, G] => return CodingResult::AminoAcid(Self::Arg),
            [C, C] => return CodingResult::AminoAcid(Self::Pro),
            [C, T] => return CodingResult::AminoAcid(Self::Leu),
            [T, C] => return CodingResult::AminoAcid(Self::Ser),
            [G, G] => return CodingResult::AminoAcid(Self::Gly),
            [G, C] => return CodingResult::AminoAcid(Self::Ala),
            [G, T] => return CodingResult::AminoAcid(Self::Val),
            [A, C] => return CodingResult::AminoAcid(Self::Thr),
            _ => (),
        }

        match codons {
            [A, T, G] => CodingResult::AminoAcid(Self::Met),
            [A, T, A] => CodingResult::AminoAcid(Self::Ile),
            [A, T, C] => CodingResult::AminoAcid(Self::Ile),
            [A, T, T] => CodingResult::AminoAcid(Self::Ile),
            [C, A, G] => CodingResult::AminoAcid(Self::Gln),
            [C, A, A] => CodingResult::AminoAcid(Self::Gln),
            [C, A, C] => CodingResult::AminoAcid(Self::His),
            [C, A, T] => CodingResult::AminoAcid(Self::His),
            [T, G, G] => CodingResult::AminoAcid(Self::Trp),
            [T, G, A] => CodingResult::StopCodon,
            [T, G, C] => CodingResult::AminoAcid(Self::Cys),
            [T, G, T] => CodingResult::AminoAcid(Self::Cys),
            [T, A, G] => CodingResult::StopCodon,
            [T, A, A] => CodingResult::StopCodon,
            [T, A, C] => CodingResult::AminoAcid(Self::Tyr),
            [T, A, T] => CodingResult::AminoAcid(Self::Tyr),
            [T, T, G] => CodingResult::AminoAcid(Self::Leu),
            [T, T, A] => CodingResult::AminoAcid(Self::Leu),
            [T, T, C] => CodingResult::AminoAcid(Self::Phe),
            [T, T, T] => CodingResult::AminoAcid(Self::Phe),
            [G, A, G] => CodingResult::AminoAcid(Self::Glu),
            [G, A, A] => CodingResult::AminoAcid(Self::Glu),
            [G, A, C] => CodingResult::AminoAcid(Self::Asp),
            [G, A, T] => CodingResult::AminoAcid(Self::Asp),
            [A, G, G] => CodingResult::AminoAcid(Self::Arg),
            [A, G, A] => CodingResult::AminoAcid(Self::Arg),
            [A, G, C] => CodingResult::AminoAcid(Self::Ser),
            [A, G, T] => CodingResult::AminoAcid(Self::Ser),
            [A, A, G] => CodingResult::AminoAcid(Self::Lys),
            [A, A, A] => CodingResult::AminoAcid(Self::Lys),
            [A, A, C] => CodingResult::AminoAcid(Self::Asn),
            [A, A, T] => CodingResult::AminoAcid(Self::Asn),
            _ => unreachable!(), // This the 2-nt pattners we handled above.
        }
    }

    pub fn category(&self) -> AaCategory {
        match self {
            Self::Arg => AaCategory::Basic,
            Self::His => AaCategory::Basic,
            Self::Lys => AaCategory::Basic,
            Self::Asp => AaCategory::Acidic,
            Self::Glu => AaCategory::Acidic,
            Self::Ser => AaCategory::Polar, // is polar equiv to hydrophilic?
            Self::Thr => AaCategory::Polar,
            Self::Asn => AaCategory::Polar,
            Self::Gln => AaCategory::Polar,
            Self::Cys => AaCategory::Polar,
            Self::Sec => AaCategory::Polar, // todo: unknown for now. placeholder
            Self::Gly => AaCategory::Hydrophobic,
            Self::Pro => AaCategory::Hydrophobic,
            Self::Ala => AaCategory::Hydrophobic,
            Self::Val => AaCategory::Hydrophobic,
            Self::Ile => AaCategory::Hydrophobic,
            Self::Leu => AaCategory::Hydrophobic,
            Self::Met => AaCategory::Hydrophobic,
            Self::Phe => AaCategory::Hydrophobic,
            Self::Tyr => AaCategory::Polar,
            Self::Trp => AaCategory::Hydrophobic, // Maybe no hydro. Mayb epolar or amphipathic?
        }
    }
}

impl FromStr for AminoAcid {
    type Err = io::Error;

    /// The N and C-prefixed variants indicate N and C terminus amino acids. They are present,
    /// for example, in Amber data files `aminoct12.lib` and `aminont12.lib`.
    fn from_str(val: &str) -> Result<Self, Self::Err> {
        Ok(match val.to_uppercase().as_str() {
            "R"   | "ARG"  | "NARG"  | "CARG"  => Self::Arg,
            "H"   | "HIS"  | "NHIS"  | "CHIS"  => Self::His,
            "K"   | "LYS"  | "NLYS"  | "CLYS"  => Self::Lys,
            "D"   | "ASP"  | "NASP"  | "CASP"  => Self::Asp,
            "E"   | "GLU"  | "NGLU"  | "CGLU"  => Self::Glu,
            "S"   | "SER"  | "NSER"  | "CSER"  => Self::Ser,
            "T"   | "THR"  | "NTHR"  | "CTHR"  => Self::Thr,
            "N"   | "ASN"  | "NASN"  | "CASN"  => Self::Asn,
            "Q"   | "GLN"  | "NGLN"  | "CGLN"  => Self::Gln,
            "C"   | "CYS"  | "NCYS"  | "CCYS"  => Self::Cys,
            "U"   | "SEC"  | "NSEC"  | "CSEC"  => Self::Sec,
            "G"   | "GLY"  | "NGLY"  | "CGLY"  => Self::Gly,
            "P"   | "PRO"  | "NPRO"  | "CPRO"  => Self::Pro,
            "A"   | "ALA"  | "NALA"  | "CALA"  => Self::Ala,
            "V"   | "VAL"  | "NVAL"  | "CVAL"  => Self::Val,
            "I"   | "ILE"  | "NILE"  | "CILE"  => Self::Ile,
            "L"   | "LEU"  | "NLEU"  | "CLEU"  => Self::Leu,
            "M"   | "MET"  | "NMET"  | "CMET"  => Self::Met,
            "F"   | "PHE"  | "NPHE"  | "CPHE"  => Self::Phe,
            "Y"   | "TYR"  | "NTYR"  | "CTYR"  => Self::Tyr,
            "W"   | "TRP"  | "NTRP"  | "CTRP"  => Self::Trp,
            _ => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "Invalid amino acid string provided",
                ));
            }
        })
    }
}

impl fmt::Display for AminoAcid {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let v = format!(
            "{} ({})",
            self.to_str(AaIdent::ThreeLetters),
            self.to_str(AaIdent::OneLetter)
        );

        write!(f, "{}", v)
    }
}

/// Representations of amino acids in non-standard tauteromic and protenation states.
/// See [Amber RM](https://ambermd.org/doc12/Amber25.pdf), section 13.2: Residue naming conventions.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug, Encode, Decode)]
pub enum AminoAcidProtenationVariant {
    /// His
    Hid,
    Hie,
    Hip,
    /// Cys
    Cym,
    Cyx,
    /// Asp
    Ash,
    /// Glu
    Glh,
    /// Lys
    Lyn,
    /// Terminals
    /// Acetyl group
    Ace,
    /// N-methylamid group
    Nhe,
    /// Neutral histidine
    Nme,
    /// Proline
    Hyp,
}

impl FromStr for AminoAcidProtenationVariant {
    type Err = io::Error;

    /// See note on N and C-prefixed on AminoAcid::from_str
    fn from_str(val: &str) -> Result<Self, Self::Err> {
        Ok(match val.to_uppercase().as_str() {
            "HID"  | "NHID"  | "CHID"  => Self::Hid,
            "HIE"  | "NHIE"  | "CHIE"  => Self::Hie,
            "HIP"  | "NHIP"  | "CHIP"  => Self::Hip,
            "CYM"  | "NCYM"  | "CCYM"  => Self::Cym,
            "CYX"  | "NCYX"  | "CCYX"  => Self::Cyx,
            "ASH"  | "NASH"  | "CASH"  => Self::Ash,
            "GLH"  | "NGLH"  | "CGLH"  => Self::Glh,
            "LYN"  | "NLYN"  | "CLYN"  => Self::Lyn,
            "ACE"  | "NACE"  | "CACE"  => Self::Ace,
            "NHE"  | "NNHE"  | "CNHE"  => Self::Nhe,
            "NME"  | "NNME"  | "CNME"  => Self::Nme,
            "HYP"  | "NHYP"  | "CHYP"  => Self::Hyp,

            _ => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "Invalid amino acid string provided",
                ));
            }
        })
    }
}

impl fmt::Display for AminoAcidProtenationVariant {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let v = match self {
            Self::Hid => "HID",
            Self::Hie => "HIE",
            Self::Hip => "HIP",
            Self::Cym => "CYM",
            Self::Cyx => "CYX",
            Self::Ash => "ASH",
            Self::Glh => "GLH",
            Self::Lyn => "LYN",
            Self::Ace => "ACE",
            Self::Nhe => "NHE",
            Self::Nme => "NME",
            Self::Hyp => "HYP",
        };

        write!(f, "{}", v)
    }
}

impl AminoAcidProtenationVariant {
    /// E.g. if Hid or Hie, get His. Returns None for the temrinal groups Ace and Nme.
    pub fn get_standard(&self) -> Option<AminoAcid> {
        match self {
            Self::Hid | Self::Hie | Self::Hip | Self::Nhe => Some(AminoAcid::His),
            Self::Cym | Self::Cyx => Some(AminoAcid::Cys),
            Self::Ash => Some(AminoAcid::Asp),
            Self::Glh => Some(AminoAcid::Glu),
            Self::Lyn => Some(AminoAcid::Lys),
            Self::Hyp => Some(AminoAcid::Pro),
            _ => None,
        }
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug, Encode, Decode)]
/// Allows either normal, or a protenation variant. Useful when parsing Amber amino acid
/// parameter files.
pub enum AminoAcidGeneral {
    Standard(AminoAcid),
    Variant(AminoAcidProtenationVariant),
}

impl FromStr for AminoAcidGeneral {
    type Err = io::Error;

    fn from_str(val: &str) -> Result<Self, Self::Err> {
        match AminoAcid::from_str(val) {
            Ok(v) => Ok(Self::Standard(v)),
            Err(_) => match AminoAcidProtenationVariant::from_str(val) {
                Ok(v) => Ok(Self::Variant(v)),
                Err(e) => Err(e),
            },
        }
    }
}
