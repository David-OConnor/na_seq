use std::{fmt, io, str::FromStr};

use bincode::{Decode, Encode};

use crate::{Nucleotide, Nucleotide::*};

#[derive(Clone, Copy, PartialEq, Debug, Encode, Decode)]
pub enum AaIdent {
    OneLetter,
    ThreeLetters,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum CodingResult {
    AminoAcid(AminoAcid),
    StopCodon,
}

impl CodingResult {
    pub fn from_codons(codons: [Nucleotide; 3]) -> Self {
        // Handle cases that are defined entirely by the first two codons.
        match codons[0..2] {
            [C, G] => return Self::AminoAcid(AminoAcid::Arg),
            [C, C] => return Self::AminoAcid(AminoAcid::Pro),
            [C, T] => return Self::AminoAcid(AminoAcid::Leu),
            [T, C] => return Self::AminoAcid(AminoAcid::Ser),
            [G, G] => return Self::AminoAcid(AminoAcid::Gly),
            [G, C] => return Self::AminoAcid(AminoAcid::Ala),
            [G, T] => return Self::AminoAcid(AminoAcid::Val),
            [A, C] => return Self::AminoAcid(AminoAcid::Thr),
            _ => (),
        }

        match codons {
            [A, T, G] => Self::AminoAcid(AminoAcid::Met),
            [A, T, A] => Self::AminoAcid(AminoAcid::Ile),
            [A, T, C] => Self::AminoAcid(AminoAcid::Ile),
            [A, T, T] => Self::AminoAcid(AminoAcid::Ile),
            [C, A, G] => Self::AminoAcid(AminoAcid::Gln),
            [C, A, A] => Self::AminoAcid(AminoAcid::Gln),
            [C, A, C] => Self::AminoAcid(AminoAcid::His),
            [C, A, T] => Self::AminoAcid(AminoAcid::His),
            [T, G, G] => Self::AminoAcid(AminoAcid::Trp),
            [T, G, A] => CodingResult::StopCodon,
            [T, G, C] => Self::AminoAcid(AminoAcid::Cys),
            [T, G, T] => Self::AminoAcid(AminoAcid::Cys),
            [T, A, G] => CodingResult::StopCodon,
            [T, A, A] => CodingResult::StopCodon,
            [T, A, C] => Self::AminoAcid(AminoAcid::Tyr),
            [T, A, T] => Self::AminoAcid(AminoAcid::Tyr),
            [T, T, G] => Self::AminoAcid(AminoAcid::Leu),
            [T, T, A] => Self::AminoAcid(AminoAcid::Leu),
            [T, T, C] => Self::AminoAcid(AminoAcid::Phe),
            [T, T, T] => Self::AminoAcid(AminoAcid::Phe),
            [G, A, G] => Self::AminoAcid(AminoAcid::Glu),
            [G, A, A] => Self::AminoAcid(AminoAcid::Glu),
            [G, A, C] => Self::AminoAcid(AminoAcid::Asp),
            [G, A, T] => Self::AminoAcid(AminoAcid::Asp),
            [A, G, G] => Self::AminoAcid(AminoAcid::Arg),
            [A, G, A] => Self::AminoAcid(AminoAcid::Arg),
            [A, G, C] => Self::AminoAcid(AminoAcid::Ser),
            [A, G, T] => Self::AminoAcid(AminoAcid::Ser),
            [A, A, G] => Self::AminoAcid(AminoAcid::Lys),
            [A, A, A] => Self::AminoAcid(AminoAcid::Lys),
            [A, A, C] => Self::AminoAcid(AminoAcid::Asn),
            [A, A, T] => Self::AminoAcid(AminoAcid::Asn),
            _ => unreachable!(), // This the 2-nt pattners we handled above.
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum AaCategory {
    Hydrophobic,
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
        use AminoAcid::*;
        
        match ident {
            AaIdent::OneLetter => match self {
                Arg => "R",
                His => "H",
                Lys => "K",
                Asp => "D",
                Glu => "E",
                Ser => "S",
                Thr => "T",
                Asn => "N",
                Gln => "Q",
                Cys => "C",
                Sec => "U",
                Gly => "G",
                Pro => "P",
                Ala => "A",
                Val => "V",
                Ile => "I",
                Leu => "L",
                Met => "M",
                Phe => "F",
                Tyr => "Y",
                Trp => "W",
            },
            AaIdent::ThreeLetters => match self {
                Arg => "Arg",
                His => "His",
                Lys => "Lys",
                Asp => "Asp",
                Glu => "Glu",
                Ser => "Ser",
                Thr => "Thr",
                Asn => "Asn",
                Gln => "Gln",
                Cys => "Cys",
                Sec => "Sec",
                Gly => "Gly",
                Pro => "Pro",
                Ala => "Ala",
                Val => "Val",
                Ile => "Ile",
                Leu => "Leu",
                Met => "Met",
                Phe => "Phe",
                Tyr => "Tyr",
                Trp => "Trp",
            },
        }
        .to_owned()
    }

    /// Convert to a byte for the associated single-letter ident.
    pub fn to_u8_upper(&self) -> u8 {
        use AminoAcid::*;
        match self {
            Arg => b'R',
            His => b'H',
            Lys => b'K',
            Asp => b'D',
            Glu => b'E',
            Ser => b'S',
            Thr => b'T',
            Asn => b'N',
            Gln => b'Q',
            Cys => b'C',
            Sec => b'U',
            Gly => b'G',
            Pro => b'P',
            Ala => b'A',
            Val => b'V',
            Ile => b'I',
            Leu => b'L',
            Met => b'M',
            Phe => b'F',
            Tyr => b'Y',
            Trp => b'W',
        }
    }

    /// Convert to a byte for the associated single-letter ident.
    pub fn to_u8_lower(&self) -> u8 {
        use AminoAcid::*;
        match self {
            Arg => b'r',
            His => b'h',
            Lys => b'k',
            Asp => b'd',
            Glu => b'e',
            Ser => b's',
            Thr => b't',
            Asn => b'n',
            Gln => b'q',
            Cys => b'c',
            Sec => b'u',
            Gly => b'g',
            Pro => b'p',
            Ala => b'a',
            Val => b'v',
            Ile => b'i',
            Leu => b'l',
            Met => b'm',
            Phe => b'f',
            Tyr => b'y',
            Trp => b'w',
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
        use AminoAcid::*;
        match self {
            Arg => 174.,
            His => 155.,
            Lys => 146.,
            Asp => 133.,
            Glu => 147.,
            Ser => 105.,
            Thr => 119.,
            Asn => 132.,
            Gln => 146.,
            Cys => 121.,
            Sec => 168.06,
            Gly => 75.,
            Pro => 115.,
            Ala => 89.,
            Val => 117.,
            Ile => 131.,
            Leu => 131.,
            Met => 149.,
            Phe => 165.,
            Tyr => 181.,
            Trp => 204.,
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

    /// Returns None if a Stop Codon.
    pub fn from_codons(codons: [Nucleotide; 3]) -> Option<Self> {
        match CodingResult::from_codons(codons) {
            CodingResult::AminoAcid(aa) => Some(aa),
            CodingResult::StopCodon => None,
        }
    }

    /// https://en.wikipedia.org/wiki/DNA_and_RNA_codon_tables#/media/File:Aminoacids_table.svg
    /// Returns a Vec of all combinations of condos which can make up the amino acid.
    ///
    /// If a resulting codon has less than 3 nucleotides, it means the third can be any; this may have both conciseness,
    /// and performance advantages.
    pub fn codons(&self) -> Vec<Vec<Nucleotide>> {
        use AminoAcid::*;
        match self {
            Ala => vec![vec![G, C]],

            Arg => vec![
                vec![C, G],    // CGN
                vec![A, G, A], // AGA
                vec![A, G, G], // AGG
            ],

            Asn => vec![vec![A, A, T], vec![A, A, C]],
            Asp => vec![vec![G, A, T], vec![G, A, C]],
            Cys => vec![vec![T, G, T], vec![T, G, C]],
            Gln => vec![vec![C, A, G], vec![C, A, A]],
            Glu => vec![vec![G, A, A], vec![G, A, G]],
            Gly => vec![vec![G, G]],

            His => vec![vec![C, A, C], vec![C, A, T]],

            Ile => vec![vec![A, T, T], vec![A, T, C], vec![A, T, A]],

            Leu => vec![
                vec![C, T],    // CTN
                vec![T, T, A], // TTA
                vec![T, T, G], // TTG
            ],

            Lys => vec![vec![A, A, A], vec![A, A, G]],
            Met => vec![vec![A, T, G]],
            Phe => vec![vec![T, T, T], vec![T, T, C]],
            Pro => vec![vec![C, C]],

            Ser => vec![
                vec![T, C],    // TCN
                vec![A, G, T], // AGT
                vec![A, G, C], // AGC
            ],

            Thr => vec![vec![A, C]],
            Trp => vec![vec![T, G, G]],
            Tyr => vec![vec![T, A, T], vec![T, A, C]],
            Val => vec![vec![G, T]],
            Sec => unimplemented!(),
        }
    }

    pub fn category(&self) -> AaCategory {
        use AminoAcid::*;
        
        match self {
            Arg => AaCategory::Basic,
            His => AaCategory::Basic,
            Lys => AaCategory::Basic,
            Asp => AaCategory::Acidic,
            Glu => AaCategory::Acidic,
            Ser => AaCategory::Polar, // is polar equiv to hydrophilic?
            Thr => AaCategory::Polar,
            Asn => AaCategory::Polar,
            Gln => AaCategory::Polar,
            Cys => AaCategory::Polar,
            Sec => AaCategory::Polar, // todo: unknown for now. placeholder
            Gly => AaCategory::Hydrophobic,
            Pro => AaCategory::Hydrophobic,
            Ala => AaCategory::Hydrophobic,
            Val => AaCategory::Hydrophobic,
            Ile => AaCategory::Hydrophobic,
            Leu => AaCategory::Hydrophobic,
            Met => AaCategory::Hydrophobic,
            Phe => AaCategory::Hydrophobic,
            Tyr => AaCategory::Polar,
            Trp => AaCategory::Hydrophobic, // Maybe no hydro. Mayb epolar or amphipathic?
        }
    }
}

impl FromStr for AminoAcid {
    type Err = io::Error;

    /// The N and C-prefixed variants indicate N and C terminus amino acids. They are present,
    /// for example, in Amber data files `aminoct12.lib` and `aminont12.lib`.
    fn from_str(val: &str) -> Result<Self, Self::Err> {
        Ok(match val.to_uppercase().as_str() {
            "R" | "ARG" | "NARG" | "CARG" => Self::Arg,
            "H" | "HIS" | "NHIS" | "CHIS" => Self::His,
            "K" | "LYS" | "NLYS" | "CLYS" => Self::Lys,
            "D" | "ASP" | "NASP" | "CASP" => Self::Asp,
            "E" | "GLU" | "NGLU" | "CGLU" => Self::Glu,
            "S" | "SER" | "NSER" | "CSER" => Self::Ser,
            "T" | "THR" | "NTHR" | "CTHR" => Self::Thr,
            "N" | "ASN" | "NASN" | "CASN" => Self::Asn,
            "Q" | "GLN" | "NGLN" | "CGLN" => Self::Gln,
            "C" | "CYS" | "NCYS" | "CCYS" => Self::Cys,
            "U" | "SEC" | "NSEC" | "CSEC" => Self::Sec,
            "G" | "GLY" | "NGLY" | "CGLY" => Self::Gly,
            "P" | "PRO" | "NPRO" | "CPRO" => Self::Pro,
            "A" | "ALA" | "NALA" | "CALA" => Self::Ala,
            "V" | "VAL" | "NVAL" | "CVAL" => Self::Val,
            "I" | "ILE" | "NILE" | "CILE" => Self::Ile,
            "L" | "LEU" | "NLEU" | "CLEU" => Self::Leu,
            "M" | "MET" | "NMET" | "CMET" => Self::Met,
            "F" | "PHE" | "NPHE" | "CPHE" => Self::Phe,
            "Y" | "TYR" | "NTYR" | "CTYR" => Self::Tyr,
            "W" | "TRP" | "NTRP" | "CTRP" => Self::Trp,
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

        write!(f, "{v}")
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
            "HID" | "NHID" | "CHID" => Self::Hid,
            "HIE" | "NHIE" | "CHIE" => Self::Hie,
            "HIP" | "NHIP" | "CHIP" => Self::Hip,
            "CYM" | "NCYM" | "CCYM" => Self::Cym,
            "CYX" | "NCYX" | "CCYX" => Self::Cyx,
            "ASH" | "NASH" | "CASH" => Self::Ash,
            "GLH" | "NGLH" | "CGLH" => Self::Glh,
            "LYN" | "NLYN" | "CLYN" => Self::Lyn,
            "ACE" | "NACE" | "CACE" => Self::Ace,
            "NHE" | "NNHE" | "CNHE" => Self::Nhe,
            "NME" | "NNME" | "CNME" => Self::Nme,
            "HYP" | "NHYP" | "CHYP" => Self::Hyp,

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
        use AminoAcidProtenationVariant::*;
        
        let v = match self {
            Hid => "HID",
            Hie => "HIE",
            Hip => "HIP",
            Cym => "CYM",
            Cyx => "CYX",
            Ash => "ASH",
            Glh => "GLH",
            Lyn => "LYN",
            Ace => "ACE",
            Nhe => "NHE",
            Nme => "NME",
            Hyp => "HYP",
        };

        write!(f, "{v}")
    }
}

impl AminoAcidProtenationVariant {
    /// E.g. if Hid or Hie, get His. Returns None for the temrinal groups Ace and Nme.
    pub fn get_standard(&self) -> Option<AminoAcid> {
        use AminoAcid::*;
        
        match self {
            Self::Hid | Self::Hie | Self::Hip | Self::Nhe => Some(His),
            Self::Cym | Self::Cyx => Some(Cys),
            Self::Ash => Some(Asp),
            Self::Glh => Some(Glu),
            Self::Lyn => Some(Lys),
            Self::Hyp => Some(Pro),
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

/// Calculates the Kyte-Doolittle value of Hydropathy index. Uses per-AA
/// experimentally-determined hydrophobicity values, averaged over a moving window.
/// 9 is a good default for window size.
///
/// https://web.expasy.org/protscale/
pub fn hydropathy_doolittle(seq: &[AminoAcid], window_size: usize) -> Vec<(usize, f32)> {
    if window_size.is_multiple_of(2) {
        eprintln!("Window size for KD must be odd");
        return Vec::new();
    }
    let mut result = Vec::new();

    let win_div_2 = window_size / 2; // Rounds down.

    if win_div_2 > seq.len() {
        eprintln!("Error with window size for hydropathy");
        return result;
    }

    // We'll center each window on `i`.
    for i in win_div_2..seq.len() - win_div_2 - 1 {
        let aas = &seq[i - win_div_2..i + win_div_2 + 1];
        let mut val_this_window = 0.;
        for aa in aas {
            val_this_window += aa.hydropathicity();
        }
        result.push((i, val_this_window / window_size as f32));
    }

    result
}
