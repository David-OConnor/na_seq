//! A general sequence type: DNA, RNA, or amino acids, with a name and metadata. This is a record
//! type, e.g. one entry of a FASTA or GenBank file; the bare nucleotide sequence is [`crate::Seq`].

use std::{
    collections::HashMap,
    fmt,
    fmt::{Display, Formatter},
    path::PathBuf,
};

use crate::{AminoAcid, Nucleotide, Seq};

/// The metadata key holding a sequence's free-text description. E.g. from the text following the
/// ID in a FASTA header, or a GenBank `DEFINITION` line.
pub const SEQ_DESCRIPTION_KEY: &str = "Description";

/// A nucleic-acid guess needs at least this fraction of its letters to be nucleotide letters
/// (including N). Protein sequences rarely come close, as they use most of the alphabet.
const NUCLEIC_LETTER_THRESH: f32 = 0.9;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum SeqType {
    AminoAcid,
    Dna,
    Rna,
}

impl Display for SeqType {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        let v = match self {
            Self::AminoAcid => "Protein",
            Self::Dna => "DNA",
            Self::Rna => "RNA",
        };

        write!(f, "{v}")
    }
}

impl SeqType {
    /// Abbreviated unit for a count of residues: "aa" for proteins, "nt" otherwise.
    pub fn residue_unit(self) -> &'static str {
        match self {
            Self::AminoAcid => "aa",
            Self::Dna | Self::Rna => "nt",
        }
    }

    /// Guess a sequence's type from its letters, for formats like FASTA that don't state it.
    /// Nucleic acid if nearly all letters are nucleotide letters; RNA if it has U and no T.
    pub fn infer(text: &str) -> Self {
        let mut letters = 0;
        let mut nucleic = 0;
        let mut has_t = false;
        let mut has_u = false;

        for b in text.bytes().filter(u8::is_ascii_alphabetic) {
            letters += 1;

            match b.to_ascii_uppercase() {
                b'A' | b'C' | b'G' | b'N' => nucleic += 1,
                b'T' => {
                    nucleic += 1;
                    has_t = true;
                }
                b'U' => {
                    nucleic += 1;
                    has_u = true;
                }
                _ => (),
            }
        }

        if letters == 0 || (nucleic as f32) < letters as f32 * NUCLEIC_LETTER_THRESH {
            return Self::AminoAcid;
        }

        if has_u && !has_t {
            Self::Rna
        } else {
            Self::Dna
        }
    }
}

/// The residues of a sequence.
///
/// [`Nucleotide`] has no uracil, so RNA is stored with `T` standing in for `U`; the RNA text
/// conversions here map between the two.
#[derive(Clone, Debug, PartialEq)]
pub enum SequenceData {
    AminoAcid(Vec<AminoAcid>),
    Dna(Seq),
    Rna(Seq),
}

impl SequenceData {
    pub fn get_type(&self) -> SeqType {
        match self {
            Self::AminoAcid(_) => SeqType::AminoAcid,
            Self::Dna(_) => SeqType::Dna,
            Self::Rna(_) => SeqType::Rna,
        }
    }

    pub fn len(&self) -> usize {
        match self {
            Self::AminoAcid(v) => v.len(),
            Self::Dna(v) | Self::Rna(v) => v.len(),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Parse single-letter residue codes, case-insensitive, as `seq_type`. Returns the data, and
    /// how many residue letters were left out because they can't be represented: ambiguity codes
    /// like N or X, and non-standard residues.
    ///
    /// Whitespace, digits, gaps (`-`, `.`), and stop codons (`*`) are ignored without being counted.
    /// For nucleic acids, T and U are both accepted.
    pub fn from_letters(text: &str, seq_type: SeqType) -> (Self, usize) {
        let mut skipped = 0;

        let letters = text
            .bytes()
            .filter(|b| !(b.is_ascii_whitespace() || b.is_ascii_digit()))
            .filter(|b| !matches!(b, b'-' | b'.' | b'*'));

        let data = match seq_type {
            SeqType::AminoAcid => {
                let mut result = Vec::with_capacity(text.len());
                for b in letters {
                    match AminoAcid::from_u8_letter(b) {
                        Ok(aa) => result.push(aa),
                        Err(_) => skipped += 1,
                    }
                }
                Self::AminoAcid(result)
            }
            SeqType::Dna | SeqType::Rna => {
                let mut result = Vec::with_capacity(text.len());
                for b in letters {
                    let b = match b {
                        b'U' => b'T',
                        b'u' => b't',
                        _ => b,
                    };

                    match Nucleotide::from_u8_letter(b) {
                        Ok(nt) => result.push(nt),
                        Err(_) => skipped += 1,
                    }
                }

                if seq_type == SeqType::Dna {
                    Self::Dna(result)
                } else {
                    Self::Rna(result)
                }
            }
        };

        (data, skipped)
    }

    /// Upper-case single-letter residue codes. RNA uses U.
    pub fn to_letters(&self) -> String {
        let bytes: Vec<u8> = match self {
            Self::AminoAcid(v) => v.iter().map(|aa| aa.to_u8_upper()).collect(),
            Self::Dna(v) => v.iter().map(|nt| nt.to_u8_upper()).collect(),
            Self::Rna(v) => v
                .iter()
                .map(|nt| match nt {
                    Nucleotide::T => b'U',
                    _ => nt.to_u8_upper(),
                })
                .collect(),
        };

        // All bytes produced above are ASCII letters.
        String::from_utf8(bytes).unwrap_or_default()
    }
}

/// A named sequence with metadata, e.g. one record of a FASTA or GenBank file.
#[derive(Clone, Debug, PartialEq)]
pub struct Sequence {
    pub data: SequenceData,
    /// E.g. a FASTA record's ID, or a GenBank `LOCUS` name. May be empty.
    pub name: String,
    /// Free-form fields, e.g. an accession or organism. The description, if present, is stored
    /// under [`SEQ_DESCRIPTION_KEY`].
    pub metadata: HashMap<String, String>,
    /// The file this was loaded from or last saved to, if any.
    pub path: Option<PathBuf>,
}

impl Sequence {
    pub fn new(data: SequenceData, name: String) -> Self {
        Self {
            data,
            name,
            metadata: HashMap::new(),
            path: None,
        }
    }

    pub fn seq_type(&self) -> SeqType {
        self.data.get_type()
    }

    pub fn description(&self) -> Option<&str> {
        self.metadata
            .get(SEQ_DESCRIPTION_KEY)
            .map(String::as_str)
            .filter(|d| !d.is_empty())
    }

    /// The name, or a placeholder if it's empty.
    pub fn display_name(&self) -> &str {
        if self.name.trim().is_empty() {
            "(Unnamed sequence)"
        } else {
            &self.name
        }
    }
}
