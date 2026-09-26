//! A general sequence type: DNA, RNA, or amino acids, with a name and metadata. This is a record
//! type, e.g. one entry of a FASTA or GenBank file; the bare nucleotide sequence is [`crate::Seq`].

use std::{
    collections::HashMap,
    fmt,
    fmt::{Display, Formatter},
    path::PathBuf,
};

use crate::{AminoAcid, Nucleotide, Seq, SeqTopology};

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

/// Characters in sequence text that aren't residues: whitespace, digits, gaps (`-`, `.`), and
/// stop codons (`*`).
fn is_ignored_letter(b: u8) -> bool {
    b.is_ascii_whitespace() || b.is_ascii_digit() || matches!(b, b'-' | b'.' | b'*')
}

/// Accepts U as well as T.
fn nucleotide_from_letter(b: u8) -> Option<Nucleotide> {
    let b = match b {
        b'U' => b'T',
        b'u' => b't',
        _ => b,
    };

    Nucleotide::from_u8_letter(b).ok()
}

/// Whether a residue letter can be represented as `seq_type`, e.g. false for N or X.
fn is_representable(b: u8, seq_type: SeqType) -> bool {
    match seq_type {
        SeqType::AminoAcid => AminoAcid::from_u8_letter(b).is_ok(),
        SeqType::Dna | SeqType::Rna => nucleotide_from_letter(b).is_some(),
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

        let letters = text.bytes().filter(|&b| !is_ignored_letter(b));

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
                    match nucleotide_from_letter(b) {
                        Some(nt) => result.push(nt),
                        None => skipped += 1,
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

/// A range of positions in a sequence: 1-based and inclusive, as in GenBank and SnapGene files.
/// On a circular sequence, `end < start` means the range wraps past the origin.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SeqRange {
    pub start: usize,
    pub end: usize,
}

impl SeqRange {
    pub fn new(start: usize, end: usize) -> Self {
        Self { start, end }
    }

    /// Shift this range from positions in a sequence's original letters to positions in the data
    /// parsed from them, which leaves out letters that can't be represented. `kept_before[i]` is the
    /// number of letters kept among the first `i`. `None` if none in the range were kept.
    fn remap(self, kept_before: &[usize]) -> Option<Self> {
        let len_orig = kept_before.len() - 1;
        let len_new = kept_before[len_orig];

        // Positions of the first and last kept letters in `start..=end`, not wrapping.
        let span = |start: usize, end: usize| {
            let start = start.clamp(1, len_orig.max(1));
            let end = end.min(len_orig);

            let new_start = kept_before[start - 1] + 1;
            let new_end = kept_before[end];

            (start <= end && new_start <= new_end).then_some((new_start, new_end))
        };

        if self.end >= self.start {
            let (start, end) = span(self.start, self.end)?;
            return Some(Self::new(start, end));
        }

        // Wraps past the origin: `start..=len`, then `1..=end`.
        match (span(self.start, len_orig), span(1, self.end)) {
            (Some((start, _)), Some((_, end))) => Some(Self::new(start, end)),
            (Some((start, _)), None) => Some(Self::new(start, len_new)),
            (None, Some((_, end))) => Some(Self::new(1, end)),
            (None, None) => None,
        }
    }
}

/// Which strand a feature is on.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Strand {
    /// Not stated, or not applicable.
    #[default]
    None,
    Forward,
    Reverse,
}

/// An annotated region of a sequence, e.g. a gene, promoter, or primer binding site. From GenBank
/// feature tables and SnapGene files.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SeqFeature {
    /// The feature key, e.g. "CDS", "promoter", or "primer_bind".
    pub kind: String,
    pub label: String,
    /// Usually one. Several for features made of segments, e.g. a GenBank `join(...)`.
    pub ranges: Vec<SeqRange>,
    pub strand: Strand,
    /// A display color, as hex, e.g. "#ff0000".
    pub color: Option<String>,
    /// Other annotations, e.g. ("note", "..."). A key may appear more than once.
    pub qualifiers: Vec<(String, String)>,
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
    /// For nucleic acids. `None` if not known, e.g. from a FASTA file.
    pub topology: Option<SeqTopology>,
    pub features: Vec<SeqFeature>,
    /// The file this was loaded from or last saved to, if any.
    pub path: Option<PathBuf>,
}

impl Sequence {
    pub fn new(data: SequenceData, name: String) -> Self {
        Self {
            data,
            name,
            metadata: HashMap::new(),
            topology: None,
            features: Vec::new(),
            path: None,
        }
    }

    /// Like [`SequenceData::from_letters`], for residue letters that `features` refer to by
    /// position, e.g. a GenBank ORIGIN section. The features' ranges are shifted to account for any
    /// letters left out, so they still cover the same residues; features covering only left-out
    /// letters are removed. Returns the number of letters left out.
    pub fn from_letters_with_features(
        letters: &str,
        seq_type: SeqType,
        name: String,
        mut features: Vec<SeqFeature>,
    ) -> (Self, usize) {
        let (data, skipped) = SequenceData::from_letters(letters, seq_type);

        if skipped > 0 && !features.is_empty() {
            let mut kept_before = vec![0];
            for b in letters.bytes().filter(|&b| !is_ignored_letter(b)) {
                let kept = kept_before[kept_before.len() - 1];
                kept_before.push(kept + is_representable(b, seq_type) as usize);
            }

            for feature in &mut features {
                feature.ranges = feature
                    .ranges
                    .iter()
                    .filter_map(|r| r.remap(&kept_before))
                    .collect();
            }
            features.retain(|f| !f.ranges.is_empty());
        }

        let mut result = Self::new(data, name);
        result.features = features;

        (result, skipped)
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
