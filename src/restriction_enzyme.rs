//! This module contains info related to Restriction Enzyme sites.
//!
//! [Wikipedia: List of RE sites](https://en.wikipedia.org/wiki/List_of_restriction_enzyme_cutting_sites:_A)
//! [NEB guide](https://www.neb.com/en-us/tools-and-resources/selection-charts/frequencies-of-restriction-sites)
//!
//! Note: This module only currently includes a selection of popular REs, and only ones that match
//! exact NTs.

use std::{
    collections::{HashMap, hash_map::Entry},
    hash::{Hash, Hasher},
};

use crate::{Nucleotide, NucleotideGeneral, Seq};

pub struct LigationProduct {
    /// 5' to 3' (both strands; they are in opposite directions.)
    pub strand_top: Seq,
    pub strand_bottom: Seq,
    /// The index, in nucleotides, from the start of the top strand, where the
    /// bottom strand starts. Example:
    /// ACTGG  (top)
    ///    CC  (bottom)   alignment=3.
    pub alignment: usize,
}

#[derive(Debug, Clone)]
pub struct ReMatch {
    pub lib_index: usize,
    /// Cuts after this index, in the "forward" direction.
    pub seq_index: usize,
    // pub direction: PrimerDirection,
    /// todo: Experimenting
    /// The number of matches found for this RE.
    pub match_count: usize,
}

#[derive(Clone, Eq)]
pub struct RestrictionEnzyme {
    pub name: String,
    /// From the 5' end.
    // pub seq: Seq, // todo: You may eventually need Vec<NucleotideGeneral>.
    pub cut_seq: Vec<NucleotideGeneral>,
    /// Index to cut after, from the 5' end. For blunt ends, this will be
    /// halfway through the seq (rounded down)
    pub cut_after: u8,
}

impl Hash for RestrictionEnzyme {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.name.hash(state);
    }
}

impl PartialEq for RestrictionEnzyme {
    fn eq(&self, other: &Self) -> bool {
        self.name == other.name
    }
}

impl RestrictionEnzyme {
    // pub fn new(name: &str, seq: Seq, cut_after: u8) -> Self {
    pub fn new(name: &str, cut_seq: Vec<NucleotideGeneral>, cut_after: u8) -> Self {
        Self {
            name: name.to_owned(),
            cut_seq,
            cut_after,
        }
    }

    pub fn makes_blunt_ends(&self) -> bool {
        self.cut_after as isize + 1 == self.cut_seq.len() as isize / 2
    }

    /// A depiction of where to cut.
    pub fn cut_depiction(&self) -> String {
        let nt_chars = seq_general_to_str(&self.cut_seq);

        let mut result = String::new();

        for (i, nt_char) in nt_chars.chars().enumerate() {
            result.push(nt_char);
            if i as u8 == self.cut_after {
                result.push_str(" | ");
            }
        }

        result
    }

    // todo: Consider replacing these with a dual-stranded model, instead of
    // todo modeling overhangs.

    /// Find the overhanging NTs 5' of a sequence's top strand.
    /// `seq_segment` must be the same size as, and aligned with the cut sequence.
    pub fn overhang_top_left(&self, seq_segment: &[Nucleotide]) -> Vec<Nucleotide> {
        let cut = self.cut_after as usize + 1;
        let len = self.cut_seq.len();

        if cut as isize - 2 >= len as isize / 2 {
            Vec::new() // No overhang on this strand.
        } else {
            if len - cut < cut {
                eprintln!("Error with cut lens. len-cut: {}, Cut: {cut}", len - cut);
                return Vec::new();
            }

            seq_segment[cut..len - cut].to_vec()
        }
    }

    pub fn overhang_top_right(&self, seq_segment: &[Nucleotide]) -> Vec<Nucleotide> {
        let cut = self.cut_after as usize + 1;
        let len = self.cut_seq.len();

        if cut as isize - 2 < len as isize / 2 {
            Vec::new() // No overhang on this strand.
        } else {
            seq_segment[len - cut..cut].to_vec()
        }
    }

    pub fn overhang_bottom_left(&self, seq_segment: &[Nucleotide]) -> Vec<Nucleotide> {
        // todo: DRY implementation of reverse without compl
        let x = self.overhang_top_right(seq_segment);
        let mut result = x.to_vec();

        for nt in &mut result {
            *nt = nt.complement();
        }

        result
    }

    pub fn overhang_bottom_right(&self, seq_segment: &[Nucleotide]) -> Vec<Nucleotide> {
        // todo: DRY implementation of reverse without compl
        let x = self.overhang_top_left(seq_segment);
        let mut result = x.to_vec();

        for nt in &mut result {
            *nt = nt.complement();
        }

        result
    }
}

/// Go through a sequence, and attempt to match each enzyme in our RE library to the sequence.
/// Note/todo: We currently only search in the forward direction; this works if all enzymes in our
/// todo library are symmetric.
pub fn find_re_matches(seq: &[Nucleotide], lib: &[RestrictionEnzyme]) -> Vec<ReMatch> {
    let mut result = Vec::new();
    let seq_len = seq.len();

    let mut match_counts = HashMap::new(); // lib index, count

    for (lib_index, re) in lib.iter().enumerate() {
        let re_seq_len = re.cut_seq.len();

        for i in 0..seq_len {
            if i + re_seq_len + 1 >= seq_len {
                continue;
            }

            let mut matches = true;
            // If the RE cut site doesn't match this sequence segment, continue.
            for (j, nt) in seq[i..i + re_seq_len].iter().enumerate() {
                if !re.cut_seq[j].matches(*nt) {
                    matches = false;
                    break;
                }
            }

            if !matches {
                continue;
            }

            result.push(ReMatch {
                lib_index,
                // direction: PrimerDirection::Forward,
                seq_index: i + 1, // +1 indexing.
                match_count: 0,   // Updated below.
            });

            if let Entry::Vacant(e) = match_counts.entry(lib_index) {
                e.insert(1);
            } else {
                *match_counts.get_mut(&lib_index).unwrap() += 1;
            }
        }
    }

    // Apply match counts.
    for re_match in &mut result {
        re_match.match_count = match_counts[&re_match.lib_index];
    }

    result
}

/// Convert a nucleotide sequence to string.
pub fn seq_general_to_str(seq: &[NucleotideGeneral]) -> String {
    let mut result = String::new();

    for nt in seq {
        result.push_str(&nt.to_str_upper());
    }

    result
}
