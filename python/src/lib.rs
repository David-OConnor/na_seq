mod amino_acids;
pub mod element;
pub mod nucleotide;

use element::*;
use na_seq_rs::{
    AaCategory as RsAaCategory, AaIdent as RsAaIdent, AminoAcid as RsAminoAcid,
    AminoAcidProtenationVariant as RsAminoAcidProtenationVariant, AtomTypeInRes as RsAtomTypeInRes,
    Element as RsElement, Nucleotide as RsNucleotide, NucleotideGeneral as RsNucleotideGeneral,
};
use nucleotide::*;
use pyo3::{
    Bound, Py, PyResult, Python, exceptions::PyValueError, prelude::*, pymodule, types::PyType,
};

use crate::amino_acids::{
    AaCategory, AaIdent, AminoAcid, AminoAcidGeneral, AminoAcidProtenationVariant, CodingResult,
};

fn map_io<T>(r: std::io::Result<T>) -> PyResult<T> {
    r.map_err(|e| PyValueError::new_err(e.to_string()))
}

macro_rules! set_variant {
    ($py:expr, $class:expr, $Wrapper:ident, $RsEnum:ident, $name:ident, $alias:expr) => {{
        let obj = Py::new(
            $py,
            $Wrapper {
                inner: $RsEnum::$name,
            },
        )?;
        $class.setattr(stringify!($name), obj.clone_ref($py))?;
        let sym: String = $alias;
        $class.setattr(&sym, obj)?;
    }};
}

#[rustfmt::skip]
#[pymodule]
fn na_seq(py: Python<'_>, m: &Bound<PyModule>) -> PyResult<()> {
    m.add_class::<Nucleotide>()?;
    m.add_class::<NucleotideGeneral>()?;
    m.add_class::<Element>()?;
    m.add_class::<AtomTypeInRes>()?;
    m.add_class::<AaIdent>()?;
    m.add_class::<AaCategory>()?;
    m.add_class::<AminoAcid>()?;
    m.add_class::<CodingResult>()?;
    m.add_class::<AminoAcidProtenationVariant>()?;
    m.add_class::<AminoAcidGeneral>()?;

    // Keep each intermediate alive to satisfy &Bound lifetimes
    let et_obj = m.getattr("Element")?;
    let et = et_obj.downcast::<PyType>()?;

    let nt_obj = m.getattr("Nucleotide")?;
    let nt_type = nt_obj.downcast::<PyType>()?;

    let ntg_obj = m.getattr("NucleotideGeneral")?;
    let ntg_type = ntg_obj.downcast::<PyType>()?;

    let atir_obj = m.getattr("AtomTypeInRes")?;
    let atir_type = atir_obj.downcast::<PyType>()?;

    let aa_ident_obj = m.getattr("AaIdent")?;
    let aa_ident = aa_ident_obj.downcast::<PyType>()?;

    let aa_cat_obj = m.getattr("AaCategory")?;
    let aa_cat = aa_cat_obj.downcast::<PyType>()?;

    let aa_obj = m.getattr("AminoAcid")?;
    let aa = aa_obj.downcast::<PyType>()?;

    let aap_obj = m.getattr("AminoAcidProtenationVariant")?;
    let aap = aap_obj.downcast::<PyType>()?;


    macro_rules! at {
        ($v:ident) => {
            set_variant!(
                py, atir_type, AtomTypeInRes, RsAtomTypeInRes, $v,
                RsAtomTypeInRes::$v.to_string()
            );
        };
    }

    // Element + symbol aliases
    set_variant!(py, et, Element, RsElement, Hydrogen,   RsElement::Hydrogen.to_letter());
    set_variant!(py, et, Element, RsElement, Carbon,     RsElement::Carbon.to_letter());
    set_variant!(py, et, Element, RsElement, Oxygen,     RsElement::Oxygen.to_letter());
    set_variant!(py, et, Element, RsElement, Nitrogen,   RsElement::Nitrogen.to_letter());
    set_variant!(py, et, Element, RsElement, Fluorine,   RsElement::Fluorine.to_letter());
    set_variant!(py, et, Element, RsElement, Sulfur,     RsElement::Sulfur.to_letter());
    set_variant!(py, et, Element, RsElement, Phosphorus, RsElement::Phosphorus.to_letter());
    set_variant!(py, et, Element, RsElement, Iron,       RsElement::Iron.to_letter());
    set_variant!(py, et, Element, RsElement, Copper,     RsElement::Copper.to_letter());
    set_variant!(py, et, Element, RsElement, Calcium,    RsElement::Calcium.to_letter());
    set_variant!(py, et, Element, RsElement, Potassium,  RsElement::Potassium.to_letter());
    set_variant!(py, et, Element, RsElement, Aluminum,   RsElement::Aluminum.to_letter());
    set_variant!(py, et, Element, RsElement, Lead,       RsElement::Lead.to_letter());
    set_variant!(py, et, Element, RsElement, Gold,       RsElement::Gold.to_letter());
    set_variant!(py, et, Element, RsElement, Silver,     RsElement::Silver.to_letter());
    set_variant!(py, et, Element, RsElement, Mercury,    RsElement::Mercury.to_letter());
    set_variant!(py, et, Element, RsElement, Tin,        RsElement::Tin.to_letter());
    set_variant!(py, et, Element, RsElement, Zinc,       RsElement::Zinc.to_letter());
    set_variant!(py, et, Element, RsElement, Magnesium,  RsElement::Magnesium.to_letter());
    set_variant!(py, et, Element, RsElement, Manganese,  RsElement::Manganese.to_letter());
    set_variant!(py, et, Element, RsElement, Iodine,     RsElement::Iodine.to_letter());
    set_variant!(py, et, Element, RsElement, Chlorine,   RsElement::Chlorine.to_letter());
    set_variant!(py, et, Element, RsElement, Tungsten,   RsElement::Tungsten.to_letter());
    set_variant!(py, et, Element, RsElement, Tellurium,  RsElement::Tellurium.to_letter());
    set_variant!(py, et, Element, RsElement, Selenium,   RsElement::Selenium.to_letter());
    set_variant!(py, et, Element, RsElement, Bromine,    RsElement::Bromine.to_letter());
    set_variant!(py, et, Element, RsElement, Rubidium,   RsElement::Rubidium.to_letter());
    set_variant!(py, et, Element, RsElement, Other,      RsElement::Other.to_letter());


   // Nucleotide + alias "A","C","G","T"
    set_variant!(py, nt_type, Nucleotide, RsNucleotide, A, RsNucleotide::A.to_str_upper());
    set_variant!(py, nt_type, Nucleotide, RsNucleotide, C, RsNucleotide::C.to_str_upper());
    set_variant!(py, nt_type, Nucleotide, RsNucleotide, G, RsNucleotide::G.to_str_upper());
    set_variant!(py, nt_type, Nucleotide, RsNucleotide, T, RsNucleotide::T.to_str_upper());

  // NucleotideGeneral + alias "A","C","G","T","N","W","S","Y","R","M","K"
    set_variant!(py, ntg_type, NucleotideGeneral, RsNucleotideGeneral, A, RsNucleotideGeneral::A.to_str_upper());
    set_variant!(py, ntg_type, NucleotideGeneral, RsNucleotideGeneral, C, RsNucleotideGeneral::C.to_str_upper());
    set_variant!(py, ntg_type, NucleotideGeneral, RsNucleotideGeneral, G, RsNucleotideGeneral::G.to_str_upper());
    set_variant!(py, ntg_type, NucleotideGeneral, RsNucleotideGeneral, T, RsNucleotideGeneral::T.to_str_upper());
    set_variant!(py, ntg_type, NucleotideGeneral, RsNucleotideGeneral, N, RsNucleotideGeneral::N.to_str_upper());
    set_variant!(py, ntg_type, NucleotideGeneral, RsNucleotideGeneral, W, RsNucleotideGeneral::W.to_str_upper());
    set_variant!(py, ntg_type, NucleotideGeneral, RsNucleotideGeneral, S, RsNucleotideGeneral::S.to_str_upper());
    set_variant!(py, ntg_type, NucleotideGeneral, RsNucleotideGeneral, Y, RsNucleotideGeneral::Y.to_str_upper());
    set_variant!(py, ntg_type, NucleotideGeneral, RsNucleotideGeneral, R, RsNucleotideGeneral::R.to_str_upper());
    set_variant!(py, ntg_type, NucleotideGeneral, RsNucleotideGeneral, M, RsNucleotideGeneral::M.to_str_upper());
    set_variant!(py, ntg_type, NucleotideGeneral, RsNucleotideGeneral, K, RsNucleotideGeneral::K.to_str_upper());

    at!(C);
    at!(CA);
    at!(CB);
    at!(CD);
    at!(CD1);
    at!(CD2);
    at!(CE);
    at!(CE1);
    at!(CE2);
    at!(CE3);
    at!(CG);
    at!(CG1);
    at!(CG2);
    at!(CH2);
    at!(CH3);
    at!(CZ);
    at!(CZ1);
    at!(CZ2);
    at!(CZ3);
    at!(O);
    at!(OD1);
    at!(OD2);
    at!(OE1);
    at!(OE2);
    at!(OG);
    at!(OH);
    at!(OXT);
    at!(N);
    at!(ND1);
    at!(ND2);
    at!(NE);
    at!(NZ);
    at!(NH1);
    at!(NH2);
    at!(NE1);
    at!(NE2);
    at!(OG1);
    at!(OG2);
    at!(SD);
    at!(SE);
    at!(SG);

    // AaIdent
    {
        let obj = Py::new(
            py,
            AaIdent {
                inner: RsAaIdent::OneLetter,
            },
        )?;
        aa_ident.setattr("OneLetter", obj.clone_ref(py))?;
        let obj2 = Py::new(
            py,
            AaIdent {
                inner: RsAaIdent::ThreeLetters,
            },
        )?;
        aa_ident.setattr("ThreeLetters", obj2)?;
    }

    // AaCategory
    macro_rules! set_cat {
        ($v:ident) => {
            aa_cat.setattr(
                stringify!($v),
                Py::new(
                    py,
                    AaCategory {
                        inner: RsAaCategory::$v,
                    },
                )?,
            )?;
        };
    }
    set_cat!(Hydrophobic);
    set_cat!(Acidic);
    set_cat!(Basic);
    set_cat!(Polar);

    // AminoAcid + one-letter aliases
    macro_rules! set_aa {
        ($v:ident) => {
            let obj = Py::new(
                py,
                AminoAcid {
                    inner: RsAminoAcid::$v,
                },
            )?;
            aa.setattr(stringify!($v), obj.clone_ref(py))?;
            let letter = RsAminoAcid::$v.to_str(RsAaIdent::OneLetter);
            aa.setattr(&letter, obj)?;
        };
    }
    set_aa!(Arg);
    set_aa!(His);
    set_aa!(Lys);
    set_aa!(Asp);
    set_aa!(Glu);
    set_aa!(Ser);
    set_aa!(Thr);
    set_aa!(Asn);
    set_aa!(Gln);
    set_aa!(Cys);
    set_aa!(Sec);
    set_aa!(Gly);
    set_aa!(Pro);
    set_aa!(Ala);
    set_aa!(Val);
    set_aa!(Ile);
    set_aa!(Leu);
    set_aa!(Met);
    set_aa!(Phe);
    set_aa!(Tyr);
    set_aa!(Trp);

    // AminoAcidProtenationVariant + 3-letter aliases
    macro_rules! set_aap {
        ($v:ident) => {
            let obj = Py::new(
                py,
                AminoAcidProtenationVariant {
                    inner: RsAminoAcidProtenationVariant::$v,
                },
            )?;
            aap.setattr(stringify!($v), obj.clone_ref(py))?;
            aap.setattr(&RsAminoAcidProtenationVariant::$v.to_string(), obj)?;
        };
    }
    set_aap!(Hid);
    set_aap!(Hie);
    set_aap!(Hip);
    set_aap!(Cym);
    set_aap!(Cyx);
    set_aap!(Ash);
    set_aap!(Glh);
    set_aap!(Lyn);
    set_aap!(Ace);
    set_aap!(Nhe);
    set_aap!(Nme);
    set_aap!(Hyp);

    Ok(())
}
