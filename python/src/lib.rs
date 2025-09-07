mod amino_acids;
pub mod element;
pub mod nucleotide;

use element::*;
use na_seq_rs::{self, AtomTypeInRes as RsAtomTypeInRes};
use nucleotide::*;
use pyo3::{
    Bound, Py, PyResult, Python, prelude::*, pymodule, types::PyType,
};

use crate::{
    amino_acids::{
        AaCategory, AaIdent, AminoAcid, AminoAcidGeneral, AminoAcidProtenationVariant, CodingResult,
    },
    nucleotide::Nucleotide,
};

/// Candidate for standalone helper lib.
#[macro_export]
macro_rules! make_enum {
    ($Py:ident, $Native:path, $( $Var:ident ),+ $(,)?) => {
        #[pyclass]
        #[derive(Clone, Copy, Debug, PartialEq, Eq)]
        pub enum $Py { $( $Var ),+ }

        impl ::core::convert::From<$Py> for $Native {
            fn from(v: $Py) -> Self { match v { $( $Py::$Var => <$Native>::$Var ),+ } }
        }

        impl ::core::convert::From<$Native> for $Py {
            fn from(v: $Native) -> Self { match v { $( <$Native>::$Var => $Py::$Var ),+ } }
        }

        impl $Py {
            pub fn to_native(self) -> $Native {
                self.into()
            }

            pub fn from_native(native: $Native) -> Self {
               native.into()
            }
        }
    };
}

/// Candidate for standalone helper lib.
macro_rules! field {
    ($name:ident, $ty:ty) => {
        #[getter]
        fn $name(&self) -> $ty {
            self.inner.$name.into()
        }

        #[setter($name)]
        fn $name##_set(&mut self, val: $ty) -> $ty {
            self.inner.$name = val.into();
            val
        }
    };
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

#[pyfunction]
fn seq_complement(seq: Vec<Nucleotide>) -> Vec<Nucleotide> {
    let seq_native: Vec<_> = seq.iter().map(|n| n.to_native()).collect();
    na_seq_rs::seq_complement(&seq_native)
        .into_iter()
        .map(|n| Nucleotide::from_native(n))
        .collect()
}

#[pyfunction]
fn seq_from_str(str: &str) -> Vec<Nucleotide> {
    na_seq_rs::seq_from_str(str)
        .into_iter()
        .map(|n| Nucleotide::from_native(n))
        .collect()
}

#[pyfunction]
fn seq_aa_from_str(str: &str) -> Vec<AminoAcid> {
    na_seq_rs::seq_aa_from_str(str)
        .into_iter()
        .map(|n| AminoAcid::from_native(n))
        .collect()
}

#[pyfunction]
fn seq_to_str_lower(seq: Vec<Nucleotide>) -> String {
    let seq_native: Vec<_> = seq.iter().map(|n| n.to_native()).collect();
    na_seq_rs::seq_to_str_lower(&seq_native)
}

#[pyfunction]
fn seq_to_str_upper(seq: Vec<Nucleotide>) -> String {
    let seq_native: Vec<_> = seq.iter().map(|n| n.to_native()).collect();
    na_seq_rs::seq_to_str_upper(&seq_native)
}

#[pyfunction]
fn seq_to_u8_lower(seq: Vec<Nucleotide>) -> Vec<u8> {
    let seq_native: Vec<_> = seq.iter().map(|n| n.to_native()).collect();
    na_seq_rs::seq_to_u8_lower(&seq_native)
}

#[pyfunction]
fn seq_to_u8_upper(seq: Vec<Nucleotide>) -> Vec<u8> {
    let seq_native: Vec<_> = seq.iter().map(|n| n.to_native()).collect();
    na_seq_rs::seq_to_u8_upper(&seq_native)
}

#[pyfunction]
fn seq_aa_to_str(seq: Vec<AminoAcid>) -> String {
    let seq_native: Vec<_> = seq.iter().map(|n| n.to_native()).collect();
    na_seq_rs::seq_aa_to_str(&seq_native)
}

#[pyfunction]
fn seq_aa_to_u8_lower(seq: Vec<AminoAcid>) -> Vec<u8> {
    let seq_native: Vec<_> = seq.iter().map(|n| n.to_native()).collect();
    na_seq_rs::seq_aa_to_u8_lower(&seq_native)
}

#[pyfunction]
fn seq_aa_to_u8_upper(seq: Vec<AminoAcid>) -> Vec<u8> {
    let seq_native: Vec<_> = seq.iter().map(|n| n.to_native()).collect();
    na_seq_rs::seq_aa_to_u8_upper(&seq_native)
}

#[pyfunction]
fn seq_weight(seq: Vec<Nucleotide>) -> f32 {
    let seq_native: Vec<_> = seq.iter().map(|n| n.to_native()).collect();
    na_seq_rs::seq_weight(&seq_native)
}

#[pyfunction]
fn calc_gc(seq: Vec<Nucleotide>) -> f32 {
    let seq_native: Vec<_> = seq.iter().map(|n| n.to_native()).collect();
    na_seq_rs::calc_gc(&seq_native)
}

#[pyfunction]
fn serialize_seq_bin(seq: Vec<Nucleotide>) -> Vec<u8> {
    let seq_native: Vec<_> = seq.iter().map(|n| n.to_native()).collect();
    na_seq_rs::serialize_seq_bin(&seq_native)
}

#[pyfunction]
fn deser_seq_bin(data: Vec<u8>) -> PyResult<Vec<Nucleotide>> {
    let result = na_seq_rs::deser_seq_bin(&data)?;
    Ok(result
        .into_iter()
        .map(|n| Nucleotide::from_native(n))
        .collect())
}

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

    m.add_function(wrap_pyfunction!(seq_complement, m)?)?;

    m.add_function(wrap_pyfunction!(seq_from_str, m)?)?;
    m.add_function(wrap_pyfunction!(seq_aa_from_str, m)?)?;

    m.add_function(wrap_pyfunction!(seq_to_str_lower, m)?)?;
    m.add_function(wrap_pyfunction!(seq_to_str_upper, m)?)?;
    m.add_function(wrap_pyfunction!(seq_to_u8_lower, m)?)?;
    m.add_function(wrap_pyfunction!(seq_to_u8_upper, m)?)?;

    m.add_function(wrap_pyfunction!(seq_aa_to_str, m)?)?;
    m.add_function(wrap_pyfunction!(seq_aa_to_u8_lower, m)?)?;
    m.add_function(wrap_pyfunction!(seq_aa_to_u8_upper, m)?)?;
    m.add_function(wrap_pyfunction!(seq_weight, m)?)?;
    m.add_function(wrap_pyfunction!(calc_gc, m)?)?;

    m.add_function(wrap_pyfunction!(serialize_seq_bin, m)?)?;
    m.add_function(wrap_pyfunction!(deser_seq_bin, m)?)?;

    // We use these for Complex enums. Better way?
    let atir_obj = m.getattr("AtomTypeInRes")?;
    let atir_type = atir_obj.downcast::<PyType>()?;

    macro_rules! at {
        ($v:ident) => {
            set_variant!(
                py,
                atir_type,
                AtomTypeInRes,
                RsAtomTypeInRes,
                $v,
                RsAtomTypeInRes::$v.to_string()
            );
        };
    }

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

    Ok(())
}
