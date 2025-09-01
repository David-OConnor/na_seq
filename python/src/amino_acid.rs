use std::str::FromStr;

use na_seq_rs::{
    AaCategory as RsAaCategory, AaIdent as RsAaIdent, AminoAcid as RsAminoAcid,
    AminoAcidGeneral as RsAminoAcidGeneral,
    AminoAcidProtenationVariant as RsAminoAcidProtenationVariant, CodingResult as RsCodingResult,
    Nucleotide as RsNucleotide,
};
use pyo3::{prelude::*, types::PyType};

use crate::{Nucleotide, map_io};

#[pyclass(module = "na_seq")]
#[derive(Clone, Copy)]
pub struct AaIdent {
    pub inner: RsAaIdent,
}

#[pymethods]
impl AaIdent {
    #[classmethod]
    fn from_str(_cls: &Bound<PyType>, s: &str) -> PyResult<Self> {
        let val = match s {
            "OneLetter" | "one" | "1" => RsAaIdent::OneLetter,
            "ThreeLetters" | "three" | "3" => RsAaIdent::ThreeLetters,
            _ => {
                return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                    "invalid AaIdent",
                ));
            }
        };
        Ok(Self { inner: val })
    }

    fn __str__(&self) -> &'static str {
        match self.inner {
            RsAaIdent::OneLetter => "OneLetter",
            RsAaIdent::ThreeLetters => "ThreeLetters",
        }
    }
    fn __repr__(&self) -> String {
        format!("AaIdent({})", self.__str__())
    }
}

#[pyclass(module = "na_seq")]
#[derive(Clone, Copy)]
pub struct AaCategory {
    pub inner: RsAaCategory,
}

#[pymethods]
impl AaCategory {
    fn __str__(&self) -> &'static str {
        match self.inner {
            RsAaCategory::Hydrophobic => "Hydrophobic",
            RsAaCategory::Acidic => "Acidic",
            RsAaCategory::Basic => "Basic",
            RsAaCategory::Polar => "Polar",
        }
    }
    fn __repr__(&self) -> String {
        format!("AaCategory({})", self.__str__())
    }
}

#[pyclass(module = "na_seq")]
#[derive(Clone, Copy)]
pub struct AminoAcid {
    pub inner: RsAminoAcid,
}

#[pymethods]
impl AminoAcid {
    #[classmethod]
    fn from_str(_cls: &Bound<PyType>, s: &str) -> PyResult<Self> {
        Ok(Self {
            inner: map_io(RsAminoAcid::from_str(s))?,
        })
    }

    fn to_str(&self, ident: &AaIdent) -> String {
        self.inner.to_str(ident.inner)
    }
    fn to_u8_upper(&self) -> u8 {
        self.inner.to_u8_upper()
    }
    fn to_u8_lower(&self) -> u8 {
        self.inner.to_u8_lower()
    }
    fn to_str_offset(&self) -> String {
        self.inner.to_str_offset()
    }
    fn weight(&self) -> f32 {
        self.inner.weight()
    }
    fn hydropathicity(&self) -> f32 {
        self.inner.hydropathicity()
    }

    fn category(&self) -> AaCategory {
        AaCategory {
            inner: self.inner.category(),
        }
    }

    fn codons(&self) -> Vec<Vec<Nucleotide>> {
        self.inner
            .codons()
            .into_iter()
            .map(|row| {
                row.into_iter()
                    .map(|nt: RsNucleotide| Nucleotide { inner: nt })
                    .collect()
            })
            .collect()
    }

    #[classmethod]
    fn from_codons(
        _cls: &Bound<PyType>,
        a: &Nucleotide,
        b: &Nucleotide,
        c: &Nucleotide,
    ) -> CodingResult {
        CodingResult {
            inner: RsAminoAcid::from_codons([a.inner, b.inner, c.inner]),
        }
    }

    fn __str__(&self) -> String {
        self.inner.to_string()
    }
    fn __repr__(&self) -> String {
        format!("AminoAcid({})", self.inner.to_string())
    }
}

#[pyclass(module = "na_seq")]
#[derive(Clone, Copy)]
pub struct CodingResult {
    pub inner: RsCodingResult,
}

#[pymethods]
impl CodingResult {
    fn is_stop(&self) -> bool {
        matches!(self.inner, RsCodingResult::StopCodon)
    }

    fn amino_acid(&self) -> Option<AminoAcid> {
        match self.inner {
            RsCodingResult::AminoAcid(aa) => Some(AminoAcid { inner: aa }),
            RsCodingResult::StopCodon => None,
        }
    }

    fn __str__(&self) -> String {
        match self.inner {
            RsCodingResult::AminoAcid(aa) => format!("{}", aa.to_string()),
            RsCodingResult::StopCodon => "StopCodon".to_string(),
        }
    }
    fn __repr__(&self) -> String {
        match self.inner {
            RsCodingResult::AminoAcid(aa) => format!("CodingResult(AminoAcid({}))", aa.to_string()),
            RsCodingResult::StopCodon => "CodingResult(StopCodon)".to_string(),
        }
    }
}

#[pyclass(module = "na_seq")]
#[derive(Clone, Copy)]
pub struct AminoAcidProtenationVariant {
    pub inner: RsAminoAcidProtenationVariant,
}

#[pymethods]
impl AminoAcidProtenationVariant {
    #[classmethod]
    fn from_str(_cls: &Bound<PyType>, s: &str) -> PyResult<Self> {
        Ok(Self {
            inner: map_io(RsAminoAcidProtenationVariant::from_str(s))?,
        })
    }

    fn get_standard(&self) -> Option<AminoAcid> {
        self.inner.get_standard().map(|aa| AminoAcid { inner: aa })
    }

    fn __str__(&self) -> String {
        self.inner.to_string()
    }
    fn __repr__(&self) -> String {
        format!("AminoAcidProtenationVariant({})", self.inner.to_string())
    }
}

#[pyclass(module = "na_seq")]
#[derive(Clone, Copy)]
pub struct AminoAcidGeneral {
    pub inner: RsAminoAcidGeneral,
}

#[pymethods]
impl AminoAcidGeneral {
    #[classmethod]
    fn from_str(_cls: &Bound<PyType>, s: &str) -> PyResult<Self> {
        Ok(Self {
            inner: map_io(RsAminoAcidGeneral::from_str(s))?,
        })
    }

    #[classmethod]
    fn from_standard(_cls: &Bound<PyType>, aa: &AminoAcid) -> Self {
        Self {
            inner: RsAminoAcidGeneral::Standard(aa.inner),
        }
    }

    #[classmethod]
    fn from_variant(_cls: &Bound<PyType>, v: &AminoAcidProtenationVariant) -> Self {
        Self {
            inner: RsAminoAcidGeneral::Variant(v.inner),
        }
    }

    fn is_standard(&self) -> bool {
        matches!(self.inner, RsAminoAcidGeneral::Standard(_))
    }

    fn is_variant(&self) -> bool {
        matches!(self.inner, RsAminoAcidGeneral::Variant(_))
    }

    fn to_standard(&self) -> Option<AminoAcid> {
        match self.inner {
            RsAminoAcidGeneral::Standard(aa) => Some(AminoAcid { inner: aa }),
            RsAminoAcidGeneral::Variant(v) => v.get_standard().map(|aa| AminoAcid { inner: aa }),
        }
    }

    fn __repr__(&self) -> String {
        match self.inner {
            RsAminoAcidGeneral::Standard(aa) => {
                format!("AminoAcidGeneral(Standard({}))", aa.to_string())
            }
            RsAminoAcidGeneral::Variant(v) => {
                format!("AminoAcidGeneral(Variant({}))", v.to_string())
            }
        }
    }
}
