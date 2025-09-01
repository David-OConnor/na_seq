use std::str::FromStr;

use na_seq_rs;
use pyo3::{prelude::*, types::PyType};

use crate::{Nucleotide, map_io};

#[pyclass(module = "na_seq")]
#[derive(Clone, Copy)]
pub struct AaIdent {
    pub inner: na_seq_rs::AaIdent,
}

#[pymethods]
impl AaIdent {
    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
}

#[pyclass(module = "na_seq")]
#[derive(Clone, Copy)]
pub struct CodingResult {
    pub inner: na_seq_rs::CodingResult,
}

#[pymethods]
impl CodingResult {
    #[classmethod]
    fn from_codons(_cls: &Bound<PyType>, codons: [Nucleotide; 3]) -> Self {
        let codons_rs = codons.map(|c| c.inner);
        Self {
            inner: na_seq_rs::CodingResult::from_codons(codons_rs),
        }
    }

    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
}

#[pyclass(module = "na_seq")]
#[derive(Clone, Copy)]
pub struct AaCategory {
    pub inner: na_seq_rs::AaCategory,
}

#[pymethods]
impl AaCategory {
    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
}

#[pyclass(module = "na_seq")]
#[derive(Clone, Copy)]
pub struct AminoAcid {
    pub inner: na_seq_rs::AminoAcid,
}

#[pymethods]
impl AminoAcid {
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

    #[classmethod]
    fn from_str(_cls: &Bound<PyType>, s: &str) -> PyResult<Self> {
        Ok(Self {
            inner: map_io(na_seq_rs::AminoAcid::from_str(s))?,
        })
    }

    fn weight(&self) -> f32 {
        self.inner.weight()
    }
    fn hydropathicity(&self) -> f32 {
        self.inner.hydropathicity()
    }

    fn codons(&self) -> Vec<Vec<Nucleotide>> {
        self.inner
            .codons()
            .into_iter()
            .map(|row| {
                row.into_iter()
                    .map(|nt: na_seq_rs::Nucleotide| Nucleotide { inner: nt })
                    .collect()
            })
            .collect()
    }

    #[classmethod]
    fn from_codons(_cls: &Bound<PyType>, codons: [Nucleotide; 3]) -> Option<Self> {
        let codons_rs = codons.map(|c| c.inner);
        match na_seq_rs::AminoAcid::from_codons(codons_rs) {
            Some(aa) => Some(Self { inner: aa }),
            None => None,
        }
    }

    fn category(&self) -> AaCategory {
        AaCategory {
            inner: self.inner.category(),
        }
    }

    fn __str__(&self) -> String {
        self.inner.to_string()
    }
    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
}

#[pyclass(module = "na_seq")]
#[derive(Clone, Copy)]
pub struct AminoAcidProtenationVariant {
    pub inner: na_seq_rs::AminoAcidProtenationVariant,
}

#[pymethods]
impl AminoAcidProtenationVariant {
    #[classmethod]
    fn from_str(_cls: &Bound<PyType>, s: &str) -> PyResult<Self> {
        Ok(Self {
            inner: map_io(na_seq_rs::AminoAcidProtenationVariant::from_str(s))?,
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
    pub inner: na_seq_rs::AminoAcidGeneral,
}

#[pymethods]
impl AminoAcidGeneral {
    #[classmethod]
    fn from_str(_cls: &Bound<PyType>, s: &str) -> PyResult<Self> {
        Ok(Self {
            inner: map_io(na_seq_rs::AminoAcidGeneral::from_str(s))?,
        })
    }

    #[classmethod]
    fn from_standard(_cls: &Bound<PyType>, aa: &AminoAcid) -> Self {
        Self {
            inner: na_seq_rs::AminoAcidGeneral::Standard(aa.inner),
        }
    }

    #[classmethod]
    fn from_variant(_cls: &Bound<PyType>, v: &AminoAcidProtenationVariant) -> Self {
        Self {
            inner: na_seq_rs::AminoAcidGeneral::Variant(v.inner),
        }
    }

    fn is_standard(&self) -> bool {
        matches!(self.inner, na_seq_rs::AminoAcidGeneral::Standard(_))
    }

    fn is_variant(&self) -> bool {
        matches!(self.inner, na_seq_rs::AminoAcidGeneral::Variant(_))
    }

    fn to_standard(&self) -> Option<AminoAcid> {
        match self.inner {
            na_seq_rs::AminoAcidGeneral::Standard(aa) => Some(AminoAcid { inner: aa }),
            na_seq_rs::AminoAcidGeneral::Variant(v) => {
                v.get_standard().map(|aa| AminoAcid { inner: aa })
            }
        }
    }

    fn __repr__(&self) -> String {
        match self.inner {
            na_seq_rs::AminoAcidGeneral::Standard(aa) => {
                format!("AminoAcidGeneral(Standard({}))", aa.to_string())
            }
            na_seq_rs::AminoAcidGeneral::Variant(v) => {
                format!("AminoAcidGeneral(Variant({}))", v.to_string())
            }
        }
    }
}
