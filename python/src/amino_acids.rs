use std::str::FromStr;

use na_seq_rs;
use pyo3::{prelude::*, types::PyType};

use crate::{Nucleotide, make_enum};

make_enum!(AaIdent, na_seq_rs::AaIdent, OneLetter, ThreeLetters);
#[pymethods]
impl AaIdent {
    fn __repr__(&self) -> String {
        format!("{:?}", self.to_native())
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
        let codons_rs = codons.map(|c| c.into());
        Self {
            inner: na_seq_rs::CodingResult::from_codons(codons_rs),
        }
    }

    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
}

make_enum!(
    AaCategory,
    na_seq_rs::AaCategory,
    Hydrophobic,
    Acidic,
    Basic,
    Polar
);

#[pymethods]
impl AaCategory {
    fn __repr__(&self) -> String {
        format!("{:?}", self.to_native())
    }
}

make_enum!(
    AminoAcid,
    na_seq_rs::AminoAcid,
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
);

#[pymethods]
impl AminoAcid {
    fn to_str(&self, ident: &AaIdent) -> String {
        self.to_native().to_str(ident.to_native())
    }

    fn to_u8_upper(&self) -> u8 {
        self.to_native().to_u8_upper()
    }
    fn to_u8_lower(&self) -> u8 {
        self.to_native().to_u8_lower()
    }
    fn to_str_offset(&self) -> String {
        self.to_native().to_str_offset()
    }

    #[classmethod]
    fn from_str(_cls: &Bound<PyType>, s: &str) -> PyResult<Self> {
        Ok(na_seq_rs::AminoAcid::from_str(s)?.into())
    }

    fn weight(&self) -> f32 {
        self.to_native().weight()
    }
    fn hydropathicity(&self) -> f32 {
        self.to_native().hydropathicity()
    }

    fn codons(&self) -> Vec<Vec<Nucleotide>> {
        self.to_native()
            .codons()
            .into_iter()
            .map(|row| row.into_iter().map(|nt| nt.into()).collect())
            .collect()
    }

    #[classmethod]
    fn from_codons(_cls: &Bound<PyType>, codons: [Nucleotide; 3]) -> Option<Self> {
        let codons_rs = codons.map(|c| c.to_native());
        match na_seq_rs::AminoAcid::from_codons(codons_rs) {
            Some(aa) => Some(aa.into()),
            None => None,
        }
    }

    fn category(&self) -> AaCategory {
        self.to_native().category().into()
    }

    fn __str__(&self) -> String {
        self.to_native().to_string()
    }
    fn __repr__(&self) -> String {
        format!("{:?}", self.to_native())
    }
}

make_enum!(
    AminoAcidProtenationVariant,
    na_seq_rs::AminoAcidProtenationVariant,
    Hid,
    Hie,
    Hip,
    Cym,
    Cyx,
    Ash,
    Glh,
    Lyn,
    Ace,
    Nhe,
    Nme,
    Hyp,
);

#[pymethods]
impl AminoAcidProtenationVariant {
    #[classmethod]
    fn from_str(_cls: &Bound<PyType>, s: &str) -> PyResult<Self> {
        Ok(na_seq_rs::AminoAcidProtenationVariant::from_str(s)?.into())
    }

    fn get_standard(&self) -> Option<AminoAcid> {
        self.to_native().get_standard().map(|aa| aa.into())
    }

    fn __str__(&self) -> String {
        self.to_native().to_string()
    }
    fn __repr__(&self) -> String {
        format!("{:?}", self.to_native().to_string())
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
            inner: na_seq_rs::AminoAcidGeneral::from_str(s)?,
        })
    }

    #[classmethod]
    fn from_standard(_cls: &Bound<PyType>, aa: &AminoAcid) -> Self {
        Self {
            inner: na_seq_rs::AminoAcidGeneral::Standard(aa.to_native()).into(),
        }
    }

    #[classmethod]
    fn from_variant(_cls: &Bound<PyType>, v: &AminoAcidProtenationVariant) -> Self {
        Self {
            inner: na_seq_rs::AminoAcidGeneral::Variant(v.to_native()),
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
            na_seq_rs::AminoAcidGeneral::Standard(aa) => Some(aa.into()),
            na_seq_rs::AminoAcidGeneral::Variant(v) => v.get_standard().map(|aa| aa.into()),
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
