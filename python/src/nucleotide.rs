use std::str::FromStr;

use na_seq_rs;
use pyo3::{prelude::*, types::PyType};

use crate::make_enum;

make_enum!(Nucleotide, na_seq_rs::Nucleotide, T, C, A, G);

#[pymethods]
impl Nucleotide {
    #[classmethod]
    fn from_str(_cls: &Bound<PyType>, s: &str) -> PyResult<Self> {
        Ok(na_seq_rs::Nucleotide::from_str(s)?.into())
    }

    #[classmethod]
    fn from_u8_letter(_cls: &Bound<PyType>, val: u8) -> PyResult<Self> {
        Ok(na_seq_rs::Nucleotide::from_u8_letter(val)?.into())
    }

    fn to_u8_upper(&self) -> u8 {
        self.to_native().to_u8_upper()
    }
    fn to_u8_lower(&self) -> u8 {
        self.to_native().to_u8_lower()
    }
    fn to_str_upper(&self) -> String {
        self.to_native().to_str_upper()
    }
    fn to_str_lower(&self) -> String {
        self.to_native().to_str_lower()
    }

    fn complement(&self) -> Self {
        self.to_native().complement().into()
    }
    fn weight(&self) -> f32 {
        self.to_native().weight()
    }
    fn a_max(&self) -> f32 {
        self.to_native().a_max()
    }
    fn molar_density(&self) -> f32 {
        self.to_native().molar_density()
    }

    #[getter]
    fn value(&self) -> u8 {
        self.to_native() as u8
    }

    fn __str__(&self) -> String {
        self.to_native().to_string()
    }
    fn __repr__(&self) -> String {
        format!("{:?}", self.to_native())
    }
}

make_enum!(
    NucleotideGeneral,
    na_seq_rs::NucleotideGeneral,
    T,
    C,
    A,
    G,
    N,
    W,
    S,
    Y,
    R,
    M,
    K,
);

#[pymethods]
impl NucleotideGeneral {
    #[classmethod]
    fn from_str(_cls: &Bound<PyType>, s: &str) -> PyResult<Self> {
        Ok(na_seq_rs::NucleotideGeneral::from_str(s)?.into())
    }

    #[classmethod]
    fn from_u8_letter(_cls: &Bound<PyType>, val: u8) -> PyResult<Self> {
        Ok(na_seq_rs::NucleotideGeneral::from_u8_letter(val)?.into())
    }

    fn matches(&self, nt: &Nucleotide) -> bool {
        self.to_native().matches(nt.to_native())
    }

    fn to_u8_lower(&self) -> u8 {
        self.to_native().to_u8_lower()
    }
    fn to_u8_upper(&self) -> u8 {
        self.to_native().to_u8_upper()
    }
    fn to_str_lower(&self) -> String {
        self.to_native().to_str_lower()
    }
    fn to_str_upper(&self) -> String {
        self.to_native().to_str_upper()
    }

    #[getter]
    fn value(&self) -> u8 {
        self.to_native() as u8
    }

    fn __str__(&self) -> String {
        self.to_native().to_string()
    }
    fn __repr__(&self) -> String {
        format!("{:?}", self.to_native())
    }
}
