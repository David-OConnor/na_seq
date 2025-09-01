use std::str::FromStr;

use na_seq_rs;
use pyo3::{prelude::*, types::PyType};

use crate::map_io;

#[pyclass(module = "na_seq")]
#[derive(Clone, Copy)]
pub struct Nucleotide {
    pub inner: na_seq_rs::Nucleotide,
}

#[pymethods]
impl Nucleotide {
    #[classmethod]
    fn from_str(_cls: &Bound<PyType>, s: &str) -> PyResult<Self> {
        Ok(Self {
            inner: map_io(na_seq_rs::Nucleotide::from_str(s))?,
        })
    }

    #[classmethod]
    fn from_u8_letter(_cls: &Bound<PyType>, val: u8) -> PyResult<Self> {
        Ok(Self {
            inner: map_io(na_seq_rs::Nucleotide::from_u8_letter(val))?,
        })
    }

    fn to_u8_upper(&self) -> u8 {
        self.inner.to_u8_upper()
    }
    fn to_u8_lower(&self) -> u8 {
        self.inner.to_u8_lower()
    }
    fn to_str_upper(&self) -> String {
        self.inner.to_str_upper()
    }
    fn to_str_lower(&self) -> String {
        self.inner.to_str_lower()
    }

    fn complement(&self) -> Self {
        Self {
            inner: self.inner.complement(),
        }
    }
    fn weight(&self) -> f32 {
        self.inner.weight()
    }
    fn a_max(&self) -> f32 {
        self.inner.a_max()
    }
    fn molar_density(&self) -> f32 {
        self.inner.molar_density()
    }

    #[getter]
    fn value(&self) -> u8 {
        self.inner as u8
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
pub struct NucleotideGeneral {
    pub inner: na_seq_rs::NucleotideGeneral,
}

#[pymethods]
impl NucleotideGeneral {
    #[classmethod]
    fn from_str(_cls: &Bound<PyType>, s: &str) -> PyResult<Self> {
        Ok(Self {
            inner: map_io(na_seq_rs::NucleotideGeneral::from_str(s))?,
        })
    }

    #[classmethod]
    fn from_u8_letter(_cls: &Bound<PyType>, val: u8) -> PyResult<Self> {
        Ok(Self {
            inner: map_io(na_seq_rs::NucleotideGeneral::from_u8_letter(val))?,
        })
    }

    fn matches(&self, nt: &Nucleotide) -> bool {
        self.inner.matches(nt.inner)
    }

    fn to_u8_lower(&self) -> u8 {
        self.inner.to_u8_lower()
    }
    fn to_u8_upper(&self) -> u8 {
        self.inner.to_u8_upper()
    }
    fn to_str_lower(&self) -> String {
        self.inner.to_str_lower()
    }
    fn to_str_upper(&self) -> String {
        self.inner.to_str_upper()
    }

    #[getter]
    fn value(&self) -> u8 {
        self.inner as u8
    }

    fn __str__(&self) -> String {
        self.inner.to_string()
    }
    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
}
