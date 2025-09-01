use std::str::FromStr;

use na_seq_rs::{Nucleotide as RsNucleotide, NucleotideGeneral as RsNucleotideGeneral};
use pyo3::{prelude::*, types::PyType};

use crate::map_io;

#[pyclass(module = "na_seq")]
#[derive(Clone, Copy)]
pub struct Nucleotide {
    pub inner: RsNucleotide,
}

#[pymethods]
impl Nucleotide {
    #[new]
    fn new_from_letter(val: u8) -> PyResult<Self> {
        Ok(Self {
            inner: map_io(RsNucleotide::from_u8_letter(val))?,
        })
    }

    #[classmethod]
    fn from_u8_letter(_cls: &Bound<PyType>, val: u8) -> PyResult<Self> {
        Ok(Self {
            inner: map_io(RsNucleotide::from_u8_letter(val))?,
        })
    }

    #[classmethod]
    fn from_str(_cls: &Bound<PyType>, s: &str) -> PyResult<Self> {
        Ok(Self {
            inner: map_io(RsNucleotide::from_str(s))?,
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
        self.inner.to_str_upper()
    }
    fn __repr__(&self) -> String {
        format!("Nucleotide({})", self.inner.to_str_upper())
    }
}

#[pyclass(module = "na_seq")]
#[derive(Clone, Copy)]
pub struct NucleotideGeneral {
    pub inner: RsNucleotideGeneral,
}

#[pymethods]
impl NucleotideGeneral {
    #[new]
    fn new_from_letter(val: u8) -> PyResult<Self> {
        Ok(Self {
            inner: map_io(RsNucleotideGeneral::from_u8_letter(val))?,
        })
    }

    #[classmethod]
    fn from_u8_letter(_cls: &Bound<PyType>, val: u8) -> PyResult<Self> {
        Ok(Self {
            inner: map_io(RsNucleotideGeneral::from_u8_letter(val))?,
        })
    }

    #[classmethod]
    fn from_str(_cls: &Bound<PyType>, s: &str) -> PyResult<Self> {
        Ok(Self {
            inner: map_io(RsNucleotideGeneral::from_str(s))?,
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
        self.inner.to_str_upper()
    }
    fn __repr__(&self) -> String {
        format!("NucleotideGeneral({})", self.inner.to_str_upper())
    }
}
