use std::str::FromStr;

use na_seq_rs;
use pyo3::{prelude::*, types::PyType};

use crate::map_io;

#[pyclass(module = "na_seq")]
#[derive(Clone, Copy)]
pub struct Element {
    pub inner: na_seq_rs::Element,
}

#[pymethods]
impl Element {
    #[classmethod]
    fn from_letter(_cls: &Bound<PyType>, s: &str) -> PyResult<Self> {
        Ok(Self {
            inner: map_io(na_seq_rs::Element::from_letter(s))?,
        })
    }

    fn to_letter(&self) -> String {
        self.inner.to_letter()
    }

    fn valence_typical(&self) -> usize {
        self.inner.valence_typical()
    }

    fn color(&self) -> (f32, f32, f32) {
        self.inner.color()
    }
    fn covalent_radius(&self) -> f64 {
        self.inner.covalent_radius()
    }
    fn vdw_radius(&self) -> f32 {
        self.inner.vdw_radius()
    }
    fn atomic_number(&self) -> u8 {
        self.inner.atomic_number()
    }
    fn atomic_weight(&self) -> f32 {
        self.inner.atomic_weight()
    }

    fn __str__(&self) -> String {
        self.inner.to_string()
    }
    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
}

#[pyclass(module = "na_seq")]
#[derive(Clone)]
pub struct AtomTypeInRes {
    pub inner: na_seq_rs::AtomTypeInRes,
}

#[pymethods]
impl AtomTypeInRes {
    #[classmethod]
    fn from_str(_cls: &Bound<PyType>, s: &str) -> PyResult<Self> {
        Ok(Self {
            inner: map_io(na_seq_rs::AtomTypeInRes::from_str(s))?,
        })
    }

    fn __str__(&self) -> String {
        self.inner.to_string()
    }
    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
}
