use std::str::FromStr;

use na_seq_rs::{AtomTypeInRes as RsAtomTypeInRes, Element as RsElement};
use pyo3::{prelude::*, types::PyType};

use crate::map_io;

#[pyclass(module = "na_seq")]
#[derive(Clone, Copy)]
pub struct Element {
    pub inner: RsElement,
}

#[pymethods]
impl Element {
    #[new]
    fn new_from_letter(s: &str) -> PyResult<Self> {
        Ok(Self {
            inner: map_io(RsElement::from_letter(s))?,
        })
    }

    #[classmethod]
    fn from_letter(_cls: &Bound<PyType>, s: &str) -> PyResult<Self> {
        Ok(Self {
            inner: map_io(RsElement::from_letter(s))?,
        })
    }

    fn to_letter(&self) -> String {
        self.inner.to_letter()
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

    fn name(&self) -> String {
        self.inner.to_string()
    }

    #[getter]
    fn letter(&self) -> String {
        self.inner.to_letter()
    }

    fn __str__(&self) -> String {
        self.inner.to_letter()
    }
    fn __repr__(&self) -> String {
        format!("Element({})", self.inner.to_letter())
    }
}

#[pyclass(module = "na_seq")]
#[derive(Clone)]
pub struct AtomTypeInRes {
    pub inner: RsAtomTypeInRes,
}

#[pymethods]
impl AtomTypeInRes {
    #[new]
    fn new_from_label(s: &str) -> PyResult<Self> {
        Ok(Self {
            inner: map_io(RsAtomTypeInRes::from_str(s))?,
        })
    }

    #[classmethod]
    fn from_str(_cls: &Bound<PyType>, s: &str) -> PyResult<Self> {
        Ok(Self {
            inner: map_io(RsAtomTypeInRes::from_str(s))?,
        })
    }

    fn label(&self) -> String {
        self.inner.to_string()
    }

    fn __str__(&self) -> String {
        self.inner.to_string()
    }
    fn __repr__(&self) -> String {
        format!("AtomTypeInRes({})", self.inner.to_string())
    }
}
