use std::str::FromStr;

use na_seq_rs;
use pyo3::{prelude::*, types::PyType};

use crate::make_enum;

make_enum!(
    Element,
    na_seq_rs::Element,
    Hydrogen,
    Carbon,
    Oxygen,
    Nitrogen,
    Fluorine,
    Sulfur,
    Phosphorus,
    Iron,
    Copper,
    Calcium,
    Potassium,
    Aluminum,
    Lead,
    Gold,
    Silver,
    Mercury,
    Tin,
    Zinc,
    Magnesium,
    Manganese,
    Iodine,
    Chlorine,
    Tungsten,
    Tellurium,
    Selenium,
    Bromine,
    Rubidium,
    Other,
);

#[pymethods]
impl Element {
    #[classmethod]
    fn from_letter(_cls: &Bound<PyType>, s: &str) -> PyResult<Self> {
        Ok(na_seq_rs::Element::from_letter(s)?.into())
    }

    fn to_letter(&self) -> String {
        self.to_native().to_letter()
    }

    fn valence_typical(&self) -> usize {
        self.to_native().valence_typical()
    }

    fn color(&self) -> (f32, f32, f32) {
        self.to_native().color()
    }
    fn covalent_radius(&self) -> f64 {
        self.to_native().covalent_radius()
    }
    fn vdw_radius(&self) -> f32 {
        self.to_native().vdw_radius()
    }
    fn atomic_number(&self) -> u8 {
        self.to_native().atomic_number()
    }
    fn atomic_weight(&self) -> f32 {
        self.to_native().atomic_weight()
    }

    fn __str__(&self) -> String {
        self.to_native().to_string()
    }
    fn __repr__(&self) -> String {
        format!("{:?}", self.to_native())
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
            inner: na_seq_rs::AtomTypeInRes::from_str(s)?,
        })
    }

    fn __str__(&self) -> String {
        self.inner.to_string()
    }
    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
}
