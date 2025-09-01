use std::{fmt, io, io::ErrorKind, str::FromStr};

use Element::*;

#[derive(Clone, Copy, PartialEq, Debug, Default, Hash, Eq)]
pub enum Element {
    Hydrogen,
    #[default]
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
}

impl Element {
    pub fn from_letter(letter: &str) -> io::Result<Self> {
        match letter.to_uppercase().as_ref() {
            "H" => Ok(Hydrogen),
            "C" => Ok(Carbon),
            "O" => Ok(Oxygen),
            "N" => Ok(Nitrogen),
            "F" => Ok(Fluorine),
            "S" => Ok(Sulfur),
            "P" => Ok(Phosphorus),
            "FE" => Ok(Iron),
            "CU" => Ok(Copper),
            "CA" => Ok(Calcium),
            "K" => Ok(Potassium),
            "AL" => Ok(Aluminum),
            "PB" => Ok(Lead),
            "AU" => Ok(Gold),
            "AG" => Ok(Silver),
            "HG" => Ok(Mercury),
            "SN" => Ok(Tin),
            "ZN" => Ok(Zinc),
            "MG" => Ok(Magnesium),
            "MN" => Ok(Manganese),
            "I" => Ok(Iodine),
            "CL" => Ok(Chlorine),
            "W" => Ok(Tungsten),
            "TE" => Ok(Tellurium),
            "SE" => Ok(Selenium),
            "BR" => Ok(Bromine),
            "RU" => Ok(Rubidium),
            // todo: Fill in if you need, or remove this fn.
            _ => Err(io::Error::new(
                ErrorKind::InvalidData,
                format!("Invalid atom letter: {letter}"),
            )),
        }
    }

    pub fn to_letter(self) -> String {
        match self {
            Hydrogen => "H".into(),
            Carbon => "C".into(),
            Oxygen => "O".into(),
            Nitrogen => "N".into(),
            Fluorine => "F".into(),
            Sulfur => "S".into(),
            Phosphorus => "P".into(),
            Iron => "Fe".into(),
            Copper => "Cu".into(),
            Calcium => "Ca".into(),
            Potassium => "K".into(),
            Aluminum => "Al".into(),
            Lead => "Pb".into(),
            Gold => "Au".into(),
            Silver => "Ag".into(),
            Mercury => "Hg".into(),
            Tin => "Sn".into(),
            Zinc => "Zz".into(),
            Magnesium => "Mg".into(),
            Manganese => "Mn".into(),
            Iodine => "I".into(),
            Chlorine => "Cl".into(),
            Tungsten => "W".into(),
            Tellurium => "Te".into(),
            Selenium => "Se".into(),
            Bromine => "Br".into(),
            Rubidium => "Ru".into(),
            Other => "X".into(),
        }
    }

    pub const fn valence_typical(&self) -> usize {
        match self {
            Hydrogen => 1,
            Carbon => 4,
            Oxygen => 2,
            Nitrogen => 3,
            Fluorine => 1,
            Sulfur => 2,     // can be 2, 4, or 6, but 2 is a common choice
            Phosphorus => 5, // can be 3 or 5, here we pick 5
            Iron => 2,       // Fe(II) is common (Fe(III) also common)
            Copper => 2,     // Cu(I) and Cu(II) both occur, pick 2 as a naive default
            Calcium => 2,
            Potassium => 1,
            Aluminum => 3,
            Lead => 2,    // Pb(II) or Pb(IV), but Pb(II) is more common/stable
            Gold => 3,    // Au(I) and Au(III) are common, pick 3
            Silver => 1,  // Ag(I) is most common
            Mercury => 2, // Hg(I) and Hg(II), pick 2
            Tin => 4,     // Sn(II) or Sn(IV), pick 4
            Zinc => 2,
            Magnesium => 2,
            Manganese => 7, // todo: Not sure
            Iodine => 1,    // can have higher, but 1 is typical in many simple compounds
            Chlorine => 1,  // can also be 3,5,7, but 1 is the simplest (e.g., HCl)
            Tungsten => 6,  // W can have multiple but 6 is a common oxidation state
            Tellurium => 2, // can also be 4 or 6, pick 2
            Selenium => 2,  // can also be 4 or 6, pick 2
            Bromine => 7,
            Rubidium => 1,
            Other => 0, // default to 0 for unknown or unhandled elements
        }
    }

    /// From [PyMol](https://pymolwiki.org/index.php/Color_Values)
    pub const fn color(&self) -> (f32, f32, f32) {
        match self {
            Hydrogen => (0.9, 0.9, 0.9),
            Carbon => (0.2, 1., 0.2),
            Oxygen => (1., 0.3, 0.3),
            Nitrogen => (0.2, 0.2, 1.0),
            Fluorine => (0.701, 1.0, 1.0),
            Sulfur => (0.9, 0.775, 0.25),
            Phosphorus => (1.0, 0.502, 0.),
            Iron => (0.878, 0.4, 0.2),
            Copper => (0.784, 0.502, 0.2),
            Calcium => (0.239, 1.0, 0.),
            Potassium => (0.561, 0.251, 0.831),
            Aluminum => (0.749, 0.651, 0.651),
            Lead => (0.341, 0.349, 0.380),
            Gold => (1., 0.820, 0.137),
            Silver => (0.753, 0.753, 0.753),
            Mercury => (0.722, 0.722, 0.816),
            Tin => (0.4, 0.502, 0.502),
            Zinc => (0.490, 0.502, 0.690),
            Magnesium => (0.541, 1., 0.),
            Manganese => (0.541, 1., 0.541),
            Iodine => (0.580, 0., 0.580),
            Chlorine => (0.121, 0.941, 0.121),
            Tungsten => (0.129, 0.580, 0.840),
            Tellurium => (0.831, 0.478, 0.),
            Selenium => (1.0, 0.631, 0.),
            Bromine => (1.0, 0.99, 0.),
            Rubidium => (0.439, 0.180, 0.690),
            Other => (5., 5., 5.),
        }
    }

    #[rustfmt::skip]
    /// Covalent radius, in angstrom.
    /// https://github.com/openbabel/openbabel/blob/master/src/elementtable.h
    /// https://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page)
    pub const fn covalent_radius(self) -> f64 {
        match self {
            Hydrogen   => 0.31,
            Carbon     => 0.76,
            Oxygen     => 0.66,
            Nitrogen   => 0.71,
            Fluorine   => 0.57,
            Sulfur     => 1.05,
            Phosphorus => 1.07,
            Iron       => 1.32,
            Copper     => 1.32,
            Calcium    => 1.76,
            Potassium  => 2.03,
            Aluminum   => 1.21,
            Lead       => 1.46,
            Gold       => 1.36,
            Silver     => 1.45,
            Mercury    => 1.32,
            Tin        => 1.39,
            Zinc       => 1.22,
            Magnesium  => 1.41, // 1.19?
            Manganese  => 1.39,
            Iodine     => 1.39,
            Chlorine   => 1.02,
            Tungsten   => 1.62,
            Tellurium  => 1.38,
            Selenium   => 1.20,
            Bromine  => 1.14, // 1.14 - 1.20
            Rubidium  => 2.20,
            Other      => 0.00,
        }
    }

    #[rustfmt::skip]
    /// Van-der-wals radius, in angstrom.
    /// https://github.com/openbabel/openbabel/blob/master/src/elementtable.h
    /// https://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page)
    pub const fn vdw_radius(&self) -> f32 {
        match self {
            Hydrogen   => 1.10, // or 120
            Carbon     => 1.70,
            Oxygen     => 1.52,
            Nitrogen   => 1.55,
            Fluorine   => 1.47,
            Sulfur     => 1.80,
            Phosphorus => 1.80,
            Iron       => 2.05,
            Copper     => 2.00,
            Calcium    => 2.31,
            Potassium  => 2.75,
            Aluminum   => 1.84,
            Lead       => 2.02,
            Gold       => 2.10,
            Silver     => 2.10,
            Mercury    => 2.05,
            Tin        => 1.93,
            Zinc       => 2.10,
            Magnesium  => 1.73,
            Manganese  => 0., // N/A?
            Iodine     => 1.98,
            Chlorine   => 1.75,
            Tungsten   => 2.10,
            Tellurium  => 2.06,
            Selenium   => 1.90,
            Bromine   => 1.85,
            Rubidium   => 3.21,
            Other      => 0.0,
        }
    }

    pub const fn atomic_number(&self) -> u8 {
        match self {
            Hydrogen => 1,
            Carbon => 6,
            Nitrogen => 7,
            Oxygen => 8,
            Fluorine => 9,
            Sulfur => 16,
            Phosphorus => 15,
            Iron => 26,
            Copper => 29,
            Calcium => 20,
            Potassium => 19,
            Aluminum => 13,
            Lead => 82,
            Gold => 79,
            Silver => 47,
            Mercury => 80,
            Tin => 50,
            Zinc => 30,
            Magnesium => 12,
            Manganese => 25,
            Iodine => 53,
            Chlorine => 17,
            Tungsten => 74,
            Tellurium => 52,
            Selenium => 34,
            Bromine => 35,
            Rubidium => 37,
            Other => 20, // fallback
        }
    }

    /// Standard atomic weight (in atomic mass units) for each element.
    pub fn atomic_weight(&self) -> f32 {
        match self {
            Hydrogen => 1.008,
            Carbon => 12.011,
            Oxygen => 15.999,
            Nitrogen => 14.007,
            Fluorine => 18.998,
            Sulfur => 32.06,
            Phosphorus => 30.974,
            Iron => 55.845,
            Copper => 63.546,
            Calcium => 40.078,
            Potassium => 39.098,
            Aluminum => 26.982,
            Lead => 207.2,
            Gold => 196.967,
            Silver => 107.8682,
            Mercury => 200.592,
            Tin => 118.71,
            Zinc => 65.38,
            Magnesium => 24.305,
            Manganese => 54.938,
            Iodine => 126.90,
            Chlorine => 35.45,
            Tungsten => 183.84,
            Tellurium => 127.60,
            Selenium => 78.971,
            Bromine => 79.904,
            Rubidium => 85.468,
            Other => 0.0, // fallback for unknowns
        }
    }
}

impl fmt::Display for Element {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let v = match self {
            Hydrogen => "Hydrogen",
            Carbon => "Carbon",
            Oxygen => "Oxygen",
            Nitrogen => "Nitrogen",
            Fluorine => "Fluorine",
            Sulfur => "Sulfur",
            Phosphorus => "Phosphorus",
            Iron => "Iron",
            Copper => "Copper",
            Calcium => "Calcium",
            Potassium => "Potassium",
            Aluminum => "Aluminum",
            Lead => "Lead",
            Gold => "Gold",
            Silver => "Silver",
            Mercury => "Mercury",
            Tin => "Tin",
            Zinc => "Zinc",
            Magnesium => "Magnesium",
            Manganese => "Manganese",
            Iodine => "Iodine",
            Chlorine => "Chlorine",
            Tungsten => "Tungsten",
            Tellurium => "Tellurium",
            Selenium => "Selenium",
            Bromine => "Bromine",
            Rubidium => "Rubidium",
            Other => "Other",
        };

        write!(f, "{v}")
    }
}

/// Identifies the atom "type" or "name", as used in a residue. This information
/// is provided, for example, in atom-coordinate mmCIF files.
#[derive(Clone, PartialEq, Debug)]
pub enum AtomTypeInRes {
    C,
    CA,
    CB,
    CD,
    CD1,
    CD2,
    CE,
    CE1,
    CE2,
    CE3,
    CG,
    CG1,
    CG2,
    CH2,
    CH3,
    CZ,
    CZ1,
    CZ2,
    CZ3,
    O,
    OD1,
    OD2,
    OE1,
    OE2,
    OG,
    OH,
    OXT,
    N,
    ND1,
    ND2,
    NE,
    NZ,
    NH1,
    NH2,
    NE1,
    NE2,
    OG1,
    OG2,
    SD,
    SE,
    SG,
    H(String),
    /// E.g. ligands and water molecules.
    Hetero(String),
}

impl FromStr for AtomTypeInRes {
    type Err = io::Error;

    /// Accepts the exact (case-sensitive) atom label.
    /// We parse hetero atoms into this manually, to prevent accidental coercion to that.
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        // We keep this caps-sensitive, as residues and
        // proteins use different capitalization conventions.

        if s.starts_with("H") {
            return Ok(Self::H(s.to_string()));
        }

        match s {
            "C" => Ok(Self::C),
            "CA" => Ok(Self::CA),
            "CB" => Ok(Self::CB),
            "CD" => Ok(Self::CD),
            "CD1" => Ok(Self::CD1),
            "CD2" => Ok(Self::CD2),
            "CE" => Ok(Self::CE),
            "CE1" => Ok(Self::CE1),
            "CE2" => Ok(Self::CE2),
            "CE3" => Ok(Self::CE3),
            "CG" => Ok(Self::CG),
            "CG1" => Ok(Self::CG1),
            "CG2" => Ok(Self::CG2),
            "CH2" => Ok(Self::CH2),
            "CH3" => Ok(Self::CH3),
            "CZ" => Ok(Self::CZ),
            "CZ1" => Ok(Self::CZ1),
            "CZ2" => Ok(Self::CZ2),
            "CZ3" => Ok(Self::CZ3),
            "O" => Ok(Self::O),
            "OD1" => Ok(Self::OD1),
            "OD2" => Ok(Self::OD2),
            "OE1" => Ok(Self::OE1),
            "OE2" => Ok(Self::OE2),
            "OG" => Ok(Self::OG),
            "OXT" => Ok(Self::OXT),
            "OH" => Ok(Self::OH),
            "N" => Ok(Self::N),
            "ND1" => Ok(Self::ND1),
            "ND2" => Ok(Self::ND2),
            "NE" => Ok(Self::NE),
            "NZ" => Ok(Self::NZ),
            "NH1" => Ok(Self::NH1),
            "NH2" => Ok(Self::NH2),
            "NE1" => Ok(Self::NE1),
            "NE2" => Ok(Self::NE2),
            "OG1" => Ok(Self::OG1),
            "OG2" => Ok(Self::OG2),
            "SD" => Ok(Self::SD),
            "SE" => Ok(Self::SE),
            "SG" => Ok(Self::SG),
            _ => Err(io::Error::new(
                ErrorKind::InvalidData,
                "Unknown atom type when parsing type in res",
            )),
        }
    }
}

impl fmt::Display for AtomTypeInRes {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let label = match self {
            Self::C => "C",
            Self::CA => "CA",
            Self::CB => "CB",
            Self::CD => "CD",
            Self::CD1 => "CD1",
            Self::CD2 => "CD2",
            Self::CE => "CE",
            Self::CE1 => "CE1",
            Self::CE2 => "CE2",
            Self::CE3 => "CE3",
            Self::CG => "CG",
            Self::CG1 => "CG1",
            Self::CG2 => "CG2",
            Self::CH2 => "CH2",
            Self::CH3 => "CH3",
            Self::CZ => "CZ",
            Self::CZ1 => "CZ1",
            Self::CZ2 => "CZ2",
            Self::CZ3 => "CZ3",
            Self::O => "O",
            Self::OD1 => "OD1",
            Self::OD2 => "OD2",
            Self::OE1 => "OE1",
            Self::OE2 => "OE2",
            Self::OG => "OG",
            Self::OH => "OH",
            Self::OXT => "OXT",
            Self::N => "N",
            Self::ND1 => "ND1",
            Self::ND2 => "ND2",
            Self::NE => "NE",
            Self::NZ => "NZ",
            Self::NH1 => "NH1",
            Self::NH2 => "NH2",
            Self::NE1 => "NE1",
            Self::NE2 => "NE2",
            Self::OG1 => "OG1",
            Self::OG2 => "OG2",
            Self::SD => "SD",
            Self::SE => "SE",
            Self::SG => "SG",
            Self::H(label) => label,
            Self::Hetero(label) => label,
        };
        f.write_str(label)
    }
}
