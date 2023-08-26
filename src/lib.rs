use std::{
    collections::{
        BTreeMap,
        HashMap,
    },
    hash::Hash,
    ops::{
        AddAssign,
        Index,
        IndexMut,
    },
};

use num::{
    Complex,
    Float,
    Num,
};

#[cfg(test)]
mod tests;

#[derive(Debug, Default, Clone, Copy, PartialEq, Hash, Eq)]
pub enum Spin {
    #[default]
    Down,
    Up,
}

impl Spin {
    pub fn is_up(&self) -> bool {
        *self == Self::Up
    }

    pub fn flip(self) -> Self {
        use Spin::*;
        match self {
            Down => Up,
            Up => Down,
        }
    }
}

impl From<Spin> for u32 {
    fn from(value: Spin) -> Self {
        match value {
            Spin::Down => 0,
            Spin::Up => 1,
        }
    }
}

impl From<u32> for Spin {
    fn from(value: u32) -> Self {
        if value == 0 {
            Spin::Down
        } else {
            Spin::Up
        }
    }
}

#[derive(Debug, Default, Clone, Copy, PartialEq, Hash, Eq)]
pub struct Orbital {
    pub j: u32,
    pub s: Spin,
}

impl Orbital {
    pub fn new(
        index: u32,
        spin: Spin,
    ) -> Self {
        Self {
            j: index, s: spin
        }
    }

    pub fn h(&self) -> u32 {
        self.j * 2 + u32::from(self.s)
    }

    pub fn g(
        &self,
        num_orbitals: u32,
    ) -> u32 {
        self.j + u32::from(self.s) * num_orbitals
    }
}

impl From<(u32, Spin)> for Orbital {
    fn from(value: (u32, Spin)) -> Self {
        Self::new(value.0, value.1)
    }
}

impl From<Orbital> for u32 {
    fn from(value: Orbital) -> Self {
        value.h()
    }
}

impl From<u32> for Orbital {
    fn from(value: u32) -> Self {
        let orbital = value / 2;
        let spin = Spin::from(value % 2);
        Self::new(orbital, spin)
    }
}

#[derive(Debug, PartialEq, Clone, Copy, Hash, Eq)]
pub enum FermiCode {
    One {
        cr: Orbital,
        an: Orbital,
    },
    Two {
        cr: (Orbital, Orbital),
        an: (Orbital, Orbital),
    },
}

impl FermiCode {
    pub fn canonical(self) -> Self {
        match self {
            FermiCode::One { cr, an } => { 
                if cr.h() <= an.h() {
                    self
                } else {
                    FermiCode::One { an, cr}
                }
            }
            FermiCode::Two { cr, an } => {
                let 
            }
        }
    }
}

impl From<(Orbital, Orbital)> for FermiCode {
    fn from(value: (Orbital, Orbital)) -> Self {
        Self::One {
            cr: value.0,
            an: value.1,
        }
    }
}

impl From<(Orbital, Orbital, Orbital, Orbital)> for FermiCode {
    fn from(value: (Orbital, Orbital, Orbital, Orbital)) -> Self {
        Self::Two {
            cr: (value.0, value.1),
            an: (value.2, value.3),
        }
    }
}

impl FermiCode {
    pub fn symmetries(self) -> impl Iterator<Item = Self> {
        ElectronSyms::new(self)
    }
}

impl IntoIterator for FermiCode {
    type IntoIter = ElectronOrbitals;
    type Item = Orbital;

    fn into_iter(self) -> Self::IntoIter {
        ElectronOrbitals::new(self)
    }
}

#[derive(Debug)]
pub struct ElectronOrbitals {
    idx:   FermiCode,
    count: u8,
}

impl ElectronOrbitals {
    pub fn new(idx: FermiCode) -> Self {
        Self {
            idx,
            count: 0,
        }
    }
}

impl Iterator for ElectronOrbitals {
    type Item = Orbital;

    fn next(&mut self) -> Option<Self::Item> {
        let item = match self.idx {
            FermiCode::One {
                cr,
                an,
            } => match self.count {
                0 => Some(cr),
                1 => Some(an),
                _ => None,
            },
            FermiCode::Two {
                cr,
                an,
            } => match self.count {
                0 => Some(cr.0),
                1 => Some(cr.1),
                2 => Some(an.0),
                3 => Some(an.1),
                _ => None,
            },
        };
        self.count += 1;
        if self.count >= 4 {
            self.count = 4;
        }
        item
    }
}

#[derive(Debug)]
pub struct ElectronSyms {
    idx:   FermiCode,
    count: u8,
}

impl ElectronSyms {
    pub fn new(integral: FermiCode) -> Self {
        Self {
            idx:   integral,
            count: 0,
        }
    }
}

impl Iterator for ElectronSyms {
    type Item = FermiCode;

    fn next(&mut self) -> Option<Self::Item> {
        let item = match self.idx {
            FermiCode::One {
                cr,
                an,
            } => match self.count {
                0 => Some((cr, an).into()),
                1 => Some((an, cr).into()),
                _ => None,
            },

            FermiCode::Two {
                cr,
                an,
            } => match self.count {
                0 => Some((cr.0, cr.1, an.0, an.1).into()),
                1 => Some((cr.0, cr.1, an.1, an.0).into()),
                2 => Some((cr.1, cr.0, an.0, an.1).into()),
                3 => Some((cr.1, cr.0, an.1, an.0).into()),
                4 => Some((an.0, an.1, cr.0, cr.1).into()),
                5 => Some((an.0, an.1, cr.1, cr.0).into()),
                6 => Some((an.1, an.0, cr.1, cr.0).into()),
                7 => Some((an.1, an.0, cr.0, cr.1).into()),
                _ => None,
            },
        };
        self.count += 1;
        if self.count >= 8 {
            self.count = 8
        }
        item
    }
}

#[derive(Debug, Default, PartialEq, Clone, Copy, Hash, Eq)]
pub enum Pauli {
    #[default]
    I,
    X,
    Y,
    Z,
}

#[derive(Debug, PartialEq, Clone, Hash, Eq)]
pub struct PauliCode(BTreeMap<u32, Pauli>);

impl Default for PauliCode {
    fn default() -> Self {
        Self::new()
    }
}

impl PauliCode {
    pub fn new() -> Self {
        Self(BTreeMap::new())
    }
}

impl From<&[Pauli]> for PauliCode {
    fn from(value: &[Pauli]) -> Self {
        let mut code = PauliCode::new();
        for (i, pauli) in value.iter().enumerate() {
            code.0.insert(i as u32, *pauli);
        }
        code
    }
}

pub trait Conj {
    fn conj(self) -> Self;
}

impl<T> Conj for Complex<T>
where
    T: Float,
{
    fn conj(self) -> Self {
        Complex::conj(&self)
    }
}

impl Conj for PauliCode {
    fn conj(self) -> Self {
        self
    }
}

impl Conj for FermiCode {
    fn conj(self) -> Self {
        match self {
            FermiCode::One {
                cr,
                an,
            } => (an, cr).into(),
            FermiCode::Two {
                cr,
                an,
            } => (an.1, an.0, cr.1, cr.0).into(),
        }
    }
}

#[derive(Debug)]
pub struct Hamiltonian<K, T> {
    terms: HashMap<K, T>,
}

impl<K, T> Hamiltonian<K, T> {
    pub fn new() -> Self {
        Self {
            terms: HashMap::new(),
        }
    }
}

type FermiHamil<T> = Hamiltonian<FermiCode, T>;

impl<T> FermiHamil<T>
where
    T: Copy + Hash + Eq,
{
    /// # Panics
    ///
    /// Panics if `idx` is out of bound, i.e. `orbital.h() >= N`
    /// for any `orbital` in `idx`
    pub fn update(
        &mut self,
        code: FermiCode,
        value: T,
    ) {
        self.terms.insert(code.canonical(), value);
    }
}

impl<T> FermiHamil<T>
where
    T: AddAssign + Copy + Hash + Eq,
{
    pub fn add(
        &mut self,
        code: FermiCode,
        value: T,
    ) {
        self.terms
            .entry(code.canonical())
            .and_modify(|coeff| *coeff += value)
            .or_insert(value);
    }
}

type PauliHamil<T> = Hamiltonian<PauliCode, T>;

impl<T> PauliHamil<T>
where
    T: Hash + Eq,
{
    pub fn update(
        &mut self,
        code: PauliCode,
        value: T,
    ) {
        self.terms.insert(code, value);
    }
}

impl<T> PauliHamil<T>
where
    T: AddAssign + Copy + Hash + Eq,
{
    pub fn add(
        &mut self,
        code: PauliCode,
        value: T,
    ) {
        self.terms
            .entry(code)
            .and_modify(|coeff| *coeff += value)
            .or_insert(value);
    }
}

impl<T> From<FermiHamil<T>> for PauliHamil<T>
where
    T: Float,
{
    fn from(hamil: FermiHamil<T>) -> Self {
        todo!()
    }
}
