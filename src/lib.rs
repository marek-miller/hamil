use std::{
    collections::HashMap,
    hash::Hash,
    ops::AddAssign,
};

use num::{
    Float,
    Num,
};

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
pub enum ElectronIntegral {
    One((Orbital, Orbital)),
    Two((Orbital, Orbital, Orbital, Orbital)),
}

impl From<(Orbital, Orbital)> for ElectronIntegral {
    fn from(value: (Orbital, Orbital)) -> Self {
        Self::One(value)
    }
}

impl From<(Orbital, Orbital, Orbital, Orbital)> for ElectronIntegral {
    fn from(value: (Orbital, Orbital, Orbital, Orbital)) -> Self {
        Self::Two(value)
    }
}

impl ElectronIntegral {
    pub fn symmetries(self) -> impl Iterator<Item = Self> {
        Symmetries::new(self)
    }
}

#[derive(Debug)]
pub struct Symmetries {
    integral: ElectronIntegral,
    index:    u8,
}

impl Symmetries {
    pub fn new(integral: ElectronIntegral) -> Self {
        Self {
            integral,
            index: 0,
        }
    }
}

impl Iterator for Symmetries {
    type Item = ElectronIntegral;

    fn next(&mut self) -> Option<Self::Item> {
        if self.index == 0 {
            self.index += 1;
            return Some(self.integral);
        }

        let item = match self.integral {
            ElectronIntegral::One(h) => {
                if self.index == 1 {
                    Some(ElectronIntegral::One((h.1, h.0)))
                } else {
                    None
                }
            }
            ElectronIntegral::Two(h) => {
                let tuple = match self.index {
                    1 => Some((h.0, h.1, h.3, h.2)),
                    2 => Some((h.1, h.0, h.2, h.3)),
                    3 => Some((h.1, h.0, h.3, h.2)),
                    4 => Some((h.2, h.3, h.0, h.1)),
                    5 => Some((h.2, h.3, h.1, h.0)),
                    6 => Some((h.3, h.2, h.1, h.0)),
                    7 => Some((h.3, h.2, h.0, h.1)),
                    _ => None,
                };

                tuple.map(ElectronIntegral::Two)
            }
        };
        self.index += 1;
        item
    }
}

#[derive(Debug, Default, PartialEq, Clone, Copy, Hash)]
pub enum Pauli {
    #[default]
    I = 0,
    X = 1,
    Y = 2,
    Z = 3,
}

#[derive(Debug)]
pub struct Hamiltonian<K, T> {
    num_qubits: u32,
    coeffs:     HashMap<K, T>,
}

impl<K, T> Hamiltonian<K, T> {
    pub fn new(num_qubits: u32) -> Self {
        Self {
            num_qubits,
            coeffs: HashMap::new(),
        }
    }

    pub fn as_hashmap(&self) -> &HashMap<K, T> {
        &self.coeffs
    }

    pub fn as_hashmap_mut(&mut self) -> &mut HashMap<K, T> {
        &mut self.coeffs
    }
}

type FermionHamiltonian<T> = Hamiltonian<ElectronIntegral, T>;

impl<T> FermionHamiltonian<T>
where
    T: Copy + Hash + Eq,
{
    pub fn update(
        &mut self,
        integral: ElectronIntegral,
        value: T,
    ) {
        for x in integral.symmetries() {
            self.coeffs.insert(x, value);
        }
    }
}

impl<T> FermionHamiltonian<T>
where
    T: AddAssign + Copy + Hash + Eq,
{
    pub fn add(
        &mut self,
        integral: ElectronIntegral,
        value: T,
    ) {
        for x in integral.symmetries() {
            self.coeffs
                .entry(x)
                .and_modify(|coeff| *coeff += value)
                .or_insert(value);
        }
    }
}
