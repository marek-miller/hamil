use std::{
    collections::{
        BTreeMap,
        HashMap,
    },
    hash::Hash,
    ops::AddAssign,
};

use num::Float;

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

/// FermiCode in canonical ordering
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
    pub fn one(
        cr: Orbital,
        an: Orbital,
    ) -> Option<Self> {
        (cr.h() <= an.h()).then(|| FermiCode::One {
            cr,
            an,
        })
    }

    pub fn two(
        cr: (Orbital, Orbital),
        an: (Orbital, Orbital),
    ) -> Option<Self> {
        (cr.0.h() < cr.1.h() && an.0.h() > an.1.h() && cr.0.h() <= an.1.h())
            .then(|| FermiCode::Two {
                cr,
                an,
            })
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

    pub fn with_pauli(
        idx: u32,
        pauli: Pauli,
    ) -> Self {
        let mut code = PauliCode::new();
        code.0.insert(idx, pauli);
        code
    }

    pub fn update(
        &mut self,
        idx: u32,
        pauli: Pauli,
    ) -> &mut Self {
        self.0.insert(idx, pauli);
        self
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

#[derive(Debug)]
pub struct Hamiltonian<K, T> {
    terms: HashMap<K, T>,
}

impl<K, T> Hamiltonian<K, T>
where
    K: Hash + Eq,
{
    pub fn new() -> Self {
        Self {
            terms: HashMap::new(),
        }
    }

    pub fn update(
        &mut self,
        code: K,
        value: T,
    ) -> Option<T> {
        self.terms.insert(code, value)
    }
}

impl<K, T> Hamiltonian<K, T>
where
    T: AddAssign + Copy,
    K: Hash + Eq,
{
    pub fn add(
        &mut self,
        code: K,
        value: T,
    ) {
        self.terms
            .entry(code)
            .and_modify(|coeff| *coeff += value)
            .or_insert(value);
    }
}

pub type FermiHamil<T> = Hamiltonian<FermiCode, T>;
pub type PauliHamil<T> = Hamiltonian<PauliCode, T>;

impl<T> From<FermiHamil<T>> for PauliHamil<T>
where
    T: Float + AddAssign,
{
    fn from(hamil: FermiHamil<T>) -> Self {
        let mut pauli_hamil = PauliHamil::new();
        for (code, value) in hamil.terms {
            match code {
                FermiCode::One {
                    cr,
                    an,
                } => {
                    let (p, q) = (cr.h(), an.h());
                    let half = T::from(0.5).unwrap();
                    if p == q {
                        pauli_hamil.add(
                            PauliCode::with_pauli(p, Pauli::I),
                            value * half,
                        );
                        pauli_hamil.add(
                            PauliCode::with_pauli(q, Pauli::Z),
                            value * -half,
                        );
                    } else {
                        // canonical ordering means p<=q
                        let mut pauli_code = PauliCode::new();
                        pauli_code.update(p, Pauli::X).update(q, Pauli::X);
                        for i in p + 1..q - 1 {
                            pauli_code.update(i, Pauli::Z);
                        }
                        pauli_hamil.add(pauli_code, value * half);

                        let mut pauli_code = PauliCode::new();
                        pauli_code.update(p, Pauli::Y).update(q, Pauli::Y);
                        for i in p + 1..q - 1 {
                            pauli_code.update(i, Pauli::Z);
                        }
                        pauli_hamil.add(pauli_code, value * half)
                    }
                }
                FermiCode::Two {
                    cr,
                    an,
                } => {
                    let (p, q, r, s) = (cr.0.h(), cr.1.h(), an.0.h(), an.1.h());
                    let quarter = T::from(0.25).unwrap();

                    if p == s && q == r {
                        // canonical ordering means q>p
                        pauli_hamil.add(
                            PauliCode::with_pauli(p, Pauli::I),
                            value * quarter,
                        );
                        pauli_hamil.add(
                            PauliCode::with_pauli(p, Pauli::Z),
                            value * -quarter,
                        );
                        pauli_hamil.add(
                            PauliCode::with_pauli(q, Pauli::Z),
                            value * -quarter,
                        );
                        let mut pauli_code = PauliCode::new();
                        pauli_code.update(p, Pauli::Z).update(q, Pauli::Z);
                        pauli_hamil.add(pauli_code, value * quarter);
                    } else if q == r {
                        let mut pauli_code = PauliCode::new();
                        for i in p + 1..s - 1 {
                            pauli_code.update(i, Pauli::Z);
                        }
                        pauli_code.update(p, Pauli::X).update(s, Pauli::X);
                        pauli_hamil.add(pauli_code, value * quarter);

                        let mut pauli_code = PauliCode::new();
                        for i in p + 1..s - 1 {
                            pauli_code.update(i, Pauli::Z);
                        }
                        pauli_code.update(p, Pauli::X).update(s, Pauli::X);
                        pauli_code.update(q, Pauli::Z);
                        pauli_hamil.add(pauli_code, value * -quarter);

                        let mut pauli_code = PauliCode::new();
                        for i in p + 1..s - 1 {
                            pauli_code.update(i, Pauli::Z);
                        }
                        pauli_code.update(p, Pauli::Y).update(s, Pauli::Y);
                        pauli_hamil.add(pauli_code, value * quarter);

                        let mut pauli_code = PauliCode::new();
                        for i in p + 1..s - 1 {
                            pauli_code.update(i, Pauli::Z);
                        }
                        pauli_code.update(p, Pauli::Y).update(s, Pauli::Y);
                        pauli_code.update(q, Pauli::Z);
                        pauli_hamil.add(pauli_code, value * -quarter);
                    } else {
                        let eighth = quarter / (T::one() + T::one());
                        let pc = {
                            let mut pc = PauliCode::new();
                            for i in p + 1..q - 1 {
                                pc.update(i, Pauli::Z);
                            }
                            for i in s + 1..r - 1 {
                                pc.update(i, Pauli::Z);
                            }
                            pc
                        };

                        let mut pauli_code = pc.clone();
                        pauli_code
                            .update(p, Pauli::X)
                            .update(q, Pauli::X)
                            .update(r, Pauli::X)
                            .update(s, Pauli::X);
                        pauli_hamil.add(pauli_code, eighth);

                        let mut pauli_code = pc.clone();
                        pauli_code
                            .update(p, Pauli::X)
                            .update(q, Pauli::X)
                            .update(r, Pauli::Y)
                            .update(s, Pauli::Y);
                        pauli_hamil.add(pauli_code, -eighth);

                        let mut pauli_code = pc.clone();
                        pauli_code
                            .update(p, Pauli::X)
                            .update(q, Pauli::Y)
                            .update(r, Pauli::X)
                            .update(s, Pauli::Y);
                        pauli_hamil.add(pauli_code, eighth);

                        let mut pauli_code = pc.clone();
                        pauli_code
                            .update(p, Pauli::Y)
                            .update(q, Pauli::X)
                            .update(r, Pauli::X)
                            .update(s, Pauli::Y);
                        pauli_hamil.add(pauli_code, eighth);

                        let mut pauli_code = pc.clone();
                        pauli_code
                            .update(p, Pauli::Y)
                            .update(q, Pauli::X)
                            .update(r, Pauli::Y)
                            .update(s, Pauli::X);
                        pauli_hamil.add(pauli_code, eighth);

                        let mut pauli_code = pc.clone();
                        pauli_code
                            .update(p, Pauli::Y)
                            .update(q, Pauli::Y)
                            .update(r, Pauli::X)
                            .update(s, Pauli::X);
                        pauli_hamil.add(pauli_code, -eighth);

                        let mut pauli_code = pc.clone();
                        pauli_code
                            .update(p, Pauli::X)
                            .update(q, Pauli::Y)
                            .update(r, Pauli::Y)
                            .update(s, Pauli::X);
                        pauli_hamil.add(pauli_code, eighth);

                        let mut pauli_code = pc;
                        pauli_code
                            .update(p, Pauli::Y)
                            .update(q, Pauli::Y)
                            .update(r, Pauli::Y)
                            .update(s, Pauli::Y);
                        pauli_hamil.add(pauli_code, eighth);
                    }
                }
            }
        }
        pauli_hamil
    }
}
