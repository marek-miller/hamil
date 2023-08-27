#![feature(step_trait)]

use std::{
    collections::{
        BTreeMap,
        HashMap,
    },
    hash::Hash,
    iter::Step,
    ops::AddAssign,
};

use num::{
    Float,
    Num,
};

#[cfg(test)]
mod tests;

pub trait Enumerate<T> {
    fn enumerate(&self) -> T;
}

#[derive(Debug, Default, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Spin {
    #[default]
    Down,
    Up,
}

impl Spin {
    #[must_use]
    pub fn is_up(&self) -> bool {
        *self == Self::Up
    }

    #[must_use]
    pub fn flip(self) -> Self {
        use Spin::{
            Down,
            Up,
        };
        match self {
            Down => Up,
            Up => Down,
        }
    }
}

impl<T> Enumerate<T> for Spin
where
    T: Num,
{
    fn enumerate(&self) -> T {
        match self {
            Self::Down => T::zero(),
            Self::Up => T::one(),
        }
    }
}

#[derive(Debug, Default, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Orbital<K> {
    pub j: K,
    pub s: Spin,
}

impl<K> Orbital<K> {
    pub fn new(
        j: K,
        s: Spin,
    ) -> Self {
        Self {
            j,
            s,
        }
    }
}
impl<K> Orbital<K>
where
    K: Num + Copy,
{
    pub fn h(&self) -> K {
        self.j * (K::one() + K::one()) + self.s.enumerate()
    }

    pub fn g(
        &self,
        num_orbitals: K,
    ) -> K {
        self.j + num_orbitals * self.s.enumerate()
    }
}

impl<K> From<(K, Spin)> for Orbital<K> {
    fn from(value: (K, Spin)) -> Self {
        Self::new(value.0, value.1)
    }
}

impl<K> Enumerate<K> for Orbital<K>
where
    K: Num + Copy,
{
    fn enumerate(&self) -> K {
        self.h()
    }
}

/// `FermiCode` in canonical ordering
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum FermiCode<K> {
    One {
        cr: Orbital<K>,
        an: Orbital<K>,
    },
    Two {
        cr: (Orbital<K>, Orbital<K>),
        an: (Orbital<K>, Orbital<K>),
    },
}

impl<K> FermiCode<K>
where
    K: Num + PartialOrd + Copy,
{
    pub fn one(
        cr: Orbital<K>,
        an: Orbital<K>,
    ) -> Option<Self> {
        (cr.enumerate() <= an.enumerate()).then_some(FermiCode::One {
            cr,
            an,
        })
    }

    pub fn two(
        cr: (Orbital<K>, Orbital<K>),
        an: (Orbital<K>, Orbital<K>),
    ) -> Option<Self> {
        (cr.0.enumerate() < cr.1.enumerate()
            && an.0.enumerate() > an.1.enumerate()
            && cr.0.enumerate() <= an.1.enumerate())
        .then_some(FermiCode::Two {
            cr,
            an,
        })
    }
}

#[derive(Debug, Default, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Pauli {
    #[default]
    X,
    Y,
    Z,
}

macro_rules! impl_enumerate_pauli {
    ($($Typ:ty)*) => {
        $(
            impl Enumerate<$Typ> for Pauli {
                fn enumerate(&self) -> $Typ {
                    use Pauli::*;
                    match self {
                        X => 0,
                        Y => 1,
                        Z => 2,
                    }
                }
            }
        )*
    };
}

impl_enumerate_pauli!(u8 u16 u32 u64 usize);

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum PauliCode<K> {
    Identity,
    Product(BTreeMap<K, Pauli>),
}

impl<K> Default for PauliCode<K> {
    fn default() -> Self {
        Self::Identity
    }
}

impl<K> PauliCode<K> {
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }
}

impl<K> PauliCode<K>
where
    K: Ord,
{
    pub fn with_pauli(
        idx: K,
        pauli: Pauli,
    ) -> Self {
        let mut tree = BTreeMap::new();
        tree.insert(idx, pauli);
        Self::Product(tree)
    }

    pub fn update(
        &mut self,
        idx: K,
        pauli: Pauli,
    ) {
        match self {
            PauliCode::Identity => {
                *self = PauliCode::with_pauli(idx, pauli);
            }
            PauliCode::Product(tree) => {
                tree.insert(idx, pauli);
            }
        }
    }
}

impl<I, K> From<I> for PauliCode<K>
where
    I: IntoIterator<Item = (K, Pauli)>,
    K: Ord,
{
    fn from(value: I) -> Self {
        let mut code = PauliCode::new();
        for (idx, pauli) in value {
            code.update(idx, pauli);
        }
        code
    }
}

macro_rules! impl_enumerate_paulicode {
    ($($Typ:ty)*) => {
        $(
            impl Enumerate<$Typ> for PauliCode<$Typ> {
                fn enumerate(&self) -> $Typ {
                    let mut count = 0;
                    if let PauliCode::Product(tree) = self {
                        for (key, value) in tree {
                            count += key * 3 + <Pauli as Enumerate<$Typ>>::enumerate(value)
                        }
                    }
                    count
                }
            }
        )*
    };
}

impl_enumerate_paulicode!(u8 u16 u32 u64 usize);

#[derive(Debug)]
pub struct Hamil<T, K> {
    terms: HashMap<K, T>,
}

impl<T, K> Default for Hamil<T, K>
where
    K: Hash + Eq,
{
    fn default() -> Self {
        Self::new()
    }
}

impl<T, K> Hamil<T, K>
where
    K: Hash + Eq,
{
    #[must_use]
    pub fn new() -> Self {
        Self {
            terms: HashMap::new(),
        }
    }
}

impl<T, K> Hamil<T, K>
where
    K: Clone + Eq + Hash,
{
    pub fn replace_or_insert(
        &mut self,
        code: &K,
        value: T,
    ) {
        if let Some(x) = self.terms.get_mut(code) {
            *x = value;
        } else {
            self.terms.insert(code.clone(), value);
        }
    }
}

impl<T, K> Hamil<T, K>
where
    T: AddAssign + Copy,
    K: Clone + Eq + Hash,
{
    pub fn add_or_insert(
        &mut self,
        code: &K,
        value: T,
    ) {
        if let Some(x) = self.terms.get_mut(code) {
            *x += value;
        } else {
            self.terms.insert(code.clone(), value);
        }
    }
}

pub type FermiHamil<T, K> = Hamil<T, FermiCode<K>>;
pub type PauliHamil<T, K> = Hamil<T, PauliCode<K>>;

impl<T, K> From<Hamil<T, FermiCode<K>>> for Hamil<T, PauliCode<K>>
where
    T: Float + AddAssign,
    K: Copy + Hash + Num + Ord + Step,
{
    fn from(fermi_hamil: FermiHamil<T, K>) -> Self {
        let mut hamil = PauliHamil::new();
        for (fermi_code, value) in fermi_hamil.terms {
            match fermi_code {
                FermiCode::One {
                    cr,
                    an,
                } => {
                    let (p, q) = (cr.enumerate(), an.enumerate());
                    update_hamil_one_pq(&mut hamil, value, p, q);
                }
                FermiCode::Two {
                    cr,
                    an,
                } => {
                    let (p, q, r, s) = (
                        cr.0.enumerate(),
                        cr.1.enumerate(),
                        an.0.enumerate(),
                        an.1.enumerate(),
                    );

                    if p == s && q == r {
                        update_hamil_two_pq(&mut hamil, value, p, q);
                    } else if q == r {
                        update_hamil_two_pqs(&mut hamil, value, p, q, s);
                    } else {
                        update_hamil_two_pqrs(&mut hamil, value, p, q, r, s);
                    }
                }
            }
        }
        hamil
    }
}

fn update_hamil_one_pq<T, K>(
    hamil: &mut Hamil<T, PauliCode<K>>,
    value: T,
    p: K,
    q: K,
) where
    T: Float + AddAssign,
    K: Copy + Hash + Num + Ord + Step,
{
    let half = T::from(0.5).unwrap();

    if p == q {
        hamil.add_or_insert(&PauliCode::Identity, value * half);
        hamil.add_or_insert(&PauliCode::with_pauli(q, Pauli::Z), value * -half);
    } else {
        // canonical ordering means p<=q
        let mut code = PauliCode::from(
            (p + K::one()..q - K::one()).map(|i| (i, Pauli::Z)),
        );
        code.update(p, Pauli::X);
        code.update(q, Pauli::X);
        hamil.add_or_insert(&code, value * half);
        code.update(p, Pauli::Y);
        code.update(q, Pauli::Y);
        hamil.add_or_insert(&code, value * half);
    }
}

fn update_hamil_two_pq<T, K>(
    hamil: &mut Hamil<T, PauliCode<K>>,
    value: T,
    p: K,
    q: K,
) where
    T: Float + AddAssign,
    K: Copy + Hash + Num + Ord + Step,
{
    let quarter = T::from(0.25).unwrap();

    // canonical ordering means q>p
    hamil.add_or_insert(&PauliCode::Identity, value * quarter);
    hamil.add_or_insert(&PauliCode::with_pauli(p, Pauli::Z), value * -quarter);
    hamil.add_or_insert(&PauliCode::with_pauli(q, Pauli::Z), value * -quarter);
    hamil.add_or_insert(
        &PauliCode::from([(p, Pauli::Z), (q, Pauli::Z)]),
        value * quarter,
    );
}

fn update_hamil_two_pqs<T, K>(
    hamil: &mut Hamil<T, PauliCode<K>>,
    value: T,
    p: K,
    q: K,
    s: K,
) where
    T: Float + AddAssign,
    K: Copy + Hash + Num + Ord + Step,
{
    let quarter = T::from(0.25).unwrap();
    let mut code =
        PauliCode::from((p + K::one()..s - K::one()).map(|i| (i, Pauli::Z)));
    code.update(p, Pauli::X);
    code.update(s, Pauli::X);
    hamil.add_or_insert(&code, value * quarter);
    code.update(p, Pauli::Y);
    code.update(s, Pauli::Y);
    hamil.add_or_insert(&code, value * quarter);

    code.update(q, Pauli::Z);
    code.update(p, Pauli::X);
    code.update(s, Pauli::X);
    hamil.add_or_insert(&code, value * -quarter);
    code.update(p, Pauli::Y);
    code.update(s, Pauli::Y);
    hamil.add_or_insert(&code, value * -quarter);
}

fn update_hamil_two_pqrs<T, K>(
    hamil: &mut Hamil<T, PauliCode<K>>,
    value: T,
    p: K,
    q: K,
    s: K,
    r: K,
) where
    T: Float + AddAssign,
    K: Copy + Hash + Num + Ord + Step,
{
    let eighth = T::from(0.125).unwrap();

    let mut code = PauliCode::new();
    for i in p + K::one()..q - K::one() {
        code.update(i, Pauli::Z);
    }
    for i in s + K::one()..r - K::one() {
        code.update(i, Pauli::Z);
    }

    code.update(p, Pauli::X);
    code.update(q, Pauli::X);
    code.update(r, Pauli::X);
    code.update(s, Pauli::X);
    hamil.add_or_insert(&code, value * eighth);

    code.update(p, Pauli::X);
    code.update(q, Pauli::X);
    code.update(r, Pauli::Y);
    code.update(s, Pauli::Y);
    hamil.add_or_insert(&code, value * -eighth);

    code.update(p, Pauli::X);
    code.update(q, Pauli::Y);
    code.update(r, Pauli::X);
    code.update(s, Pauli::Y);
    hamil.add_or_insert(&code, value * eighth);

    code.update(p, Pauli::Y);
    code.update(q, Pauli::X);
    code.update(r, Pauli::X);
    code.update(s, Pauli::Y);
    hamil.add_or_insert(&code, value * eighth);

    code.update(p, Pauli::Y);
    code.update(q, Pauli::X);
    code.update(r, Pauli::Y);
    code.update(s, Pauli::X);
    hamil.add_or_insert(&code, value * eighth);

    code.update(p, Pauli::Y);
    code.update(q, Pauli::Y);
    code.update(r, Pauli::X);
    code.update(s, Pauli::X);
    hamil.add_or_insert(&code, value * -eighth);

    code.update(p, Pauli::X);
    code.update(q, Pauli::Y);
    code.update(r, Pauli::Y);
    code.update(s, Pauli::X);
    hamil.add_or_insert(&code, value * eighth);

    code.update(p, Pauli::Y);
    code.update(q, Pauli::Y);
    code.update(r, Pauli::Y);
    code.update(s, Pauli::Y);
    hamil.add_or_insert(&code, value * eighth);
}
