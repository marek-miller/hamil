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

    #[must_use]
    pub fn variants() -> std::array::IntoIter<Self, 2> {
        [Self::Down, Self::Up].into_iter()
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

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum FermiCode<K> {
    Offset,
    One(Orbital<K>, Orbital<K>),
    Two((Orbital<K>, Orbital<K>), (Orbital<K>, Orbital<K>)),
}

impl<K> FermiCode<K>
where
    K: Num + PartialOrd + Copy,
{
    pub fn one(
        cr: Orbital<K>,
        an: Orbital<K>,
    ) -> Option<Self> {
        (cr.enumerate() <= an.enumerate()).then_some(FermiCode::One(cr, an))
    }

    pub fn two(
        cr: (Orbital<K>, Orbital<K>),
        an: (Orbital<K>, Orbital<K>),
    ) -> Option<Self> {
        (cr.0.enumerate() < cr.1.enumerate()
            && an.0.enumerate() > an.1.enumerate()
            && cr.0.enumerate() <= an.1.enumerate())
        .then_some(FermiCode::Two(cr, an))
    }
}

#[derive(Debug, Default, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Pauli {
    #[default]
    X,
    Y,
    Z,
}

impl Pauli {
    #[must_use]
    pub fn variants() -> std::array::IntoIter<Self, 3> {
        [Self::X, Self::Y, Self::Z].into_iter()
    }
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
    pub fn replace(
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
    pub fn add_to(
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
                FermiCode::Offset => {
                    hamil.add_to(&PauliCode::Identity, value);
                }
                FermiCode::One(cr, an) => {
                    let (p, q) = (cr.enumerate(), an.enumerate());
                    update_hamil_one_pq(&mut hamil, value, p, q);
                }
                FermiCode::Two(cr, an) => {
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
        hamil.add_to(&PauliCode::Identity, value * half);
        hamil.add_to(&PauliCode::with_pauli(q, Pauli::Z), value * -half);
    } else {
        // canonical ordering means p<=q
        let mut code = PauliCode::from(
            (p + K::one()..q - K::one()).map(|i| (i, Pauli::Z)),
        );
        code.update(p, Pauli::X);
        code.update(q, Pauli::X);
        hamil.add_to(&code, value * half);
        code.update(p, Pauli::Y);
        code.update(q, Pauli::Y);
        hamil.add_to(&code, value * half);
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
    hamil.add_to(&PauliCode::Identity, value * quarter);
    hamil.add_to(&PauliCode::with_pauli(p, Pauli::Z), value * -quarter);
    hamil.add_to(&PauliCode::with_pauli(q, Pauli::Z), value * -quarter);
    hamil.add_to(
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
    hamil.add_to(&code, value * quarter);
    code.update(p, Pauli::Y);
    code.update(s, Pauli::Y);
    hamil.add_to(&code, value * quarter);

    code.update(q, Pauli::Z);
    code.update(p, Pauli::X);
    code.update(s, Pauli::X);
    hamil.add_to(&code, value * -quarter);
    code.update(p, Pauli::Y);
    code.update(s, Pauli::Y);
    hamil.add_to(&code, value * -quarter);
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
    hamil.add_to(&code, value * eighth);

    code.update(p, Pauli::X);
    code.update(q, Pauli::X);
    code.update(r, Pauli::Y);
    code.update(s, Pauli::Y);
    hamil.add_to(&code, value * -eighth);

    code.update(p, Pauli::X);
    code.update(q, Pauli::Y);
    code.update(r, Pauli::X);
    code.update(s, Pauli::Y);
    hamil.add_to(&code, value * eighth);

    code.update(p, Pauli::Y);
    code.update(q, Pauli::X);
    code.update(r, Pauli::X);
    code.update(s, Pauli::Y);
    hamil.add_to(&code, value * eighth);

    code.update(p, Pauli::Y);
    code.update(q, Pauli::X);
    code.update(r, Pauli::Y);
    code.update(s, Pauli::X);
    hamil.add_to(&code, value * eighth);

    code.update(p, Pauli::Y);
    code.update(q, Pauli::Y);
    code.update(r, Pauli::X);
    code.update(s, Pauli::X);
    hamil.add_to(&code, value * -eighth);

    code.update(p, Pauli::X);
    code.update(q, Pauli::Y);
    code.update(r, Pauli::Y);
    code.update(s, Pauli::X);
    hamil.add_to(&code, value * eighth);

    code.update(p, Pauli::Y);
    code.update(q, Pauli::Y);
    code.update(r, Pauli::Y);
    code.update(s, Pauli::Y);
    hamil.add_to(&code, value * eighth);
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum MullikenCode<K> {
    Coulomb,
    One(K, K),
    Two(K, K, K, K),
}

impl<K> MullikenCode<K>
where
    K: PartialOrd,
{
    pub fn one(
        i: K,
        j: K,
    ) -> Option<Self> {
        (i >= j).then_some(MullikenCode::One(i, j))
    }

    pub fn two(
        i: K,
        j: K,
        k: K,
        l: K,
    ) -> Option<Self> {
        Some(MullikenCode::Two(i, j, k, l))
    }
}

impl<T, K> From<Hamil<T, MullikenCode<K>>> for Hamil<T, FermiCode<K>>
where
    T: AddAssign + Copy,
    K: Copy + Eq + Hash + Num + PartialOrd,
{
    fn from(mull_hamil: Hamil<T, MullikenCode<K>>) -> Self {
        let mut hamil = FermiHamil::new();
        for (mull_code, value) in mull_hamil.terms {
            match mull_code {
                MullikenCode::Coulomb => {
                    hamil.add_to(&FermiCode::Offset, value);
                }
                MullikenCode::One(i, j) => {
                    for spin in Spin::variants() {
                        let p = Orbital::new(j, spin);
                        let q = Orbital::new(i, spin);
                        hamil.add_to(&FermiCode::one(p, q).unwrap(), value);
                    }
                }
                // we change from chemists' to Dirac convention here:
                MullikenCode::Two(i, k, l, j) => {
                    use Spin::{
                        Down,
                        Up,
                    };
                    for (s1, s2) in
                        [(Down, Down), (Down, Up), (Up, Down), (Up, Up)]
                    {
                        for (is, js, ks, ls) in [
                            (i, j, l, k),
                            (j, i, k, l),
                            (j, i, l, k),
                            (k, l, i, j),
                            (k, l, j, i),
                            (l, k, j, i),
                            (l, k, i, j),
                        ] {
                            let p = Orbital::new(is - K::one(), s1);
                            let q = Orbital::new(js - K::one(), s2);
                            let r = Orbital::new(ks - K::one(), s2);
                            let s = Orbital::new(ls - K::one(), s1);

                            if let Some(code) = FermiCode::two((p, q), (r, s)) {
                                hamil.add_to(&code, value);
                            }
                        }
                    }
                }
            }
        }
        hamil
    }
}
