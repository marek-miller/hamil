use std::{
    collections::HashMap,
    hash::Hash,
    ops::AddAssign,
};

use num::Num;

pub trait Idx: Clone + Copy + Eq + Hash + Num + Ord {}
impl Idx for u8 {}
impl Idx for u16 {}
impl Idx for u32 {}
impl Idx for u64 {}
impl Idx for usize {}

pub trait Float: num::Float + AddAssign {}
impl Float for f32 {}
impl Float for f64 {}

pub trait Code: Clone + Eq + Hash {
    const IDENTITY: Self;
}

pub trait Terms<T, K>
where
    T: Float,
    K: Code,
{
    fn add_to(
        &mut self,
        repr: &mut HashMap<K, T>,
    );
}

pub enum Hamil<T, K> {
    Offset(T),
    Terms(Box<dyn Terms<T, K>>),
    Sum(Box<Self>, Box<Self>),
}

impl<T, K> Default for Hamil<T, K>
where
    T: Default,
{
    fn default() -> Self {
        Self::Offset(T::default())
    }
}

impl<T, K> Terms<T, K> for Hamil<T, K>
where
    T: Float + AddAssign,
    K: Code,
{
    fn add_to(
        &mut self,
        repr: &mut HashMap<K, T>,
    ) {
        match self {
            Self::Offset(t) => {
                repr.entry(K::IDENTITY)
                    .and_modify(|coeff| *coeff += *t)
                    .or_insert(*t);
            }
            Self::Terms(terms) => terms.add_to(repr),
            Self::Sum(h1, h2) => {
                h1.add_to(repr);
                h2.add_to(repr);
            }
        }
    }
}

impl<T, K> From<Hamil<T, K>> for HashMap<K, T>
where
    T: Float + AddAssign,
    K: Code,
{
    fn from(value: Hamil<T, K>) -> Self {
        let mut hamil = value;
        let mut repr = HashMap::new();
        hamil.add_to(&mut repr);
        repr
    }
}

pub use pauli::{
    Pauli,
    PauliCode,
    PauliProduct,
};

mod pauli {
    use std::collections::BTreeMap;

    use crate::{
        Code,
        Idx,
    };

    #[derive(Debug, Default, Clone, Copy, PartialEq, Eq, Hash)]
    pub enum Pauli {
        #[default]
        X,
        Y,
        Z,
    }

    #[derive(Debug, Clone, PartialEq, Eq, Hash)]
    pub struct PauliProduct<U> {
        tree: BTreeMap<U, Pauli>,
    }

    #[derive(Debug, Clone, PartialEq, Eq, Hash)]
    pub enum PauliCode<U> {
        Identity,
        Product(PauliProduct<U>),
    }

    impl<U> Default for PauliCode<U> {
        fn default() -> Self {
            Self::new()
        }
    }

    impl<U> PauliCode<U> {
        #[must_use]
        pub fn new() -> Self {
            Self::Identity
        }
    }

    impl<U> PauliCode<U>
    where
        U: Idx,
    {
        pub fn with_pauli(
            index: U,
            pauli: Pauli,
        ) -> Self {
            let mut tree = BTreeMap::new();
            tree.insert(index, pauli);
            Self::Product(PauliProduct {
                tree,
            })
        }

        pub fn update(
            &mut self,
            index: U,
            pauli: Pauli,
        ) -> Option<Pauli> {
            match self {
                Self::Identity => {
                    *self = PauliCode::with_pauli(index, pauli);
                    None
                }
                Self::Product(prod) => prod.tree.insert(index, pauli),
            }
        }
    }

    impl<I, U> From<I> for PauliCode<U>
    where
        I: IntoIterator<Item = (U, Pauli)>,
        U: Idx,
    {
        fn from(value: I) -> Self {
            let mut code = PauliCode::new();
            for (index, pauli) in value {
                code.update(index, pauli);
            }
            code
        }
    }

    impl<U> Code for PauliCode<U>
    where
        U: Idx,
    {
        const IDENTITY: Self = Self::Identity;
    }
}

pub trait Enumerate<U>
where
    U: Idx,
{
    fn enumerate(&self) -> U;
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

impl<U> Enumerate<U> for Spin
where
    U: Idx,
{
    fn enumerate(&self) -> U {
        match self {
            Self::Down => U::zero(),
            Self::Up => U::one(),
        }
    }
}

#[derive(Debug, Default, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Orbital<U> {
    j: U,
    s: Spin,
}

impl<U> Orbital<U> {
    pub fn new(
        j: U,
        s: Spin,
    ) -> Self {
        Self {
            j,
            s,
        }
    }
}
impl<U> Orbital<U>
where
    U: Idx,
{
    pub fn h(&self) -> U {
        self.j * (U::one() + U::one()) + self.s.enumerate()
    }

    pub fn g(
        &self,
        num_orbitals: U,
    ) -> U {
        self.j + num_orbitals * self.s.enumerate()
    }
}

impl<U> From<(U, Spin)> for Orbital<U> {
    fn from(value: (U, Spin)) -> Self {
        Self::new(value.0, value.1)
    }
}

impl<U> Enumerate<U> for Orbital<U>
where
    U: Idx,
{
    fn enumerate(&self) -> U {
        self.h()
    }
}

pub struct OneElectron<U>(Orbital<U>, Orbital<U>);

impl<U> OneElectron<U>
where
    U: Idx,
{
    pub fn new(
        cr: Orbital<U>,
        an: Orbital<U>,
    ) -> Option<Self> {
        (cr.enumerate() <= an.enumerate()).then_some(Self(cr, an))
    }
}

impl<I, T, U> Terms<T, PauliCode<U>> for I
where
    I: Iterator<Item = (OneElectron<U>, T)>,
    T: Float,
    U: Idx,
{
    fn add_to(
        &mut self,
        repr: &mut HashMap<PauliCode<U>, T>,
    ) {
        let half = T::from(0.5).unwrap();
        for (elec, t) in self {
            let (p, q) = (elec.0.enumerate(), elec.1.enumerate());
            if p == q {
                repr.entry(PauliCode::Identity)
                    .and_modify(|coeff| *coeff += t * half)
                    .or_insert(t * half);
                todo!()
            } else {
                todo!()
            }
        }
    }
}
