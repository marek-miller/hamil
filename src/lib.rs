#![feature(step_trait)]

use std::{
    collections::HashMap,
    hash::Hash,
    iter::Step,
    ops::{
        Add,
        AddAssign,
    },
};

use num::Num;

pub trait Idx: Clone + Copy + Eq + Hash + Num + Ord + Step {}
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

#[derive(Debug)]
pub struct SumRepr<T, K>
where
    T: Float,
    K: Code,
{
    repr: HashMap<K, T>,
}

impl<T, K> Default for SumRepr<T, K>
where
    T: Float,
    K: Code,
{
    fn default() -> Self {
        Self::new()
    }
}

impl<T, K> SumRepr<T, K>
where
    T: Float,
    K: Code,
{
    #[must_use]
    pub fn new() -> Self {
        Self {
            repr: HashMap::new(),
        }
    }

    pub fn add(
        &mut self,
        code: &K,
        value: T,
    ) {
        if let Some(t) = self.repr.get_mut(code) {
            *t += value;
        } else {
            self.repr.insert(code.clone(), value);
        }
    }

    pub fn insert(
        &mut self,
        code: K,
        value: T,
    ) {
        self.repr.insert(code, value);
    }

    pub fn remove(
        &mut self,
        code: &K,
    ) -> Option<T> {
        self.repr.remove(code)
    }

    #[must_use]
    pub fn as_map(&self) -> &HashMap<K, T> {
        &self.repr
    }

    pub fn as_map_mut(&mut self) -> &mut HashMap<K, T> {
        &mut self.repr
    }
}

pub trait Terms<T, K>
where
    T: Float,
    K: Code,
{
    fn add_to(
        &mut self,
        repr: &mut SumRepr<T, K>,
    );
}

impl<T, K> Terms<T, K> for SumRepr<T, K>
where
    T: Float,
    K: Code,
{
    fn add_to(
        &mut self,
        repr: &mut SumRepr<T, K>,
    ) {
        for (code, value) in self.as_map() {
            repr.add(code, *value);
        }
    }
}

#[derive(Debug, Clone)]
pub struct TermsIter<I> {
    iter: I,
}

impl<I> TermsIter<I>
where
    I: Iterator + Clone,
{
    pub fn new(iter: I) -> Self {
        Self {
            iter,
        }
    }
}

impl<I, T, K> Terms<T, K> for TermsIter<I>
where
    I: Iterator + Clone,
    I::Item: Terms<T, K>,
    T: Float,
    K: Code,
{
    fn add_to(
        &mut self,
        repr: &mut SumRepr<T, K>,
    ) {
        self.iter.clone().for_each(|mut s| s.add_to(repr));
    }
}

pub trait TermsIterator<T, K>: Sized + Iterator + Clone
where
    <Self as Iterator>::Item: Terms<T, K>,
    T: Float,
    K: Code,
{
    fn into_terms(self) -> TermsIter<Self> {
        TermsIter::new(self)
    }
}

impl<T, K, I> TermsIterator<T, K> for I
where
    I: Iterator + Clone + Sized,
    I::Item: Terms<T, K>,
    T: Float,
    K: Code,
{
}

pub enum Hamil<T, K> {
    Offset(T),
    Sum(Box<Self>, Box<Self>),
    Terms(Box<dyn Terms<T, K>>),
}

impl<T, K> Default for Hamil<T, K>
where
    T: Default,
{
    fn default() -> Self {
        Self::Offset(T::default())
    }
}

impl<T, K> Add for Hamil<T, K> {
    type Output = Self;

    fn add(
        self,
        rhs: Self,
    ) -> Self::Output {
        Self::Sum(Box::new(self), Box::new(rhs))
    }
}

impl<T, K> Hamil<T, K> {
    pub fn add_hamil(
        self,
        other: Self,
    ) -> Self {
        self + other
    }

    pub fn add_offset(
        self,
        value: T,
    ) -> Self {
        self + Self::Offset(value)
    }

    pub fn add_terms(
        self,
        terms: Box<dyn Terms<T, K>>,
    ) -> Self
    where
        T: Float,
        K: Code,
    {
        self + Self::Terms(terms)
    }
}

impl<T, K> Terms<T, K> for Hamil<T, K>
where
    T: Float,
    K: Code,
{
    fn add_to(
        &mut self,
        repr: &mut SumRepr<T, K>,
    ) {
        match self {
            Self::Offset(t) => {
                repr.add(&K::IDENTITY, *t);
            }
            Self::Terms(terms) => terms.add_to(repr),
            Self::Sum(h1, h2) => {
                h1.add_to(repr);
                h2.add_to(repr);
            }
        }
    }
}

impl<T, K> From<Hamil<T, K>> for SumRepr<T, K>
where
    T: Float,
    K: Code,
{
    fn from(value: Hamil<T, K>) -> Self {
        let mut hamil = value;
        let mut repr = SumRepr::new();
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

pub trait Enum<U>
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

impl<U> Enum<U> for Spin
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

impl<U> Enum<U> for Orbital<U>
where
    U: Idx,
{
    fn enumerate(&self) -> U {
        self.h()
    }
}

pub struct OneElectron<T, U> {
    orbs:  (Orbital<U>, Orbital<U>),
    value: T,
}

impl<T, U> OneElectron<T, U>
where
    U: Idx,
{
    pub fn new(
        cr: Orbital<U>,
        an: Orbital<U>,
        value: T,
    ) -> Option<Self> {
        (cr.enumerate() <= an.enumerate()).then_some(Self {
            orbs: (cr, an),
            value,
        })
    }
}

impl<T, U> Terms<T, PauliCode<U>> for OneElectron<T, U>
where
    T: Float,
    U: Idx,
{
    fn add_to(
        &mut self,
        repr: &mut SumRepr<T, PauliCode<U>>,
    ) {
        let half = T::from(0.5).unwrap();

        let (p, q) = (self.orbs.0.enumerate(), self.orbs.1.enumerate());
        if p == q {
            repr.add(&PauliCode::Identity, self.value * half);
            repr.add(&PauliCode::with_pauli(q, Pauli::Z), self.value * -half);
        } else {
            let mut code = PauliCode::from(
                (p + U::one()..q - U::one()).map(|i| (i, Pauli::Z)),
            );
            code.update(p, Pauli::X);
            code.update(q, Pauli::X);
            repr.add(&code, self.value * half);
            code.update(p, Pauli::Y);
            code.update(q, Pauli::Y);
            repr.add(&code, self.value * half);
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct TwoElectron<T, U> {
    orbs:  ((Orbital<U>, Orbital<U>), (Orbital<U>, Orbital<U>)),
    value: T,
}

impl<T, U> TwoElectron<T, U>
where
    U: Idx,
{
    pub fn new(
        cr: (Orbital<U>, Orbital<U>),
        an: (Orbital<U>, Orbital<U>),
        value: T,
    ) -> Option<Self> {
        (cr.0.enumerate() < cr.1.enumerate()
            && an.0.enumerate() > an.1.enumerate()
            && cr.0.enumerate() <= an.1.enumerate())
        .then_some(Self {
            orbs: (cr, an),
            value,
        })
    }
}

impl<T, U> Terms<T, PauliCode<U>> for TwoElectron<T, U>
where
    T: Float,
    U: Idx,
{
    fn add_to(
        &mut self,
        repr: &mut SumRepr<T, PauliCode<U>>,
    ) {
        let (p, q, r, s) = (
            self.orbs.0 .0.enumerate(),
            self.orbs.0 .1.enumerate(),
            self.orbs.1 .0.enumerate(),
            self.orbs.1 .1.enumerate(),
        );

        if p == s && q == r {
            let quarter = T::from(0.25).unwrap();

            // canonical ordering means q>p
            repr.add(&PauliCode::IDENTITY, self.value * quarter);
            repr.add(
                &PauliCode::with_pauli(p, Pauli::Z),
                self.value * -quarter,
            );
            repr.add(
                &PauliCode::with_pauli(q, Pauli::Z),
                self.value * -quarter,
            );
            repr.add(
                &PauliCode::from([(p, Pauli::Z), (q, Pauli::Z)]),
                self.value * quarter,
            );
        } else if q == r {
            let quarter = T::from(0.25).unwrap();
            let mut code = PauliCode::from(
                (p + U::one()..s - U::one()).map(|i| (i, Pauli::Z)),
            );
            code.update(p, Pauli::X);
            code.update(s, Pauli::X);
            repr.add(&code, self.value * quarter);
            code.update(p, Pauli::Y);
            code.update(s, Pauli::Y);
            repr.add(&code, self.value * quarter);

            code.update(q, Pauli::Z);
            code.update(p, Pauli::X);
            code.update(s, Pauli::X);
            repr.add(&code, self.value * -quarter);
            code.update(p, Pauli::Y);
            code.update(s, Pauli::Y);
            repr.add(&code, self.value * -quarter);
        } else {
            let eighth = T::from(0.125).unwrap();

            let mut code = PauliCode::new();
            for i in p + U::one()..q - U::one() {
                code.update(i, Pauli::Z);
            }
            for i in s + U::one()..r - U::one() {
                code.update(i, Pauli::Z);
            }

            code.update(p, Pauli::X);
            code.update(q, Pauli::X);
            code.update(r, Pauli::X);
            code.update(s, Pauli::X);
            repr.add(&code, self.value * eighth);

            code.update(p, Pauli::X);
            code.update(q, Pauli::X);
            code.update(r, Pauli::Y);
            code.update(s, Pauli::Y);
            repr.add(&code, self.value * -eighth);

            code.update(p, Pauli::X);
            code.update(q, Pauli::Y);
            code.update(r, Pauli::X);
            code.update(s, Pauli::Y);
            repr.add(&code, self.value * eighth);

            code.update(p, Pauli::Y);
            code.update(q, Pauli::X);
            code.update(r, Pauli::X);
            code.update(s, Pauli::Y);
            repr.add(&code, self.value * eighth);

            code.update(p, Pauli::Y);
            code.update(q, Pauli::X);
            code.update(r, Pauli::Y);
            code.update(s, Pauli::X);
            repr.add(&code, self.value * eighth);

            code.update(p, Pauli::Y);
            code.update(q, Pauli::Y);
            code.update(r, Pauli::X);
            code.update(s, Pauli::X);
            repr.add(&code, self.value * -eighth);

            code.update(p, Pauli::X);
            code.update(q, Pauli::Y);
            code.update(r, Pauli::Y);
            code.update(s, Pauli::X);
            repr.add(&code, self.value * eighth);

            code.update(p, Pauli::Y);
            code.update(q, Pauli::Y);
            code.update(r, Pauli::Y);
            code.update(s, Pauli::Y);
            repr.add(&code, self.value * eighth);
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum MullikenCode<U> {
    One(U, U),
    Two(U, U, U, U),
}

pub struct Mulliken<T, U> {
    code:  MullikenCode<U>,
    value: T,
}

impl<T, U> Terms<T, PauliCode<U>> for Mulliken<T, U>
where
    T: Float,
    U: Idx,
{
    fn add_to(
        &mut self,
        repr: &mut SumRepr<T, PauliCode<U>>,
    ) {
        match self.code {
            MullikenCode::One(i, j) => {
                for spin in [Spin::Down, Spin::Up] {
                    let p = Orbital::new(j, spin);
                    let q = Orbital::new(i, spin);
                    let mut code = OneElectron::new(p, q, self.value).unwrap();
                    code.add_to(repr);
                }
            }
            // we change from chemists' to Dirac convention here:
            MullikenCode::Two(i, k, l, j) => {
                use Spin::{
                    Down,
                    Up,
                };
                for (s1, s2) in [(Down, Down), (Down, Up), (Up, Down), (Up, Up)]
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
                        let op = Orbital::new(is - U::one(), s1);
                        let oq = Orbital::new(js - U::one(), s2);
                        let or = Orbital::new(ks - U::one(), s2);
                        let os = Orbital::new(ls - U::one(), s1);

                        if let Some(mut code) =
                            TwoElectron::new((op, oq), (or, os), self.value)
                        {
                            code.add_to(repr);
                        }
                    }
                }
            }
        }
    }
}
