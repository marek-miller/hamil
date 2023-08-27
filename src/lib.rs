#![feature(step_trait)]

use std::{
    collections::{
        BTreeMap,
        HashMap,
    },
    hash::Hash,
    iter::Step,
    ops::{
        Add,
        AddAssign,
    },
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
    j: K,
    s: Spin,
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
pub struct Fermi<K> {
    kind: FermiCode<K>,
}

impl<K> Default for Fermi<K> {
    fn default() -> Self {
        Self {
            kind: FermiCode::default(),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
enum FermiCode<K> {
    #[default]
    Offset,
    One(Orbital<K>, Orbital<K>),
    Two((Orbital<K>, Orbital<K>), (Orbital<K>, Orbital<K>)),
}

impl<K> Fermi<K>
where
    K: Num + PartialOrd + Copy,
{
    #[must_use]
    pub fn offset() -> Self {
        Self {
            kind: FermiCode::Offset,
        }
    }

    pub fn is_offset(&self) -> bool {
        matches!(self.kind, FermiCode::Offset)
    }

    pub fn one(
        cr: Orbital<K>,
        an: Orbital<K>,
    ) -> Option<Self> {
        (cr.enumerate() <= an.enumerate()).then_some(Self {
            kind: FermiCode::One(cr, an),
        })
    }

    pub fn is_one(&self) -> bool {
        matches!(self.kind, FermiCode::One(..))
    }

    pub fn two(
        cr: (Orbital<K>, Orbital<K>),
        an: (Orbital<K>, Orbital<K>),
    ) -> Option<Self> {
        (cr.0.enumerate() < cr.1.enumerate()
            && an.0.enumerate() > an.1.enumerate()
            && cr.0.enumerate() <= an.1.enumerate())
        .then_some(Self {
            kind: FermiCode::Two(cr, an),
        })
    }

    pub fn is_two(&self) -> bool {
        matches!(self.kind, FermiCode::Two(..))
    }

    fn kind(&self) -> &FermiCode<K> {
        &self.kind
    }
}

#[derive(Debug, Default, Clone, Copy, PartialEq, Eq, Hash)]
pub enum PauliOp {
    #[default]
    X,
    Y,
    Z,
}

impl PauliOp {
    #[must_use]
    pub fn variants() -> std::array::IntoIter<Self, 3> {
        [Self::X, Self::Y, Self::Z].into_iter()
    }
}

macro_rules! impl_enumerate_pauli {
    ($($Typ:ty)*) => {
        $(
            impl Enumerate<$Typ> for PauliOp {
                fn enumerate(&self) -> $Typ {
                    use PauliOp::*;
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
pub struct Pauli<K> {
    kind: PauliCode<K>,
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
enum PauliCode<K> {
    Identity,
    Product(BTreeMap<K, PauliOp>),
}

impl<K> Default for Pauli<K> {
    fn default() -> Self {
        Self::new()
    }
}

impl<K> Pauli<K> {
    #[must_use]
    pub fn new() -> Self {
        Self {
            kind: PauliCode::Identity,
        }
    }

    fn kind(&self) -> &PauliCode<K> {
        &self.kind
    }

    fn kind_mut(&mut self) -> &mut PauliCode<K> {
        &mut self.kind
    }

    #[must_use]
    pub fn identity() -> Self {
        Self {
            kind: PauliCode::Identity,
        }
    }

    #[must_use]
    pub fn is_identity(&self) -> bool {
        matches!(self.kind, PauliCode::Identity)
    }
}

impl<K> Pauli<K>
where
    K: Ord,
{
    pub fn with_pauli(
        index: K,
        pauli: PauliOp,
    ) -> Self {
        let mut tree = BTreeMap::new();
        tree.insert(index, pauli);
        Self {
            kind: PauliCode::Product(tree),
        }
    }

    pub fn update(
        &mut self,
        index: K,
        pauli: PauliOp,
    ) {
        match self.kind_mut() {
            PauliCode::Identity => {
                *self = Pauli::with_pauli(index, pauli);
            }
            PauliCode::Product(tree) => {
                tree.insert(index, pauli);
            }
        }
    }
}

impl<I, K> From<I> for Pauli<K>
where
    I: IntoIterator<Item = (K, PauliOp)>,
    K: Ord,
{
    fn from(value: I) -> Self {
        let mut code = Pauli::new();
        for (index, pauli) in value {
            code.update(index, pauli);
        }
        code
    }
}

macro_rules! impl_enumerate_paulicode {
    ($($Typ:ty)*) => {
        $(
            impl Enumerate<$Typ> for Pauli<$Typ> {
                fn enumerate(&self) -> $Typ {
                    let mut count = 0;
                    if let PauliCode::Product(tree) = self.kind() {
                        for (key, value) in tree {
                            count += key * 3 + <PauliOp as Enumerate<$Typ>>::enumerate(value)
                        }
                    }
                    count
                }
            }
        )*
    };
}

impl_enumerate_paulicode!(u8 u16 u32 u64 usize);

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Mulliken<K> {
    kind: MullikenCode<K>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
enum MullikenCode<K> {
    #[default]
    Coulomb,
    One(K, K),
    Two(K, K, K, K),
}

impl<K> Default for Mulliken<K> {
    fn default() -> Self {
        Self {
            kind: MullikenCode::Coulomb,
        }
    }
}

impl<K> Mulliken<K>
where
    K: PartialOrd,
{
    #[must_use]
    pub fn coulomb() -> Self {
        Self {
            kind: MullikenCode::Coulomb,
        }
    }

    pub fn one(
        i: K,
        j: K,
    ) -> Option<Self> {
        (i >= j).then_some(Self {
            kind: MullikenCode::One(i, j),
        })
    }

    pub fn two(
        i: K,
        j: K,
        k: K,
        l: K,
    ) -> Option<Self> {
        Some(Self {
            kind: MullikenCode::Two(i, j, k, l),
        })
    }
}

pub trait Code: Clone + Default + Eq + Hash {}

impl<K> Code for Fermi<K> where K: Clone + Eq + Hash {}

impl<K> Code for Pauli<K> where K: Clone + Eq + Hash {}

impl<K> Code for Mulliken<K> where K: Clone + Eq + Hash {}

#[derive(Debug)]
pub struct Terms<T, K> {
    terms: HashMap<K, T>,
}

impl<T, K> Default for Terms<T, K>
where
    K: Code,
{
    fn default() -> Self {
        Self::new()
    }
}

impl<T, K> Terms<T, K>
where
    K: Code,
{
    #[must_use]
    pub fn new() -> Self {
        Self {
            terms: HashMap::new(),
        }
    }
}

impl<T, K> Terms<T, K>
where
    K: Code,
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

impl<T, K> Terms<T, K>
where
    T: AddAssign + Copy,
    K: Code,
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

impl<T, K> From<Terms<T, Fermi<K>>> for Terms<T, Pauli<K>>
where
    T: Float + AddAssign,
    K: Copy + Hash + Num + Ord + Step,
{
    fn from(fermi_hamil: Terms<T, Fermi<K>>) -> Self {
        let mut hamil = Terms::new();
        for (fermi_code, value) in fermi_hamil.terms {
            match fermi_code.kind() {
                FermiCode::Offset => {
                    hamil.add_to(&Pauli::identity(), value);
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
    hamil: &mut Terms<T, Pauli<K>>,
    value: T,
    p: K,
    q: K,
) where
    T: Float + AddAssign,
    K: Copy + Hash + Num + Ord + Step,
{
    let half = T::from(0.5).unwrap();

    if p == q {
        hamil.add_to(&Pauli::identity(), value * half);
        hamil.add_to(&Pauli::with_pauli(q, PauliOp::Z), value * -half);
    } else {
        // canonical ordering means p<=q
        let mut code =
            Pauli::from((p + K::one()..q - K::one()).map(|i| (i, PauliOp::Z)));
        code.update(p, PauliOp::X);
        code.update(q, PauliOp::X);
        hamil.add_to(&code, value * half);
        code.update(p, PauliOp::Y);
        code.update(q, PauliOp::Y);
        hamil.add_to(&code, value * half);
    }
}

fn update_hamil_two_pq<T, K>(
    hamil: &mut Terms<T, Pauli<K>>,
    value: T,
    p: K,
    q: K,
) where
    T: Float + AddAssign,
    K: Copy + Hash + Num + Ord + Step,
{
    let quarter = T::from(0.25).unwrap();

    // canonical ordering means q>p
    hamil.add_to(&Pauli::identity(), value * quarter);
    hamil.add_to(&Pauli::with_pauli(p, PauliOp::Z), value * -quarter);
    hamil.add_to(&Pauli::with_pauli(q, PauliOp::Z), value * -quarter);
    hamil.add_to(
        &Pauli::from([(p, PauliOp::Z), (q, PauliOp::Z)]),
        value * quarter,
    );
}

fn update_hamil_two_pqs<T, K>(
    hamil: &mut Terms<T, Pauli<K>>,
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
        Pauli::from((p + K::one()..s - K::one()).map(|i| (i, PauliOp::Z)));
    code.update(p, PauliOp::X);
    code.update(s, PauliOp::X);
    hamil.add_to(&code, value * quarter);
    code.update(p, PauliOp::Y);
    code.update(s, PauliOp::Y);
    hamil.add_to(&code, value * quarter);

    code.update(q, PauliOp::Z);
    code.update(p, PauliOp::X);
    code.update(s, PauliOp::X);
    hamil.add_to(&code, value * -quarter);
    code.update(p, PauliOp::Y);
    code.update(s, PauliOp::Y);
    hamil.add_to(&code, value * -quarter);
}

fn update_hamil_two_pqrs<T, K>(
    hamil: &mut Terms<T, Pauli<K>>,
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

    let mut code = Pauli::new();
    for i in p + K::one()..q - K::one() {
        code.update(i, PauliOp::Z);
    }
    for i in s + K::one()..r - K::one() {
        code.update(i, PauliOp::Z);
    }

    code.update(p, PauliOp::X);
    code.update(q, PauliOp::X);
    code.update(r, PauliOp::X);
    code.update(s, PauliOp::X);
    hamil.add_to(&code, value * eighth);

    code.update(p, PauliOp::X);
    code.update(q, PauliOp::X);
    code.update(r, PauliOp::Y);
    code.update(s, PauliOp::Y);
    hamil.add_to(&code, value * -eighth);

    code.update(p, PauliOp::X);
    code.update(q, PauliOp::Y);
    code.update(r, PauliOp::X);
    code.update(s, PauliOp::Y);
    hamil.add_to(&code, value * eighth);

    code.update(p, PauliOp::Y);
    code.update(q, PauliOp::X);
    code.update(r, PauliOp::X);
    code.update(s, PauliOp::Y);
    hamil.add_to(&code, value * eighth);

    code.update(p, PauliOp::Y);
    code.update(q, PauliOp::X);
    code.update(r, PauliOp::Y);
    code.update(s, PauliOp::X);
    hamil.add_to(&code, value * eighth);

    code.update(p, PauliOp::Y);
    code.update(q, PauliOp::Y);
    code.update(r, PauliOp::X);
    code.update(s, PauliOp::X);
    hamil.add_to(&code, value * -eighth);

    code.update(p, PauliOp::X);
    code.update(q, PauliOp::Y);
    code.update(r, PauliOp::Y);
    code.update(s, PauliOp::X);
    hamil.add_to(&code, value * eighth);

    code.update(p, PauliOp::Y);
    code.update(q, PauliOp::Y);
    code.update(r, PauliOp::Y);
    code.update(s, PauliOp::Y);
    hamil.add_to(&code, value * eighth);
}

impl<T, K> From<Terms<T, Mulliken<K>>> for Terms<T, Fermi<K>>
where
    T: AddAssign + Copy,
    K: Copy + Eq + Hash + Num + PartialOrd,
{
    fn from(mull_hamil: Terms<T, Mulliken<K>>) -> Self {
        let mut hamil = Terms::new();
        for (mull_code, value) in mull_hamil.terms {
            match mull_code.kind {
                MullikenCode::Coulomb => {
                    hamil.add_to(&Fermi::offset(), value);
                }
                MullikenCode::One(i, j) => {
                    for spin in [Spin::Down, Spin::Up] {
                        let p = Orbital::new(j, spin);
                        let q = Orbital::new(i, spin);
                        hamil.add_to(&Fermi::one(p, q).unwrap(), value);
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

                            if let Some(code) = Fermi::two((p, q), (r, s)) {
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

#[derive(Debug)]
pub enum Hamil<T, K> {
    Offset(T),
    Terms(Terms<T, K>),
    Sum(Box<Self>, Box<Self>),
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

impl<T, K> From<Hamil<T, Fermi<K>>> for Hamil<T, Pauli<K>>
where
    T: Float + AddAssign,
    K: Copy + Hash + Num + Ord + Step,
{
    fn from(value: Hamil<T, Fermi<K>>) -> Self {
        match value {
            Hamil::Offset(x) => Hamil::Offset(x),
            Hamil::Terms(terms) => Hamil::Terms(terms.into()),
            Hamil::Sum(h1, h2) => Self::from(*h1) + Self::from(*h2),
        }
    }
}

impl<T, K> From<Hamil<T, Mulliken<K>>> for Hamil<T, Fermi<K>>
where
    T: AddAssign + Copy,
    K: Copy + Eq + Hash + Num + PartialOrd,
{
    fn from(value: Hamil<T, Mulliken<K>>) -> Self {
        match value {
            Hamil::Offset(x) => Hamil::Offset(x),
            Hamil::Terms(terms) => Hamil::Terms(terms.into()),
            Hamil::Sum(h1, h2) => Self::from(*h1) + Self::from(*h2),
        }
    }
}
