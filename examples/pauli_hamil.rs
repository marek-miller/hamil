use hamil::{
    Fermi,
    Hamil,
    Mulliken,
    Pauli,
    Terms,
};

fn main() {
    let mut terms = Terms::new();
    terms.replace(&Mulliken::coulomb(), 0.1);
    terms.replace(&Mulliken::one(1, 0).unwrap(), 0.2);
    let hamil = Hamil::Terms(terms);
    let mut terms = Terms::new();
    terms.replace(&Mulliken::two(2, 1, 1, 1).unwrap(), 0.3);
    let hamil = hamil + Hamil::Terms(terms);

    println!("{hamil:?}");
    let fermi = Hamil::<_, Fermi<_>>::from(hamil);
    println!("{fermi:?}");
    let pauli = Hamil::<_, Pauli<_>>::from(fermi);
    println!("{pauli:?}");
}
