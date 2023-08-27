use hamil::{
    FermiHamil,
    Hamil,
    Mulliken,
    PauliHamil,
};

fn main() {
    let mut hamil = Hamil::new();
    hamil.replace(&Mulliken::coulomb(), 0.1);
    hamil.replace(&Mulliken::one(1, 0).unwrap(), 0.2);
    hamil.replace(&Mulliken::two(2, 1, 1, 1).unwrap(), 0.3);

    println!("{hamil:?}");
    let fermi = FermiHamil::from(hamil);
    println!("{fermi:?}");
    let pauli = PauliHamil::from(fermi);
    println!("{pauli:?}");
}
