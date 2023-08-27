use hamil::{
    FermiHamil,
    Hamil,
    MullikenCode,
    PauliHamil,
};

fn main() {
    let mut hamil = Hamil::new();
    hamil.replace(&MullikenCode::Coulomb, 0.1);
    hamil.replace(&MullikenCode::one(1, 0).unwrap(), 0.2);
    hamil.replace(&MullikenCode::two(2, 1, 1, 1).unwrap(), 0.3);

    println!("{hamil:?}");
    let fermi = FermiHamil::from(hamil);
    println!("{fermi:?}");
    let pauli = PauliHamil::from(fermi);
    println!("{pauli:?}");
}
