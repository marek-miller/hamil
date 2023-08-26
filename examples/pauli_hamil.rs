use hamil::{
    FermiCode,
    FermiHamil,
    Orbital,
    PauliHamil,
    Spin,
};

fn main() {
    let mut fermi = FermiHamil::new();

    let o1 = Orbital::new(0, Spin::Down);
    let o2 = Orbital::new(1, Spin::Down);
    fermi.update(FermiCode::one(o1, o2).unwrap(), 0.5);

    let o1 = Orbital::new(0, Spin::Up);
    let o2 = Orbital::new(1, Spin::Up);
    fermi.update(FermiCode::one(o1, o2).unwrap(), 0.5);

    let o1 = Orbital::new(0, Spin::Down);
    let o2 = Orbital::new(1, Spin::Down);
    let o3 = Orbital::new(2, Spin::Down);
    let o4 = Orbital::new(1, Spin::Down);
    fermi.update(FermiCode::two((o1, o2), (o3, o4)).unwrap(), 0.5);

    let o1 = Orbital::new(0, Spin::Up);
    let o2 = Orbital::new(1, Spin::Up);
    let o3 = Orbital::new(2, Spin::Up);
    let o4 = Orbital::new(1, Spin::Up);
    fermi.update(FermiCode::two((o1, o2), (o3, o4)).unwrap(), 0.5);

    println!("{:?}", fermi);
    let pauli = PauliHamil::from(fermi);
    println!("{:?}", pauli);
}
