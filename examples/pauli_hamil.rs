// use hamil::{
//     Fermi,
//     Hamil,
//     Mulliken,
//     Pauli,
//     Terms,
// };

// fn main() {
//     let mut terms = Terms::new();
//     terms.replace(&Mulliken::coulomb(), 0.1);
//     terms.replace(&Mulliken::one(1, 0).unwrap(), 0.2);
//     let hamil = Hamil::Terms(terms);
//     let mut terms = Terms::new();
//     terms.replace(&Mulliken::two(2, 1, 1, 1).unwrap(), 0.3);
//     let hamil = hamil + Hamil::Terms(terms);

//     println!("{hamil:?}");
//     let fermi = Hamil::<_, Fermi<_>>::from(hamil);
//     println!("{fermi:?}");
//     let pauli = Hamil::<_, Pauli<_>>::from(fermi);
//     println!("{pauli:?}");
// }

use hamil::{
    Hamil,
    OneElectron,
    Spin,
    SumRepr,
    Terms,
    TwoElectron,
};

fn main() {
    let codes = (0..2u32)
        .map(|x| {
            OneElectron::new(
                (x, Spin::Down).into(),
                (x + 1, Spin::Up).into(),
                0.5,
            )
            .unwrap()
        })
        .collect::<Box<_>>();

    let twoelec = TwoElectron::new(
        ((0u32, Spin::Down).into(), (1, Spin::Up).into()),
        ((1, Spin::Down).into(), (0, Spin::Up).into()),
        0.11,
    )
    .unwrap();
    let hamil = Hamil::Terms(Box::new(codes)) + Hamil::Offset(1.1);
    let mut hamil =
        hamil + Hamil::Offset(0.12) + Hamil::Terms(Box::new(twoelec));
    let mut repr = SumRepr::new();
    hamil.add_to(&mut repr);

    println!("{repr:?}");
}
