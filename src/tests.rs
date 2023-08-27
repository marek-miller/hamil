use crate::Spin;

#[test]
fn spin_intoiter() {
    for spin in Spin::variants() {
        println!("{spin:?}");
    }
}
