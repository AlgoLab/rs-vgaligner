// Reverse-complement adapted from: https://github.com/delagoya/rusty-bio
// This highlights how powerful Rust's match is!
pub fn reverse_complement(sequence : &str) -> String {
    let mut rc_seq = String::with_capacity(sequence.len());

    // Iterate over the bases of the sequence in reverse
    for base in sequence.chars().rev() {
        match is_dna(base) {
            false => panic!("Input sequence base is not DNA: {}", base),
            true => rc_seq.push(switch_base(base))
        }
    }

    rc_seq
}

// By using match bases can be switched or...
fn switch_base(c:char) -> char {
    match c {
        'a' => 't',
        'c' => 'g',
        't' => 'a',
        'g' => 'c',
        'u' => 'a',
        'A' => 'T',
        'C' => 'G',
        'T' => 'A',
        'G' => 'C',
        'U' => 'A',
        _ => 'N'
    }
}

// ...one can easily check if some condition holds or not
fn is_dna(base: char) -> bool {
    match base {
        'A' | 'a' | 'C' | 'c' | 'G' | 'g' | 'T' | 't' | 'U'| 'u'  => true,
        _ => false
    }
}

#[test]
fn test_is_dna() {
    assert!(is_dna('A'))
}

#[test]
#[should_panic]
fn test_is_dna_false() {
    assert!(is_dna('z'))
}

#[test]
fn test_revcomp() {
    // ATGC =>  CGTA => GCAT
    assert_eq!("GCAT".to_string(), reverse_complement("ATGC"))
}