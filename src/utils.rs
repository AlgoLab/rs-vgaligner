// Adapted from: https://github.com/delagoya/rusty-bio

pub fn reverse_complement(sequence : &str) -> String {
    let mut rc_seq = String::with_capacity(sequence.len());

    for base in sequence.chars().rev() {
        match is_dna(base) {
            false => panic!("Input sequence base is not DNA: {}", base),
            true => rc_seq.push(switch_base(base))
        }
    }

    rc_seq
}

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


