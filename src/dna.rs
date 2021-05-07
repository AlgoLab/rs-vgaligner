// Reverse-complement adapted from: https://github.com/delagoya/rusty-bio

/// This function returns the reverse-complement of a sequence
/// given as input.
pub fn reverse_complement(sequence: &str) -> String {
    let mut rc_seq = String::with_capacity(sequence.len());

    // Iterate over the bases of the sequence in reverse
    for base in sequence.chars().rev() {
        match is_dna(base) {
            false => panic!("Input sequence base is not DNA: {}", base),
            true => rc_seq.push(switch_base(base)),
        }
    }

    rc_seq
}

fn switch_base(c: char) -> char {
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
        _ => 'N',
    }
}

fn is_dna(base: char) -> bool {
    match base {
        'A' | 'a' | 'C' | 'c' | 'G' | 'g' | 'T' | 't' | 'U' | 'u' | 'N' => true,
        _ => false,
    }
}

#[test]
fn test_is_dna() {
    assert!(is_dna('A'))
}

#[test]
fn test_revcomp() {
    // ATGC =>  CGTA => GCAT
    assert_eq!("GCAT".to_string(), reverse_complement("ATGC"))
}
