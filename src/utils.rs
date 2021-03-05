use bv::*;

// Returns the rank vector of a given bitvector bv
// The rank for an element bv[i] is defined as the number of 1s in [bv[0]...bv[i]]
// We compute the rank for every element in bv, and put these ranks in a Vec
pub fn get_bv_rank(bv : &BitVec) -> Vec<u32> {
    let mut bv_rank = Vec::with_capacity(bv.len() as usize);

    let mut curr_rank = 0;
    for i in 0..bv.len() {
        if bv.get(i) == true {
            curr_rank += 1;
        }
        bv_rank.push(curr_rank);
    }

    bv_rank
}

#[test]
fn test_simple_rank_1_a() {
    assert_eq!(get_bv_rank(&bit_vec![false]), vec![0])
}

#[test]
fn test_simple_rank_1_b() {
    assert_eq!(get_bv_rank(&bit_vec![true]), vec![1])
}

#[test]
fn test_simple_rank_2() {
    assert_eq!(get_bv_rank(&bit_vec![false, true, false]), vec![0,1,1])
}

#[test]
fn test_simple_rank_3() {
    assert_eq!(get_bv_rank(&bit_vec![true, false, true, false]), vec![1,1,2,2])
}


