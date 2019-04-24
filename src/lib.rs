use na::{DMatrix, DVector, Schur};
use num::Complex;

pub fn roots(p: &[f64]) -> Vec<Complex<f64>> {
    // short-circuit if p is empty or all zero
    if p.len() == 0 || p.iter().all(|&x| x == 0.0) {
        return vec![];
    }

    // split off trailing zeros
    let trailing_zeros = p.iter().rev().take_while(|x| **x == 0.0).count();
    let trailing_zeros_idx = if trailing_zeros == 0 {
        None
    } else {
        Some(p.len() - trailing_zeros)
    };

    let (non_zero, trailing_zeros) = match trailing_zeros_idx {
        Some(idx) => p.split_at(idx),
        None => (p, &[][..]),
    };

    // now strip leading zeros
    let non_zero_start = p.iter()
        .take_while(|x| **x == 0.0)
        .count();
    let non_zero = &non_zero[non_zero_start..];

    let matrix_size = p.len() - 1;

    let roots = if matrix_size > 1 {
        // Now create companion matrix
        let mut companion_matrix = DMatrix::from_fn(matrix_size, matrix_size,
            |i,j| if i == j+1 { 1.0 } else { 0.0 }
        );
        let divisor = non_zero[0];
        let p_col = non_zero.into_iter()
            .skip(1)
            .rev()
            .cloned()
            .map(|x| -1.0 * (x / divisor));

        companion_matrix.set_column(
            matrix_size-1,
            &DVector::from_iterator(matrix_size, p_col),
        );

        Some(Schur::new(companion_matrix).complex_eigenvalues())
    } else {
        None
    };

    // now append zero roots
    let trailing_zeros = trailing_zeros.into_iter()
        .cloned()
        .map(|x| Complex::new(x, 0.0));

    match roots {
        Some(rs) => {
            let mut rs = rs.as_slice().to_vec();
            rs.extend(trailing_zeros);
            rs
        },
        None => trailing_zeros.collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        let rs = roots(&[1.0, -10.0, 31.0, -30.0]);
        println!("{:?}", rs);

        // doesn't do complex
        let rs = roots(&[3.2, 2.0, 1.0]);
        println!("{:?}", rs);
        //panic!()
    }
}
