use std::io;
use std::cmp;

//two-opt
//two-opt-move
//three-opt-move
//greedy-tour
//node-shift 
//k-opt



fn fill_distance_matrix(n: i32, dm: &mut Vec<Vec<i32>>) {
    let mut points : Vec<f64> = vec![];
    let stdin = io::stdin();

    //Read in std
    for line in stdin.lock().lines().map(|l| l.unwrap()) {
        let nums: Vec<f64> = line.split_whitespace()
            .map(|num| num.parse().unwrap())
            .collect();
        points.push(nums[0]);
        points.push(nums[1]);
    }


    
}

fn greedy_tour(n: i32, dm: &Vec<Vec<i32>>) -> Vec<i32> {
    let mut tour : Vec<i32> = vec![0; n as usize];
    let mut used : Vec<bool> = vec![false; n as usize];
    used[0] = true;

    for i in 1..n {
        let mut best = -1i32;

        for j in 0..n {
            if !used[j] && (best == -1 || dm[tour[i-1]][j] < dm[tour[i - 1]][best]) {
                best = j;
            }
        }

        tour[i] = best;
        used[best] = true;
    }

    return tour;
}

fn n_mod_m <T: std::ops::Rem<Output = T> + std::ops::Add<Output = T> + Copy>
(n: T, m: T) -> T {
    ((n % m) + m) % m
}

fn two_opt(tour: &Vec<i32>, n: i32, dm: &Vec<Vec<i32>>) {
    let mut improved = true;
    let mut capital_delta = 0;

    while improved {
        improved = false;

        for i in 0..n {
            for j in i+1..n-1 {
                let i_m = n_mod_m((i - 1), n);
                let j_p = (j + 1) % n;

                let new_cost = dm[tour[i]][tour[j_p]] + dm[tour[i_m]][tour[j]];
                let old_cost = dm[tour[i_m]][tour[i]] + dm[tour[j]][tour[j_p]];

                let delta = new_cost - old_cost;

                if delta < 0 {
                    &tour[i..j+1].reverse();
                    improved = true;

                    capital_delta += delta;
                }
            }
        }
    }
}

fn two_opt_move(tour: &Vec<i32>, n: i32, dm: &Vec<Vec<i32>>, i: i32, j: i32) {
    let i_m = n_mod_m((i - 1), n);
    let j_p = (j + 1) % n;

    let new_cost = dm[tour[i]][tour[j_p]] + dm[tour[i_m]][tour[j]];
    let old_cost = dm[tour[i_m]][tour[i]] + dm[tour[j]][tour[j_p]];

    &tour[i..j+1].reverse();

    return new_cost - old_cost;
}

fn three_opt_move(tour: &Vec<i32>, n: i32, dm: &Vec<Vec<i32>>, i: i32, j: i32, k: i32) {
    let a = tour[i - 1];
    let b = tour[i];
    let c = tour[j - 1];
    let d = tour[j];
    let e = tour[k - 1];
    let f = tour[k];

    let d0 = dm[a][b] + dm[c][d] + dm[e][f];
    let d1 = dm[a][c] + dm[b][d] + dm[e][f];
    let d2 = dm[a][b] + dm[c][e] + dm[d][f];
    let d3 = dm[a][d] + dm[e][b] + dm[c][f];
    let d4 = dm[f][b] + dm[c][d] + dm[e][a];

    let case_1 = -d0 + d1;
    let case_2 = -d0 + d2;
    let case_3 = -d0 + d3;
    let case_4 = -d0 + d4;
    let case_min = cmp::min(case_2, cmp::min(case_3, case_4));

    if case_1 <= case_min {

        &tour[i..j].reverse();
        return case_1;

    } else if case_2 == case_min {

        &tour[j..k].reverse();
        return case_2;

    } else if case_4 == case_min {

        &tour[i..k].reverse();
        return case_3;

    } else {

        // let mut tmp = vec![0i32; (k - j) as usize];


    }
}



fn main() {
    let n = 10;

    let mut dm = vec![vec![n]; n];
    // fill_distance_matrix(dm);
    
    println!("Hello, world!");
}
