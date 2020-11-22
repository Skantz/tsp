use std::io::{self, BufRead};

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

fn main() {
    let n = 10;

    let mut dm = vec![vec![n]; n];
    fill_distance_matrix(dm);
    
    println!("Hello, world!");
}
