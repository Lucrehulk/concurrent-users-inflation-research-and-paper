// See paper for code derivation.

use rand::seq::SliceRandom;
use rand::Rng;
use std::collections::HashSet;
use std::f32;
use std::fs::File;
use std::io;
use std::io::{BufWriter, Write};

// Set constants.
const K: f32 = 0.2;
const MU: f32 = 0.0;
const SIG: f32 = 4.5;
const P_0: i32 = 470;
const P_MIN: i32 = 350;
const P_MAX: i32 = 630;
const MIN_T: i32 = 1;
const MAX_T: i32 = 500;
const T_JUMP: usize = 1;
const CONTROL_FRACT: f32 = 1.0;

// Generate Norm(x) for a given x.
fn norm(x: f32) -> f32 {
    let amplitude = 1.0 / (SIG * (2.0 * f32::consts::PI).sqrt());
    let exponent = -0.5 * ((x - MU) / SIG).powi(2);
    amplitude * exponent.exp()
}

// Monte Carlo hit & miss implementation.
fn monte_carlo(p: i32) -> i32 {
    let mut rng = rand::thread_rng();

    let x_bounds: [f32; 2] = [MU - 3.0 * SIG, MU + 3.0 * SIG];
    let y_bounds: [f32; 2] = [0.0, 1.0 / (SIG * (2.0 * f32::consts::PI).sqrt())];

    // Loop selections until a selected y <= norm(x).
    loop {
        let x_continuous = rng.gen_range(x_bounds[0]..=x_bounds[1]);
        let x_rounded = x_continuous.round();

        // Ensure the guess will not fluctuate players below or above bounds.
        if p + (x_rounded as i32) < P_MIN || p + (x_rounded as i32) > P_MAX {
            continue;
        }

        let y = rng.gen_range(y_bounds[0]..=y_bounds[1]);

        if y <= norm(x_rounded) {
            return x_rounded as i32;
        }
    }
}

// Generate a vector of player fluctuations per tick.
fn generate_player_fluctuations(len: usize) -> Vec<i32> {
    let mut result = Vec::new();
    let mut p = P_0;
    // Calculate each guess and track the p value with fluctuations to ensure it will stay within bounds;
    // the Monte Carlo function itself checks bounds to filter out guesses.
    for _ in 0..len {
        let monte_carlo_result = monte_carlo(p);
        result.push(monte_carlo_result);
        p += monte_carlo_result;
    }
    result
}

// Build a single trial.
fn build_trial_vector(len: usize, t_0: i32, fluctuations: Vec<i32>, inflate: bool) -> (Vec<i32>, f32, usize, usize) {
    let mut rng = rand::thread_rng();

    let mut p = P_0;
    let mut result = Vec::new();
    // The latest valid start is len - t_0 - 1, so the inflation reaches its full value before the vector ends.
    let inflation_period_start = rng.gen_range(0..(len - (t_0 as usize)));
    let inflation_period_end = inflation_period_start + (t_0 as usize);

    // Fluctuate p, calculate f, push total CCU as p + f.
    for t in 0..len {
        p += fluctuations[t];

        let f = if inflate && t >= inflation_period_start {
            let time_fract = f32::min(1.0, ((t - inflation_period_start) as f32) / (t_0 as f32));

            f32::floor(K * (p as f32) * time_fract) as i32
        } else {
            0
        };

        result.push(p + f);
    }

    if inflate {
        (
            result,
            K / (t_0 as f32),
            inflation_period_start,
            inflation_period_end,
        )
    } else {
        (result, 0.0, 0, 0)
    }
}

// Generate random u16 value to then map to hexadecimal IDs.
fn generate_ids(count: usize) -> Vec<String> {
    let mut rng = rand::thread_rng();
    let mut seen: HashSet<u16> = HashSet::new();
    let mut ids = Vec::new();

    // Generate a random u16, check if the value has already been inserted.
    // If it is fresh, then insert it.
    while ids.len() < count {
        let value: u16 = rng.gen();
        if seen.insert(value) {
            ids.push(format!("{:04x}", value));
        }
    }

    ids
}

// Generate a full set of individual trials from MIN_T to MAX_T
// (in steps of T_JUMP), plus a set of control trials.
fn generate_trials(len: usize) -> (Vec<(String, f32, usize, usize)>, Vec<(String, Vec<i32>)>) {
    let mut trials: Vec<(Vec<i32>, f32, usize, usize)> = Vec::new();

    // Generate trials from MIN_T to MAX_T, advancing t_0 by T_JUMP each time.
    for t_0 in (MIN_T..=MAX_T).step_by(T_JUMP) {
        let fluctuations = generate_player_fluctuations(len);
        trials.push(build_trial_vector(len, t_0, fluctuations, true));
    }

    // The number of inflated trials built depends on T_JUMP, so count them directly.
    let inflated_count = trials.len() as f32;
    let control_count = (inflated_count * CONTROL_FRACT).round() as i32;

    // Build control trials.
    for _ in 0..control_count {
        let fluctuations = generate_player_fluctuations(len);
        trials.push(build_trial_vector(len, MIN_T, fluctuations, false));
    }

    let ids = generate_ids(trials.len());

    let mut data = Vec::new();
    let mut vectors = Vec::new();

    // Combine (zip) the IDs with the trial data, then iterate over each
    // and push the analysis data to one, and the trial data to the other,
    // with IDs passed to both.
    for ((trial, d_t, start, end), id) in trials.into_iter().zip(ids.into_iter()) {
        data.push((id.clone(), d_t, start, end));
        vectors.push((id, trial));
    }

    // Randomly order the trials to reduce bias.
    let mut rng = rand::thread_rng();
    vectors.shuffle(&mut rng);

    (data, vectors)
}

// Write each trial in form:
// [number] (id). [...data]
// to ./trials.txt
fn write_trials(trials: &Vec<(String, Vec<i32>)>) -> io::Result<()> {
    let file = File::create("./trials.txt")?;
    let mut writer = BufWriter::new(file);

    for (i, (id, trial)) in trials.iter().enumerate() {
        let values: Vec<String> = trial.iter().map(|v| v.to_string()).collect();
        writeln!(writer, "{} ({}). [{}]", i + 1, id, values.join(", "))?;
    }

    writer.flush()?;
    Ok(())
}

// Write each trial's analysis data in form:
// [number] (id). [d_t] [start] [end]
// to ./trial_data.txt
fn write_trial_data(data: &Vec<(String, f32, usize, usize)>) -> io::Result<()> {
    let file = File::create("./trial_data.txt")?;
    let mut writer = BufWriter::new(file);

    for (i, (id, d_t, start, end)) in data.iter().enumerate() {
        writeln!(writer, "{} ({}). {} {} {}", i + 1, id, d_t, start, end)?;
    }

    writer.flush()?;
    Ok(())
}

// Generate trials and write info to files.
fn main() -> io::Result<()> {
    let (data, vectors) = generate_trials((2 * MAX_T) as usize);
    write_trials(&vectors)?;
    write_trial_data(&data)
}