use configparser::ini::Ini;
use sprs::*;
use csv::ReaderBuilder;
use std::collections::HashMap;

mod cluster;
use cluster::*;
mod generation;
use generation::*;
mod genetic;
use genetic::*;

fn read_vars(config_file: &str) -> (usize, usize, usize, f64) {
    let mut config = Ini::new();
    config.load(config_file).expect("File not found.");

    (
        config.getint("vars", "budget").unwrap().unwrap().try_into().unwrap(),
        config.getint("vars", "streams").unwrap().unwrap().try_into().unwrap(),
        config.getint("vars", "num_of_ind").unwrap().unwrap().try_into().unwrap(),
        config.getfloat("vars", "pclo").unwrap().unwrap(),
    )
}

fn read_csv(path_to_file: &str) -> CsVecBase<Vec<usize>, Vec<Vec<f32>>, Vec<f32>> {
    let mut rdr = ReaderBuilder::new().has_headers(true).delimiter(b',').from_path(path_to_file).unwrap();
    let mut data = CsVec::empty(6);
    for (i, result) in rdr.records().enumerate() {
        let row: Vec<f32> = result.into_iter().map(|row| -> f32 {row.deserialize(None).unwrap()}).collect();
        data.append(i, row)
    }

    return data;
}

fn main() {
    let config_file = "config.ini";
    let _data = read_csv("basket_1h.csv");
    let _dim: usize = 6;
    let _generation_count: usize = 0;
    let (_budget, _streams, _num_of_ind, _pclo) = read_vars(config_file);

    let mut _initial = Generation { num_of_ind: _num_of_ind, deterministic: HashMap::new(), random: Vec::new(), streamChromosomes: Vec::new(), chromosomes: HashMap::new(), generationCount: _generation_count, data: _data.clone(), k: Vec::new(), streams: _streams, dim: _dim };
    _initial.generate_chromosomes();

    let _kmax = 14;

    let mut clustering = Clustering {generation: _initial.clone(), data: _data.clone(), dim: _dim, kmax: _kmax};

    let mut _generation = clustering.calcChromosomesFit();
    _generation.sortChromosomes();

    let mut ga = Genetic {
        num_of_ind: _num_of_ind,
        pclo: _pclo,
        budget: _budget,
        data: _data.clone(),
        generationCount: _generation_count,
        kmax: _kmax,
        dim: _dim,
        prevGeneration: HashMap::new(),
    };

    let (mut _generation, mut _generation_count) = ga.geneticProcess(&_initial, &_initial.deterministic);
    while _generation_count <= _budget {
        let (_generation, _generation_count) = ga.geneticProcess(&_generation, &_generation.deterministic);
    }

    let iBest = _generation.getBestChromosome();
    clustering.printIBest(&iBest);
}
