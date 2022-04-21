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

fn read_vars(config_file: &str) -> (i64, i64, i64, f64) {
    let mut config = Ini::new();
    config.load(config_file).expect("File not found.");

    (
        config.getint("vars", "budget").unwrap().unwrap(),
        config.getint("vars", "streams").unwrap().unwrap(),
        config.getint("vars", "num_of_ind").unwrap().unwrap(),
        config.getfloat("vars", "pclo").unwrap().unwrap(),
    )
}

fn read_csv(path_to_file: &str) -> CsVecBase<Vec<usize>, Vec<Vec<f32>>, Vec<f32>> {
    let mut rdr = ReaderBuilder::new().has_headers(true).delimiter(b',').from_path(path_to_file).unwrap();
    let data = CsVec::empty(6);
    for (i, result) in rdr.records().enumerate() {
        let val: Vec<f32> = result.into_iter().map(|dim| -> f32 {dim.deserialize(None).unwrap()}).collect();
        data.append(i, val)
    }

    return data
}

fn main() {
    let config_file = "config.ini";
    let _data = read_csv("basket_1h.csv");
    let _dim: usize = 6;
    let _generation_count = 0;
    let (_budget, _streams, _num_of_ind, _pclo) = read_vars(config_file);

    let initial = Generation { num_of_ind: _num_of_ind, deterministic: HashMap::new(), random: Vec::new(), streamChromosomes: Vec::new(), chromosomes: HashMap::new(), generationCount: _generation_count, data: _data, k: Vec::new(), streams: _streams, dim: _dim };
    initial.generate_chromosomes();

    let _kmax = 14;

    let clustering = Clustering {generation: initial, data: _data, dim: _dim, kmax: _kmax};

    let generation = clustering.calcChromosomesFit();
    generation.sortChromosomes();

    let ga = Genetic {
        num_of_ind: _num_of_ind,
        pclo: _pclo,
        budget: _budget,
        data: _data,
        generationCount: _generation_count,
        kmax: _kmax,
        dim: _dim,
        prevGeneration: HashMap::new(),
    };

    let (generation, _generation_count) = ga.geneticProcess(initial, initial.deterministic);
    while _generation_count <= _budget {
        let (generation, _generation_count) = ga.geneticProcess(generation, generation.deterministic);
    }

    let iBest = generation.getBestChromosome();
    clustering.printIBest(iBest);
}
