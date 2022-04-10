use std::collections::HashMap;
use crate::Error;
use ndarray::{ArrayBase, OwnedRepr, Dim};
use rand::prelude::*;

pub struct Chromosome {
    pub genes: Vec<f32>,
    pub length: u16,
    pub fitness: f64,
    pub mj: f32,
}

impl Chromosome {
    pub fn randomGenerateChromosome(&mut self) {
        for i in 0..self.length {
            let gene: f32 = rand::thread_rng().gen_range(0.0..1.0);
            self.genes.push(gene);
        }
    }
}


pub struct Generation<> {
    pub num_of_ind: i64,
    pub deterministic: HashMap<i64, Vec<Chromosome>>,
    pub random: Vec<Chromosome>,
    pub streamChromosomes: Vec<Chromosome>,
    pub chromosomes: HashMap<i64, Vec<Chromosome>>,
    pub generationCount: i64,
    pub data: Result<ArrayBase<OwnedRepr<u32>, Dim<[usize; 2]>>, Box<dyn Error>>,
    pub k: Vec<u16>,
    pub streams: i64,
    pub dim: usize,
}

impl Generation {
    pub fn sortChromosomes(&mut self) -> HashMap<i64, Vec<Chromosome>> {
        for streamOfChromosome in self.chromosomes.values() {
            streamOfChromosome.sort_by(|a, b| b.fitness.partial_cmp(&a.fitness).unwrap());
        }

        return self.chromosomes;
    }

    pub fn generate_chromosomes(&mut self) {
        for s in 0..self.streams {
            let mut deterministic = Vec::new();
            let mut random = Vec::new();
            let mut stream_chromosomes = Vec::new();

            let chromosomes_from_individual_phase = (self.num_of_ind / 2) as usize;

            for k in self.dim..chromosomes_from_individual_phase {
                for _ in 0..5 {
                    let mut chromosome = Chromosome {
                        genes: Vec::new(),
                        length: k as u16,
                        fitness: 0.0,
                        mj: 0.0,
                    };
                    chromosome.randomGenerateChromosome();
                    deterministic.push(chromosome);
                }
            }

            for i in 0..deterministic.len() {
                let n: f32 = 640992.0;
                let k = rand::thread_rng().gen_range(self.dim..(n.sqrt() as usize));
                self.k.push((k as u8).into());
            }

            for i in 0..deterministic.len() {
                let k = self.k.choose(&mut rand::thread_rng());
                let mut chromosome = Chromosome {
                    genes: Vec::new(),
                    length: *k.unwrap(),
                    fitness: 0.0,
                    mj: 0.0,
                };
                chromosome.randomGenerateChromosome();
                random.push(chromosome);
            }

            self.deterministic.insert(s, deterministic);

            for i in 0..chromosomes_from_individual_phase {
                stream_chromosomes.push(deterministic[i]);
                stream_chromosomes.push(random[i]);
            }

            self.chromosomes.insert(s, stream_chromosomes);
        }
    }

    pub fn getBestChromosome(&mut self) -> Chromosome {
        let bestChromosomeInStream = Vec::new();

        for stream in self.chromosomes.values() {
            bestChromosomeInStream.push(stream[0]);
        }

        let bestChromosome;
        for i in 0..bestChromosomeInStream.len() {
            if bestChromosomeInStream[i].fitness > bestChromosomeInStream[i - 1].fitness {
                bestChromosome = bestChromosomeInStream[i];
            }
        }

        return bestChromosome;
    }
}
