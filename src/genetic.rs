use std::collections::{HashMap, HashSet};
use ordered_float::OrderedFloat;
use either::*;
use ndarray::{ArrayBase, OwnedRepr, Dim};
use rand::prelude::*;
use crate::{generation::Generation, generation::Chromosome, cluster::Clustering};

pub struct Genetic {
    pub num_of_ind: i64,
    pub pclo: f64,
    pub budget: i64,
    pub data: ArrayBase<OwnedRepr<u32>, Dim<[usize; 2]>>,
    pub generationCount: i64,
    pub kmax: usize,
    pub dim: usize,
    pub prevGeneration: HashMap<i64, Vec<Chromosome>>,
}

impl Genetic {
    pub fn geneticProcess(&mut self, generation: Generation, deterministic: HashMap<i64, Vec<Chromosome>>) -> (Generation, i64) {
        println!("------------Generation: {} -----------------", self.generationCount);

        //  ----------------------------------Information Sharing-----------------------

        if self.generationCount % 10 == 0 {
            generation = self.informationSharing(generation);
        }

        let s = 0;

        for stream in generation.chromosomes.values() {
            //  ------------------------------Simple noise based selection--------------

            if self.generationCount > 0 {
                stream = self.selection(stream, self.prevGeneration[&s]);
            }
            else {
                continue;
            }

            //  ------------------------------Crossover---------------------------------

            stream = self.crossover(stream, generation);

            //  ------------------------------Twin Removal------------------------------

            stream = self.twinRemoval(stream);

            //  ------------------------------Mutation----------------------------------

            stream = self.mutation(stream);

            //  ------------------------------Health Improvement------------------------

            stream = self.healthImprovement(stream, generation, deterministic[&s]);
            
            //  ------------------------------Cleansing---------------------------------

            stream = self.cleansing(stream, deterministic[&s]);

            //  ------------------------------Elitist-----------------------------------

            if self.generationCount > 0 {
                stream = self.elitist(stream, self.prevGeneration[&s]);
            }
            else {
                continue;
            }

            s += 1;
        }

        generation.sortChromosomes();
        self.prevGeneration = generation.chromosomes;
        self.generationCount += 1;

        return (generation, self.generationCount);
        
    }

    fn informationSharing(&mut self, generation: Generation) -> Generation {
        let neighbors = Vec::new();

        for i in 0..generation.chromosomes.len() as i64 {
            let chromo = generation.chromosomes[&i];
            let bestChromo = chromo[i as usize];
            neighbors.push(bestChromo);
        }

        for i in 0..generation.chromosomes.len() as i64 {
            let chromo = generation.chromosomes[&i].last().unwrap();
            chromo = &self.maxOfNeighbors(neighbors, i);
        }

        return generation;
    }

    fn maxOfNeighbors(&mut self, neighbors: Vec<Chromosome>, i: i64) -> Chromosome {
        let d = neighbors.len() - i as usize;
        let p1: usize;
        if d >= 2 {
            p1 = i as usize;
        } else {
            p1 = neighbors.len() - 3;
        }

        let p2 = p1 + 2;
        let max_chromosome;
        for i in p1..p2 {
            let neighbors_set = Vec::new();
            neighbors_set.push(neighbors[i]);
            neighbors_set.sort_by(|a, b| b.fitness.partial_cmp(&a.fitness).unwrap());
            let max_chromosome = neighbors_set[0];
        }

        return max_chromosome;
    }

    fn selection<'a>(&mut self, stream: &'a Vec<Chromosome>, prevStream: Vec<Chromosome>) -> &'a Vec<Chromosome> {
        let num_of_ind = self.num_of_ind as usize;
        let noise = rand::thread_rng().gen_range(0.0..1.0);

        for j in 0..num_of_ind {
            if prevStream[j].fitness > stream[j].fitness {
                if prevStream[j].fitness > (stream[j].fitness + noise) {
                    stream[j] = prevStream[j]
                }
            }
        }

        return stream;
    }

    fn crossover<'a>(&mut self, stream: &'a Vec<Chromosome>, generation: Generation) -> &'a Vec<Chromosome> {
        let dad = stream[0];
        let Fsum = 0.0;

        for i in 0..stream.len() {
            Fsum += stream[i].fitness;
        }

        for j in 1..stream.len() {
            let tj = stream[j].fitness / Fsum;
            let die = rand::thread_rng().gen_range(0.0..1.0);

            if tj >= die {
                let mom = stream[j];

                let cut: usize;
                if dad.length < mom.length {
                    cut = rand::thread_rng().gen_range(1..(dad.length - 1)) as usize;
                } else {
                    cut = rand::thread_rng().gen_range(1..(mom.length - 1)) as usize;
                }

                let genesChild1 = Vec::new();
                for i in 0..cut {
                    genesChild1.push(dad.genes[i]);
                }
                for i in cut..((mom.length - 1) as usize) {
                    genesChild1.push(mom.genes[i]);
                }

                while genesChild1.len() < self.dim as usize {
                    genesChild1.push(OrderedFloat(rand::thread_rng().gen_range(0.0..1.0)));
                }

                let genesChild2 = Vec::new();
                for i in 0..cut {
                    genesChild2.push(mom.genes[i]);
                }
                for i in cut..((dad.length - 1) as usize) {
                    genesChild2.push(dad.genes[i]);
                }

                while genesChild2.len() < self.dim as usize {
                    genesChild2.push(OrderedFloat(rand::thread_rng().gen_range(0.0..1.0)));
                }

                let child1 = Chromosome {
                    genes: genesChild1,
                    length: genesChild1.len() as u16,
                    fitness: 0.0,
                    mj: 0.0,
                };

                let child2 = Chromosome {
                    genes: genesChild2,
                    length: genesChild2.len() as u16,
                    fitness: 0.0,
                    mj: 0.0,
                };

                let clustering = Clustering{ generation: generation, data: self.data, dim: self.dim, kmax: self.kmax };
                child1 = clustering.calcChildFit(child1);
                child2 = clustering.calcChildFit(child2);

                let family = Vec::new();
                family.push(dad);
                family.push(mom);
                family.push(child1);
                family.push(child2);

                family.sort_by(|a, b| b.fitness.partial_cmp(&a.fitness).unwrap());
                stream[0] = family[0];
                stream[j] = family[1];
           }
        }
        return stream;
    }

    fn twinRemoval<'a>(&mut self, stream: &'a Vec<Chromosome>) -> &'a Vec<Chromosome> {
        let gene;
        for i in 0..stream.len() {
            gene = self.checkForIdenticalGenes(stream[i]);

            if stream[i].length > self.dim as u16 || gene != either::Right(false) {
                let idx = stream[i].genes.iter().position(|r| r == &gene.unwrap_left()).unwrap();
                stream[i].genes.remove(idx);
            } else if stream[i].length == self.dim as u16|| gene != either::Right(false) {
                let idx = stream[i].genes.iter().position(|r| r == &gene.unwrap_left()).unwrap();
                stream[i].genes[idx] = OrderedFloat(rand::thread_rng().gen_range(0.0..1.0));
            }
        }

        return stream;

    }

    fn checkForIdenticalGenes(&mut self, chromosome: Chromosome) -> Either<OrderedFloat<f32>, bool> {
        let mut setOfGenes: HashSet<OrderedFloat<f32>> = HashSet::new();
        for gene in chromosome.genes {
            if setOfGenes.contains(&gene) {
                return Left(gene);
            } else {
                setOfGenes.insert(gene);
            }
        }
        return Right(false);
    }

    fn mutation<'a>(&mut self, stream: &'a Vec<Chromosome>) -> &'a Vec<Chromosome> {
        let k1 = 0.5;
        let k2 = 0.5;

        let fmax;
        for i in 0..stream.len() {
            if stream[i].fitness > stream[i - 1].fitness {
                fmax = stream[i].fitness;
            }
        }

        let totalFitness = 0.0;
        for k in 0..stream.len() {
            totalFitness += stream[k].fitness;
        }

        let fave = totalFitness / stream.len() as f64;

        let fj;
        for j in 0..stream.len() {
            fj = stream[j].fitness;

            if fj > fave {
                stream[j].mj = k1 * ((fmax - fj) / (fmax - fave)) as f32;
            } else {
                stream[j].mj = k2;
            }
        }

        for i in 0..stream.len() {
            let die = rand::thread_rng().gen_range(0.0..1.0);

            if stream[i].mj >= die {

            }
        }

        return stream;
    }

    fn healthImprovement<'a>(&mut self, stream: &'a Vec<Chromosome>, generation: Generation, deterministic: Vec<Chromosome>) -> &'a Vec<Chromosome> {
        let num_of_ind = self.num_of_ind;
        let twentyPercent = &Vec::new();

        for i in 0..(num_of_ind as f32 * 0.2) as usize {
            twentyPercent.push(stream[i]);
        }

        twentyPercent = self.crossover(twentyPercent, generation);
        let thirtyPercent = Vec::new();

        for i in 0..(num_of_ind as f64 * 0.3) as usize {
            thirtyPercent.push(deterministic[i]);
        }

        let idx;
        for i in 0..thirtyPercent.len() {
            idx = thirtyPercent[i].genes.iter().position(|r| r == thirtyPercent[i].genes.choose(&mut rand::thread_rng()).unwrap()).unwrap();
            thirtyPercent[i].genes[idx] = OrderedFloat(rand::thread_rng().gen_range(0.0..1.0));
        }

        let i: usize = (num_of_ind - 1).try_into().unwrap();

        while i >= (num_of_ind as f64 * 0.5) as usize {
            stream.remove(i);
            i -= 1;
        }

        for i in 0..twentyPercent.len() {
            stream.push(twentyPercent[i]);
        }

        for i in 0..thirtyPercent.len() {
            stream.push(thirtyPercent[i]);
        }
        return stream;
    }

    fn cleansing<'a>(&mut self, stream: &'a Vec<Chromosome>, deterministic: Vec<Chromosome>) -> &'a Vec<Chromosome> {
        let num_of_ind = self.num_of_ind as usize;

        let mr: f32;
        for i in 0..num_of_ind {
            if stream[i].length < stream[i - 1].length {
                mr = stream[i].length.into();
            }
        }

        for i in 0..num_of_ind {
            let mn = stream[i].genes.iter().fold(f32::INFINITY, |a, &b| a.min(*b));
            let mx = stream[i].genes.iter().fold(f32::INFINITY, |a, &b| a.min(*b));
            let t = 0.1;

            mx = mx + (mx * t);
            mn = mn - (mn * t);
            mr = mr - (mr * t);

            stream[i] = *self.doCleansing(stream[i], mn, mx, mr, deterministic);
        }

        return stream;
    }

    fn doCleansing(
        &mut self,
        chromosome: Chromosome,
        mn: f32,
        mx: f32,
        mr: f32,
        deterministic: Vec<Chromosome>,
    ) -> &Chromosome {
        if !(chromosome.length >= mn as u16) || !(chromosome.length <= mx as u16) || !(chromosome.length > mr as u16) {
            let new_chromosome = self.cloning(&chromosome, deterministic);
            return new_chromosome;
        } else {
            return &chromosome;
        }
    }

    fn cloning<'a>(&mut self, sickChromosome: &'a Chromosome, deterministic: Vec<Chromosome>) -> &'a Chromosome {
        let pclo = self.pclo;
        let die = rand::thread_rng().gen_range(0.0..1.0);

        if pclo >= die {
            sickChromosome = deterministic.choose(&mut rand::thread_rng()).unwrap();
            let idx = sickChromosome.genes.iter().position(|r| r == sickChromosome.genes.choose(&mut rand::thread_rng()).unwrap()).unwrap();
            sickChromosome.genes[idx] = OrderedFloat(rand::thread_rng().gen_range(0.0..1.0));
            let chromosome = sickChromosome;
            return chromosome;
        } else {
            return sickChromosome;
        }
    }

    fn elitist<'a>(&mut self, stream: &'a Vec<Chromosome>, prevStream: Vec<Chromosome>) -> &'a Vec<Chromosome> {
        let lastChromo = stream.last().unwrap();
        let worstFitnessCurrent = lastChromo.fitness;
        let bestFitnessCurrent = stream[0].fitness;
        let bestFitnessAll = prevStream[0].fitness;

        if worstFitnessCurrent < bestFitnessAll {
            worstFitnessCurrent = bestFitnessAll;
        }

        if bestFitnessCurrent > bestFitnessAll {
            bestFitnessCurrent = bestFitnessAll;
        }

        return stream;
    }
}
