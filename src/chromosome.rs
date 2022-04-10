use rand::Rng;

pub struct Chromosome {
    pub genes: Vec<f32>,
    pub length: u16,
    pub fitness: f32,
    pub mj: f32,
}

impl Chromosome {
    pub fn randomGenerateChromosome(&mut self) -> Self {
        for i in 0..self.length {
            let gene: f32 = rng.gen_range(0.0..1.0);
            self.genes.push(gene);
        }
    }
}
