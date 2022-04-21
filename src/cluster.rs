use ndarray::Array;
use rand::Rng;
use sprs::CsVecBase;

use crate::generation::{Generation, Chromosome};

#[derive(Debug)]
struct Point {
    length: u8,
    pattern_id: Vec<f32>,
    z: usize,
}

impl Point {
    fn __str__(&mut self) {
        let s = format!("{:?}", self.pattern_id);
    }

    fn toJSON(&mut self) {
        let json = r#"
            {
                "pattern_id": self.pattern_id,
            }
        "#;
    }
}

struct Cluster {
    dim: u8,
    centroid: Point,
    points: Vec<Point>,
    distances: Vec<f32>,
}

impl Cluster {
    fn dispersion(&mut self) -> f32 {
        let t = self.points.len();

        let s = 0.0;

        for x in self.distances {
            s += x;
        }
        s = s / t as f32;

        return s;
    }
}

pub struct Clustering {
    pub generation: Generation,
    pub data: CsVecBase<Vec<usize>, Vec<Vec<f32>>, Vec<f32>>,
    pub dim: usize,
    pub kmax: usize,
}

impl Clustering {
    fn daviesBouldin(&mut self, clusters: Vec<Cluster>) -> f64 {
        let sigmaR = 0.0;
        let k = clusters.len();

        for i in 0..k {
            sigmaR += self.computeR(clusters);
        }

        let dbindex = sigmaR / k as f64;

        return dbindex;
    }

    fn computeR(mut self, clusters: Vec<Cluster>) -> f64 {
        let listR = Vec::new();

        for (i, icluster) in clusters.iter().enumerate() {
            for (j, jcluster) in clusters.iter().enumerate() {
                if i != j {
                    let temp = self.computeRij(icluster, jcluster);
                    listR.push(temp);
                }
            }
        }

        let max = listR.iter().fold(f64::INFINITY, |a, &b| a.max(b.into()));
        return max;
    }

    fn computeRij(&mut self, icluster: &Cluster, jcluster: &Cluster) -> f64 {
        let mij = self.minkowskiDistance(icluster.centroid, jcluster.centroid) as f64;
        if mij == 0.0 {
            mij = 1e-10;
        }

        let rij= (icluster.dispersion() + jcluster.dispersion()) as f64;
        rij = rij / mij;

        return rij;
    }

    pub fn calcChromosomesFit(&mut self) -> Generation {
        let kmax = self.kmax;
        let generation = self.generation;
        let num_of_ind = self.generation.num_of_ind as usize;

        for stream in generation.chromosomes.values() {
            for i in 0..num_of_ind {
                let dim = self.generation.dim;
                let clusters = Vec::new();

                let start_idx: usize;
                for _ in 0..kmax {
                    if stream[i].length > dim as u16 {
                        start_idx = rand::thread_rng().gen_range(0..(stream[i].length as usize - dim));
                    } else {
                        start_idx = 0;
                    }

                    let genes = Vec::new();
                    for j in start_idx..(start_idx + dim as usize) {
                        genes.push(stream[i].genes[j]);
                    }

                    let centroid = Point {
                        length: genes.len() as u8,
                        pattern_id: genes,
                        z: usize::MAX,
                    };
                    let cluster = Cluster {
                        dim: dim as u8,
                        centroid: centroid,
                        points: Vec::new(),
                        distances: Vec::new(),
                    };
                    clusters.push(cluster);
                }

                clusters = self.calcDistance(clusters);
                let dbindex = self.daviesBouldin(clusters);

                stream[i].fitness = 1.0 / dbindex;
            }
        }

        return self.generation;
    }

    pub fn calcChildFit(&mut self, childChromosome: Chromosome) -> Chromosome {
        let kmax = self.kmax;
        let dim = self.dim;
        let clusters = Vec::new();

        let start_idx: usize;
        for _ in 0..kmax {
            if childChromosome.length > dim as u16 {
                start_idx = rand::thread_rng().gen_range(0..(childChromosome.length as usize - dim));
            } else {
                start_idx = 0;
            }

            let genes = Vec::new();
            for j in start_idx..(start_idx + dim) {
                genes.push(childChromosome.genes[j]);
            }

            let centroid = Point {
                length: genes.len() as u8,
                pattern_id: genes,
                z: usize::MAX,
            };
            let cluster = Cluster {
                dim: dim as u8,
                centroid: centroid,
                points: Vec::new(),
                distances: Vec::new(),
            };
            clusters.push(cluster);
        }

        clusters = self.calcDistance(clusters);
        let dbindex = self.daviesBouldin(clusters);

        childChromosome.fitness = 1.0 / dbindex;

        return childChromosome;
    }

    fn calcDistance(&mut self, clusters: Vec<Cluster>) -> Vec<Cluster>{
        let kmax = self.kmax;
        let dim = self.dim;
        let data = self.data;
        let dis = 0;

        for z in data.indices() {
            let point = Point {
                length: dim as u8,
                pattern_id: *data.get(*z).unwrap(),
                z: *z,
            };

            let disSet = Vec::new();

            for i in 0..kmax {
                let dis = self.minkowskiDistance(point, clusters[i].centroid);
                disSet.push(dis);
            }

            clusters = self.findMin(disSet, clusters, point);
        }

        return clusters;
    }

    fn minkowskiDistance(&mut self, point1: Point, point2: Point) -> f32 {
        let dim = self.dim;
        let p = 2.0;
        let d = 0.0;

        for i in 0..dim {
            d += (point1.pattern_id[i] - point2.pattern_id[i]).abs().powf(p);
        }

        let mij = d.powf(1.0/p);

        return mij;
    }

    fn findMin(&mut self, dis_set: Vec<f32>, clusters: Vec<Cluster>, point: Point) -> Vec<Cluster> {
        let idx = dis_set.iter().position(|r| r == dis_set.iter().fold(&f32::INFINITY, |a, b| &a.min(*b))).unwrap();
        clusters[idx].points.push(point);
        clusters[idx].distances.push(dis_set.iter().fold(f32::INFINITY, |a, b| a.min(*b)));

        return clusters;
    }

    pub fn printIBest(&mut self, iBest: Chromosome) {
        let kmax = self.kmax;
        let dim = self.dim;
        let data = self.data;
        let clusters = Vec::new();

        let start_idx: usize;
        for _ in 0..kmax {
            if iBest.length > dim as u16 {
                start_idx = rand::thread_rng().gen_range(0..(iBest.length as usize - dim));
            } else {
                start_idx = 0;
            }

            let genes = Vec::new();
            for j in start_idx..(start_idx + dim) {
                genes.push(iBest.genes[j]);
            }

            let centroid = Point {
                length: genes.len() as u8,
                pattern_id: genes,
                z: usize::MAX,
            };
            let cluster = Cluster {
                dim: dim as u8,
                centroid: centroid,
                points: Vec::new(),
                distances: Vec::new(),
            };
            clusters.push(cluster);
        }

        clusters = self.calcDistance(clusters);
        let dbindex = self.daviesBouldin(clusters);
        let mut z = Array::zeros(640992);

        for (i, cluster) in clusters.iter().enumerate() {
            for j in cluster.points {
                z[j.z] = i;
            }
        }

        let correct_answer = 0;

        for i in 0..(640992 * 0.33 as usize) {
            if z[i] == 2 {
                correct_answer += 1;
            }
        }
        for i in (640992 * 0.33 as usize)..(640992 * 0.67 as usize) {
            if z[i] == 1 {
                correct_answer += 1;
            }
        }
        for i in (640992 * 0.67 as usize)..640992 {
            if z[i] == 0 {
                correct_answer += 1;
            }
        }

        let accuracy = (correct_answer / 640992) * 100;
        println!("accuracy: {}", accuracy);
        println!("iBest Fitness: {}", (1.0 / dbindex));
        println!("Clusters centroid:");

        for (i, cluster) in clusters.iter().enumerate() {
            println!("centroid {}: {:?}", i, cluster.centroid);
        }
    }
}
