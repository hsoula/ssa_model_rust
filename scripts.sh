#Generate a random chemistry:
cargo run --release --bin gen_chemistry -- -n 20 -f 4 -s 4 --seed 42 -o data/my_chem.json

#Run a simulation on it:
cargo run --release --bin simulate -- \
  -r data/my_chem.json \
  -a 0x0100:300 -a 0x0200:300 -a 0x0300:300 \
  --tmax 10000 --seed 42 -o data/output/run1

#Stack multiple reaction files exactly like the old code:
cargo run --release --bin simulate -- \
  -r data/replicator.json -r data/prod_degra.json -r data/detergent.json \
  -a 0x0100:300 -a 0x0200:300 -a 0x0300:300 \
  --tmax 10000 --seed 42 -o data/output/run1

#Evolve a chemistry with a GA:
cargo run --release --bin evolve -- \
  -a 0x0100:200 -a 0x0200:200 \
  --population 30 --generations 100 --tmax 500 \
  --seed 42 -o data/evolve/run1


# Single target
cargo run --release --bin evolve_target -- \
  -m molecules/json/C01101.json \
  --n-atoms 100 --genome-size 300 --population 30 \
  --generations 200 --tmax 2000 --seed 42 \
  -o data/evolve/L-Ribulose

# Multiple targets at once
cargo run --release --bin evolve_target -- \
  -m molecules/json/C01101.json \
  -m molecules/json/C01802.json \
  --population 50 --generations 500 --tmax 5000 \
  -o data/evolve/multi_target