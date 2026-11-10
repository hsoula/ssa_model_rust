#!/bin/sh
cargo run --release --bin simulate -- -r data/my_chem.json -a 0x0403:20 -a 0x0203:20 -a 0x0101:20 -a 0x0102:20 -a 0x0401:20 -a 0x0201:20 -a 0x0300:20 --mol data/env/initial_bonds.mol --tmax 10000.0 --seed 42 -o data/env/run
