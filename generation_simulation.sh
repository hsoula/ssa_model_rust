# generate random chemistry
seed=2
cargo run --release  --bin gen_chemistry \
-- -n 10 -f 4 -s 6 --seed $seed \
--output chemistries/chem_gen_$seed.json
cat chemistries/chem_gen_$seed.md
# and run it
seed2=42
cargo run --release --bin simulate -- \
  -r chemistries/chem_gen_$seed.json \
  -a 0x0100:3 -a 0x0101:3 -a 0x0102:3 -a 0x0103:3 -a 0x0104:3 -a 0x1305:3\
  -a 0x0200:3 -a 0x0201:3 -a 0x0202:3 -a 0x0203:3 -a 0x0204:3 -a 0x0205:3 \
  -a 0x0300:3 -a 0x0301:3 -a 0x0302:3 -a 0x0303:3 -a 0x0304:3 -a 0x0305:3\
  -a 0x0400:3 -a 0x0401:3 -a 0x0402:3 -a 0x0403:3 -a 0x0404:3 -a 0x0405:3\
  --tmax 10000 --seed $seed2 -o results/random_chems/random_$seed
