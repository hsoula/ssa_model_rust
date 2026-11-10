# Gene8 Encoding

Gene8 is a linear genome of u8 tokens, each in `[0x00, 0x0F]` (one hex digit).
Reactions are encoded as sub-sequences delimited by **E** (start) and **F** (end).
The genome is created by sampling uniformly from 16 symbols and reducing mod 16.

## Token table

| Token | Role | Field |
|-------|------|-------|
| `E (0x0E)` | START marker | opens a reaction pattern |
| `F (0x0F)` | END marker | closes a reaction pattern |
| `0,1` | sf1++ | substrate 1 form  (x2 probability) |
| `2,3` | ss1++ | substrate 1 state (x2 probability) |
| `4,5` | sf2++ | substrate 2 form  (x2 probability) |
| `6,7` | ss2++ | substrate 2 state (x2 probability) |
| `8` | contact++ | reaction-type bits (value mod 4) |
| `9,A` | ps1++ | product 1 state   (x2 probability) |
| `B,C` | ps2++ | product 2 state   (x2 probability) |
| `D` | speed += 0.1 | reaction speed increment |

## Decoding rules

For each `E…F` pattern the payload tokens are tallied into counters:

```
sf1, ss1, sf2, ss2   ← substrate forms/states
ps1, ps2             ← product states  (token 5 → ps1, token 6 → ps2)
contact              ← reaction-type bits (mod 4)
speed                ← 0.1 × (count of token 9)
```

Atom keys are 16-bit words: `(form << 8) | state`.

```
s1 = (sf1 << 8) | ss1      s2 = (sf2 << 8) | ss2
p1 = (sf1 << 8) | ps1      p2 = (sf2 << 8) | ps2
```

Product forms **always** equal substrate forms — only states can change.
A reaction is **invalid** (discarded) if `sf1 == 0` or `sf2 == 0`.

## Reaction-type bits

| contact mod 4 | Meaning |
|---------------|---------|
| `0b00` (0) | contact (no bond change) |
| `0b01` (1) | product bond |
| `0b10` (2) | substrate bond |
| `0b11` (3) | bond exchange |

## Worked examples

### Minimal bond-forming reaction

Pattern: `E024F`

| Counter | Value |
|---------|-------|
| sf1 | 1 |
| ss1 | 1 |
| sf2 | 1 |
| ss2 | 0 |
| contact (mod 4) | 0 → contact (no bond change) |
| ps1 | 0 |
| ps2 | 0 |
| speed | 0.0 |

**Reaction:** `a1.a0 → a0.a0`  *(speed=0.0)*

### State-change reaction with speed

Pattern: `E01259F`

| Counter | Value |
|---------|-------|
| sf1 | 2 |
| ss1 | 1 |
| sf2 | 1 |
| ss2 | 0 |
| contact (mod 4) | 0 → contact (no bond change) |
| ps1 | 1 |
| ps2 | 0 |
| speed | 0.0 |

**Reaction:** `b1.a0 → b1.a0`  *(speed=0.0)*

### Two-substrate bond-exchange with product state changes

Pattern: `E024445799F`

| Counter | Value |
|---------|-------|
| sf1 | 1 |
| ss1 | 1 |
| sf2 | 4 |
| ss2 | 1 |
| contact (mod 4) | 0 → contact (no bond change) |
| ps1 | 2 |
| ps2 | 0 |
| speed | 0.0 |

**Reaction:** `a1.d1 → a2.d0`  *(speed=0.0)*

## Random genome example (seed=42, len=60)

Genome: `247237C6D1DF3A0A141E10F2206D8D97AF18F0A881FE3CAAB889887C77EF`

0 valid reaction(s) decoded:


