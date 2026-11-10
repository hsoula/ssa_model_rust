# Binary rules

## Description 

Only one atom type 'a' - and two states: 0, 1 (binary)

Rules description : b1b2b3->b4b5b6 

Total : 6 binary numbers -> 64 rules for complete chemistry 

Identical rules : 
eg: 000->000 a0.a0->a0.a0 no changes
we  remove identical rules 
the 3 strong bits equal to the 3 weak bits 
r such as r&0b111000 >> 3 == r&0b000111 

56 rules in total 

010->xxx needed for starters but with random soup of a0 and a1 any rules will work 

Simple test :  9-57-43 
- a0 + a0 -> a1.a2  010-001 ->  9d
- a1 + a1 -> a1.a2  111-001 ->  57d 
- a2.a2 -> a1 + a2  101-011 ->  43d 
with only a0 