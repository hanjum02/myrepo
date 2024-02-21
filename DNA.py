#!/usr/bin/python3
s = open('rosalind_dna.txt','r').read()
a = s.count("A")
c = s.count("C")
g = s.count("G")
t = s.count("T")
print(a, c, g, t)
s = open('rosalind_dna.txt','r').read()
for n in ["A","C","G","T"]:
    print(s.count(n), end=' ')
