# Hamiltonian of a Majorana Wire

Majorana particles can be found in condensed matter systems as zero energy bound states. Here we develop two different codes for implement the well known model [1,2] of Majoranas. One code in python, and other in Mathematica. 

## Python implementation

The code found in `majorana_wire.py` was initally developed during Jul/Aug 2018 to study Majorana Bound States in a semiconducting nanowire following references [1,2]. We obtain the spectra and wavefunctions. Then, we use it to obtain the local density of states in the following way:

```python
ham = Hamiltonian(N)
ham = Set_Hamiltonian(ham,N)
eigenvalues,eigenvectors = eigenv(ham)
ldos = surfLDOS(eigenvalues,eigenvectors,N,q)
```


## Mathematica implementation

## References

  [1] Paweł Szumniak, Denis Chevallier, Daniel Loss, and Jelena Klinovaja. Phys. Rev. B 96, 041401(R) – Published 5 July 2017.
  
  [2] Denis Chevallier and Jelena Klinovaja. Phys. Rev. B 94, 035417 – Published 12 July 2016.
