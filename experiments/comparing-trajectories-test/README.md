Let $\varphi(t, x)$ denote the position in the PCR3BP system after time $t$ and with initial position $x$.

This program compares the values of $\varphi$ obtained by two different approaches:
- explicit solution in the normal form - implemented in [the `normal-forms` repository](https://github.com/5arcoramphu5/normal-forms)
- numerical integration - implemented in [the CAPD library](https://github.com/CAPDgroup/CAPD)

The comparison is used to validate the correctness the implementation of the algorithm in the normal form approach.