Let $\varphi(t, x)$ denote the position in the PCR3BP system after time $t$ and with initial position $x$.

This program compares the derivatives of $\varphi$ obtained by two different approaches:
- calculating the derivative of the explicit solution in the normal form - implemented in [the `normal-forms` repository](https://github.com/5arcoramphu5/normal-forms)
- numerical integration - implemented in [the CAPD library](https://github.com/CAPDgroup/CAPD)

The comparison is used to validate the correctness of the derivative implementation in the normal form approach.