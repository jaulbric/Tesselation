# Tesselation

## Gauge Invariant Operator
Lowest dimension gauge-invariant operators in SU(N) Yang-Mills theories.

We look at operators of the form

$$\mathcal{L}_{\mathrm{int}} \supset \frac{1}{\Lambda^{k}} X Q_{i} G^{\otimes k}_{N},$$

where $X \cong \left(p_{1}, p_{2}, \ldots, p_{N-1} \right)$ is some arbitrary representation of $\mathrm{SU}\left(N\right)$, $Q_{i}$ are generalized matter fields in the i<sup>th</sup> fundamental representation, and $G_{N} \cong \left(1, 0, \ldots, 0, 1\right)$ is the adjoint representation. For a given $X$ we then seek the minimum number of copies $k_{\mathrm{min}}$ of the adjoint such that this operator is invariant under a global gauge transformation.

## Installation

The code is very simple to install. Simply clone this repository and navigate to the root directory and do

```bash
pip install .
```

## Code

The code is a simple collection of functions, all located in `tesselation/functions.py`. An example of their use is in `Tesselation.ipynb`. There are four primary functions of interest to the user:

```python
l_pair(p, i)
```
> Pair of l that results in the optimal irreducible representation in tensor product of X times Q_i.
>
> Parameters:<br>
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;p : array_like<br>
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Dynkin label of X.<br>
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;i : int<br>
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Index of only nonzero Dynkin coefficient in Q_i.
>
> Returns:<br>
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;l1, l2 : int, int<br>
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Pair of l such that p[l1] !=0, p[l2] != 0, l1 is largest l such that l <= N - i, and is smallest l such that l > l1.

```python
Nality(p)
```
> N-ality of irreducible representation X. The number of boxes in the Young tableau for X mod N.
>
> Parameters:<br>
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;p : array_like<br>
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Dynkin label of X.<br>
>
> Returns:<br>
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;t : int<br>
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;The number of boxes in the Young tableau for X mod N.

```python
fj(p)
```
> Optimum irreducible representation in the tensor product of X times Q_i.
>
> Parameters:<br>
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;p  : array_like<br>
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Dynkin label of X.<br>
>
> Returns:<br>
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;fj : array_like<br>
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Resulting Dynkin label for optimal irreducible representation.<br>

```python
kmin(p)
```
> Minimum number of copies of the adjoint such that tensor product X Q_i G_N ... G_N  contains a trivial subspace.
>
> Parameters:<br>
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;p : array_like<br>
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Dynkin label of irreducible representation X.<br>
>
> Returns:<br>
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;kmin : int<br>
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Minimum number of copies of the adjoint representation such that the operator is gauge invariant.
