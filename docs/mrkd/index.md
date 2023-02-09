# Moment Tools

The package provide tools for moment optimization problems on Positive Moment Sequences (PMS).
    
A PMS is a sequence of moments $\mu=(\mu_{\alpha})$ or equivalently a linear functional
$\mu: p \in \mathbb{R}[\mathbf{x}] \mapsto \langle \mu, p \rangle = \sum_{\alpha} p_{\alpha}
\mu_{\alpha}$, which is positive on the square of the polynomials:
 $\langle \mu, p^2 \rangle \geq 0$ for all $p\in \mathbb{R}[\mathbf{x}]$.


## Optimization

Optimization problems of the following form are considered:

```math
\begin{array}{rl}
\mathrm{inf}_{\mu_i \in PMS} & \sum_i \langle f_i\star \mu_i, 1 \rangle \\
s.t. &  \sum_i g_{i,j}\star \mu_i \succeq 0, \quad j=1,\ldots, n_1 \\ 
     &  \sum_i h_{i,j}\star \mu_i = 0, \quad j=1,\ldots, n_2\\
     &  \sum_i \langle p_{i,j}\star \mu_i, 1 \rangle \ge 0 , \quad j=1,\ldots, n_3\\
     &  \sum_i \langle q_{i,j}\star \mu_i, 1 \rangle = 0, \quad j=1,\ldots, n_4 \\
\end{array}
```
where
    
-  $\mu_i$ are Positive Moment Sequences,
-  $f_i, g_{i,j}, h_{i,j}, p_{i,j}, q_{i,j} \in \mathbb{R}[\mathbf{x}]$ are multivariate polynomials.

The solution of such optimization problem is approximated by the
solution of a truncated relaxation of the problem, which is a convex
optimization problem on Positive SemiDefinite matrices. Tools to
construct such moment relaxation of a given order are available in the package.

## Decomposition

Decomposition tools are available to decompose or approximate a PMS by
a weighted sum of Dirac measures:

```math
\mu \approx \sum_k \omega_k \, \delta_{\xi_k}
```        

where $\omega_k\in \mathbb{R}$ (resp. $\mathbb{C}$), $\xi_k \in \mathbb{R}^n$ (resp. $\mathbb{C}^n$) and $\delta_{\xi}$ is
the Dirac measure at the point $\xi$. 
