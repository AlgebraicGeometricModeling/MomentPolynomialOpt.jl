# Moment Tools

The package provide tools for moment optimization problems of the form

```math
\begin{array}{rl}
\mathrm{inf} & \sum_i \langle f_i\star \mu_i, 1 \rangle \\
s.t. &  \sum_i \langle g_{j,i}\star \mu_i \succeq 0 \\ 
        \sum_i \langle h_{j,i}\star \mu_i = 0 \\
        \sum_i \langle p_{j,i}\star \mu_i, 1 \rangle \ge 0 \\
        \sum_i \langle q_{j,i}\star \mu_i, 1 \rangle = 0 \\
\end{array}
```
