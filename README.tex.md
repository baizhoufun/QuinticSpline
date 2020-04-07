# QuinticSpline

## Problem statement

Given a set of k-dimensional knots in $\mathbb{R}^k$

$$\{\boldsymbol{x}_0, \boldsymbol{x}_1,\boldsymbol{x}_2, ...\, , \boldsymbol{x}_N\}~,
$$

and a set of intrinsic coordinates

$$
\{0, l_1,l_2, ...\, , l_N\}~,
$$

we compute a piecewise quintic function (i.e. a quintic spline)

$$\boldsymbol{s}(l)=\left\{
\begin{aligned}
\boldsymbol{s}_0(l) && \textrm{if}\,\, l\in&[0,l_1]~,\\
\boldsymbol{s}_1(l) && \textrm{if}\,\,l\in&[l_1,l_2]~,\\
\vdots&&&\\
\boldsymbol{s}_{N-1}(l) &&\textrm{if}\,\, l\in&[l_{N-1},l_N]~,\\
\end{aligned}\right.
$$

where each spline segment $\boldsymbol{s}_i(l)$
is in the form of a fifth-degree polynomial,

$$
\boldsymbol{s}_i(l)=
\sum_{j=0}^5
\boldsymbol{c}^{(j)}_i\Big(\frac{l-l_i}{l_{i+1}-l_i}\Big)^j~,
$$

over interval $t\in[0,1]$ such that not only
the resulting quintic spline goes through all knots,

$$
\boldsymbol{s}(l_i) =\boldsymbol{s}_i(l_i) =\boldsymbol{x}_i~,
$$

but also enforces continuity up to the fourth order,

$$\left\{
\begin{aligned}
\boldsymbol{s}_i(l_{i+1})&=\boldsymbol{s}_{i+1}(l_{i+1})~,\\
\frac{\mathrm{d}\boldsymbol{s}_i}{\mathrm{d}l}(l_{i+1}) &=
\frac{\mathrm{d}\boldsymbol{s}_{i+1}}{\mathrm{d}l}(l_{i+1})~,\\
\frac{\mathrm{d}^2\boldsymbol{s}_i}{\mathrm{d}l^2}(l_{i+1}) &=
\frac{\mathrm{d}^2\boldsymbol{s}_{i+1}}{\mathrm{d}l^2}(l_{i+1})~,\\
\frac{\mathrm{d}^3\boldsymbol{s}_i}{\mathrm{d}l^3}(l_{i+1}) &=
\frac{\mathrm{d}^3\boldsymbol{s}_{i+1}}{\mathrm{d}l^3}(l_{i+1})~,\\
\frac{\mathrm{d}^4\boldsymbol{s}_i}{\mathrm{d}l^4}(l_{i+1}) &=
\frac{\mathrm{d}^4\boldsymbol{s}_{i+1}}{\mathrm{d}l^4}(l_{i+1})~.\\
\end{aligned}\right.
$$

## Features

1. Re-parametrize spline intrinsic coordinates according to arc length 
2. Support sampling of non-knot points along spline curve
3. Can impose 1st and 2nd order boundary condition at both ends
