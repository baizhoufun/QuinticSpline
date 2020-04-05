# QuinticSpline

## Problem statement

Given a set of k-dimensional points    

$$
      \{\boldsymbol{x}_0, \boldsymbol{x}_1,\boldsymbol{x}_2, ...\, , \boldsymbol{x}_N\}~, 
$$

we compute a set of polynomials (i.e. quintic spline)
        
$$\{\boldsymbol{s}_0(t), \boldsymbol{s}_1(t), \boldsymbol{s}_2(t), ...\, , \boldsymbol{s}_{N-1}(t)\} $$

where each spline segment is in the form
a fifth-degree polynomial
,
$$\boldsymbol{s}_i(t)=
\boldsymbol{a}_i
+\boldsymbol{b}_it
+\boldsymbol{c}_it^2
+\boldsymbol{d}_it^3
+\boldsymbol{e}_it^4
+\boldsymbol{f}_it^5~,
$$

over interval $t\in[0,1]$ such that
 not only $\boldsymbol{s}_i(0) = \boldsymbol{x}_i$ but also
        
$$\left\{
\begin{aligned}
\boldsymbol{s}_i(1) =\boldsymbol{s}_{i+1}(0) \\
\boldsymbol{s}_i'(1) =\boldsymbol{s}_{i+1}'(0) \\
\boldsymbol{s}_i''(1) =\boldsymbol{s}_{i+1}''(0) \\
\boldsymbol{s}_i'''(1) =\boldsymbol{s}_{i+1}'''(0) \\
\boldsymbol{s}_i'''(1) =\boldsymbol{s}_{i+1}''''(0) \\
\end{aligned}\right.
$$

hold for $i = 0, 1, 2, ... \,, N-1$ 















