# QuinticSpline

## Problem statement

Given a set of k-dimensional knots in <img src="/tex/76b11a20d53ed4d10c9d38e8b4ecd46a.svg?invert_in_darkmode&sanitize=true" align=middle width=19.13820809999999pt height=27.91243950000002pt/>

<p align="center"><img src="/tex/ad29cc071e37ebf7a97f637df0bd36ed.svg?invert_in_darkmode&sanitize=true" align=middle width=150.0706713pt height=16.438356pt/></p>

and a set of intrinsic coordinates

<p align="center"><img src="/tex/4c70ced8b69102238cf037f31909e76b.svg?invert_in_darkmode&sanitize=true" align=middle width=122.29675755pt height=16.438356pt/></p>

we compute a piecewise quintic function (i.e. a quintic spline)

<p align="center"><img src="/tex/5dafbb30d09a9768c6ccb96a5169543b.svg?invert_in_darkmode&sanitize=true" align=middle width=249.63393674999998pt height=102.9954717pt/></p>

where each spline segment <img src="/tex/a5e6e63dee617a780bd4c040a604b665.svg?invert_in_darkmode&sanitize=true" align=middle width=32.21945099999999pt height=24.65753399999998pt/>
is in the form of a fifth-degree polynomial,

<p align="center"><img src="/tex/042cbda0c9011a36a557e7056ed8f12c.svg?invert_in_darkmode&sanitize=true" align=middle width=204.3037326pt height=49.59602339999999pt/></p>

over interval <img src="/tex/489dfef0eefc2611fce620116590ce89.svg?invert_in_darkmode&sanitize=true" align=middle width=58.903985249999984pt height=24.65753399999998pt/> such that not only
the resulting quintic spline goes through all knots,

<p align="center"><img src="/tex/c66608ece98d1e21f12a1d508fe3b60f.svg?invert_in_darkmode&sanitize=true" align=middle width=139.45180755pt height=16.438356pt/></p>

but also enforces continuity up to the fourth order,

<p align="center"><img src="/tex/13048de7793bfa44585e3e507d77d5c8.svg?invert_in_darkmode&sanitize=true" align=middle width=212.14867409999997pt height=187.3991658pt/></p>

## Features

1. Re-parametrize spline intrinsic coordinates according to arc length 
2. Support sampling of non-knot points along spline curve
3. Can impose 1st and 2nd order boundary condition at both ends
