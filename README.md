# QuinticSpline

## Problem statement

Given a set of k-dimensional points  
        
        {x0, x1, x2, ... , xN} 

<p align="center"><img src="/tex/14dc45ce284c7354c0d8ca7534e97ef4.svg?invert_in_darkmode&sanitize=true" align=middle width=137.28534434999997pt height=16.438356pt/></p>

Compute a set of quintic polynomials over interval `[0,1]` (i.e. quitic spline)
        
        {s0(t), s1(t), s2(t), ... , sN-1(t)} 

such that not only `si(0) = xi` but also
        
        si(1) = si+1(0), 
        si'(1) = si+1'(0), 
        si''(1) = si+1''(0), 
        si'''(1) = si+1'''(0), 
        si''''(1) = si+1''''(0)

hold for `i = 0, 1, 2, ... , N-1` where each spline segment 
is a fifth-degree polynomial,

        si()²











