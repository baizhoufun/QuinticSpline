# QuinticSpline

## Problem statement

Given a set of k-dimensional points  
        
        {x0, x1, x2, ... , xN} 

Compute a set of quintic polynomials over interval `[0,1]` (i.e. quitic spline)
        
        {s0(t), s1(t), s2(t), ... , sN-1(t)} 

such that for each `i` not only `si(0) = xi` but also
        
        si(1) = si+1(0), 
        si'(1) = si+1'(0), 
        si''(1) = si+1''(0), 
        si'''(1) = si+1'''(0), 
        si''''(1) = si+1''''(0)

hold.




