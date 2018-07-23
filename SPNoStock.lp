\ENCODING=ISO-8859-1
\Problem name: IloCplex

Maximize
 obj: 3 x2 + 8 x3 + 10.5 x4 - x5
Subject To
 c1: x1 <= 0
 c2: 2 x2 + 4 x3 + 6 x4 - x6  = 0
Bounds
 0 <= x1 <= 1
 0 <= x2 <= 1
 0 <= x3 <= 1
 0 <= x4 <= 1
      x5 Free
      x6 Free
Binaries
 x1  x2  x3  x4 
Pwl
 p1: x5 = x6 2 (5, 10) 1.5
End
