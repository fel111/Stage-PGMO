\ENCODING=ISO-8859-1
\Problem name: IloCplex

Maximize
 obj: - 10 x5 - 3000 x6 - 10 x8
Subject To
 c1: x1 <= 0
 c2: x2 <= 0
 c3: x3 <= 0
 c4: x4 <= 0
 c5: 5 x5 + 5 x6 <= 20
 c6: - 1000 x5 - 1000 x6 - 5 x7 - 5 x8 <= 20
 c7: x5 + x6 + x7 + x8  = 1
Bounds
 0 <= x1 <= 1
 0 <= x2 <= 1
 0 <= x3 <= 1
 0 <= x4 <= 1
 0 <= x5 <= 1
 0 <= x6 <= 1
 0 <= x7 <= 1
 0 <= x8 <= 1
Binaries
 x1  x2  x3  x4  x5  x6  x7  x8 
End
