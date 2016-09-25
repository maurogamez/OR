\ENCODING=ISO-8859-1
\Problem name: ex1

Maximize
 obj: x + y + z
Subject To
 c1: x + y  = 1
 c2: x + 5 y + 2 z <= 10
 c3: x <= 5
 c4: z >= 2
 q1: x + y + [ x ^2 - 2 x * y + 3 y ^2 ] <= 5
Bounds
      x >= 0
      y >= 0
      z >= 0
Generals
 x  y  z 
End
