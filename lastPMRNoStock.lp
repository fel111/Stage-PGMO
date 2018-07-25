\ SCIP STATISTICS
\   Problem name     : t_Problem_to_Solve_for_Test
\   Variables        : 65 (31 binary, 0 integer, 0 implicit integer, 34 continuous)
\   Constraints      : 49
\   Obj. scale       : 1
\   Obj. offset      : 0
Minimize
 Obj: +20.5 t_y_lt0,0 +20.5 t_y_lt0,1 +23.5 t_y_lt1,1 +23.5 t_y_lt1,2 +23.5 t_y_lt2,2 +23.5 t_y_lt2,3
      +17.5 t_y_lt3,2 +17.5 t_y_lt3,3 +17.5 t_y_lt3,4 +17.5 t_y_lt3,5 +17.5 t_y_lt3,6 +17.5 t_y_lt3,7
      +16 t_y_lt4,4 +16 t_y_lt4,5 +16 t_y_lt4,6 +16 t_y_lt4,7 +16 t_y_lt4,8 +16 t_y_lt4,9 +6 t_y_lt5,4
      +6 t_y_lt5,5 +6 t_y_lt5,6 +6 t_y_lt5,7 +6 t_y_lt5,8 +6 t_y_lt5,9 +17.5 y_lt6,0 +26.5 y_lt7,2 +17.5 y_lt8,1
      +26.5 y_lt9,1 +32.5 y_lt10,2 +22 y_lt11,4 +22 y_lt11,5 +26.5 y_lt12,2 +23.5 y_lt13,2 +32.5 y_lt14,2
Subject to
 cons_1,0,0: +1 t_x_it0,0 -1 t_y_lt0,0 = +0
 cons_1,0,1: +1 t_x_it0,1 -1 t_y_lt0,1 -1 y_lt8,1 -1 y_lt9,1 = +0
 cons_1,1,1: +1 t_x_it1,1 -1 t_y_lt1,1 -1 y_lt8,1 -1 y_lt9,1 = +0
 cons_1,1,2: +1 t_x_it1,2 -1 t_y_lt1,2 -1 y_lt10,2 = +0
 cons_1,2,0: +1 t_x_it2,0 -1 t_y_lt0,0 -1 y_lt6,0 = +0
 cons_1,2,1: +1 t_x_it2,1 -1 t_y_lt0,1 -1 t_y_lt1,1 -1 y_lt9,1 = +0
 cons_1,2,2: +1 t_x_it2,2 -1 t_y_lt1,2 -1 y_lt7,2 -1 y_lt10,2 -1 y_lt12,2 -1 y_lt13,2 -1 y_lt14,2 = +0
 cons_1,3,2: +1 t_x_it3,2 -1 t_y_lt2,2 -1 y_lt12,2 -1 y_lt13,2 -1 y_lt14,2 = +0
 cons_1,3,3: +1 t_x_it3,3 -1 t_y_lt2,3 = +0
 cons_1,4,0: +1 t_x_it4,0 -1 t_y_lt0,0 -1 y_lt6,0 = +0
 cons_1,4,1: +1 t_x_it4,1 -1 t_y_lt0,1 -1 t_y_lt1,1 -1 y_lt8,1 -1 y_lt9,1 = +0
 cons_1,4,2: +1 t_x_it4,2 -1 t_y_lt1,2 -1 t_y_lt2,2 -1 t_y_lt3,2 -1 y_lt7,2 -1 y_lt10,2 -1 y_lt13,2 -1 y_lt14,2
  = +0
 cons_1,4,3: +1 t_x_it4,3 -1 t_y_lt2,3 -1 t_y_lt3,3 = +0
 cons_1,4,4: +1 t_x_it4,4 -1 t_y_lt3,4 -1 y_lt11,4 = +0
 cons_1,4,5: +1 t_x_it4,5 -1 t_y_lt3,5 -1 y_lt11,5 = +0
 cons_1,4,6: +1 t_x_it4,6 -1 t_y_lt3,6 = +0
 cons_1,4,7: +1 t_x_it4,7 -1 t_y_lt3,7 = +0
 cons_1,5,2: +1 t_x_it5,2 -1 t_y_lt2,2 -1 t_y_lt3,2 -1 y_lt7,2 -1 y_lt10,2 -1 y_lt12,2 -1 y_lt14,2 = +0
 cons_1,5,3: +1 t_x_it5,3 -1 t_y_lt2,3 -1 t_y_lt3,3 = +0
 cons_1,5,4: +1 t_x_it5,4 -1 t_y_lt3,4 -1 t_y_lt4,4 -1 y_lt11,4 = +0
 cons_1,5,5: +1 t_x_it5,5 -1 t_y_lt3,5 -1 t_y_lt4,5 -1 y_lt11,5 = +0
 cons_1,5,6: +1 t_x_it5,6 -1 t_y_lt3,6 -1 t_y_lt4,6 = +0
 cons_1,5,7: +1 t_x_it5,7 -1 t_y_lt3,7 -1 t_y_lt4,7 = +0
 cons_1,5,8: +1 t_x_it5,8 -1 t_y_lt4,8 = +0
 cons_1,5,9: +1 t_x_it5,9 -1 t_y_lt4,9 = +0
 cons_1,6,4: +1 t_x_it6,4 -1 t_y_lt4,4 -1 t_y_lt5,4 -1 y_lt11,4 = +0
 cons_1,6,5: +1 t_x_it6,5 -1 t_y_lt4,5 -1 t_y_lt5,5 -1 y_lt11,5 = +0
 cons_1,6,6: +1 t_x_it6,6 -1 t_y_lt4,6 -1 t_y_lt5,6 = +0
 cons_1,6,7: +1 t_x_it6,7 -1 t_y_lt4,7 -1 t_y_lt5,7 = +0
 cons_1,6,8: +1 t_x_it6,8 -1 t_y_lt4,8 -1 t_y_lt5,8 = +0
 cons_1,6,9: +1 t_x_it6,9 -1 t_y_lt4,9 -1 t_y_lt5,9 = +0
 cons_2,0: +1 t_y_lt0,0 +1 t_y_lt0,1 +1 y_lt8,1 +1 y_lt9,1 >= +1
 cons_2,1: +1 t_y_lt1,1 +1 t_y_lt1,2 +1 y_lt8,1 +1 y_lt9,1 +1 y_lt10,2 >= +1
 cons_2,2: +1 t_y_lt0,0 +1 t_y_lt0,1 +1 t_y_lt1,1 +1 t_y_lt1,2 +1 y_lt6,0 +1 y_lt7,2 +1 y_lt9,1 +1 y_lt10,2
  +1 y_lt12,2 +1 y_lt13,2 +1 y_lt14,2 >= +2
 cons_2,3: +1 t_y_lt2,2 +1 t_y_lt2,3 +1 y_lt12,2 +1 y_lt13,2 +1 y_lt14,2 >= +1
 cons_2,4: +1 t_y_lt0,0 +1 t_y_lt0,1 +1 t_y_lt1,1 +1 t_y_lt1,2 +1 t_y_lt2,2 +1 t_y_lt2,3 +1 t_y_lt3,2
  +1 t_y_lt3,3 +1 t_y_lt3,4 +1 t_y_lt3,5 +1 t_y_lt3,6 +1 t_y_lt3,7 +1 y_lt6,0 +1 y_lt7,2 +1 y_lt8,1 +1 y_lt9,1
  +1 y_lt10,2 +1 y_lt11,4 +1 y_lt11,5 +1 y_lt13,2 +1 y_lt14,2 >= +4
 cons_2,5: +1 t_y_lt2,2 +1 t_y_lt2,3 +1 t_y_lt3,2 +1 t_y_lt3,3 +1 t_y_lt3,4 +1 t_y_lt3,5 +1 t_y_lt3,6
  +1 t_y_lt3,7 +1 t_y_lt4,4 +1 t_y_lt4,5 +1 t_y_lt4,6 +1 t_y_lt4,7 +1 t_y_lt4,8 +1 t_y_lt4,9 +1 y_lt7,2
  +1 y_lt10,2 +1 y_lt11,4 +1 y_lt11,5 +1 y_lt12,2 +1 y_lt14,2 >= +3
 cons_2,6: +1 t_y_lt4,4 +1 t_y_lt4,5 +1 t_y_lt4,6 +1 t_y_lt4,7 +1 t_y_lt4,8 +1 t_y_lt4,9 +1 t_y_lt5,4
  +1 t_y_lt5,5 +1 t_y_lt5,6 +1 t_y_lt5,7 +1 t_y_lt5,8 +1 t_y_lt5,9 +1 y_lt11,4 +1 y_lt11,5 >= +5
 cons_3,0: -1 t_y_lt0,0 -1 y_lt6,0 >= -1
 cons_3,1: -1 t_y_lt0,1 -1 t_y_lt1,1 -1 y_lt8,1 -1 y_lt9,1 >= -1
 cons_3,2: -1 t_y_lt1,2 -1 t_y_lt2,2 -1 t_y_lt3,2 -1 y_lt7,2 -1 y_lt10,2 -1 y_lt12,2 -1 y_lt13,2 -1 y_lt14,2
  >= -1
 cons_3,3: -1 t_y_lt2,3 -1 t_y_lt3,3 >= -1
 cons_3,4: -1 t_y_lt3,4 -1 t_y_lt4,4 -1 t_y_lt5,4 -1 y_lt11,4 >= -1
 cons_3,5: -1 t_y_lt3,5 -1 t_y_lt4,5 -1 t_y_lt5,5 -1 y_lt11,5 >= -1
 cons_3,6: -1 t_y_lt3,6 -1 t_y_lt4,6 -1 t_y_lt5,6 >= -1
 cons_3,7: -1 t_y_lt3,7 -1 t_y_lt4,7 -1 t_y_lt5,7 >= -1
 cons_3,8: -1 t_y_lt4,8 -1 t_y_lt5,8 >= -1
 cons_3,9: -1 t_y_lt4,9 -1 t_y_lt5,9 >= -1
 cons_3,10: >= -1
Bounds
 0 <= t_x_it0,0 <= 1
 0 <= t_x_it0,1 <= 1
 0 <= t_x_it1,1 <= 1
 0 <= t_x_it1,2 <= 1
 0 <= t_x_it2,0 <= 1
 0 <= t_x_it2,1 <= 1
 0 <= t_x_it2,2 <= 1
 0 <= t_x_it3,2 <= 1
 0 <= t_x_it3,3 <= 1
 0 <= t_x_it4,0 <= 1
 0 <= t_x_it4,1 <= 1
 0 <= t_x_it4,2 <= 1
 0 <= t_x_it4,3 <= 1
 0 <= t_x_it4,4 <= 1
 0 <= t_x_it4,5 <= 1
 0 <= t_x_it4,6 <= 1
 0 <= t_x_it4,7 <= 1
 0 <= t_x_it5,2 <= 1
 0 <= t_x_it5,3 <= 1
 0 <= t_x_it5,4 <= 1
 0 <= t_x_it5,5 <= 1
 0 <= t_x_it5,6 <= 1
 0 <= t_x_it5,7 <= 1
 0 <= t_x_it5,8 <= 1
 0 <= t_x_it5,9 <= 1
 0 <= t_x_it6,4 <= 1
 0 <= t_x_it6,5 <= 1
 0 <= t_x_it6,6 <= 1
 0 <= t_x_it6,7 <= 1
 0 <= t_x_it6,8 <= 1
 0 <= t_x_it6,9 <= 1
Binaries
 t_x_it0,0 t_x_it0,1 t_x_it1,1 t_x_it1,2 t_x_it2,0 t_x_it2,1 t_x_it2,2 t_x_it3,2 t_x_it3,3 t_x_it4,0 t_x_it4,1
 t_x_it4,2 t_x_it4,3 t_x_it4,4 t_x_it4,5 t_x_it4,6 t_x_it4,7 t_x_it5,2 t_x_it5,3 t_x_it5,4 t_x_it5,5 t_x_it5,6
 t_x_it5,7 t_x_it5,8 t_x_it5,9 t_x_it6,4 t_x_it6,5 t_x_it6,6 t_x_it6,7 t_x_it6,8 t_x_it6,9
End
