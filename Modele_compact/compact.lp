\ SCIP STATISTICS
\   Problem name     : Compact
\   Variables        : 191 (130 binary, 0 integer, 0 implicit integer, 61 continuous)
\   Constraints      : 266
\   Obj. scale       : 1
\   Obj. offset      : 0
Minimize
 Obj: +1 ct0 +1 ct1 +1 ct2 +1 ct3 +1 ct4 +1 ct5 +1 ct6 +1 ct7 +1 ct8 +1 ct9
Subject to
 cons_ct0: +1 ct0 -3 xtij00 -2 xtij01 -10 pwd01 -1 xtij02 -30 pwd02 = +0
 cons_ct1: +1 ct1 -3 xtij10 -2 xtij11 -10 pwd11 -1 xtij12 -30 pwd12 = +0
 cons_ct2: +1 ct2 -3 xtij20 -2 xtij21 -10 pwd21 -1 xtij22 -30 pwd22 = +0
 cons_ct3: +1 ct3 -3 xtij30 -2 xtij31 -10 pwd31 -1 xtij32 -30 pwd32 = +0
 cons_ct4: +1 ct4 -3 xtij40 -2 xtij41 -10 pwd41 -1 xtij42 -30 pwd42 = +0
 cons_ct5: +1 ct5 -3 xtij50 -2 xtij51 -10 pwd51 -1 xtij52 -30 pwd52 = +0
 cons_ct6: +1 ct6 -3 xtij60 -2 xtij61 -10 pwd61 -1 xtij62 -30 pwd62 = +0
 cons_ct7: +1 ct7 -3 xtij70 -2 xtij71 -10 pwd71 -1 xtij72 -30 pwd72 = +0
 cons_ct8: +1 ct8 -3 xtij80 -2 xtij81 -10 pwd81 -1 xtij82 -30 pwd82 = +0
 cons_ct9: +1 ct9 -3 xtij90 -2 xtij91 -10 pwd91 -1 xtij92 -30 pwd92 = +0
 cons_pwd_binf00: +1 xtij00 >= +0
 cons_pwd_binf01: +1 xtij01 -10 pwd01 >= +0
 cons_pwd_binf02: +1 xtij02 -20 pwd02 >= +0
 cons_pwd_binf10: +1 xtij10 >= +0
 cons_pwd_binf11: +1 xtij11 -10 pwd11 >= +0
 cons_pwd_binf12: +1 xtij12 -20 pwd12 >= +0
 cons_pwd_binf20: +1 xtij20 >= +0
 cons_pwd_binf21: +1 xtij21 -10 pwd21 >= +0
 cons_pwd_binf22: +1 xtij22 -20 pwd22 >= +0
 cons_pwd_binf30: +1 xtij30 >= +0
 cons_pwd_binf31: +1 xtij31 -10 pwd31 >= +0
 cons_pwd_binf32: +1 xtij32 -20 pwd32 >= +0
 cons_pwd_binf40: +1 xtij40 >= +0
 cons_pwd_binf41: +1 xtij41 -10 pwd41 >= +0
 cons_pwd_binf42: +1 xtij42 -20 pwd42 >= +0
 cons_pwd_binf50: +1 xtij50 >= +0
 cons_pwd_binf51: +1 xtij51 -10 pwd51 >= +0
 cons_pwd_binf52: +1 xtij52 -20 pwd52 >= +0
 cons_pwd_binf60: +1 xtij60 >= +0
 cons_pwd_binf61: +1 xtij61 -10 pwd61 >= +0
 cons_pwd_binf62: +1 xtij62 -20 pwd62 >= +0
 cons_pwd_binf70: +1 xtij70 >= +0
 cons_pwd_binf71: +1 xtij71 -10 pwd71 >= +0
 cons_pwd_binf72: +1 xtij72 -20 pwd72 >= +0
 cons_pwd_binf80: +1 xtij80 >= +0
 cons_pwd_binf81: +1 xtij81 -10 pwd81 >= +0
 cons_pwd_binf82: +1 xtij82 -20 pwd82 >= +0
 cons_pwd_binf90: +1 xtij90 >= +0
 cons_pwd_binf91: +1 xtij91 -10 pwd91 >= +0
 cons_pwd_binf92: +1 xtij92 -20 pwd92 >= +0
 cons_pwd_bsup00: +1 xtij00 -9 pwd00 <= +0
 cons_pwd_bsup01: +1 xtij01 -19 pwd01 <= +0
 cons_pwd_bsup02: +1 xtij02 -2147483646 pwd02 <= +0
 cons_pwd_bsup10: +1 xtij10 -9 pwd10 <= +0
 cons_pwd_bsup11: +1 xtij11 -19 pwd11 <= +0
 cons_pwd_bsup12: +1 xtij12 -2147483646 pwd12 <= +0
 cons_pwd_bsup20: +1 xtij20 -9 pwd20 <= +0
 cons_pwd_bsup21: +1 xtij21 -19 pwd21 <= +0
 cons_pwd_bsup22: +1 xtij22 -2147483646 pwd22 <= +0
 cons_pwd_bsup30: +1 xtij30 -9 pwd30 <= +0
 cons_pwd_bsup31: +1 xtij31 -19 pwd31 <= +0
 cons_pwd_bsup32: +1 xtij32 -2147483646 pwd32 <= +0
 cons_pwd_bsup40: +1 xtij40 -9 pwd40 <= +0
 cons_pwd_bsup41: +1 xtij41 -19 pwd41 <= +0
 cons_pwd_bsup42: +1 xtij42 -2147483646 pwd42 <= +0
 cons_pwd_bsup50: +1 xtij50 -9 pwd50 <= +0
 cons_pwd_bsup51: +1 xtij51 -19 pwd51 <= +0
 cons_pwd_bsup52: +1 xtij52 -2147483646 pwd52 <= +0
 cons_pwd_bsup60: +1 xtij60 -9 pwd60 <= +0
 cons_pwd_bsup61: +1 xtij61 -19 pwd61 <= +0
 cons_pwd_bsup62: +1 xtij62 -2147483646 pwd62 <= +0
 cons_pwd_bsup70: +1 xtij70 -9 pwd70 <= +0
 cons_pwd_bsup71: +1 xtij71 -19 pwd71 <= +0
 cons_pwd_bsup72: +1 xtij72 -2147483646 pwd72 <= +0
 cons_pwd_bsup80: +1 xtij80 -9 pwd80 <= +0
 cons_pwd_bsup81: +1 xtij81 -19 pwd81 <= +0
 cons_pwd_bsup82: +1 xtij82 -2147483646 pwd82 <= +0
 cons_pwd_bsup90: +1 xtij90 -9 pwd90 <= +0
 cons_pwd_bsup91: +1 xtij91 -19 pwd91 <= +0
 cons_pwd_bsup92: +1 xtij92 -2147483646 pwd92 <= +0
 cons_pwdsum0: +1 pwd00 +1 pwd01 +1 pwd02 = +1
 cons_pwdsum1: +1 pwd10 +1 pwd11 +1 pwd12 = +1
 cons_pwdsum2: +1 pwd20 +1 pwd21 +1 pwd22 = +1
 cons_pwdsum3: +1 pwd30 +1 pwd31 +1 pwd32 = +1
 cons_pwdsum4: +1 pwd40 +1 pwd41 +1 pwd42 = +1
 cons_pwdsum5: +1 pwd50 +1 pwd51 +1 pwd52 = +1
 cons_pwdsum6: +1 pwd60 +1 pwd61 +1 pwd62 = +1
 cons_pwdsum7: +1 pwd70 +1 pwd71 +1 pwd72 = +1
 cons_pwdsum8: +1 pwd80 +1 pwd81 +1 pwd82 = +1
 cons_pwdsum9: +1 pwd90 +1 pwd91 +1 pwd92 = +1
 cons_xt0: +1 xt0 -1 xtij00 -1 xtij01 -1 xtij02 = +0
 cons_xt1: +1 xt1 -1 xtij10 -1 xtij11 -1 xtij12 = +0
 cons_xt2: +1 xt2 -1 xtij20 -1 xtij21 -1 xtij22 = +0
 cons_xt3: +1 xt3 -1 xtij30 -1 xtij31 -1 xtij32 = +0
 cons_xt4: +1 xt4 -1 xtij40 -1 xtij41 -1 xtij42 = +0
 cons_xt5: +1 xt5 -1 xtij50 -1 xtij51 -1 xtij52 = +0
 cons_xt6: +1 xt6 -1 xtij60 -1 xtij61 -1 xtij62 = +0
 cons_xt7: +1 xt7 -1 xtij70 -1 xtij71 -1 xtij72 = +0
 cons_xt8: +1 xt8 -1 xtij80 -1 xtij81 -1 xtij82 = +0
 cons_xt9: +1 xt9 -1 xtij90 -1 xtij91 -1 xtij92 = +0
 cons_y_jkt00: +1 y_jkt000 +1 y_jkt010 <= +1
 cons_y_jkt01: +1 y_jkt001 +1 y_jkt011 <= +1
 cons_y_jkt02: +1 y_jkt002 +1 y_jkt012 <= +1
 cons_y_jkt03: +1 y_jkt003 +1 y_jkt013 <= +1
 cons_y_jkt04: +1 y_jkt004 +1 y_jkt014 <= +1
 cons_y_jkt05: +1 y_jkt005 +1 y_jkt015 <= +1
 cons_y_jkt06: +1 y_jkt006 +1 y_jkt016 <= +1
 cons_y_jkt07: +1 y_jkt007 +1 y_jkt017 <= +1
 cons_y_jkt08: +1 y_jkt008 +1 y_jkt018 <= +1
 cons_y_jkt09: +1 y_jkt009 +1 y_jkt019 <= +1
 cons_y_jkt10: +1 y_jkt100 +1 y_jkt110 <= +1
 cons_y_jkt11: +1 y_jkt101 +1 y_jkt111 <= +1
 cons_y_jkt12: +1 y_jkt102 +1 y_jkt112 <= +1
 cons_y_jkt13: +1 y_jkt103 +1 y_jkt113 <= +1
 cons_y_jkt14: +1 y_jkt104 +1 y_jkt114 <= +1
 cons_y_jkt15: +1 y_jkt105 +1 y_jkt115 <= +1
 cons_y_jkt16: +1 y_jkt106 +1 y_jkt116 <= +1
 cons_y_jkt17: +1 y_jkt107 +1 y_jkt117 <= +1
 cons_y_jkt18: +1 y_jkt108 +1 y_jkt118 <= +1
 cons_y_jkt19: +1 y_jkt109 +1 y_jkt119 <= +1
 cons_y_jkt20: +1 y_jkt200 +1 y_jkt210 <= +1
 cons_y_jkt21: +1 y_jkt201 +1 y_jkt211 <= +1
 cons_y_jkt22: +1 y_jkt202 +1 y_jkt212 <= +1
 cons_y_jkt23: +1 y_jkt203 +1 y_jkt213 <= +1
 cons_y_jkt24: +1 y_jkt204 +1 y_jkt214 <= +1
 cons_y_jkt25: +1 y_jkt205 +1 y_jkt215 <= +1
 cons_y_jkt26: +1 y_jkt206 +1 y_jkt216 <= +1
 cons_y_jkt27: +1 y_jkt207 +1 y_jkt217 <= +1
 cons_y_jkt28: +1 y_jkt208 +1 y_jkt218 <= +1
 cons_y_jkt29: +1 y_jkt209 +1 y_jkt219 <= +1
 cons_y_jkt30: +1 y_jkt300 +1 y_jkt310 <= +1
 cons_y_jkt31: +1 y_jkt301 +1 y_jkt311 <= +1
 cons_y_jkt32: +1 y_jkt302 +1 y_jkt312 <= +1
 cons_y_jkt33: +1 y_jkt303 +1 y_jkt313 <= +1
 cons_y_jkt34: +1 y_jkt304 +1 y_jkt314 <= +1
 cons_y_jkt35: +1 y_jkt305 +1 y_jkt315 <= +1
 cons_y_jkt36: +1 y_jkt306 +1 y_jkt316 <= +1
 cons_y_jkt37: +1 y_jkt307 +1 y_jkt317 <= +1
 cons_y_jkt38: +1 y_jkt308 +1 y_jkt318 <= +1
 cons_y_jkt39: +1 y_jkt309 +1 y_jkt319 <= +1
 cons_y_jkt20: +1 y_jkt000 +1 y_jkt001 +1 y_jkt002 +1 y_jkt003 +1 y_jkt004 +1 y_jkt005 +1 y_jkt006 +1 y_jkt007
  +1 y_jkt008 +1 y_jkt009 +1 y_jkt010 +1 y_jkt011 +1 y_jkt012 +1 y_jkt013 +1 y_jkt014 +1 y_jkt015 +1 y_jkt016
  +1 y_jkt017 +1 y_jkt018 +1 y_jkt019 = +4
 cons_y_jkt21: +1 y_jkt100 +1 y_jkt101 +1 y_jkt102 +1 y_jkt103 +1 y_jkt104 +1 y_jkt105 +1 y_jkt106 +1 y_jkt107
  +1 y_jkt108 +1 y_jkt109 +1 y_jkt110 +1 y_jkt111 +1 y_jkt112 +1 y_jkt113 +1 y_jkt114 +1 y_jkt115 +1 y_jkt116
  +1 y_jkt117 +1 y_jkt118 +1 y_jkt119 = +2
 cons_y_jkt22: +1 y_jkt200 +1 y_jkt201 +1 y_jkt202 +1 y_jkt203 +1 y_jkt204 +1 y_jkt205 +1 y_jkt206 +1 y_jkt207
  +1 y_jkt208 +1 y_jkt209 +1 y_jkt210 +1 y_jkt211 +1 y_jkt212 +1 y_jkt213 +1 y_jkt214 +1 y_jkt215 +1 y_jkt216
  +1 y_jkt217 +1 y_jkt218 +1 y_jkt219 = +8
 cons_y_jkt23: +1 y_jkt300 +1 y_jkt301 +1 y_jkt302 +1 y_jkt303 +1 y_jkt304 +1 y_jkt305 +1 y_jkt306 +1 y_jkt307
  +1 y_jkt308 +1 y_jkt309 +1 y_jkt310 +1 y_jkt311 +1 y_jkt312 +1 y_jkt313 +1 y_jkt314 +1 y_jkt315 +1 y_jkt316
  +1 y_jkt317 +1 y_jkt318 +1 y_jkt319 = +5
 cons_ckr000: +20 y_jkt000 +15 y_jkt100 +35 y_jkt200 +10 y_jkt300 <= +30
 cons_ckr001: +20 y_jkt001 +15 y_jkt101 +35 y_jkt201 +10 y_jkt301 <= +30
 cons_ckr002: +20 y_jkt002 +15 y_jkt102 +35 y_jkt202 +10 y_jkt302 <= +30
 cons_ckr003: +20 y_jkt003 +15 y_jkt103 +35 y_jkt203 +10 y_jkt303 <= +30
 cons_ckr004: +20 y_jkt004 +15 y_jkt104 +35 y_jkt204 +10 y_jkt304 <= +30
 cons_ckr005: +20 y_jkt005 +15 y_jkt105 +35 y_jkt205 +10 y_jkt305 <= +30
 cons_ckr006: +20 y_jkt006 +15 y_jkt106 +35 y_jkt206 +10 y_jkt306 <= +30
 cons_ckr007: +20 y_jkt007 +15 y_jkt107 +35 y_jkt207 +10 y_jkt307 <= +30
 cons_ckr008: +20 y_jkt008 +15 y_jkt108 +35 y_jkt208 +10 y_jkt308 <= +30
 cons_ckr009: +20 y_jkt009 +15 y_jkt109 +35 y_jkt209 +10 y_jkt309 <= +30
 cons_ckr010: +0.400000005960464 y_jkt000 +0.200000002980232 y_jkt100 +1.5 y_jkt200 +1.70000004768372 y_jkt300
  <= +1.5
 cons_ckr011: +0.400000005960464 y_jkt001 +0.200000002980232 y_jkt101 +1.5 y_jkt201 +1.70000004768372 y_jkt301
  <= +1.5
 cons_ckr012: +0.400000005960464 y_jkt002 +0.200000002980232 y_jkt102 +1.5 y_jkt202 +1.70000004768372 y_jkt302
  <= +1.5
 cons_ckr013: +0.400000005960464 y_jkt003 +0.200000002980232 y_jkt103 +1.5 y_jkt203 +1.70000004768372 y_jkt303
  <= +1.5
 cons_ckr014: +0.400000005960464 y_jkt004 +0.200000002980232 y_jkt104 +1.5 y_jkt204 +1.70000004768372 y_jkt304
  <= +1.5
 cons_ckr015: +0.400000005960464 y_jkt005 +0.200000002980232 y_jkt105 +1.5 y_jkt205 +1.70000004768372 y_jkt305
  <= +1.5
 cons_ckr016: +0.400000005960464 y_jkt006 +0.200000002980232 y_jkt106 +1.5 y_jkt206 +1.70000004768372 y_jkt306
  <= +1.5
 cons_ckr017: +0.400000005960464 y_jkt007 +0.200000002980232 y_jkt107 +1.5 y_jkt207 +1.70000004768372 y_jkt307
  <= +1.5
 cons_ckr018: +0.400000005960464 y_jkt008 +0.200000002980232 y_jkt108 +1.5 y_jkt208 +1.70000004768372 y_jkt308
  <= +1.5
 cons_ckr019: +0.400000005960464 y_jkt009 +0.200000002980232 y_jkt109 +1.5 y_jkt209 +1.70000004768372 y_jkt309
  <= +1.5
 cons_ckr100: +20 y_jkt010 +15 y_jkt110 +35 y_jkt210 +10 y_jkt310 <= +50
 cons_ckr101: +20 y_jkt011 +15 y_jkt111 +35 y_jkt211 +10 y_jkt311 <= +50
 cons_ckr102: +20 y_jkt012 +15 y_jkt112 +35 y_jkt212 +10 y_jkt312 <= +50
 cons_ckr103: +20 y_jkt013 +15 y_jkt113 +35 y_jkt213 +10 y_jkt313 <= +50
 cons_ckr104: +20 y_jkt014 +15 y_jkt114 +35 y_jkt214 +10 y_jkt314 <= +50
 cons_ckr105: +20 y_jkt015 +15 y_jkt115 +35 y_jkt215 +10 y_jkt315 <= +50
 cons_ckr106: +20 y_jkt016 +15 y_jkt116 +35 y_jkt216 +10 y_jkt316 <= +50
 cons_ckr107: +20 y_jkt017 +15 y_jkt117 +35 y_jkt217 +10 y_jkt317 <= +50
 cons_ckr108: +20 y_jkt018 +15 y_jkt118 +35 y_jkt218 +10 y_jkt318 <= +50
 cons_ckr109: +20 y_jkt019 +15 y_jkt119 +35 y_jkt219 +10 y_jkt319 <= +50
 cons_ckr110: +0.400000005960464 y_jkt010 +0.200000002980232 y_jkt110 +1.5 y_jkt210 +1.70000004768372 y_jkt310
  <= +1.39999997615814
 cons_ckr111: +0.400000005960464 y_jkt011 +0.200000002980232 y_jkt111 +1.5 y_jkt211 +1.70000004768372 y_jkt311
  <= +1.39999997615814
 cons_ckr112: +0.400000005960464 y_jkt012 +0.200000002980232 y_jkt112 +1.5 y_jkt212 +1.70000004768372 y_jkt312
  <= +1.39999997615814
 cons_ckr113: +0.400000005960464 y_jkt013 +0.200000002980232 y_jkt113 +1.5 y_jkt213 +1.70000004768372 y_jkt313
  <= +1.39999997615814
 cons_ckr114: +0.400000005960464 y_jkt014 +0.200000002980232 y_jkt114 +1.5 y_jkt214 +1.70000004768372 y_jkt314
  <= +1.39999997615814
 cons_ckr115: +0.400000005960464 y_jkt015 +0.200000002980232 y_jkt115 +1.5 y_jkt215 +1.70000004768372 y_jkt315
  <= +1.39999997615814
 cons_ckr116: +0.400000005960464 y_jkt016 +0.200000002980232 y_jkt116 +1.5 y_jkt216 +1.70000004768372 y_jkt316
  <= +1.39999997615814
 cons_ckr117: +0.400000005960464 y_jkt017 +0.200000002980232 y_jkt117 +1.5 y_jkt217 +1.70000004768372 y_jkt317
  <= +1.39999997615814
 cons_ckr118: +0.400000005960464 y_jkt018 +0.200000002980232 y_jkt118 +1.5 y_jkt218 +1.70000004768372 y_jkt318
  <= +1.39999997615814
 cons_ckr119: +0.400000005960464 y_jkt019 +0.200000002980232 y_jkt119 +1.5 y_jkt219 +1.70000004768372 y_jkt319
  <= +1.39999997615814
 cons_zkt000: -1 y_jkt000 +1 z_kt00 >= +0
 cons_zkt001: -1 y_jkt001 +1 z_kt01 >= +0
 cons_zkt002: -1 y_jkt002 +1 z_kt02 >= +0
 cons_zkt003: -1 y_jkt003 +1 z_kt03 >= +0
 cons_zkt004: -1 y_jkt004 +1 z_kt04 >= +0
 cons_zkt005: -1 y_jkt005 +1 z_kt05 >= +0
 cons_zkt006: -1 y_jkt006 +1 z_kt06 >= +0
 cons_zkt007: -1 y_jkt007 +1 z_kt07 >= +0
 cons_zkt008: -1 y_jkt008 +1 z_kt08 >= +0
 cons_zkt009: -1 y_jkt009 +1 z_kt09 >= +0
 cons_zkt010: -1 y_jkt100 +1 z_kt00 >= +0
 cons_zkt011: -1 y_jkt101 +1 z_kt01 >= +0
 cons_zkt012: -1 y_jkt102 +1 z_kt02 >= +0
 cons_zkt013: -1 y_jkt103 +1 z_kt03 >= +0
 cons_zkt014: -1 y_jkt104 +1 z_kt04 >= +0
 cons_zkt015: -1 y_jkt105 +1 z_kt05 >= +0
 cons_zkt016: -1 y_jkt106 +1 z_kt06 >= +0
 cons_zkt017: -1 y_jkt107 +1 z_kt07 >= +0
 cons_zkt018: -1 y_jkt108 +1 z_kt08 >= +0
 cons_zkt019: -1 y_jkt109 +1 z_kt09 >= +0
 cons_zkt020: -1 y_jkt200 +1 z_kt00 >= +0
 cons_zkt021: -1 y_jkt201 +1 z_kt01 >= +0
 cons_zkt022: -1 y_jkt202 +1 z_kt02 >= +0
 cons_zkt023: -1 y_jkt203 +1 z_kt03 >= +0
 cons_zkt024: -1 y_jkt204 +1 z_kt04 >= +0
 cons_zkt025: -1 y_jkt205 +1 z_kt05 >= +0
 cons_zkt026: -1 y_jkt206 +1 z_kt06 >= +0
 cons_zkt027: -1 y_jkt207 +1 z_kt07 >= +0
 cons_zkt028: -1 y_jkt208 +1 z_kt08 >= +0
 cons_zkt029: -1 y_jkt209 +1 z_kt09 >= +0
 cons_zkt030: -1 y_jkt300 +1 z_kt00 >= +0
 cons_zkt031: -1 y_jkt301 +1 z_kt01 >= +0
 cons_zkt032: -1 y_jkt302 +1 z_kt02 >= +0
 cons_zkt033: -1 y_jkt303 +1 z_kt03 >= +0
 cons_zkt034: -1 y_jkt304 +1 z_kt04 >= +0
 cons_zkt035: -1 y_jkt305 +1 z_kt05 >= +0
 cons_zkt036: -1 y_jkt306 +1 z_kt06 >= +0
 cons_zkt037: -1 y_jkt307 +1 z_kt07 >= +0
 cons_zkt038: -1 y_jkt308 +1 z_kt08 >= +0
 cons_zkt039: -1 y_jkt309 +1 z_kt09 >= +0
 cons_zkt100: -1 y_jkt010 +1 z_kt10 >= +0
 cons_zkt101: -1 y_jkt011 +1 z_kt11 >= +0
 cons_zkt102: -1 y_jkt012 +1 z_kt12 >= +0
 cons_zkt103: -1 y_jkt013 +1 z_kt13 >= +0
 cons_zkt104: -1 y_jkt014 +1 z_kt14 >= +0
 cons_zkt105: -1 y_jkt015 +1 z_kt15 >= +0
 cons_zkt106: -1 y_jkt016 +1 z_kt16 >= +0
 cons_zkt107: -1 y_jkt017 +1 z_kt17 >= +0
 cons_zkt108: -1 y_jkt018 +1 z_kt18 >= +0
 cons_zkt109: -1 y_jkt019 +1 z_kt19 >= +0
 cons_zkt110: -1 y_jkt110 +1 z_kt10 >= +0
 cons_zkt111: -1 y_jkt111 +1 z_kt11 >= +0
 cons_zkt112: -1 y_jkt112 +1 z_kt12 >= +0
 cons_zkt113: -1 y_jkt113 +1 z_kt13 >= +0
 cons_zkt114: -1 y_jkt114 +1 z_kt14 >= +0
 cons_zkt115: -1 y_jkt115 +1 z_kt15 >= +0
 cons_zkt116: -1 y_jkt116 +1 z_kt16 >= +0
 cons_zkt117: -1 y_jkt117 +1 z_kt17 >= +0
 cons_zkt118: -1 y_jkt118 +1 z_kt18 >= +0
 cons_zkt119: -1 y_jkt119 +1 z_kt19 >= +0
 cons_zkt120: -1 y_jkt210 +1 z_kt10 >= +0
 cons_zkt121: -1 y_jkt211 +1 z_kt11 >= +0
 cons_zkt122: -1 y_jkt212 +1 z_kt12 >= +0
 cons_zkt123: -1 y_jkt213 +1 z_kt13 >= +0
 cons_zkt124: -1 y_jkt214 +1 z_kt14 >= +0
 cons_zkt125: -1 y_jkt215 +1 z_kt15 >= +0
 cons_zkt126: -1 y_jkt216 +1 z_kt16 >= +0
 cons_zkt127: -1 y_jkt217 +1 z_kt17 >= +0
 cons_zkt128: -1 y_jkt218 +1 z_kt18 >= +0
 cons_zkt129: -1 y_jkt219 +1 z_kt19 >= +0
 cons_zkt130: -1 y_jkt310 +1 z_kt10 >= +0
 cons_zkt131: -1 y_jkt311 +1 z_kt11 >= +0
 cons_zkt132: -1 y_jkt312 +1 z_kt12 >= +0
 cons_zkt133: -1 y_jkt313 +1 z_kt13 >= +0
 cons_zkt134: -1 y_jkt314 +1 z_kt14 >= +0
 cons_zkt135: -1 y_jkt315 +1 z_kt15 >= +0
 cons_zkt136: -1 y_jkt316 +1 z_kt16 >= +0
 cons_zkt137: -1 y_jkt317 +1 z_kt17 >= +0
 cons_zkt138: -1 y_jkt318 +1 z_kt18 >= +0
 cons_zkt139: -1 y_jkt319 +1 z_kt19 >= +0
 cons_equil0: +1 xt0 +1 st0 -0.200000002980232 y_jkt000 -0.300000011920929 y_jkt010 -0.100000001490116 y_jkt100
  -0.200000002980232 y_jkt110 -0.400000005960464 y_jkt200 -0.200000002980232 y_jkt210 -0.100000001490116 y_jkt300
  -0.600000023841858 y_jkt310 -10 z_kt00 -15 z_kt10 = +0
 cons_equil1: +1 xt1 +1 st1 -1 st0 -0.200000002980232 y_jkt001 -0.300000011920929 y_jkt011 -0.100000001490116 y_jkt101
  -0.200000002980232 y_jkt111 -0.400000005960464 y_jkt201 -0.200000002980232 y_jkt211 -0.100000001490116 y_jkt301
  -0.600000023841858 y_jkt311 -10 z_kt01 -15 z_kt11 = +0
 cons_equil2: +1 xt2 +1 st2 -1 st1 -0.200000002980232 y_jkt002 -0.300000011920929 y_jkt012 -0.100000001490116 y_jkt102
  -0.200000002980232 y_jkt112 -0.400000005960464 y_jkt202 -0.200000002980232 y_jkt212 -0.100000001490116 y_jkt302
  -0.600000023841858 y_jkt312 -10 z_kt02 -15 z_kt12 = +0
 cons_equil3: +1 xt3 +1 st3 -1 st2 -0.200000002980232 y_jkt003 -0.300000011920929 y_jkt013 -0.100000001490116 y_jkt103
  -0.200000002980232 y_jkt113 -0.400000005960464 y_jkt203 -0.200000002980232 y_jkt213 -0.100000001490116 y_jkt303
  -0.600000023841858 y_jkt313 -10 z_kt03 -15 z_kt13 = +0
 cons_equil4: +1 xt4 +1 st4 -1 st3 -0.200000002980232 y_jkt004 -0.300000011920929 y_jkt014 -0.100000001490116 y_jkt104
  -0.200000002980232 y_jkt114 -0.400000005960464 y_jkt204 -0.200000002980232 y_jkt214 -0.100000001490116 y_jkt304
  -0.600000023841858 y_jkt314 -10 z_kt04 -15 z_kt14 = +0
 cons_equil5: +1 xt5 +1 st5 -1 st4 -0.200000002980232 y_jkt005 -0.300000011920929 y_jkt015 -0.100000001490116 y_jkt105
  -0.200000002980232 y_jkt115 -0.400000005960464 y_jkt205 -0.200000002980232 y_jkt215 -0.100000001490116 y_jkt305
  -0.600000023841858 y_jkt315 -10 z_kt05 -15 z_kt15 = +0
 cons_equil6: +1 xt6 +1 st6 -1 st5 -0.200000002980232 y_jkt006 -0.300000011920929 y_jkt016 -0.100000001490116 y_jkt106
  -0.200000002980232 y_jkt116 -0.400000005960464 y_jkt206 -0.200000002980232 y_jkt216 -0.100000001490116 y_jkt306
  -0.600000023841858 y_jkt316 -10 z_kt06 -15 z_kt16 = +0
 cons_equil7: +1 xt7 +1 st7 -1 st6 -0.200000002980232 y_jkt007 -0.300000011920929 y_jkt017 -0.100000001490116 y_jkt107
  -0.200000002980232 y_jkt117 -0.400000005960464 y_jkt207 -0.200000002980232 y_jkt217 -0.100000001490116 y_jkt307
  -0.600000023841858 y_jkt317 -10 z_kt07 -15 z_kt17 = +0
 cons_equil8: +1 xt8 +1 st8 -1 st7 -0.200000002980232 y_jkt008 -0.300000011920929 y_jkt018 -0.100000001490116 y_jkt108
  -0.200000002980232 y_jkt118 -0.400000005960464 y_jkt208 -0.200000002980232 y_jkt218 -0.100000001490116 y_jkt308
  -0.600000023841858 y_jkt318 -10 z_kt08 -15 z_kt18 = +0
 cons_equil9: +1 xt9 +1 st9 -1 st8 -0.200000002980232 y_jkt009 -0.300000011920929 y_jkt019 -0.100000001490116 y_jkt109
  -0.200000002980232 y_jkt119 -0.400000005960464 y_jkt209 -0.200000002980232 y_jkt219 -0.100000001490116 y_jkt309
  -0.600000023841858 y_jkt319 -10 z_kt09 -15 z_kt19 = +0
 cons_st: +1 st9 -1 st0 >= +0
 cons_st0: +1 st0 = +0
Bounds
 0 <= pwd00 <= 1
 0 <= pwd01 <= 1
 0 <= pwd02 <= 1
 0 <= pwd10 <= 1
 0 <= pwd11 <= 1
 0 <= pwd12 <= 1
 0 <= pwd20 <= 1
 0 <= pwd21 <= 1
 0 <= pwd22 <= 1
 0 <= pwd30 <= 1
 0 <= pwd31 <= 1
 0 <= pwd32 <= 1
 0 <= pwd40 <= 1
 0 <= pwd41 <= 1
 0 <= pwd42 <= 1
 0 <= pwd50 <= 1
 0 <= pwd51 <= 1
 0 <= pwd52 <= 1
 0 <= pwd60 <= 1
 0 <= pwd61 <= 1
 0 <= pwd62 <= 1
 0 <= pwd70 <= 1
 0 <= pwd71 <= 1
 0 <= pwd72 <= 1
 0 <= pwd80 <= 1
 0 <= pwd81 <= 1
 0 <= pwd82 <= 1
 0 <= pwd90 <= 1
 0 <= pwd91 <= 1
 0 <= pwd92 <= 1
 0 <= y_jkt000 <= 1
 0 <= y_jkt001 <= 1
 0 <= y_jkt002 <= 1
 0 <= y_jkt003 <= 1
 0 <= y_jkt004 <= 1
 0 <= y_jkt005 <= 1
 0 <= y_jkt006 <= 1
 0 <= y_jkt007 <= 1
 0 <= y_jkt008 <= 1
 0 <= y_jkt009 <= 1
 0 <= y_jkt010 <= 1
 0 <= y_jkt011 <= 1
 0 <= y_jkt012 <= 1
 0 <= y_jkt013 <= 1
 0 <= y_jkt014 <= 1
 0 <= y_jkt015 <= 1
 0 <= y_jkt016 <= 1
 0 <= y_jkt017 <= 1
 0 <= y_jkt018 <= 1
 0 <= y_jkt019 <= 1
 0 <= y_jkt100 <= 1
 0 <= y_jkt101 <= 1
 0 <= y_jkt102 <= 1
 0 <= y_jkt103 <= 1
 0 <= y_jkt104 <= 1
 0 <= y_jkt105 <= 1
 0 <= y_jkt106 <= 1
 0 <= y_jkt107 <= 1
 0 <= y_jkt108 <= 1
 0 <= y_jkt109 <= 1
 0 <= y_jkt110 <= 1
 0 <= y_jkt111 <= 1
 0 <= y_jkt112 <= 1
 0 <= y_jkt113 <= 1
 0 <= y_jkt114 <= 1
 0 <= y_jkt115 <= 1
 0 <= y_jkt116 <= 1
 0 <= y_jkt117 <= 1
 0 <= y_jkt118 <= 1
 0 <= y_jkt119 <= 1
 0 <= y_jkt200 <= 1
 0 <= y_jkt201 <= 1
 0 <= y_jkt202 <= 1
 0 <= y_jkt203 <= 1
 0 <= y_jkt204 <= 1
 0 <= y_jkt205 <= 1
 0 <= y_jkt206 <= 1
 0 <= y_jkt207 <= 1
 0 <= y_jkt208 <= 1
 0 <= y_jkt209 <= 1
 0 <= y_jkt210 <= 1
 0 <= y_jkt211 <= 1
 0 <= y_jkt212 <= 1
 0 <= y_jkt213 <= 1
 0 <= y_jkt214 <= 1
 0 <= y_jkt215 <= 1
 0 <= y_jkt216 <= 1
 0 <= y_jkt217 <= 1
 0 <= y_jkt218 <= 1
 0 <= y_jkt219 <= 1
 0 <= y_jkt300 <= 1
 0 <= y_jkt301 <= 1
 0 <= y_jkt302 <= 1
 0 <= y_jkt303 <= 1
 0 <= y_jkt304 <= 1
 0 <= y_jkt305 <= 1
 0 <= y_jkt306 <= 1
 0 <= y_jkt307 <= 1
 0 <= y_jkt308 <= 1
 0 <= y_jkt309 <= 1
 0 <= y_jkt310 <= 1
 0 <= y_jkt311 <= 1
 0 <= y_jkt312 <= 1
 0 <= y_jkt313 <= 1
 0 <= y_jkt314 <= 1
 0 <= y_jkt315 <= 1
 0 <= y_jkt316 <= 1
 0 <= y_jkt317 <= 1
 0 <= y_jkt318 <= 1
 0 <= y_jkt319 <= 1
 0 <= z_kt00 <= 1
 0 <= z_kt01 <= 1
 0 <= z_kt02 <= 1
 0 <= z_kt03 <= 1
 0 <= z_kt04 <= 1
 0 <= z_kt05 <= 1
 0 <= z_kt06 <= 1
 0 <= z_kt07 <= 1
 0 <= z_kt08 <= 1
 0 <= z_kt09 <= 1
 0 <= z_kt10 <= 1
 0 <= z_kt11 <= 1
 0 <= z_kt12 <= 1
 0 <= z_kt13 <= 1
 0 <= z_kt14 <= 1
 0 <= z_kt15 <= 1
 0 <= z_kt16 <= 1
 0 <= z_kt17 <= 1
 0 <= z_kt18 <= 1
 0 <= z_kt19 <= 1
 0 <= st0 <= 20
 0 <= st1 <= 20
 0 <= st2 <= 20
 0 <= st3 <= 20
 0 <= st4 <= 20
 0 <= st5 <= 20
 0 <= st6 <= 20
 0 <= st7 <= 20
 0 <= st8 <= 20
 0 <= st9 <= 20
 0 <= st10 <= 20
Binaries
 pwd00 pwd01 pwd02 pwd10 pwd11 pwd12 pwd20 pwd21 pwd22 pwd30 pwd31 pwd32 pwd40 pwd41 pwd42 pwd50 pwd51
 pwd52 pwd60 pwd61 pwd62 pwd70 pwd71 pwd72 pwd80 pwd81 pwd82 pwd90 pwd91 pwd92 y_jkt000 y_jkt001 y_jkt002
 y_jkt003 y_jkt004 y_jkt005 y_jkt006 y_jkt007 y_jkt008 y_jkt009 y_jkt010 y_jkt011 y_jkt012 y_jkt013 y_jkt014
 y_jkt015 y_jkt016 y_jkt017 y_jkt018 y_jkt019 y_jkt100 y_jkt101 y_jkt102 y_jkt103 y_jkt104 y_jkt105 y_jkt106
 y_jkt107 y_jkt108 y_jkt109 y_jkt110 y_jkt111 y_jkt112 y_jkt113 y_jkt114 y_jkt115 y_jkt116 y_jkt117 y_jkt118
 y_jkt119 y_jkt200 y_jkt201 y_jkt202 y_jkt203 y_jkt204 y_jkt205 y_jkt206 y_jkt207 y_jkt208 y_jkt209 y_jkt210
 y_jkt211 y_jkt212 y_jkt213 y_jkt214 y_jkt215 y_jkt216 y_jkt217 y_jkt218 y_jkt219 y_jkt300 y_jkt301 y_jkt302
 y_jkt303 y_jkt304 y_jkt305 y_jkt306 y_jkt307 y_jkt308 y_jkt309 y_jkt310 y_jkt311 y_jkt312 y_jkt313 y_jkt314
 y_jkt315 y_jkt316 y_jkt317 y_jkt318 y_jkt319 z_kt00 z_kt01 z_kt02 z_kt03 z_kt04 z_kt05 z_kt06 z_kt07
 z_kt08 z_kt09 z_kt10 z_kt11 z_kt12 z_kt13 z_kt14 z_kt15 z_kt16 z_kt17 z_kt18 z_kt19
End
\ SCIP STATISTICS
\   Problem name     : Compact
\   Variables        : 191 (130 binary, 0 integer, 0 implicit integer, 61 continuous)
\   Constraints      : 266
\   Obj. scale       : 1
\   Obj. offset      : 0
Minimize
 Obj: +1 ct0 +1 ct1 +1 ct2 +1 ct3 +1 ct4 +1 ct5 +1 ct6 +1 ct7 +1 ct8 +1 ct9
Subject to
 cons_ct0: +1 ct0 -3 xtij00 -2 xtij01 -10 pwd01 -1 xtij02 -30 pwd02 = +0
 cons_ct1: +1 ct1 -3 xtij10 -2 xtij11 -10 pwd11 -1 xtij12 -30 pwd12 = +0
 cons_ct2: +1 ct2 -3 xtij20 -2 xtij21 -10 pwd21 -1 xtij22 -30 pwd22 = +0
 cons_ct3: +1 ct3 -3 xtij30 -2 xtij31 -10 pwd31 -1 xtij32 -30 pwd32 = +0
 cons_ct4: +1 ct4 -3 xtij40 -2 xtij41 -10 pwd41 -1 xtij42 -30 pwd42 = +0
 cons_ct5: +1 ct5 -3 xtij50 -2 xtij51 -10 pwd51 -1 xtij52 -30 pwd52 = +0
 cons_ct6: +1 ct6 -3 xtij60 -2 xtij61 -10 pwd61 -1 xtij62 -30 pwd62 = +0
 cons_ct7: +1 ct7 -3 xtij70 -2 xtij71 -10 pwd71 -1 xtij72 -30 pwd72 = +0
 cons_ct8: +1 ct8 -3 xtij80 -2 xtij81 -10 pwd81 -1 xtij82 -30 pwd82 = +0
 cons_ct9: +1 ct9 -3 xtij90 -2 xtij91 -10 pwd91 -1 xtij92 -30 pwd92 = +0
 cons_pwd_binf00: +1 xtij00 >= +0
 cons_pwd_binf01: +1 xtij01 -10 pwd01 >= +0
 cons_pwd_binf02: +1 xtij02 -20 pwd02 >= +0
 cons_pwd_binf10: +1 xtij10 >= +0
 cons_pwd_binf11: +1 xtij11 -10 pwd11 >= +0
 cons_pwd_binf12: +1 xtij12 -20 pwd12 >= +0
 cons_pwd_binf20: +1 xtij20 >= +0
 cons_pwd_binf21: +1 xtij21 -10 pwd21 >= +0
 cons_pwd_binf22: +1 xtij22 -20 pwd22 >= +0
 cons_pwd_binf30: +1 xtij30 >= +0
 cons_pwd_binf31: +1 xtij31 -10 pwd31 >= +0
 cons_pwd_binf32: +1 xtij32 -20 pwd32 >= +0
 cons_pwd_binf40: +1 xtij40 >= +0
 cons_pwd_binf41: +1 xtij41 -10 pwd41 >= +0
 cons_pwd_binf42: +1 xtij42 -20 pwd42 >= +0
 cons_pwd_binf50: +1 xtij50 >= +0
 cons_pwd_binf51: +1 xtij51 -10 pwd51 >= +0
 cons_pwd_binf52: +1 xtij52 -20 pwd52 >= +0
 cons_pwd_binf60: +1 xtij60 >= +0
 cons_pwd_binf61: +1 xtij61 -10 pwd61 >= +0
 cons_pwd_binf62: +1 xtij62 -20 pwd62 >= +0
 cons_pwd_binf70: +1 xtij70 >= +0
 cons_pwd_binf71: +1 xtij71 -10 pwd71 >= +0
 cons_pwd_binf72: +1 xtij72 -20 pwd72 >= +0
 cons_pwd_binf80: +1 xtij80 >= +0
 cons_pwd_binf81: +1 xtij81 -10 pwd81 >= +0
 cons_pwd_binf82: +1 xtij82 -20 pwd82 >= +0
 cons_pwd_binf90: +1 xtij90 >= +0
 cons_pwd_binf91: +1 xtij91 -10 pwd91 >= +0
 cons_pwd_binf92: +1 xtij92 -20 pwd92 >= +0
 cons_pwd_bsup00: +1 xtij00 -9 pwd00 <= +0
 cons_pwd_bsup01: +1 xtij01 -19 pwd01 <= +0
 cons_pwd_bsup02: +1 xtij02 -2147483646 pwd02 <= +0
 cons_pwd_bsup10: +1 xtij10 -9 pwd10 <= +0
 cons_pwd_bsup11: +1 xtij11 -19 pwd11 <= +0
 cons_pwd_bsup12: +1 xtij12 -2147483646 pwd12 <= +0
 cons_pwd_bsup20: +1 xtij20 -9 pwd20 <= +0
 cons_pwd_bsup21: +1 xtij21 -19 pwd21 <= +0
 cons_pwd_bsup22: +1 xtij22 -2147483646 pwd22 <= +0
 cons_pwd_bsup30: +1 xtij30 -9 pwd30 <= +0
 cons_pwd_bsup31: +1 xtij31 -19 pwd31 <= +0
 cons_pwd_bsup32: +1 xtij32 -2147483646 pwd32 <= +0
 cons_pwd_bsup40: +1 xtij40 -9 pwd40 <= +0
 cons_pwd_bsup41: +1 xtij41 -19 pwd41 <= +0
 cons_pwd_bsup42: +1 xtij42 -2147483646 pwd42 <= +0
 cons_pwd_bsup50: +1 xtij50 -9 pwd50 <= +0
 cons_pwd_bsup51: +1 xtij51 -19 pwd51 <= +0
 cons_pwd_bsup52: +1 xtij52 -2147483646 pwd52 <= +0
 cons_pwd_bsup60: +1 xtij60 -9 pwd60 <= +0
 cons_pwd_bsup61: +1 xtij61 -19 pwd61 <= +0
 cons_pwd_bsup62: +1 xtij62 -2147483646 pwd62 <= +0
 cons_pwd_bsup70: +1 xtij70 -9 pwd70 <= +0
 cons_pwd_bsup71: +1 xtij71 -19 pwd71 <= +0
 cons_pwd_bsup72: +1 xtij72 -2147483646 pwd72 <= +0
 cons_pwd_bsup80: +1 xtij80 -9 pwd80 <= +0
 cons_pwd_bsup81: +1 xtij81 -19 pwd81 <= +0
 cons_pwd_bsup82: +1 xtij82 -2147483646 pwd82 <= +0
 cons_pwd_bsup90: +1 xtij90 -9 pwd90 <= +0
 cons_pwd_bsup91: +1 xtij91 -19 pwd91 <= +0
 cons_pwd_bsup92: +1 xtij92 -2147483646 pwd92 <= +0
 cons_pwdsum0: +1 pwd00 +1 pwd01 +1 pwd02 = +1
 cons_pwdsum1: +1 pwd10 +1 pwd11 +1 pwd12 = +1
 cons_pwdsum2: +1 pwd20 +1 pwd21 +1 pwd22 = +1
 cons_pwdsum3: +1 pwd30 +1 pwd31 +1 pwd32 = +1
 cons_pwdsum4: +1 pwd40 +1 pwd41 +1 pwd42 = +1
 cons_pwdsum5: +1 pwd50 +1 pwd51 +1 pwd52 = +1
 cons_pwdsum6: +1 pwd60 +1 pwd61 +1 pwd62 = +1
 cons_pwdsum7: +1 pwd70 +1 pwd71 +1 pwd72 = +1
 cons_pwdsum8: +1 pwd80 +1 pwd81 +1 pwd82 = +1
 cons_pwdsum9: +1 pwd90 +1 pwd91 +1 pwd92 = +1
 cons_xt0: +1 xt0 -1 xtij00 -1 xtij01 -1 xtij02 = +0
 cons_xt1: +1 xt1 -1 xtij10 -1 xtij11 -1 xtij12 = +0
 cons_xt2: +1 xt2 -1 xtij20 -1 xtij21 -1 xtij22 = +0
 cons_xt3: +1 xt3 -1 xtij30 -1 xtij31 -1 xtij32 = +0
 cons_xt4: +1 xt4 -1 xtij40 -1 xtij41 -1 xtij42 = +0
 cons_xt5: +1 xt5 -1 xtij50 -1 xtij51 -1 xtij52 = +0
 cons_xt6: +1 xt6 -1 xtij60 -1 xtij61 -1 xtij62 = +0
 cons_xt7: +1 xt7 -1 xtij70 -1 xtij71 -1 xtij72 = +0
 cons_xt8: +1 xt8 -1 xtij80 -1 xtij81 -1 xtij82 = +0
 cons_xt9: +1 xt9 -1 xtij90 -1 xtij91 -1 xtij92 = +0
 cons_y_jkt00: +1 y_jkt000 +1 y_jkt010 <= +1
 cons_y_jkt01: +1 y_jkt001 +1 y_jkt011 <= +1
 cons_y_jkt02: +1 y_jkt002 +1 y_jkt012 <= +1
 cons_y_jkt03: +1 y_jkt003 +1 y_jkt013 <= +1
 cons_y_jkt04: +1 y_jkt004 +1 y_jkt014 <= +1
 cons_y_jkt05: +1 y_jkt005 +1 y_jkt015 <= +1
 cons_y_jkt06: +1 y_jkt006 +1 y_jkt016 <= +1
 cons_y_jkt07: +1 y_jkt007 +1 y_jkt017 <= +1
 cons_y_jkt08: +1 y_jkt008 +1 y_jkt018 <= +1
 cons_y_jkt09: +1 y_jkt009 +1 y_jkt019 <= +1
 cons_y_jkt10: +1 y_jkt100 +1 y_jkt110 <= +1
 cons_y_jkt11: +1 y_jkt101 +1 y_jkt111 <= +1
 cons_y_jkt12: +1 y_jkt102 +1 y_jkt112 <= +1
 cons_y_jkt13: +1 y_jkt103 +1 y_jkt113 <= +1
 cons_y_jkt14: +1 y_jkt104 +1 y_jkt114 <= +1
 cons_y_jkt15: +1 y_jkt105 +1 y_jkt115 <= +1
 cons_y_jkt16: +1 y_jkt106 +1 y_jkt116 <= +1
 cons_y_jkt17: +1 y_jkt107 +1 y_jkt117 <= +1
 cons_y_jkt18: +1 y_jkt108 +1 y_jkt118 <= +1
 cons_y_jkt19: +1 y_jkt109 +1 y_jkt119 <= +1
 cons_y_jkt20: +1 y_jkt200 +1 y_jkt210 <= +1
 cons_y_jkt21: +1 y_jkt201 +1 y_jkt211 <= +1
 cons_y_jkt22: +1 y_jkt202 +1 y_jkt212 <= +1
 cons_y_jkt23: +1 y_jkt203 +1 y_jkt213 <= +1
 cons_y_jkt24: +1 y_jkt204 +1 y_jkt214 <= +1
 cons_y_jkt25: +1 y_jkt205 +1 y_jkt215 <= +1
 cons_y_jkt26: +1 y_jkt206 +1 y_jkt216 <= +1
 cons_y_jkt27: +1 y_jkt207 +1 y_jkt217 <= +1
 cons_y_jkt28: +1 y_jkt208 +1 y_jkt218 <= +1
 cons_y_jkt29: +1 y_jkt209 +1 y_jkt219 <= +1
 cons_y_jkt30: +1 y_jkt300 +1 y_jkt310 <= +1
 cons_y_jkt31: +1 y_jkt301 +1 y_jkt311 <= +1
 cons_y_jkt32: +1 y_jkt302 +1 y_jkt312 <= +1
 cons_y_jkt33: +1 y_jkt303 +1 y_jkt313 <= +1
 cons_y_jkt34: +1 y_jkt304 +1 y_jkt314 <= +1
 cons_y_jkt35: +1 y_jkt305 +1 y_jkt315 <= +1
 cons_y_jkt36: +1 y_jkt306 +1 y_jkt316 <= +1
 cons_y_jkt37: +1 y_jkt307 +1 y_jkt317 <= +1
 cons_y_jkt38: +1 y_jkt308 +1 y_jkt318 <= +1
 cons_y_jkt39: +1 y_jkt309 +1 y_jkt319 <= +1
 cons_y_jkt20: +1 y_jkt000 +1 y_jkt001 +1 y_jkt002 +1 y_jkt003 +1 y_jkt004 +1 y_jkt005 +1 y_jkt006 +1 y_jkt007
  +1 y_jkt008 +1 y_jkt009 +1 y_jkt010 +1 y_jkt011 +1 y_jkt012 +1 y_jkt013 +1 y_jkt014 +1 y_jkt015 +1 y_jkt016
  +1 y_jkt017 +1 y_jkt018 +1 y_jkt019 = +4
 cons_y_jkt21: +1 y_jkt100 +1 y_jkt101 +1 y_jkt102 +1 y_jkt103 +1 y_jkt104 +1 y_jkt105 +1 y_jkt106 +1 y_jkt107
  +1 y_jkt108 +1 y_jkt109 +1 y_jkt110 +1 y_jkt111 +1 y_jkt112 +1 y_jkt113 +1 y_jkt114 +1 y_jkt115 +1 y_jkt116
  +1 y_jkt117 +1 y_jkt118 +1 y_jkt119 = +2
 cons_y_jkt22: +1 y_jkt200 +1 y_jkt201 +1 y_jkt202 +1 y_jkt203 +1 y_jkt204 +1 y_jkt205 +1 y_jkt206 +1 y_jkt207
  +1 y_jkt208 +1 y_jkt209 +1 y_jkt210 +1 y_jkt211 +1 y_jkt212 +1 y_jkt213 +1 y_jkt214 +1 y_jkt215 +1 y_jkt216
  +1 y_jkt217 +1 y_jkt218 +1 y_jkt219 = +8
 cons_y_jkt23: +1 y_jkt300 +1 y_jkt301 +1 y_jkt302 +1 y_jkt303 +1 y_jkt304 +1 y_jkt305 +1 y_jkt306 +1 y_jkt307
  +1 y_jkt308 +1 y_jkt309 +1 y_jkt310 +1 y_jkt311 +1 y_jkt312 +1 y_jkt313 +1 y_jkt314 +1 y_jkt315 +1 y_jkt316
  +1 y_jkt317 +1 y_jkt318 +1 y_jkt319 = +5
 cons_ckr000: +20 y_jkt000 +15 y_jkt100 +35 y_jkt200 +10 y_jkt300 <= +300
 cons_ckr001: +20 y_jkt001 +15 y_jkt101 +35 y_jkt201 +10 y_jkt301 <= +300
 cons_ckr002: +20 y_jkt002 +15 y_jkt102 +35 y_jkt202 +10 y_jkt302 <= +300
 cons_ckr003: +20 y_jkt003 +15 y_jkt103 +35 y_jkt203 +10 y_jkt303 <= +300
 cons_ckr004: +20 y_jkt004 +15 y_jkt104 +35 y_jkt204 +10 y_jkt304 <= +300
 cons_ckr005: +20 y_jkt005 +15 y_jkt105 +35 y_jkt205 +10 y_jkt305 <= +300
 cons_ckr006: +20 y_jkt006 +15 y_jkt106 +35 y_jkt206 +10 y_jkt306 <= +300
 cons_ckr007: +20 y_jkt007 +15 y_jkt107 +35 y_jkt207 +10 y_jkt307 <= +300
 cons_ckr008: +20 y_jkt008 +15 y_jkt108 +35 y_jkt208 +10 y_jkt308 <= +300
 cons_ckr009: +20 y_jkt009 +15 y_jkt109 +35 y_jkt209 +10 y_jkt309 <= +300
 cons_ckr010: +0.400000005960464 y_jkt000 +0.200000002980232 y_jkt100 +1.5 y_jkt200 +1.70000004768372 y_jkt300
  <= +10.5
 cons_ckr011: +0.400000005960464 y_jkt001 +0.200000002980232 y_jkt101 +1.5 y_jkt201 +1.70000004768372 y_jkt301
  <= +10.5
 cons_ckr012: +0.400000005960464 y_jkt002 +0.200000002980232 y_jkt102 +1.5 y_jkt202 +1.70000004768372 y_jkt302
  <= +10.5
 cons_ckr013: +0.400000005960464 y_jkt003 +0.200000002980232 y_jkt103 +1.5 y_jkt203 +1.70000004768372 y_jkt303
  <= +10.5
 cons_ckr014: +0.400000005960464 y_jkt004 +0.200000002980232 y_jkt104 +1.5 y_jkt204 +1.70000004768372 y_jkt304
  <= +10.5
 cons_ckr015: +0.400000005960464 y_jkt005 +0.200000002980232 y_jkt105 +1.5 y_jkt205 +1.70000004768372 y_jkt305
  <= +10.5
 cons_ckr016: +0.400000005960464 y_jkt006 +0.200000002980232 y_jkt106 +1.5 y_jkt206 +1.70000004768372 y_jkt306
  <= +10.5
 cons_ckr017: +0.400000005960464 y_jkt007 +0.200000002980232 y_jkt107 +1.5 y_jkt207 +1.70000004768372 y_jkt307
  <= +10.5
 cons_ckr018: +0.400000005960464 y_jkt008 +0.200000002980232 y_jkt108 +1.5 y_jkt208 +1.70000004768372 y_jkt308
  <= +10.5
 cons_ckr019: +0.400000005960464 y_jkt009 +0.200000002980232 y_jkt109 +1.5 y_jkt209 +1.70000004768372 y_jkt309
  <= +10.5
 cons_ckr100: +20 y_jkt010 +15 y_jkt110 +35 y_jkt210 +10 y_jkt310 <= +500
 cons_ckr101: +20 y_jkt011 +15 y_jkt111 +35 y_jkt211 +10 y_jkt311 <= +500
 cons_ckr102: +20 y_jkt012 +15 y_jkt112 +35 y_jkt212 +10 y_jkt312 <= +500
 cons_ckr103: +20 y_jkt013 +15 y_jkt113 +35 y_jkt213 +10 y_jkt313 <= +500
 cons_ckr104: +20 y_jkt014 +15 y_jkt114 +35 y_jkt214 +10 y_jkt314 <= +500
 cons_ckr105: +20 y_jkt015 +15 y_jkt115 +35 y_jkt215 +10 y_jkt315 <= +500
 cons_ckr106: +20 y_jkt016 +15 y_jkt116 +35 y_jkt216 +10 y_jkt316 <= +500
 cons_ckr107: +20 y_jkt017 +15 y_jkt117 +35 y_jkt217 +10 y_jkt317 <= +500
 cons_ckr108: +20 y_jkt018 +15 y_jkt118 +35 y_jkt218 +10 y_jkt318 <= +500
 cons_ckr109: +20 y_jkt019 +15 y_jkt119 +35 y_jkt219 +10 y_jkt319 <= +500
 cons_ckr110: +0.400000005960464 y_jkt010 +0.200000002980232 y_jkt110 +1.5 y_jkt210 +1.70000004768372 y_jkt310
  <= +10.3999996185303
 cons_ckr111: +0.400000005960464 y_jkt011 +0.200000002980232 y_jkt111 +1.5 y_jkt211 +1.70000004768372 y_jkt311
  <= +10.3999996185303
 cons_ckr112: +0.400000005960464 y_jkt012 +0.200000002980232 y_jkt112 +1.5 y_jkt212 +1.70000004768372 y_jkt312
  <= +10.3999996185303
 cons_ckr113: +0.400000005960464 y_jkt013 +0.200000002980232 y_jkt113 +1.5 y_jkt213 +1.70000004768372 y_jkt313
  <= +10.3999996185303
 cons_ckr114: +0.400000005960464 y_jkt014 +0.200000002980232 y_jkt114 +1.5 y_jkt214 +1.70000004768372 y_jkt314
  <= +10.3999996185303
 cons_ckr115: +0.400000005960464 y_jkt015 +0.200000002980232 y_jkt115 +1.5 y_jkt215 +1.70000004768372 y_jkt315
  <= +10.3999996185303
 cons_ckr116: +0.400000005960464 y_jkt016 +0.200000002980232 y_jkt116 +1.5 y_jkt216 +1.70000004768372 y_jkt316
  <= +10.3999996185303
 cons_ckr117: +0.400000005960464 y_jkt017 +0.200000002980232 y_jkt117 +1.5 y_jkt217 +1.70000004768372 y_jkt317
  <= +10.3999996185303
 cons_ckr118: +0.400000005960464 y_jkt018 +0.200000002980232 y_jkt118 +1.5 y_jkt218 +1.70000004768372 y_jkt318
  <= +10.3999996185303
 cons_ckr119: +0.400000005960464 y_jkt019 +0.200000002980232 y_jkt119 +1.5 y_jkt219 +1.70000004768372 y_jkt319
  <= +10.3999996185303
 cons_zkt000: -1 y_jkt000 +1 z_kt00 >= +0
 cons_zkt001: -1 y_jkt001 +1 z_kt01 >= +0
 cons_zkt002: -1 y_jkt002 +1 z_kt02 >= +0
 cons_zkt003: -1 y_jkt003 +1 z_kt03 >= +0
 cons_zkt004: -1 y_jkt004 +1 z_kt04 >= +0
 cons_zkt005: -1 y_jkt005 +1 z_kt05 >= +0
 cons_zkt006: -1 y_jkt006 +1 z_kt06 >= +0
 cons_zkt007: -1 y_jkt007 +1 z_kt07 >= +0
 cons_zkt008: -1 y_jkt008 +1 z_kt08 >= +0
 cons_zkt009: -1 y_jkt009 +1 z_kt09 >= +0
 cons_zkt010: -1 y_jkt100 +1 z_kt00 >= +0
 cons_zkt011: -1 y_jkt101 +1 z_kt01 >= +0
 cons_zkt012: -1 y_jkt102 +1 z_kt02 >= +0
 cons_zkt013: -1 y_jkt103 +1 z_kt03 >= +0
 cons_zkt014: -1 y_jkt104 +1 z_kt04 >= +0
 cons_zkt015: -1 y_jkt105 +1 z_kt05 >= +0
 cons_zkt016: -1 y_jkt106 +1 z_kt06 >= +0
 cons_zkt017: -1 y_jkt107 +1 z_kt07 >= +0
 cons_zkt018: -1 y_jkt108 +1 z_kt08 >= +0
 cons_zkt019: -1 y_jkt109 +1 z_kt09 >= +0
 cons_zkt020: -1 y_jkt200 +1 z_kt00 >= +0
 cons_zkt021: -1 y_jkt201 +1 z_kt01 >= +0
 cons_zkt022: -1 y_jkt202 +1 z_kt02 >= +0
 cons_zkt023: -1 y_jkt203 +1 z_kt03 >= +0
 cons_zkt024: -1 y_jkt204 +1 z_kt04 >= +0
 cons_zkt025: -1 y_jkt205 +1 z_kt05 >= +0
 cons_zkt026: -1 y_jkt206 +1 z_kt06 >= +0
 cons_zkt027: -1 y_jkt207 +1 z_kt07 >= +0
 cons_zkt028: -1 y_jkt208 +1 z_kt08 >= +0
 cons_zkt029: -1 y_jkt209 +1 z_kt09 >= +0
 cons_zkt030: -1 y_jkt300 +1 z_kt00 >= +0
 cons_zkt031: -1 y_jkt301 +1 z_kt01 >= +0
 cons_zkt032: -1 y_jkt302 +1 z_kt02 >= +0
 cons_zkt033: -1 y_jkt303 +1 z_kt03 >= +0
 cons_zkt034: -1 y_jkt304 +1 z_kt04 >= +0
 cons_zkt035: -1 y_jkt305 +1 z_kt05 >= +0
 cons_zkt036: -1 y_jkt306 +1 z_kt06 >= +0
 cons_zkt037: -1 y_jkt307 +1 z_kt07 >= +0
 cons_zkt038: -1 y_jkt308 +1 z_kt08 >= +0
 cons_zkt039: -1 y_jkt309 +1 z_kt09 >= +0
 cons_zkt100: -1 y_jkt010 +1 z_kt10 >= +0
 cons_zkt101: -1 y_jkt011 +1 z_kt11 >= +0
 cons_zkt102: -1 y_jkt012 +1 z_kt12 >= +0
 cons_zkt103: -1 y_jkt013 +1 z_kt13 >= +0
 cons_zkt104: -1 y_jkt014 +1 z_kt14 >= +0
 cons_zkt105: -1 y_jkt015 +1 z_kt15 >= +0
 cons_zkt106: -1 y_jkt016 +1 z_kt16 >= +0
 cons_zkt107: -1 y_jkt017 +1 z_kt17 >= +0
 cons_zkt108: -1 y_jkt018 +1 z_kt18 >= +0
 cons_zkt109: -1 y_jkt019 +1 z_kt19 >= +0
 cons_zkt110: -1 y_jkt110 +1 z_kt10 >= +0
 cons_zkt111: -1 y_jkt111 +1 z_kt11 >= +0
 cons_zkt112: -1 y_jkt112 +1 z_kt12 >= +0
 cons_zkt113: -1 y_jkt113 +1 z_kt13 >= +0
 cons_zkt114: -1 y_jkt114 +1 z_kt14 >= +0
 cons_zkt115: -1 y_jkt115 +1 z_kt15 >= +0
 cons_zkt116: -1 y_jkt116 +1 z_kt16 >= +0
 cons_zkt117: -1 y_jkt117 +1 z_kt17 >= +0
 cons_zkt118: -1 y_jkt118 +1 z_kt18 >= +0
 cons_zkt119: -1 y_jkt119 +1 z_kt19 >= +0
 cons_zkt120: -1 y_jkt210 +1 z_kt10 >= +0
 cons_zkt121: -1 y_jkt211 +1 z_kt11 >= +0
 cons_zkt122: -1 y_jkt212 +1 z_kt12 >= +0
 cons_zkt123: -1 y_jkt213 +1 z_kt13 >= +0
 cons_zkt124: -1 y_jkt214 +1 z_kt14 >= +0
 cons_zkt125: -1 y_jkt215 +1 z_kt15 >= +0
 cons_zkt126: -1 y_jkt216 +1 z_kt16 >= +0
 cons_zkt127: -1 y_jkt217 +1 z_kt17 >= +0
 cons_zkt128: -1 y_jkt218 +1 z_kt18 >= +0
 cons_zkt129: -1 y_jkt219 +1 z_kt19 >= +0
 cons_zkt130: -1 y_jkt310 +1 z_kt10 >= +0
 cons_zkt131: -1 y_jkt311 +1 z_kt11 >= +0
 cons_zkt132: -1 y_jkt312 +1 z_kt12 >= +0
 cons_zkt133: -1 y_jkt313 +1 z_kt13 >= +0
 cons_zkt134: -1 y_jkt314 +1 z_kt14 >= +0
 cons_zkt135: -1 y_jkt315 +1 z_kt15 >= +0
 cons_zkt136: -1 y_jkt316 +1 z_kt16 >= +0
 cons_zkt137: -1 y_jkt317 +1 z_kt17 >= +0
 cons_zkt138: -1 y_jkt318 +1 z_kt18 >= +0
 cons_zkt139: -1 y_jkt319 +1 z_kt19 >= +0
 cons_equil0: +1 xt0 +1 st0 -10 y_jkt000 -20.2999992370605 y_jkt010 -12.1000003814697 y_jkt100 -30.2000007629395 y_jkt110
  -12.3999996185303 y_jkt200 -5.19999980926514 y_jkt210 -17.1000003814697 y_jkt300 -54.5999984741211 y_jkt310
  -10 z_kt00 -15 z_kt10 = +0
 cons_equil1: +1 xt1 +1 st1 -1 st0 -10 y_jkt001 -20.2999992370605 y_jkt011 -12.1000003814697 y_jkt101
  -30.2000007629395 y_jkt111 -12.3999996185303 y_jkt201 -5.19999980926514 y_jkt211 -17.1000003814697 y_jkt301
  -54.5999984741211 y_jkt311 -10 z_kt01 -15 z_kt11 = +0
 cons_equil2: +1 xt2 +1 st2 -1 st1 -10 y_jkt002 -20.2999992370605 y_jkt012 -12.1000003814697 y_jkt102
  -30.2000007629395 y_jkt112 -12.3999996185303 y_jkt202 -5.19999980926514 y_jkt212 -17.1000003814697 y_jkt302
  -54.5999984741211 y_jkt312 -10 z_kt02 -15 z_kt12 = +0
 cons_equil3: +1 xt3 +1 st3 -1 st2 -10 y_jkt003 -20.2999992370605 y_jkt013 -12.1000003814697 y_jkt103
  -30.2000007629395 y_jkt113 -12.3999996185303 y_jkt203 -5.19999980926514 y_jkt213 -17.1000003814697 y_jkt303
  -54.5999984741211 y_jkt313 -10 z_kt03 -15 z_kt13 = +0
 cons_equil4: +1 xt4 +1 st4 -1 st3 -10 y_jkt004 -20.2999992370605 y_jkt014 -12.1000003814697 y_jkt104
  -30.2000007629395 y_jkt114 -12.3999996185303 y_jkt204 -5.19999980926514 y_jkt214 -17.1000003814697 y_jkt304
  -54.5999984741211 y_jkt314 -10 z_kt04 -15 z_kt14 = +0
 cons_equil5: +1 xt5 +1 st5 -1 st4 -10 y_jkt005 -20.2999992370605 y_jkt015 -12.1000003814697 y_jkt105
  -30.2000007629395 y_jkt115 -12.3999996185303 y_jkt205 -5.19999980926514 y_jkt215 -17.1000003814697 y_jkt305
  -54.5999984741211 y_jkt315 -10 z_kt05 -15 z_kt15 = +0
 cons_equil6: +1 xt6 +1 st6 -1 st5 -10 y_jkt006 -20.2999992370605 y_jkt016 -12.1000003814697 y_jkt106
  -30.2000007629395 y_jkt116 -12.3999996185303 y_jkt206 -5.19999980926514 y_jkt216 -17.1000003814697 y_jkt306
  -54.5999984741211 y_jkt316 -10 z_kt06 -15 z_kt16 = +0
 cons_equil7: +1 xt7 +1 st7 -1 st6 -10 y_jkt007 -20.2999992370605 y_jkt017 -12.1000003814697 y_jkt107
  -30.2000007629395 y_jkt117 -12.3999996185303 y_jkt207 -5.19999980926514 y_jkt217 -17.1000003814697 y_jkt307
  -54.5999984741211 y_jkt317 -10 z_kt07 -15 z_kt17 = +0
 cons_equil8: +1 xt8 +1 st8 -1 st7 -10 y_jkt008 -20.2999992370605 y_jkt018 -12.1000003814697 y_jkt108
  -30.2000007629395 y_jkt118 -12.3999996185303 y_jkt208 -5.19999980926514 y_jkt218 -17.1000003814697 y_jkt308
  -54.5999984741211 y_jkt318 -10 z_kt08 -15 z_kt18 = +0
 cons_equil9: +1 xt9 +1 st9 -1 st8 -10 y_jkt009 -20.2999992370605 y_jkt019 -12.1000003814697 y_jkt109
  -30.2000007629395 y_jkt119 -12.3999996185303 y_jkt209 -5.19999980926514 y_jkt219 -17.1000003814697 y_jkt309
  -54.5999984741211 y_jkt319 -10 z_kt09 -15 z_kt19 = +0
 cons_st: +1 st9 -1 st0 >= +0
 cons_st0: +1 st0 = +0
Bounds
 0 <= pwd00 <= 1
 0 <= pwd01 <= 1
 0 <= pwd02 <= 1
 0 <= pwd10 <= 1
 0 <= pwd11 <= 1
 0 <= pwd12 <= 1
 0 <= pwd20 <= 1
 0 <= pwd21 <= 1
 0 <= pwd22 <= 1
 0 <= pwd30 <= 1
 0 <= pwd31 <= 1
 0 <= pwd32 <= 1
 0 <= pwd40 <= 1
 0 <= pwd41 <= 1
 0 <= pwd42 <= 1
 0 <= pwd50 <= 1
 0 <= pwd51 <= 1
 0 <= pwd52 <= 1
 0 <= pwd60 <= 1
 0 <= pwd61 <= 1
 0 <= pwd62 <= 1
 0 <= pwd70 <= 1
 0 <= pwd71 <= 1
 0 <= pwd72 <= 1
 0 <= pwd80 <= 1
 0 <= pwd81 <= 1
 0 <= pwd82 <= 1
 0 <= pwd90 <= 1
 0 <= pwd91 <= 1
 0 <= pwd92 <= 1
 0 <= y_jkt000 <= 1
 0 <= y_jkt001 <= 1
 0 <= y_jkt002 <= 1
 0 <= y_jkt003 <= 1
 0 <= y_jkt004 <= 1
 0 <= y_jkt005 <= 1
 0 <= y_jkt006 <= 1
 0 <= y_jkt007 <= 1
 0 <= y_jkt008 <= 1
 0 <= y_jkt009 <= 1
 0 <= y_jkt010 <= 1
 0 <= y_jkt011 <= 1
 0 <= y_jkt012 <= 1
 0 <= y_jkt013 <= 1
 0 <= y_jkt014 <= 1
 0 <= y_jkt015 <= 1
 0 <= y_jkt016 <= 1
 0 <= y_jkt017 <= 1
 0 <= y_jkt018 <= 1
 0 <= y_jkt019 <= 1
 0 <= y_jkt100 <= 1
 0 <= y_jkt101 <= 1
 0 <= y_jkt102 <= 1
 0 <= y_jkt103 <= 1
 0 <= y_jkt104 <= 1
 0 <= y_jkt105 <= 1
 0 <= y_jkt106 <= 1
 0 <= y_jkt107 <= 1
 0 <= y_jkt108 <= 1
 0 <= y_jkt109 <= 1
 0 <= y_jkt110 <= 1
 0 <= y_jkt111 <= 1
 0 <= y_jkt112 <= 1
 0 <= y_jkt113 <= 1
 0 <= y_jkt114 <= 1
 0 <= y_jkt115 <= 1
 0 <= y_jkt116 <= 1
 0 <= y_jkt117 <= 1
 0 <= y_jkt118 <= 1
 0 <= y_jkt119 <= 1
 0 <= y_jkt200 <= 1
 0 <= y_jkt201 <= 1
 0 <= y_jkt202 <= 1
 0 <= y_jkt203 <= 1
 0 <= y_jkt204 <= 1
 0 <= y_jkt205 <= 1
 0 <= y_jkt206 <= 1
 0 <= y_jkt207 <= 1
 0 <= y_jkt208 <= 1
 0 <= y_jkt209 <= 1
 0 <= y_jkt210 <= 1
 0 <= y_jkt211 <= 1
 0 <= y_jkt212 <= 1
 0 <= y_jkt213 <= 1
 0 <= y_jkt214 <= 1
 0 <= y_jkt215 <= 1
 0 <= y_jkt216 <= 1
 0 <= y_jkt217 <= 1
 0 <= y_jkt218 <= 1
 0 <= y_jkt219 <= 1
 0 <= y_jkt300 <= 1
 0 <= y_jkt301 <= 1
 0 <= y_jkt302 <= 1
 0 <= y_jkt303 <= 1
 0 <= y_jkt304 <= 1
 0 <= y_jkt305 <= 1
 0 <= y_jkt306 <= 1
 0 <= y_jkt307 <= 1
 0 <= y_jkt308 <= 1
 0 <= y_jkt309 <= 1
 0 <= y_jkt310 <= 1
 0 <= y_jkt311 <= 1
 0 <= y_jkt312 <= 1
 0 <= y_jkt313 <= 1
 0 <= y_jkt314 <= 1
 0 <= y_jkt315 <= 1
 0 <= y_jkt316 <= 1
 0 <= y_jkt317 <= 1
 0 <= y_jkt318 <= 1
 0 <= y_jkt319 <= 1
 0 <= z_kt00 <= 1
 0 <= z_kt01 <= 1
 0 <= z_kt02 <= 1
 0 <= z_kt03 <= 1
 0 <= z_kt04 <= 1
 0 <= z_kt05 <= 1
 0 <= z_kt06 <= 1
 0 <= z_kt07 <= 1
 0 <= z_kt08 <= 1
 0 <= z_kt09 <= 1
 0 <= z_kt10 <= 1
 0 <= z_kt11 <= 1
 0 <= z_kt12 <= 1
 0 <= z_kt13 <= 1
 0 <= z_kt14 <= 1
 0 <= z_kt15 <= 1
 0 <= z_kt16 <= 1
 0 <= z_kt17 <= 1
 0 <= z_kt18 <= 1
 0 <= z_kt19 <= 1
 0 <= st0 <= 20
 0 <= st1 <= 20
 0 <= st2 <= 20
 0 <= st3 <= 20
 0 <= st4 <= 20
 0 <= st5 <= 20
 0 <= st6 <= 20
 0 <= st7 <= 20
 0 <= st8 <= 20
 0 <= st9 <= 20
 0 <= st10 <= 20
Binaries
 pwd00 pwd01 pwd02 pwd10 pwd11 pwd12 pwd20 pwd21 pwd22 pwd30 pwd31 pwd32 pwd40 pwd41 pwd42 pwd50 pwd51
 pwd52 pwd60 pwd61 pwd62 pwd70 pwd71 pwd72 pwd80 pwd81 pwd82 pwd90 pwd91 pwd92 y_jkt000 y_jkt001 y_jkt002
 y_jkt003 y_jkt004 y_jkt005 y_jkt006 y_jkt007 y_jkt008 y_jkt009 y_jkt010 y_jkt011 y_jkt012 y_jkt013 y_jkt014
 y_jkt015 y_jkt016 y_jkt017 y_jkt018 y_jkt019 y_jkt100 y_jkt101 y_jkt102 y_jkt103 y_jkt104 y_jkt105 y_jkt106
 y_jkt107 y_jkt108 y_jkt109 y_jkt110 y_jkt111 y_jkt112 y_jkt113 y_jkt114 y_jkt115 y_jkt116 y_jkt117 y_jkt118
 y_jkt119 y_jkt200 y_jkt201 y_jkt202 y_jkt203 y_jkt204 y_jkt205 y_jkt206 y_jkt207 y_jkt208 y_jkt209 y_jkt210
 y_jkt211 y_jkt212 y_jkt213 y_jkt214 y_jkt215 y_jkt216 y_jkt217 y_jkt218 y_jkt219 y_jkt300 y_jkt301 y_jkt302
 y_jkt303 y_jkt304 y_jkt305 y_jkt306 y_jkt307 y_jkt308 y_jkt309 y_jkt310 y_jkt311 y_jkt312 y_jkt313 y_jkt314
 y_jkt315 y_jkt316 y_jkt317 y_jkt318 y_jkt319 z_kt00 z_kt01 z_kt02 z_kt03 z_kt04 z_kt05 z_kt06 z_kt07
 z_kt08 z_kt09 z_kt10 z_kt11 z_kt12 z_kt13 z_kt14 z_kt15 z_kt16 z_kt17 z_kt18 z_kt19
End
