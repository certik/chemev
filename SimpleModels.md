# Introduction #

First install chemev using the installation instructions (see the main page).

# Simple models #

## non-gaussian ##

Execute
```
cd chemev/calc/halo_nongaussian
python fit.py
```

Currently, the simplex method works much better than bfgs. Both for ExpLogistics and ReflectLogistics the bfgs is slower.

Change the parameter `iter` at line 88 from 5 to something bigger, like 50 or 100 to get a better fit. The initial values of all parameters are randomized, so you can rerun the fit many times to see how it converges.

At the end it will print "metallicity vs age" and "sfr vs age" for both the (random) initial state and the converged fit.

I also tried simulated annealing, with no more success. It's slower than simplex. And simplex by itself is terribly slow and gets stuck at a different "fit" every time. But the annealing and bfgs behaves even worse.... Well, ill defined problem.

## brute force ##

By trying random values of all parameters in every iteration, the best fit is:

```
iteration: likelihood

1: -530468.700881
3: -533037.766045
14: -536549.618184
22: -540911.195635
81: -541088.061026
241: -544192.5148
447: -546756.018438
28092: -546783.139839
```

As you can see, it converges worse than exponentially, so it's completely inneficient.

## bfgs ##

The best fit was -5.46308E+05 (the last example):
```
At iterate   12    f= -5.49110E+05    |proj g|=  7.97344E+01
--grad 17: 3.84748005867
--grad 18: 3.83931493759

At iterate   13    f= -5.49122E+05    |proj g|=  8.76562E+01
--grad 19: 3.90191888809

At iterate   14    f= -5.49128E+05    |proj g|=  9.49688E+01

.....

At iterate   26    f= -5.49200E+05    |proj g|=  1.48344E+02
--grad 36: 3.83281207085

At iterate   27    f= -5.49205E+05    |proj g|=  1.71609E+02
--grad 37: 3.84526515007

At iterate   28    f= -5.49211E+05    |proj g|=  7.54453E+01

```

And some other run:

```
At iterate    0    f= -5.07055E+05    |proj g|=  2.28011E+03
--grad 2: 2.35150790215

At iterate    1    f= -5.10541E+05    |proj g|=  2.13931E+03
--grad 3: 2.35605692863

At iterate    2    f= -5.18564E+05    |proj g|=  1.01841E+03
--grad 4: 2.35718607903

At iterate    3    f= -5.20395E+05    |proj g|=  2.83711E+02
--grad 5: 2.3593480587

....

At iterate   54    f= -5.21368E+05    |proj g|=  2.65625E-01

```

another run:
```
At iterate    0    f= -5.24782E+05    |proj g|=  4.98662E+03
--grad 2: 3.83779811859
--grad 3: 3.84446811676

At iterate    1    f= -5.43362E+05    |proj g|=  4.57234E+02
--grad 4: 3.85416007042

At iterate    2    f= -5.43630E+05    |proj g|=  3.54672E+02
--grad 5: 3.82967615128
--grad 6: 3.82982206345

At iterate    3    f= -5.44465E+05    |proj g|=  8.15812E+02
--grad 7: 3.83754587173

At iterate    4    f= -5.45047E+05    |proj g|=  2.27594E+02
--grad 8: 3.8793900013

....

At iterate   38    f= -5.46308E+05    |proj g|=  3.76484E+01
```

As you can see, it's very important where the fit starts.

## annealing ##

```
4: -532924.880542
25: -544021.743115
53: -544573.530498
Warning: Cooled to -544059.454113 at [ 378.89861582 -356.3312686   -33.49153929 -252.41633047  -25.61446198
 -152.23255953  115.83448233 -272.57470662  -10.20756802 -251.94127281
   12.80704446 -193.57570497   59.31342182   -2.66832333  190.43410323
  -79.50789665] but this is not the smallest point found.
```

or

```
1: -533338.814798
16: -541814.354498
91: -542843.361448
92: -551104.615066
Warning: Cooled to -535401.817951 at [-114.68920893  -31.75837635  104.31319836 -285.82926233   -7.5852325
 -153.96706124  152.86936021  109.16249443   93.69621906  391.64218616
 -196.24933993   22.88764212  -10.92760903  142.75433682  -84.79176758
 -281.30257051] but this is not the smallest point found.
```

## simplex ##

```
henry: -534981.173797 tom: 54193.2791883 iter: 1
henry: -537423.892667 tom: 49307.8414476 iter: 13
henry: -544302.382697 tom: 35550.8613883 iter: 14
henry: -547265.684768 tom: 29624.2572457 iter: 17
henry: -547421.766356 tom: 29312.0940701 iter: 19
henry: -548733.514434 tom: 26688.5979142 iter: 20
henry: -549006.209884 tom: 26143.2070146 iter: 27
henry: -549012.652118 tom: 26130.3225472 iter: 33
henry: -549173.553862 tom: 25808.519058 iter: 35
henry: -549856.397806 tom: 24442.8311699 iter: 37
henry: -549938.406718 tom: 24278.8133469 iter: 39
henry: -550314.086327 tom: 23527.4541276 iter: 46
henry: -550331.251725 tom: 23493.1233332 iter: 49
henry: -550403.638676 tom: 23348.3494299 iter: 50
henry: -550893.948665 tom: 22367.7294526 iter: 61
henry: -551115.179602 tom: 21925.2675781 iter: 67
henry: -551124.583559 tom: 21906.4596641 iter: 85
henry: -551148.320722 tom: 21858.9853381 iter: 91
henry: -551150.262225 tom: 21855.1023329 iter: 98
henry: -551150.641107 tom: 21854.3445692 iter: 99
```

So simplex is the best method, but it needs to be restarted:

```
cd chemev/calc/halo_nongaussian-simplex
./fit.py
```

and wait a long time. This should in theory converge to the best result. Example:

```
ondra@pc232:~/chemev/calc/halo_nongaussian-simplex$ ./fit.py
start
### doing 50 iterations ###
henry: -544065.970737 tom: 36023.6853091 iter: 1
henry: -544403.078988 tom: 35349.4688061 iter: 5
henry: -544443.504214 tom: 35268.6183546 iter: 6
henry: -545170.267874 tom: 33815.0910338 iter: 7
henry: -545187.669986 tom: 33780.2868111 iter: 8
henry: -545194.184236 tom: 33767.2583096 iter: 9
henry: -545195.242143 tom: 33765.1424964 iter: 13
henry: -545331.446369 tom: 33492.7340444 iter: 14
henry: -545523.729272 tom: 33108.1682389 iter: 15

...... (after 30 min) 

henry: -546018.33447 tom: 32118.9578426 iter: 7
henry: -546018.33447 tom: 32118.9578426 iter: 9
henry: -546018.33447 tom: 32118.9578426 iter: 10
henry: -546018.33447 tom: 32118.9578426 iter: 11
henry: -546018.373784 tom: 32118.8792141 iter: 31
### doing 50 iterations ###
henry: -546018.373784 tom: 32118.8792144 iter: 1
henry: -546018.373784 tom: 32118.8792141 iter: 3
```

another run:

```
ondra@pc232:~/chemev/calc/halo_nongaussian-simplex$ ./fit.py
start
### doing 50 iterations ###
henry: -540008.864984 tom: 44137.8968152 iter: 1
henry: -543091.530887 tom: 37972.5650075 iter: 4
henry: -543970.689458 tom: 36214.2478655 iter: 6
henry: -544133.629978 tom: 35888.3668255 iter: 8
henry: -544627.462267 tom: 34900.7022488 iter: 10
henry: -544747.295892 tom: 34661.0349991 iter: 12
henry: -544900.060329 tom: 34355.5061237 iter: 14
henry: -544910.560695 tom: 34334.5053924 iter: 33
henry: -544910.821713 tom: 34333.9833556 iter: 41
henry: -544913.202733 tom: 34329.2213162 iter: 43
### doing 50 iterations ###
henry: -545748.426824 tom: 32658.7731338 iter: 1
henry: -546296.075175 tom: 31563.476432 iter: 2
henry: -546603.851813 tom: 30947.9231555 iter: 3
henry: -546606.089162 tom: 30943.4484587 iter: 7
henry: -546639.651301 tom: 30876.3241812 iter: 11

...

henry: -547649.067181 tom: 28857.4924201 iter: 262
henry: -547649.067181 tom: 28857.4924201 iter: 266
henry: -547649.067181 tom: 28857.4924201 iter: 267
henry: -547649.067181 tom: 28857.4924201 iter: 276
henry: -547649.067181 tom: 28857.4924201 iter: 277
```

(it doesn't change any more, but clearly this is _not_ a global minimum)

another run:
```
ondra@pc232:~/chemev/calc/halo_nongaussian-simplex$ ./fit.py
start
### doing 50 iterations ###
henry: -541942.298402 tom: 40271.0299775 iter: 1
henry: -543877.000958 tom: 36401.624867 iter: 4
henry: -544395.56861 tom: 35364.4895632 iter: 5
henry: -544787.913063 tom: 34579.800657 iter: 6
henry: -545002.844498 tom: 34149.9377857 iter: 8
henry: -546033.451388 tom: 32088.7240062 iter: 9
henry: -546622.48638 tom: 30910.6540216 iter: 12
henry: -546631.548711 tom: 30892.5293608 iter: 15
henry: -546638.120732 tom: 30879.3853179 iter: 16
henry: -546662.538626 tom: 30830.54953 iter: 19
henry: -546708.563466 tom: 30738.4998506 iter: 20

...
henry: -547935.832462 tom: 28283.9618585 iter: 260
henry: -547935.832462 tom: 28283.9618585 iter: 261
henry: -547935.832462 tom: 28283.9618585 iter: 266
henry: -547935.832462 tom: 28283.9618585 iter: 276
```

and get's stuck again.

## Differential evolution ##

run 1:
```
ondra@pc232:~/chemev/calc/halo_nongaussian-de$ ./fit.py
start
henry: 1e+20 tom: 2e+20 iter: 1
henry: -498839.401363 tom: 126476.824056 iter: 2
henry: -511176.711402 tom: 101802.203978 iter: 4
henry: -528132.817883 tom: 67889.9910153 iter: 5
henry: -544040.779415 tom: 36074.0679517 iter: 8
henry: -547537.431215 tom: 29080.7643532 iter: 36
henry: -548005.107823 tom: 28145.4111363 iter: 48
henry: -548696.808909 tom: 26762.0089642 iter: 203
henry: -548776.309903 tom: 26603.0069772 iter: 471
henry: -548877.588222 tom: 26400.4503378 iter: 649
henry: -549544.188828 tom: 25067.2491258 iter: 1618

...
```

run 2:
```
ondra@pc232:~/chemev/calc/halo_nongaussian-de$ ./fit.py
start
henry: 1e+20 tom: 2e+20 iter: 1
henry: -525756.279912 tom: 72643.066958 iter: 2
henry: -538642.076633 tom: 46871.4735154 iter: 3
henry: -546073.015773 tom: 32009.5952364 iter: 47
henry: -547267.375993 tom: 29620.8747971 iter: 145
henry: -547584.670017 tom: 28986.2867492 iter: 481
henry: -549162.128186 tom: 25831.3704102 iter: 497

...

henry: -550163.392339 tom: 23828.8421047 iter: 2350
henry: -550165.190275 tom: 23825.2462313 iter: 2941
henry: -550183.428008 tom: 23788.7707664 iter: 3259

...

henry: -551437.27611 tom: 21281.0745622 iter: 11325
henry: -551439.376299 tom: 21276.8741849 iter: 11350

...

henry: -552068.943265 tom: 20017.7402529 iter: 54949
henry: -552069.396697 tom: 20016.8333883 iter: 55101
```
(this run is a little bit slower than the other ones, but important is that it converges by itself)

run 3:
```
start
henry: 1e+20 tom: 2e+20 iter: 1
henry: -534871.469286 tom: 54412.6882108 iter: 2
henry: -534921.377024 tom: 54312.8727349 iter: 5
henry: -539595.65334 tom: 44964.3201018 iter: 6
henry: -540539.619734 tom: 43076.3873138 iter: 8
henry: -541108.397883 tom: 41938.8310162 iter: 9

...

henry: -550045.195602 tom: 24065.2355773 iter: 2001
henry: -550069.315391 tom: 24016.9960008 iter: 2011
henry: -550158.23985 tom: 23839.147082 iter: 2305

...

henry: -551943.591426 tom: 20268.4439306 iter: 6593
henry: -552347.346376 tom: 19460.9340311 iter: 6622

...

henry: -552969.649553 tom: 18216.3276771 iter: 11714
henry: -552978.195681 tom: 18199.2354195 iter: 12965
```

run 4:
```
ondra@pc232:~/chemev/calc/halo_nongaussian-de$ ./fit.py
start
henry: 1e+20 tom: 2e+20 iter: 1
henry: -510979.936703 tom: 102195.753377 iter: 2
henry: -541364.039692 tom: 41427.5473974 iter: 6
henry: -541750.975982 tom: 40653.674818 iter: 18
henry: -544609.631485 tom: 34936.3638115 iter: 28
henry: -545714.583977 tom: 32726.4588278 iter: 30

...

henry: -552164.674221 tom: 19826.2783401 iter: 8574
henry: -552260.915492 tom: 19633.7957986 iter: 8986
```

# Conclusion #

All the methods fail except simplex and DE. Even simplex get's stuck however quite soon, so only DE seems to produce likelihoods around -550000 reliably in each run.

My experience with all the other algorithms (besides DE) so far is that the random initial parameters have some likelihood and then the optimization algorithm just improves a very little bit. And the next time I run the optimization, the initial random values can easily have much bigger likelihood than the previous run ever achieved.

The only explanation is that there are many local minima and all the algorithms just get stuck in the nearest.

The problem (simple models) seems to be ill defined  - which means, that sophisticated algorithms like bfgs totally fail. The simpler the algorithm, the better.