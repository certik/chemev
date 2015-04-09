# Introduction #

First install chemev using the installation instructions (see the main page).

# 117-isochrone fit #

This fits the CMD with a restricted set of 117 isochrones and we are interested in the best fit.

The best fit for halo, that Ondrej was able to achieve is: -558047.711479 (~8060.2 in Tom's units). Tom's best value is 8063.4943.

To calculate the best fit using the **simplex** method:

```
cd chemev/calc/halosimplex/
python 117isofit.py
```

and wait couple of hours. The best fit weights will be stored in bestfit117. This way you should be able to achieve the best fit, but I was never patient enough:
in the first 15 minutes, you'll get the likelihood of around -557938 (8279 in Tom's units), in 30 minutes: -558000.460779 tom: 8154.70522397, then in the following hours it progresses very slowly up to the best fit.

A superior and recommended method to simplex is a **BFGS** (http://en.wikipedia.org/wiki/BFGS_method) method. Implemented by http://www.ece.northwestern.edu/~nocedal/lbfgsb.html. With this method, we can achieve the best fit

henry: -558047.711479 tom: 8060.20382375

in less than 7 minutes:

```
cd chemev/calc/halobfgs/
python 117isofit.py
```

and wait 7 minutes (on intel core duo 2.4GHz). The result will be written to the bestfit117. Then you can execute
```
./prot.py
```
to generate the pictures below.


Here are the CMDs of the best fit:

![http://chemev.googlecode.com/svn/trunk/calc/halobfgs/graph.png](http://chemev.googlecode.com/svn/trunk/calc/halobfgs/graph.png)

And the isochrones weights:

![http://chemev.googlecode.com/svn/trunk/calc/halobfgs/graph-weights.png](http://chemev.googlecode.com/svn/trunk/calc/halobfgs/graph-weights.png)

You can also check the bestfit117 with mine:

[2949.41959028, 33.6824649377, 83.7193680479, 3308.50079189, 19649.5611506, 23700.6191434, 2300.13406547, 93474.9033253, 662555.634906, 14.5725882805, 28.2140115067, 154.290521416, 830.579588499, 307.661132528, 218.064090104, 1200.00030536, 3050.36699796, 26603.8351979, 3171.68652241, 34.8214584914, 12904.2119536, 133.490058355, 230.169407526, 390.401749082, 1250.63555334, 11333.1146815, 1053508.19063, 4272.77236201, 18036.4877756, 312.163860758, 702.326060352, 311.799449719, 1237.09243838, 5622.12487506, 563162.134817, 617416.903335, 10812.6718256, 20466.974502, 59042.7657227, 989.420451917, 12390.3309263, 11383.5504313, 169604.154022, 1154480.07607, 9950.37677349, 168.15392544, 20717.9035002, 20177.6129301, 83551.5103979, 7937.34049668, 7150.30850094, 213036.372334, 7094.69949367, 368.212925949, 7883.5755782, 18045.9792811, 43845.2899473, 60785.5695913, 5897.69246158, 47364.8721468, 747435.359797, 796506.304169, 3681.09904594, 7723.04657859, 1961.75127526, 7777.54990892, 30572.0430169, 114466.782006, 183779.308465, 320834.233165, 373710.445173, 3743.31281465, 2264.76256707, 14810.332248, 385.20293551, 11714.4453513, 45743.3889923, 46891.0518305, 7156.15168657, 6037.36099605, 5127.79931504, 269.784377356, 2380.77907529, 449.938373876, 2077.22581518, 310639.955855, 173309.044279, 17776.0691518, 3332.84003583, 21944.4805507, 3350.27009189, 17021.1794497, 253.872500762, 2121.67983136, 434516.021427, 326498.474806, 6177.67359521, 13550.2462819, 104351.199581, 4072.61347857, 615.36467784, 94.7923196453, 1763.9840412, 5281.85165737, 102225.72878, 24670.7758983, 2649.86771051, 93817.9638891, 81.8081125319, 94.9507955793, 63.9040685632, 884.600527732, 83890.7103162, 8174.73150571, 1915.55318393, 1691.77723794, 9312.45518624]


In the same way you can fit the disk and stream as well:

```
cd /chemev/calc/diskbfgs/
python 117isofit.py
```

producing:

henry: -379747.404264 tom: 7027.64661303

and
```
cd /chemev/calc/streambfgs/
python 117isofit.py
```

producing:

henry: -220554.820055 tom: 6767.81391931