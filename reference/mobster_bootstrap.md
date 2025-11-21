# Bootstrap a MOBSTER fit.

Parametric and non-parametric implementation of bootstrap estimates for
MOBSTER fits. This computation is parallel and uses `?easypar`.

## Usage

``` r
mobster_bootstrap(
  x,
  n.resamples = 100,
  bootstrap = "nonparametric",
  cores.ratio = 0.8,
  cache = NULL,
  save_data = NULL,
  ...
)
```

## Arguments

- x:

  An object of class `"dbpmm"`.

- n.resamples:

  Number of boostrap resamples.

- bootstrap:

  Type of boostrap: `"parametric"` or `"nonparametric"`

- cores.ratio:

  Ratio of cores to use for the parallel; see `?easypar`.

- cache:

  Cache for the computation; see `?easypar`.

- ...:

  fit parameters for `mobster_fit`

## Value

Data from the fits, resamples and a plottable figure.

## Examples

``` r
# Random small dataset
dataset = random_dataset(N = 200, seed = 123, Beta_variance_scaling = 100)
x = mobster_fit(dataset$data, auto_setup = 'FAST')
#>  [ MOBSTER fit ] 
#> 
#> ✔ Loaded input data, n = 200.
#> ❯ n = 200. Mixture with k = 1,2 Beta(s). Pareto tail: TRUE and FALSE. Output
#> clusters with π > 0.02 and n > 10.
#> ! mobster automatic setup FAST for the analysis.
#> ❯ Scoring (without parallel) 2 x 2 x 2 = 8 models by reICL.
#> 
#> [easypar] 2025-11-21 08:50:30.220264 - Overriding parallel execution setup [FALSE] with global option : FALSE
#> 
#> 
#> ℹ MOBSTER fits completed in 3.4s.
#> 
#> ── [ MOBSTER ] My MOBSTER model n = 200 with k = 2 Beta(s) without tail ────────
#> ● Clusters: π = 52% [C1] and 48% [C2], with π > 0.
#> ✖ No tail fit.
#> 
#> ● Beta C1 [n = 102, 52%] with mean = 0.58.
#> ● Beta C2 [n = 98, 48%] with mean = 0.11.
#> ℹ Score(s): NLL = -108.12; ICL = -182.58 (-182.58), H = 1.88 (1.88). Fit
#> converged by MM in 70 steps.

# Just 5 resamples of a nonparametric bootstrap run, disabling the parallel engine
options(easypar.parallel = FALSE)
boot_results = mobster_bootstrap(x$best, n.resamples = 5, auto_setup = 'FAST')
#>  [ MOBSTER bootstrap ~ 5 resamples from nonparametric bootstrap ] 
#> 
#> ── [ MOBSTER ] My MOBSTER model n = 200 with k = 2 Beta(s) without tail ────────
#> ● Clusters: π = 52% [C1] and 48% [C2], with π > 0.
#> ✖ No tail fit.
#> 
#> ● Beta C1 [n = 102, 52%] with mean = 0.58.
#> ● Beta C2 [n = 98, 48%] with mean = 0.11.
#> ℹ Score(s): NLL = -108.12; ICL = -182.58 (-182.58), H = 1.88 (1.88). Fit
#> converged by MM in 70 steps.
#> 
#> ℹ Creating nonparametric bootstrap resamples
#> ✔ Creating nonparametric bootstrap resamples ... done
#> 
#> 
#> ── Running fits ─────────────────────────────────── Might take some time ...  ──
#> [easypar] 2025-11-21 08:50:33.717127 - Overriding parallel execution setup [TRUE] with global option : FALSE
#>  [ MOBSTER fit ] 
#> 
#> ✔ Loaded input data, n = 200.
#> ❯ n = 200. Mixture with k = 1,2 Beta(s). Pareto tail: TRUE and FALSE. Output
#> clusters with π > 0.02 and n > 10.
#> ! mobster automatic setup FAST for the analysis.
#> ❯ Scoring (without parallel) 2 x 2 x 2 = 8 models by reICL.
#> 
#> [easypar] 2025-11-21 08:50:33.762079 - Overriding parallel execution setup [FALSE] with global option : FALSE
#> 
#> 
#> ℹ MOBSTER fits completed in 3.9s.
#> 
#> ── [ MOBSTER ] My MOBSTER model n = 200 with k = 2 Beta(s) without tail ────────
#> ● Clusters: π = 50% [C1] and 50% [C2], with π > 0.
#> ✖ No tail fit.
#> 
#> ● Beta C1 [n = 98, 50%] with mean = 0.59.
#> ● Beta C2 [n = 102, 50%] with mean = 0.11.
#> ℹ Score(s): NLL = -101.04; ICL = -168.07 (-168.07), H = 2.23 (2.23). Fit
#> interrupted by MM in 100 steps.
#>  [ MOBSTER fit ] 
#> 
#> ✔ Loaded input data, n = 200.
#> ❯ n = 200. Mixture with k = 1,2 Beta(s). Pareto tail: TRUE and FALSE. Output
#> clusters with π > 0.02 and n > 10.
#> ! mobster automatic setup FAST for the analysis.
#> ❯ Scoring (without parallel) 2 x 2 x 2 = 8 models by reICL.
#> 
#> [easypar] 2025-11-21 08:50:37.739139 - Overriding parallel execution setup [FALSE] with global option : FALSE
#> 
#> 
#> ℹ MOBSTER fits completed in 3.3s.
#> 
#> ── [ MOBSTER ] My MOBSTER model n = 200 with k = 2 Beta(s) without tail ────────
#> ● Clusters: π = 51% [C1] and 49% [C2], with π > 0.
#> ✖ No tail fit.
#> 
#> ● Beta C1 [n = 101, 51%] with mean = 0.59.
#> ● Beta C2 [n = 99, 49%] with mean = 0.11.
#> ℹ Score(s): NLL = -115.21; ICL = -197.42 (-197.42), H = 1.21 (1.21). Fit
#> converged by MM in 62 steps.
#>  [ MOBSTER fit ] 
#> 
#> ✔ Loaded input data, n = 200.
#> ❯ n = 200. Mixture with k = 1,2 Beta(s). Pareto tail: TRUE and FALSE. Output
#> clusters with π > 0.02 and n > 10.
#> ! mobster automatic setup FAST for the analysis.
#> ❯ Scoring (without parallel) 2 x 2 x 2 = 8 models by reICL.
#> 
#> [easypar] 2025-11-21 08:50:41.090482 - Overriding parallel execution setup [FALSE] with global option : FALSE
#> 
#> 
#> ℹ MOBSTER fits completed in 3.4s.
#> 
#> ── [ MOBSTER ] My MOBSTER model n = 200 with k = 2 Beta(s) without tail ────────
#> ● Clusters: π = 50% [C2] and 50% [C1], with π > 0.
#> ✖ No tail fit.
#> 
#> ● Beta C1 [n = 99, 50%] with mean = 0.59.
#> ● Beta C2 [n = 101, 50%] with mean = 0.12.
#> ℹ Score(s): NLL = -108.28; ICL = -182.9 (-182.9), H = 1.87 (1.87). Fit
#> converged by MM in 70 steps.
#>  [ MOBSTER fit ] 
#> 
#> ✔ Loaded input data, n = 200.
#> ❯ n = 200. Mixture with k = 1,2 Beta(s). Pareto tail: TRUE and FALSE. Output
#> clusters with π > 0.02 and n > 10.
#> ! mobster automatic setup FAST for the analysis.
#> ❯ Scoring (without parallel) 2 x 2 x 2 = 8 models by reICL.
#> 
#> [easypar] 2025-11-21 08:50:44.554851 - Overriding parallel execution setup [FALSE] with global option : FALSE
#> 
#> 
#> ℹ MOBSTER fits completed in 3.6s.
#> 
#> ── [ MOBSTER ] My MOBSTER model n = 200 with k = 2 Beta(s) without tail ────────
#> ● Clusters: π = 57% [C1] and 43% [C2], with π > 0.
#> ✖ No tail fit.
#> 
#> ● Beta C1 [n = 113, 57%] with mean = 0.59.
#> ● Beta C2 [n = 87, 43%] with mean = 0.11.
#> ℹ Score(s): NLL = -115.77; ICL = -199.22 (-199.22), H = 0.52 (0.52). Fit
#> converged by MM in 25 steps.
#>  [ MOBSTER fit ] 
#> 
#> ✔ Loaded input data, n = 200.
#> ❯ n = 200. Mixture with k = 1,2 Beta(s). Pareto tail: TRUE and FALSE. Output
#> clusters with π > 0.02 and n > 10.
#> ! mobster automatic setup FAST for the analysis.
#> ❯ Scoring (without parallel) 2 x 2 x 2 = 8 models by reICL.
#> 
#> [easypar] 2025-11-21 08:50:48.279657 - Overriding parallel execution setup [FALSE] with global option : FALSE
#> 
#> 
#> ℹ MOBSTER fits completed in 3.7s.
#> 
#> ── [ MOBSTER ] My MOBSTER model n = 200 with k = 2 Beta(s) without tail ────────
#> ● Clusters: π = 55% [C1] and 45% [C2], with π > 0.
#> ✖ No tail fit.
#> 
#> ● Beta C1 [n = 109, 55%] with mean = 0.59.
#> ● Beta C2 [n = 91, 45%] with mean = 0.13.
#> ℹ Score(s): NLL = -102.08; ICL = -171.04 (-171.04), H = 1.34 (1.34). Fit
#> converged by MM in 29 steps.

# The resample data is available in a list
print(boot_results$resamples[[1]])
#> [[1]]
#>      id        VAF original.id
#> 1     1 0.52355677          12
#> 2     2 0.84646937           3
#> 3     3 0.05393297         137
#> 4     4 0.43401627          14
#> 5     5 0.08026471         141
#> 6     6 0.08394884         192
#> 7     7 0.04428176         148
#> 8     8 0.08810669         106
#> 9     9 0.04328793         119
#> 10   10 0.47680339          16
#> 11   11 0.74474824          80
#> 12   12 0.46164409          62
#> 13   13 0.03345377          98
#> 14   14 0.06208767         123
#> 15   15 0.72430352          60
#> 16   16 0.04140095         105
#> 17   17 0.53245520          32
#> 18   18 0.63063057          25
#> 19   19 0.60339029          36
#> 20   20 0.59300401         166
#> 21   21 0.05393297         137
#> 22   22 0.14138817         133
#> 23   23 0.16810276         145
#> 24   24 0.02120704         127
#> 25   25 0.16810276         145
#> 26   26 0.16570645         165
#> 27   27 0.55147229          58
#> 28   28 0.08810669         106
#> 29   29 0.67301958         195
#> 30   30 0.71187995          13
#> 31   31 0.62264994          67
#> 32   32 0.11033293         177
#> 33   33 0.54932251          74
#> 34   34 0.71166697          56
#> 35   35 0.05714139         107
#> 36   36 0.30370976         121
#> 37   37 0.51231087          91
#> 38   38 0.14894925         116
#> 39   39 0.65434252          34
#> 40   40 0.08293964         172
#> 41   41 0.09575708         163
#> 42   42 0.17316129         109
#> 43   43 0.06036449         179
#> 44   44 0.64917596          55
#> 45   45 0.64982215          23
#> 46   46 0.52103494          95
#> 47   47 0.14648825         154
#> 48   48 0.08318635         183
#> 49   49 0.10144527         135
#> 50   50 0.13620251         143
#> 51   51 0.58490241          46
#> 52   52 0.20522784         140
#> 53   53 0.58780098          90
#> 54   54 0.20998868         124
#> 55   55 0.16312791          99
#> 56   56 0.50179798           5
#> 57   57 0.05714139         107
#> 58   58 0.53701189          68
#> 59   59 0.30370976         121
#> 60   60 0.65404812          11
#> 61   61 0.48903433          49
#> 62   62 0.62264994          67
#> 63   63 0.54932251          74
#> 64   64 0.71007760          15
#> 65   65 0.53021870           7
#> 66   66 0.62684659          42
#> 67   67 0.58780098          90
#> 68   68 0.64917596          55
#> 69   69 0.09192788         132
#> 70   70 0.09192788         132
#> 71   71 0.17316129         109
#> 72   72 0.53245520          32
#> 73   73 0.69144367          57
#> 74   74 0.09462104         193
#> 75   75 0.25434749         176
#> 76   76 0.51683769          86
#> 77   77 0.64917596          55
#> 78   78 0.16810276         145
#> 79   79 0.84646937           3
#> 80   80 0.57785081          18
#> 81   81 0.15346219         197
#> 82   82 0.10144527         135
#> 83   83 0.52103494          95
#> 84   84 0.04428176         148
#> 85   85 0.20522784         140
#> 86   86 0.03345377          98
#> 87   87 0.53701189          68
#> 88   88 0.09575708         163
#> 89   89 0.02233042         134
#> 90   90 0.07211605         120
#> 91   91 0.64982215          23
#> 92   92 0.02979764         111
#> 93   93 0.53701189          68
#> 94   94 0.12286390         187
#> 95   95 0.46040413          30
#> 96   96 0.14806814         178
#> 97   97 0.57466625          31
#> 98   98 0.16312791          99
#> 99   99 0.18577837         153
#> 100 100 0.25434749         176
#> 101 101 0.63063057          25
#> 102 102 0.27925367         156
#> 103 103 0.51106345          76
#> 104 104 0.47006049          92
#> 105 105 0.08415006         100
#> 106 106 0.69142420           8
#> 107 107 0.52259202          78
#> 108 108 0.17617266         146
#> 109 109 0.04428176         148
#> 110 110 0.63865136          85
#> 111 111 0.74760698          21
#> 112 112 0.66034681          47
#> 113 113 0.65404812          11
#> 114 114 0.63495385          51
#> 115 115 0.54101481          71
#> 116 116 0.08293964         172
#> 117 117 0.07584066         152
#> 118 118 0.57319452          63
#> 119 119 0.46164409          62
#> 120 120 0.04140095         105
#> 121 121 0.08391171         110
#> 122 122 0.07533878         200
#> 123 123 0.71166697          56
#> 124 124 0.53788908          37
#> 125 125 0.65404812          11
#> 126 126 0.05393297         137
#> 127 127 0.11370181         164
#> 128 128 0.60126291          10
#> 129 129 0.11806906         199
#> 130 130 0.55147229          58
#> 131 131 0.68651983          77
#> 132 132 0.03345377          98
#> 133 133 0.49730321          82
#> 134 134 0.29749078         174
#> 135 135 0.12928620         129
#> 136 136 0.14806814         178
#> 137 137 0.04140095         105
#> 138 138 0.08318635         183
#> 139 139 0.49730321          82
#> 140 140 0.51796731          93
#> 141 141 0.10662581         189
#> 142 142 0.15346219         197
#> 143 143 0.71187995          13
#> 144 144 0.05902402         175
#> 145 145 0.05902402         175
#> 146 146 0.61440798          79
#> 147 147 0.12286390         187
#> 148 148 0.46164409          62
#> 149 149 0.05278961         161
#> 150 150 0.04328793         119
#> 151 151 0.23287406         190
#> 152 152 0.08668096         122
#> 153 153 0.15346219         197
#> 154 154 0.43401627          14
#> 155 155 0.67301958         195
#> 156 156 0.57466625          31
#> 157 157 0.08434598         136
#> 158 158 0.02979764         111
#> 159 159 0.29749078         174
#> 160 160 0.51683769          86
#> 161 161 0.72430352          60
#> 162 162 0.46164409          62
#> 163 163 0.53161319          66
#> 164 164 0.61513016           1
#> 165 165 0.07562712         113
#> 166 166 0.51683769          86
#> 167 167 0.67912870          38
#> 168 168 0.46627543          88
#> 169 169 0.71166697          56
#> 170 170 0.08810669         106
#> 171 171 0.13620251         143
#> 172 172 0.11612617         138
#> 173 173 0.49730321          82
#> 174 174 0.06208767         123
#> 175 175 0.70284180          94
#> 176 176 0.06208767         123
#> 177 177 0.42751004          81
#> 178 178 0.67862073           4
#> 179 179 0.08747705         168
#> 180 180 0.08668096         122
#> 181 181 0.65434252          34
#> 182 182 0.08415006         100
#> 183 183 0.67301958         195
#> 184 184 0.08391171         110
#> 185 185 0.57319452          63
#> 186 186 0.62264994          67
#> 187 187 0.58635506          59
#> 188 188 0.61645161          45
#> 189 189 0.06731439         185
#> 190 190 0.52103494          95
#> 191 191 0.20522784         140
#> 192 192 0.18409048         159
#> 193 193 0.45979987          40
#> 194 194 0.08434598         136
#> 195 195 0.44827642          24
#> 196 196 0.06010055         191
#> 197 197 0.05393297         137
#> 198 198 0.15942611         160
#> 199 199 0.71007760          15
#> 200 200 0.05310553         170
#> 

# The best fits are returned
print(boot_results$fits)
#> $`1`
#> ── [ MOBSTER ] My MOBSTER model n = 200 with k = 2 Beta(s) without tail ────────
#> ● Clusters: π = 50% [C1] and 50% [C2], with π > 0.
#> ✖ No tail fit.
#> 
#> ● Beta C1 [n = 98, 50%] with mean = 0.59.
#> ● Beta C2 [n = 102, 50%] with mean = 0.11.
#> ℹ Score(s): NLL = -101.04; ICL = -168.07 (-168.07), H = 2.23 (2.23). Fit
#> interrupted by MM in 100 steps.
#> 
#> $`2`
#> ── [ MOBSTER ] My MOBSTER model n = 200 with k = 2 Beta(s) without tail ────────
#> ● Clusters: π = 51% [C1] and 49% [C2], with π > 0.
#> ✖ No tail fit.
#> 
#> ● Beta C1 [n = 101, 51%] with mean = 0.59.
#> ● Beta C2 [n = 99, 49%] with mean = 0.11.
#> ℹ Score(s): NLL = -115.21; ICL = -197.42 (-197.42), H = 1.21 (1.21). Fit
#> converged by MM in 62 steps.
#> 
#> $`3`
#> ── [ MOBSTER ] My MOBSTER model n = 200 with k = 2 Beta(s) without tail ────────
#> ● Clusters: π = 50% [C2] and 50% [C1], with π > 0.
#> ✖ No tail fit.
#> 
#> ● Beta C1 [n = 99, 50%] with mean = 0.59.
#> ● Beta C2 [n = 101, 50%] with mean = 0.12.
#> ℹ Score(s): NLL = -108.28; ICL = -182.9 (-182.9), H = 1.87 (1.87). Fit
#> converged by MM in 70 steps.
#> 
#> $`4`
#> ── [ MOBSTER ] My MOBSTER model n = 200 with k = 2 Beta(s) without tail ────────
#> ● Clusters: π = 57% [C1] and 43% [C2], with π > 0.
#> ✖ No tail fit.
#> 
#> ● Beta C1 [n = 113, 57%] with mean = 0.59.
#> ● Beta C2 [n = 87, 43%] with mean = 0.11.
#> ℹ Score(s): NLL = -115.77; ICL = -199.22 (-199.22), H = 0.52 (0.52). Fit
#> converged by MM in 25 steps.
#> 
#> $`5`
#> ── [ MOBSTER ] My MOBSTER model n = 200 with k = 2 Beta(s) without tail ────────
#> ● Clusters: π = 55% [C1] and 45% [C2], with π > 0.
#> ✖ No tail fit.
#> 
#> ● Beta C1 [n = 109, 55%] with mean = 0.59.
#> ● Beta C2 [n = 91, 45%] with mean = 0.13.
#> ℹ Score(s): NLL = -102.08; ICL = -171.04 (-171.04), H = 1.34 (1.34). Fit
#> converged by MM in 29 steps.
#> 
```
