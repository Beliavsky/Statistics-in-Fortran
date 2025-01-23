# Statistics in Fortran
Statistical programs in Fortran

### RANSAC program for robust multiple linear regression

Compile with

`gfortran -std=f2018 kind.f90 random.f90 loss_functions.f90 linear_algebra.f90 linear_regressor.f90 ransac.f90 xransac.f90 `

It also compiles and runs with ifort. In the output below, one sees that the RANSAC estimates of the regression coefficients are much closer to the true ones than the OLS estimates.

```
                #obs:       1000
            #bad_obs:         10
                #var:          2
            sd_noise:     1.0000
              sd_bad:  1000.0000

          true param:    10.8427     2.0371    -0.5356
        RANSAC param:    10.9294     1.9477    -0.5918
           OLS param:    11.2765     6.9411    -3.4181
          Best error:     0.2899

          true param:    11.2390     1.7272    -2.6399
        RANSAC param:    11.3630     1.4328    -2.4856
           OLS param:    13.0491     7.5301     1.3012
          Best error:     0.2844

          true param:    10.7725    -0.1381     0.7915
        RANSAC param:    10.8311    -0.0749     0.7074
           OLS param:    10.2791    -2.2170    -0.8155
          Best error:     0.2879

          true param:    10.3266    -0.4218    -0.2343
        RANSAC param:    10.5467    -0.5626    -0.2460
           OLS param:    11.3553    -4.8739    -0.4523
          Best error:     0.2945

          true param:    11.5727     0.2366     0.1696
        RANSAC param:    11.3434     0.1319     0.4588
           OLS param:    10.8339    -4.4357     1.3800
          Best error:     0.3456
```

### Compare the theoretical and empirical autocorrelations of ARMA processes

Compile with

`gfortran -std=f2018 kind.f90 stats.f90 random.f90 arma.f90 xarma_sim.f90`

Sample results:

```
                #obs:  1000000
               AR(1):   0.4000
            true ACF:   0.4000   0.1600   0.0640   0.0256   0.0102
            est. ACF:   0.3999   0.1587   0.0627   0.0228   0.0071
            est. ACF:   0.3992   0.1592   0.0632   0.0242   0.0088
            est. ACF:   0.3995   0.1609   0.0644   0.0258   0.0099
            est. ACF:   0.3985   0.1586   0.0631   0.0263   0.0117
            est. ACF:   0.4008   0.1610   0.0652   0.0269   0.0115

               MA(1):   0.5000
            true ACF:   0.4000   0.0000   0.0000   0.0000   0.0000
            true ACF:   0.4000
            est. ACF:   0.3996   0.0001   0.0009  -0.0005  -0.0005
            est. ACF:   0.4005   0.0002  -0.0007  -0.0008  -0.0020
            est. ACF:   0.3989   0.0000   0.0009   0.0009   0.0002
            est. ACF:   0.4003   0.0001   0.0004   0.0003  -0.0002
            est. ACF:   0.3998  -0.0016  -0.0019  -0.0010  -0.0002

               AR(1):   0.4000
               MA(1):   0.5000
            true ACF:   0.6545   0.2618   0.1047   0.0419   0.0168
            est. ACF:   0.6544   0.2608   0.1030   0.0403   0.0159
            est. ACF:   0.6544   0.2611   0.1039   0.0410   0.0168
            est. ACF:   0.6537   0.2605   0.1045   0.0426   0.0176
            est. ACF:   0.6551   0.2624   0.1040   0.0401   0.0156
            est. ACF:   0.6540   0.2619   0.1062   0.0438   0.0184

               MA(1):   0.5000
               MA(2):   0.3000
            true ACF:   0.4851   0.2239
            est. ACF:   0.4851   0.2239   0.0012   0.0017   0.0020
            est. ACF:   0.4862   0.2254   0.0026   0.0028   0.0036
            est. ACF:   0.4840   0.2240  -0.0005  -0.0002  -0.0003
            est. ACF:   0.4871   0.2268   0.0021   0.0018   0.0006
            est. ACF:   0.4850   0.2247   0.0015   0.0022   0.0027

               AR(1):   0.4000
               AR(2):  -0.7000
            true ACF:   0.2353  -0.6059  -0.4071   0.2613   0.3895
            est. ACF:   0.2350  -0.6074  -0.4081   0.2627   0.3915
            est. ACF:   0.2354  -0.6059  -0.4063   0.2619   0.3889
            est. ACF:   0.2351  -0.6056  -0.4068   0.2606   0.3887
            est. ACF:   0.2351  -0.6075  -0.4082   0.2634   0.3918
            est. ACF:   0.2354  -0.6053  -0.4075   0.2595   0.3889
```

### Solve an nxn set of linear equations

Compile with 

`gfortran kind.f90 linear_algebra.f90 xsolve.f90`

Sample output:

```
 n:        1000
 min, max b:   7.4492505928791530E-004  0.99843222000461973     
 min, max Ax-b:  -2.2204460492503131E-015   2.5535129566378600E-015
 sum of absolute errors:   3.6690052038212961E-013
```
