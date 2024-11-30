# Statistics in Fortran
Statistical procedures in Fortran

Compile the RANSAC program for robust multiple linear regression with

`gfortran kind.f90 random.f90 loss_functions.f90 linear_algebra.f90 linear_regressor.f90 ransac.f90 xransac.f90 `

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
