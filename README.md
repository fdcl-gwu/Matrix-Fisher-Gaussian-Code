# Matrix Fisher-Gaussian Distribution
This repository contains Matlab implementations for using matrix Fisher-Gaussian distribution in estimating the attitude and gyro-bias concurrently.

The mathematical fomulation of the presented algorithms are available at the following paper:

- W. Wang and T. Lee , ["*Matrix Fisher-Gaussian Distribution on SO(3)×R^n for Attitude Estimation with a Gyro Bias*"](https://arxiv.org/abs/2003.02180) 	arXiv:2003.02180, 2020

## Major functions
```
MFGAnalytic.m : filtering with MFG using analytical moment propagation
MFGUnscented.m: filtering with MFG using unscented moment propagation
MEKF.m : filtering using the standard multiplicative extended Kalman filter
```

## First step to use
Try the following script
```
test.m : generate a reference motion and implement the above three filters
```
