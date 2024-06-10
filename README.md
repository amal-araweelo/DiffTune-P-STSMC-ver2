# Auto-tuning of a P-STSMC controller applied to a machine tool drive train (via DiffTune+)

To read more on the DiffTune+ tool set visit the following repository: [DiffTune: Auto-Tuning through Auto-Differentiation](https://github.com/Sheng-Cheng/DiffTuneOpenSource/tree/main)

## Pre-requisites
* Install [casADI](https://web.casadi.org/get/) and add casADi's directory to your MATLAB's path
* Get a C/C++ compiler add-on in MATLAB in for building `mex`-files

## Dependencies
Make sure to add the `/Common`- and the generated `/mex` folders to the path of the root folder in MATLAB.

## Run

1. Run `AutoGeneration.m`
2. Run `runDiffTune.m`

You will be able to see the tracking performance and the RMSE reduction at run time. Also, an `mp4` file will be generated as well.
