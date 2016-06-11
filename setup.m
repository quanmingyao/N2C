clear; clc;

cd tools
mex proximalRegC.cpp
cd ..

addpath('FISTA');
addpath('SCP');
addpath('N2C');
addpath('GDPAN');
addpath('GIST');
addpath(genpath('tools'));