clear all
clear globals
close all
clc

% --- Algorithm parameters
N                       = 2^6;    % --- Number of particles
maxNumPointsPerNode     = 1;      % --- Maximum number of particles per node
maxNumLevels            = 20;     % --- Maximum tree depth

% --- Particle coordinates
particleCoordinates     = rand(2, N); 

% --- Particle masses
particleMasses = rand(1, N) / N;

nBodyAssessment(maxNumPointsPerNode, maxNumLevels, particleCoordinates, particleMasses);





