clear all
clear globals
close all
clc

% --- Algorithm parameters
N                       = 2^6;    % --- Number of particles
maxNumPointsPerNode     = 1;      % --- Maximum number of particleCoordinates per node
maxNumLevels            = 20;     % --- Maximum tree depth

% --- Particle coordinates
particleCoordinates     = rand(2, N); 

% --- Particle masses
particleMasses = rand(1, N) / N;

[potential, tree] = nbody(particleCoordinates, particleMasses, maxNumPointsPerNode, maxNumLevels);


figure(2),tree.plotTree(4); axis off; hold on;
plot(particleCoordinates(1,:),particleCoordinates(2,:),'or','MarkerSize',4);




