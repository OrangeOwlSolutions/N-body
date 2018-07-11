clear all
clear globals
close all
clc

% --- Note: the particles are assumed to reside within the [0, 1] x [0, 1] square

% --- Algorithm parameters
N                       = 2^4;    % --- Number of particles
maxNumPointsPerNode     = 1;      % --- Maximum number of points per node
maxNumLevels            = 20;     % --- Maximum tree depth
verbose                 = true;

% --- Particle coordinates
particleCoordinates     = rand(2, N); 

% --- Running the algorithm
globalIds               = 1 : N;  % --- Global particle IDs
quadTreeObject          = qtree;  % --- Initializing the quad tree object  
quadTreeObject.insertPoints(globalIds, particleCoordinates, maxNumPointsPerNode, maxNumLevels);
quadTreeObject.plottree;
axis off; hold on;
plot(particleCoordinates(1,:),particleCoordinates(2,:),'or','MarkerSize',5);
axis off;

if verbose
   % print morton ids, all nodes
   disp(' all nodes');
   quadTreeObject.print_mids;
   % print morton ids, leaves only
   disp('  leaves only');
   quadTreeObject.print_mids(true);
end
depth=find_depth(quadTreeObject);
fprintf('tree depth is %d\n', depth);


