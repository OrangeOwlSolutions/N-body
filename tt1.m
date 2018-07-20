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
globalIDs               = 1 : N;  % --- Global particle IDs
quadTreeObject          = qtree;  % --- Initializing the quad tree object  
% --- Inserting the particles
quadTreeObject.insertPoints(globalIDs, particleCoordinates, maxNumPointsPerNode, maxNumLevels);
% --- Plotting the tree
quadTreeObject.plotTree;
axis off; hold on;
plot(particleCoordinates(1, :),particleCoordinates(2, :), 'or', 'MarkerSize', 5);
axis off;

if verbose
   % --- Prints Morton IDs for all nodes
   disp('All nodes');
   quadTreeObject.printMortonIDs;
   % --- Prints Morton IDs for leaves only
   disp('Leaves only');
   quadTreeObject.printMortonIDs(true);
end

depth = findDepth(quadTreeObject);
fprintf('Tree depth is %d\n', depth);


