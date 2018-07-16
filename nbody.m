function [potential, tree] = nbody(particleCoordinates, particleMasses, maxNumPointsPerNode, maxNumLevels)
% function [potential, tree] = nbody(particleCoordinates, particleMasses, maxNumPointsPerNode, maxNumLevels)
%
% --- particleCoordinates           : 2 x N array of the point coordinates (N is the number of particles)
% --- particleMasses                : particle masses - all particle masses must be positive
% --- maxNumPointsPerNode           : used to decide wether to split the node or not  
%
% --- potential                     : computed potential at each point
% --- tree                          : pointer to the tree of the tree  
%  
% --- Assumptions:
%     the algorithm applies to arbitrary particleMasses, but the simple
%         averaging implemented here works only for positive particleMasses.
% here we assume that we only compute self interactions, that is the sourceNode and targetNode 
%  particleCoordinates coincide.

% --- Default parameters
if nargin < 3 
    maxNumPointsPerNode = 1; 
end
if nargin < 4 
    maxNumLevels = 20;      
end
assert(all(particleMasses >= 0), 'All particle masses must be positive');

% --- Builds the quad-tree
globalIDs                       = 1 : size(particleCoordinates, 2);
tree                            = qtree;
tree.insertPoints(globalIDs, particleCoordinates, maxNumPointsPerNode, maxNumLevels);

userData.particleCoordinates    = particleCoordinates;
userData.particleMasses         = particleMasses;

% --- Collects all the leaf nodes into a single array. This will be needed to expose parallelism
leaves = tree.leaves(); 

% --- Computes center of mass and total mass for each leaf
% --- Alternatively, the following line could be used:
%     tree.postorder(@averageLeaves, [], userData);
%     However, in the case of many particles per leaf, the loop version can be parallelized
for l = 1 : length(leaves)
    leaf = leaves{l};
    averageLeaves(leaf, userData);
end

% --- Computes averages for the internal nodes
tree.postorderTraversal(@averageInternal, [], []);

% --- Evaluate potential
potential = zeros(size(particleCoordinates, 2), 1);
% --- Loops over the leaves. After having selected the leave, considers all
% the particles in that leave. They are the targetNode particles, namely, the
% particles where we want to compute the potential.
for l = 1:length(leaves)

    leaf = leaves{l};
    if ~isempty(leaf.globalIDs)
        potential(leaf.globalIDs) = evaluate(leaf, tree, userData);
    end
  
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EVALUATE POTENTIAL FUNCTION %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function potential = evaluate(leaf, tree, userData)

%   if userData.plotTree, leaf.plotnode(4); end  % optional plot
% --- If leaf has no particle, then exit
if isempty(leaf.globalIDs) 
    return; 
end

targetCoordinates   = userData.particleCoordinates(:, leaf.globalIDs);          % --- Get targetNode point coordinates
potential           = zeros(length(leaf.globalIDs), 1);                         % --- Initialize potential to zero
  
% prunning function,  a handle function to use for traversals
prune = @(node,dummy) checkIfWellSeparated(node, leaf); 

% evaluation function, again to be used in the traversals
  function visit(node, userData)                           
%     if userData.plotTree, node.plotnode(); pause(0.01); end  % optional plot
    if node.isleaf & ~isempty(node.globalIDs)  % if sourceNode box has not sources, nothing to do
      src = userData.particleCoordinates(:,node.globalIDs);       % get sourceNode positions
      den = userData.particleMasses(:,node.globalIDs);    % get sourceNode particleMasses values
      potential = potential + direct(targetCoordinates,src,den);        % direct evaluation between boxes
      return;
    end
  % if we prune, we have to add the nodes contribution to the leaf
    if prune(node, []) 
      potential = potential + direct(targetCoordinates, node.data.centerOfMass, node.data.totalMass);  % direct evaluat
    end
  end

% preorderTraversal traversal for evaluation
  tree.preorderTraversal(@visit,prune,userData);   % if couldn't prune, continue with traversal
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION TO COMPUTE THE AVERAGE MASS ON EACH LEAF %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function averageLeaves(node, userData)
 
% --- If the node is not a leaf, then exits
if ~node.isleaf
    return; 
end

particleCoordinates     = userData.particleCoordinates;
particleMasses          = userData.particleMasses;
globalIDs               = node.globalIDs;

if isempty(globalIDs) % --- Empty leaf
    centerOfMass    = zeros(2, 1); 
    totalMass       = 0; 
else
    totalMass           = sum(particleMasses(globalIDs));  % --- Total mass
    assert(totalMass > 0, 'Only positive particle masses allowed');
    centerOfMass(1)     = sum(particleCoordinates(1, globalIDs) .* particleMasses(globalIDs)) / totalMass;  % --- x-coordinate of the center of mass
    centerOfMass(2)     = sum(particleCoordinates(2, globalIDs) .* particleMasses(globalIDs)) / totalMass;  % --- y-coordinate of the center of mass
  end
  node.data.centerOfMass    = centerOfMass(:);
  node.data.totalMass       = totalMass;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION TO COMPUTE TOTAL MASS AND CENTER OF MASS FOR EACH NODE %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function averageInternal(node, userData)

if node.isleaf
    return; 
end

centerOfMass = zeros(2, 1);
totalMass = 0;
for k = 1 : 4
    totalMass       = totalMass + node.children{k}.data.totalMass;
    centerOfMass    = centerOfMass + node.children{k}.data.centerOfMass * node.children{k}.data.totalMass;
end
node.data.centerOfMass  = centerOfMass / totalMass;
node.data.totalMass     = totalMass;

% --- If we have empty leaves totalMass = 0
if ~(totalMass > 0) 
    centerOfMass = zeros(2,1); 
end  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION TO CHECK IF TWO NODES ARE WELL SEPARATED %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function itis = checkIfWellSeparated(sourceNode, targetNode)

% --- Checks whether two boxes are well separated  
  
[t_xmin, t_xmax, t_ymin, t_ymax] = targetNode.getCornerCoordinates();
[s_xmin, s_xmax, s_ymin, s_ymax] = sourceNode.getCornerCoordinates();
  
sourceWidth = sourceNode.getNodeWidth();
  
% --- Neighbor region of the source node
s_xmin = s_xmin - sourceWidth;   
s_ymin = s_ymin - sourceWidth;
s_xmax = s_xmax + sourceWidth;   
s_ymax = s_ymax + sourceWidth;

% check overlap of sourceNode box neigbhorhood region with targetNode
  flagx = t_xmax > s_xmin & t_xmin < s_xmax;
  flagy = t_ymax > s_ymin & t_ymin < s_ymax;

% if both flags  are true the boxes are not well separated
  itis = ~ (flagx & flagy);
end


%/* ************************************************** */
function potential = direct(targetCoordinates,src,particleMasses)
% this is a direct N^2 evaluation between particles  
  N = size(targetCoordinates,2);
  potential = zeros(N,1);

% MAIN LOOP OVER POINTS
  for k=1:N
  % compute distance  
    rx = targetCoordinates(1,k) - src(1,:);  
    ry = targetCoordinates(2,k) - src(2,:);
    r =  sqrt(rx.*rx + ry.*ry);

  % compute greens function (natural logarithm in 2D)
    g = -1/2/pi * log (r);

  % correct for singularities (r=0->log=inf)
    idx = find (g==inf | g==-inf);
    g(idx)=0;

    potential(k) = sum ( g.*particleMasses );
  end
  
end


%/* ************************************************** */
% function selftest()
% 
% N              = 2^6;    % number of particles
% maxNumPointsPerNode = 10;      % particleCoordinates per box
% maxNumLevels       = 20;     % maximum tree depth
% verbose        = false;  % debugging
% particleCoordinates = rand(2,N);  
% particleMasses = rand(1,N)/N;
% 
% [potential,tree] = nbody(particleCoordinates,particleMasses,maxNumPointsPerNode,maxNumLevels,false);
% uex=direct(particleCoordinates(:,1:10), particleCoordinates, particleMasses);
% error = norm(potential(1:10)-uex)/norm(uex)
% tree.plotTree();
% 
% end
% 

