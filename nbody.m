function [potential, tree] = nbody(particleCoordinates, particleCharges, maxNumPointsPerNode, maxNumLevels)
% function [potential, tree] = nbody(particleCoordinates, particleCharges, maxNumPointsPerNode, maxNumLevels, plotTree)
%
% --- particleCoordinates           : 2 x N array of the point coordinates (N is the number of particles)
% --- particleCharges                : particle masses - all particle masses must be positive
% --- maxNumPointsPerNode           : used to decide wether to split the node or not  
%
% --- potential                     : computed potential at each point
% --- tree                          : pointer to the root of the tree  
%  
% --- Assumptions:
%     the algorithm applies to arbitrary particleCharges, but the simple
%         averaging implemented here works only for positive particleCharges.
% here we assume that we only compute self interactions, that is the source and target 
%  particleCoordinates coincide.

% --- Default parameters
if nargin < 3 
    maxNumPointsPerNode = 1; 
end
if nargin < 4 
    maxNumLevels = 20;      
end
assert(all(particleCharges >= 0), 'All particle masses must be positive');

% BUILD TREE ---------------------------------------
globalIds = 1:size(particleCoordinates,2);
root = qtree;
root.insertPoints(globalIds,particleCoordinates,maxNumPointsPerNode,maxNumLevels);

user_data.particleCoordinates    = particleCoordinates;
user_data.particleCharges = particleCharges;
% user_data.plotTree  = plotTree; % to plot tree as you evaluate()

% collect all the leaves on a single lits
% used in averaging and evaluation to expose parallelism  
leaves = root.leaves(); 

% AVERAGE  ---------------------------------------
% we could use the commented line below, but for many particleCoordinates per box
% the loop version can be parallelized for the leaf nodes
%root.postorder(@average_leaves,[], user_data);
for l=1:length(leaves)
  leaf = leaves{l};
  average_leaves(leaf,user_data);
end
% average internal nodes
root.postorder(@average_internal,[],[]);

% EVALUATE  ---------------------------------------
u = zeros(size(particleCoordinates,2),1);
for l = 1:length(leaves)
  if(plotTree) pause(0.05); clf; hold on; end
  leaf = leaves{l};
  if ~isempty(leaf.globalIds)
    u(leaf.globalIds) = evaluate( leaf, root, user_data);
  end
end

potential = u;
tree      = root;
end


%/* ************************************************** */
function u = evaluate(leaf,root,ud)
% the main evaluation routine for the NBODY   

  if ud.plotTree, leaf.plotnode(4); end  % optional plot
  if isempty(leaf.globalIds) return; end;     % if leaf has no particleCoordinates nothing to do
  trg = ud.particleCoordinates(:,leaf.globalIds);          % get target point coordinates
  u   = zeros(length(leaf.globalIds),1);      % initialize potential to zero
  
% prunning function,  a handle function to use for traversals
  prune=@(node,dummy) well_separated(node,leaf); 

% evaluation function, again to be used in the traversals
  function visit(node, ud)                           
    if ud.plotTree, node.plotnode(); pause(0.01); end  % optional plot
    if node.isleaf & ~isempty(node.globalIds)  % if source box has not sources, nothing to do
      src = ud.particleCoordinates(:,node.globalIds);       % get source positions
      den = ud.particleCharges(:,node.globalIds);    % get source particleCharges values
      u = u + direct(trg,src,den);        % direct evaluation between boxes
      return;
    end
  % if we prune, we have to add the nodes contribution to the leaf
    if prune(node,[]) 
      u = u + direct(trg, node.data.x_a, node.data.d_a);  % direct evaluat
    end
  end

% preorderTraversal traversal for evaluation
  root.preorderTraversal(@visit,prune,ud);   % if couldn't prune, continue with traversal
end


%/* ************************************************** */
function average_leaves(node,user_data)
  if ~node.isleaf, return; end;
  particleCoordinates    = user_data.particleCoordinates;
  particleCharges = user_data.particleCharges;
  globalIds = node.globalIds;

  if isempty(globalIds) % empty leaf
    x_a=zeros(2,1); d_a=0; 
  else
    d_a = sum( particleCharges(globalIds) );  % average 
    assert(d_a>0,'Only positive particleCharges allowed');
    x_a(1) = sum( particleCoordinates   (1,globalIds).*particleCharges(globalIds) )/d_a;  % average position
    x_a(2) = sum( particleCoordinates   (2,globalIds).*particleCharges(globalIds) )/d_a;  % average position
  end
  node.data.x_a = x_a(:);
  node.data.d_a = d_a;
end

%/* ************************************************** */
function average_internal(node,user_data)
if node.isleaf, return; end;
  x_a = zeros(2,1); 
  d_a = 0;
  for k=1:4
    d_a = d_a + node.children{k}.data.d_a;
    x_a = x_a + node.children{k}.data.x_a * node.children{k}.data.d_a;
  end
  node.data.x_a = x_a/d_a;
  node.data.d_a = d_a;
  if ~(d_a > 0) x_a = zeros(2,1); end  % if we have empty leaves d_a=0
end


%/* ************************************************** */
function itis=well_separated(source,target)
% check whether two boxes are well separated  
  [t_xmin,t_xmax,t_ymin,t_ymax] = target.getCornerCoordinates();
  [s_xmin,s_xmax,s_ymin,s_ymax] = source.getCornerCoordinates();
  h = source.getNodeWidth();
  
% neighbor region of target
  s_xmin = s_xmin -h;   s_ymin = s_ymin -h;
  s_xmax = s_xmax +h;   s_ymax = s_ymax +h;

% check overlap of source box neigbhorhood region with target
  flagx = t_xmax > s_xmin & t_xmin < s_xmax;
  flagy = t_ymax > s_ymin & t_ymin < s_ymax;

% if both flags  are true the boxes are not well separated
  itis = ~ (flagx & flagy);
end


%/* ************************************************** */
function u = direct(trg,src,particleCharges)
% this is a direct N^2 evaluation between particles  
  N = size(trg,2);
  u = zeros(N,1);

% MAIN LOOP OVER POINTS
  for k=1:N
  % compute distance  
    rx = trg(1,k) - src(1,:);  
    ry = trg(2,k) - src(2,:);
    r =  sqrt(rx.*rx + ry.*ry);

  % compute greens function (natural logarithm in 2D)
    g = -1/2/pi * log (r);

  % correct for singularities (r=0->log=inf)
    idx = find (g==inf | g==-inf);
    g(idx)=0;

    u(k) = sum ( g.*particleCharges );
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
% particleCharges = rand(1,N)/N;
% 
% [u,tree] = nbody(particleCoordinates,particleCharges,maxNumPointsPerNode,maxNumLevels,false);
% uex=direct(particleCoordinates(:,1:10), particleCoordinates, particleCharges);
% error = norm(u(1:10)-uex)/norm(uex)
% tree.plotTree();
% 
% end
% 

