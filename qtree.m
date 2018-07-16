classdef qtree < handle
% --- Implements a quatree data structure for point-like particles
% --- Preorder traversal and postorder traversal implemented.
%      
% --- See "doc qtree" for documentation on methods and classes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QUAD-TREE CLASS PROPERTIES %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
properties
    globalIDs           % --- Global particle IDs 
    children            % --- children{i}, i = 1 : 4, are the children of the current node
    parent              % --- Node parent
    level               % --- Node level
    data                % --- Averaging-related data
    isleaf              % --- True if node is a leaf or false otherwise
    lowerLeftCorner     % --- Coordinates of lower left point of a node
end

methods

%%%%%%%%%%%%%%%%%%%%%%%%%
% QUAD-TREE CONSTRUCTOR %
%%%%%%%%%%%%%%%%%%%%%%%%%
function this = qtree(parent, level, lowerLeftCorner)
% function this = qtree(parent, level, lowerLeftCorner)

% --- Quad-tree constructor 

% --- Defaults parameters
%     parent          = []
%     level           = 0
%     lowerLeftCorner = [0; 0]
if nargin < 1 
    parent = [];           
end
if nargin < 2 
    level  = 0;            
end
if nargin < 3 
    lowerLeftCorner = [0, 0];        
end

this.globalIDs          = [];
this.children           = [];
this.parent             = parent;
this.level              = level;
this.data               = [];
this.isleaf             = true;
this.lowerLeftCorner    = lowerLeftCorner;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET THE WIDTH OF THE NODE METHOD %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = getNodeWidth(this)
% function h = getNodeWidth(this)
%
% --- Gets the width of the square that corresponds to this node
  h = (1 / 2)^this.level;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET CORNER COORDINATES OF A NODE METHOD %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xmin, xmax, ymin, ymax] = getCornerCoordinates(this)
% function [xmin, xmax, ymin, ymax] = getCornerCoordinates(this)  
%
% --- Gets the coordinates of lower-left and upper-right corners
% --- Lower-left point: [xmin; ymin], upper-left point: [xmin; ymax] 
  
h       = getNodeWidth(this);
xmin    = this.lowerLeftCorner(1);
ymin    = this.lowerLeftCorner(2);

xmax    = xmin + h;
ymax    = ymin + h;

% --- To capture particles on the domain boundaries, enlarge xmax and ymax of a tiny quantity.
if xmax == 1 
    xmax = 1 + eps; 
end 
if ymax == 1 
    ymax = 1 + eps; 
end
end
  
%%%%%%%%%%%%%%%%%%%%%
% SPLIT NODE METHOD %
%%%%%%%%%%%%%%%%%%%%%
function splitNode(this)  
% function splitNode(this)    
% 
% --- Creates the four children of a parent node

assert(isempty(this.children), 'Node already split.');

childWidth    = this.getNodeWidth / 2;

% --- The four children are created.
% --- Lower-left child
this.children{1} = qtree(this, this.level + 1, this.lowerLeftCorner                     );
% --- Upper-left child
this.children{2} = qtree(this, this.level + 1, this.lowerLeftCorner + childWidth * [1 0]);
% --- Lower-right child
this.children{3} = qtree(this, this.level + 1, this.lowerLeftCorner + childWidth * [0 1]);
% --- Upper-right child
this.children{4} = qtree(this, this.level + 1, this.lowerLeftCorner + childWidth * [1 1]);

% --- The parent is not anymore a leaf.
this.isleaf = false;
end

%%%%%%%%%%%%%%%%%%%%%%%%
% INSERT POINTS METHOD %
%%%%%%%%%%%%%%%%%%%%%%%%
function insertPoints(this, globalIDs, particleCoordinates, maxNumPointsPerNode, maxNumLevels)
% function insertPoints(this, globalIDs, particleCoordinates, maxNumPointsPerNode, maxNumLevels)
%
% --- Function adding points to the quad tree.
% --- particleCoordinates           : 2 x N array of the point coordinates (N is the number of particles)
% --- globalIDs                     : ids of particleCoordinates (unsigned int)
% --- maxNumPointsPerNode           : used to decide wether to split the node or not  
%
%     Note that this function will break if globalIDs is not a permuation of 1 : N
  
if isempty(globalIDs), return; end;

% --- In the node is a leaf...
if this.isleaf
    numParticlesToBeInserted      = length(globalIDs);
    numAlreadyPresentParticles    = length(this.globalIDs);
    % ...the particle IDs are merged...
    this.globalIDs                = unique([this.globalIDs(:); globalIDs(:)]);

    % ...if the node can host the already contained particles plus the new particles, or if the node level is the maximum possible, 
    %    then there is nothing to do and the routine returns...
    if (numAlreadyPresentParticles + numParticlesToBeInserted <= maxNumPointsPerNode || this.level == maxNumLevels) 
        return;  
    end

    % ...if the node cannot host the already contained particles plus the new particles, or if the node level is less than maximum 
    %    possible, then the node is split.
    splitNode(this);
        
    % --- The particle IDs are temporary saved and canceled from the current node that has become a parent.
    globalIDs       = this.globalIDs;
    this.globalIDs  = [];
end

% --- Now we have to insert the particle coordinates to the children.
for k = 1 : 4
    % --- First, we have to check which points belong to which child...
    idx = this.children{k}.getPointIndicesInNode(particleCoordinates(:, globalIDs));
    % ... and then the function recursively calls itself.
    this.children{k}.insertPoints(globalIDs(idx), particleCoordinates, maxNumPointsPerNode,maxNumLevels  );
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET INDICES OF THE PARTICLES BELONGING TO THE CURRENT NODE METHOD %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function idx = getPointIndicesInNode(this, particleCoordinates)
% function idx = getPointIndicesInNode(this, particleCoordinates)  
%
% --- Gets the indices of the particles belonging to the current node

[xmin, xmax, ymin, ymax] = getCornerCoordinates(this);

idx_x =  find(xmin <= particleCoordinates(1, :) & particleCoordinates(1,:) < xmax);
idx_y =  find(ymin <= particleCoordinates(2, :) & particleCoordinates(2,:) < ymax);
  
idx = intersect(idx_x, idx_y);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREORDER TRAVERSAL METHOD %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function preorderTraversal(this, visitFunction, prune, userData)
% function preorderTraversal(this, visitFunction, prune, userData)  
%
% --- Performs preorder traversal with prunning  
%     visitFunction represents the operations that must be performed on each not upon traversal
% prune(node,data) can be empty
% userData: user data structure for calculations

doPrune = false;
if ~isempty(prune)
    doPrune = prune(this, userData);  
end
visitFunction(this, userData);
% --- If the node is not a leaf and pruning must not be performed, then
% visit the four children
if ~this.isleaf && ~doPrune
   for k = 1 : 4
     this.children{k}.preorderTraversal(visitFunction, prune, userData);
   end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POSTORDER TRAVERSAL METHOD %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function postorderTraversal(this, visit, prune, userData)
% function  postorder(this, visit, prune, userData)  
% preorderTraversal traversal with prunning  
% visit(node,data)
% prune(node,data) can be empty
% userData: user data structure for calculations

doPrune = false;
if ~isempty(prune),   doPrune= prune(this, userData);  end
if ~this.isleaf & ~doPrune
   for k=1:4
     this.children{k}.postorderTraversal( visit, prune, userData);
   end
end
visit(this,userData);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LEAVES COLLECTION METHOD %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function leavesArray = leaves(this)
% function leavesArray = leaves(this)  
%
% --- Collects all the leaves in a single array  

leavesArray     = {};
counter         = 0;
function visit(this, dummyParameter)
    if this.isleaf
        counter = counter + 1;
        leavesArray{counter} = this;
    end
end
this.preorderTraversal(@visit, [], []);
end
  

%/* ************************************************** */
function print(this)
% function print(this)  
% print information about a node  
  printNode = @(this,dummyParameter)...
      fprintf('node  at level %d: lowerLeftCorner:[%1.4f %1.4f]\N',...
              this.level,this.lowerLeftCorner(1),this.lowerLeftCorner(2));
  this.preorderTraversal(printNode,[],[]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRINT NODE MORTON IDs METHOD %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function printMortonIDs(this, leavesOnly)
% function printMortonIDs(this, leavesOnly)  
%
% --- Prints the Morton IDs of the nodes (using the class mortonID)
%     If leavesOnly is true, it only prints the leaves.

mID = mortonID;
if nargin < 2
    leavesOnly = false; 
end
function printNode(this, dummyParameter)
    if leavesOnly && ~this.isleaf
        return; 
    end
    fprintf('Morton ID:');
    ID = mID.id(this.level, this.lowerLeftCorner);
    mID.print(ID);
    fprintf(': %20u at level %2d: lowerLeftCorner: [%1.4f %1.4f]\n', ID, this.level, this.lowerLeftCorner(1), this.lowerLeftCorner(2));
end
this.preorderTraversal(@printNode, [], []);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIND TREE DEPTH METHOD %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function depth = findDepth(this)
% function depth = findDepth(this)  
%
% --- Returns the tree depth
   depth = 0;
   function v = fdt(this, dummyParameter)
      depth = max(this.level, depth);
   end
   this.preorderTraversal(@fdt, [], []);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING A SINGLE NODE FUNCTION %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotNode(this, linewidth) 
% function plotNode(this, width)   
%
% --- Plots the rectangle corresponding to the node  
% --- Sets the linewidth to 2 by default
if nargin < 2
    linewidth = 2;
end
rectangle('Curvature', [0.00, 0.00], ...
          'Position', [this.lowerLeftCorner(1), this.lowerLeftCorner(2), 1/2^this.level, 1/2^this.level], ...
          'LineWidth', linewidth);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING THE ENTIRE TREE METHOD %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotTree(this, linewidth)
% function plotTree(this, linewidth)  
%
% --- Plots the entire tree.  
hold on;
% --- Sets the linewidth to 2 by default
if nargin < 2
    linewidth = 2;
end
% --- Defines an inline function plotting a node
plotNode = @(node, dummyParameter)node.plotNode(linewidth);
this.preorderTraversal(plotNode, [], []);
hold off;
end

end % methods
end % classdef
