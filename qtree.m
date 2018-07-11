classdef qtree < handle
% implements a quatree data structure
%    point based construction in 2D.
%    preorder and postorder traversals
%      
% see "doc qtree" for documentation on methods and classes
% see "nbody.m" for the 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QUAD-TREE CLASS PROPERTIES %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
properties
    globalIds           % --- Global particle IDs 
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

this.globalIds          = [];
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
function insertPoints(this, globalIds, particleCoordinates, maxNumPointsPerNode, maxNumLevels)
% function insertPoints(this, globalIds, particleCoordinates, maxNumPointsPerNode, maxNumLevels)
%
% --- Function adding points to the quad tree.
% --- particleCoordinates           : 2 x N array of the point coordinates (N is the number of particles)
% --- globalIds                     : ids of particleCoordinates (unsigned int)
% --- maxNumPointsPerNode           : used to decide wether to split the node or not  
%
%     Note that this function will break if globalIds is not a permuation of 1 : N
  
if isempty(globalIds), return; end;

% --- In the node is a leaf...
if this.isleaf
    numParticlesToBeInserted      = length(globalIds);
    numAlreadyPresentParticles    = length(this.globalIds);
    % ...the particle IDs are merged...
    this.globalIds                = unique([this.globalIds(:); globalIds(:)]);

    % ...if the node can host the already contained particles plus the new particles, or if the node level is the maximum possible, 
    %    then there is nothing to do and the routine returns...
    if (numAlreadyPresentParticles + numParticlesToBeInserted <= maxNumPointsPerNode || this.level == maxNumLevels) 
        return;  
    end

    % ...if the node cannot host the already contained particles plus the new particles, or if the node level is less than maximum 
    %    possible, then the node is split.
    splitNode(this);
        
    % --- The particle IDs are temporary saved and canceled from the current node that has become a parent.
    globalIds       = this.globalIds;
    this.globalIds  = [];
end

% --- Now we have to insert the particle coordinates to the children.
for k = 1 : 4
    % --- First, we have to check which points belong to which child...
    idx = this.children{k}.getPointIndicesInNode(particleCoordinates(:, globalIds));
    % ... and then the function recursively calls itself.
    this.children{k}.insertPoints(globalIds(idx), particleCoordinates, maxNumPointsPerNode,maxNumLevels  );
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

%/* ************************************************** */
function  preorder(this, visit, prune, user_data)
% function  preorder(this, visit, prune, user_data)  
% preorder traversal with prunning  
% visit(node,data)
% prune(node,data) can be empty
% user_data: user data structure for calculations

doPrune = false;
if ~isempty(prune),   doPrune= prune(this, user_data);  end
visit(this,user_data);
if ~this.isleaf & ~doPrune
   for k=1:4
     this.children{k}.preorder( visit, prune, user_data);
   end
end
end


%/* ************************************************** */
function  postorder(this, visit, prune, user_data)
% function  postorder(this, visit, prune, user_data)  
% preorder traversal with prunning  
% visit(node,data)
% prune(node,data) can be empty
% user_data: user data structure for calculations

doPrune = false;
if ~isempty(prune),   doPrune= prune(this, user_data);  end
if ~this.isleaf & ~doPrune
   for k=1:4
     this.children{k}.postorder( visit, prune, user_data);
   end
end
visit(this,user_data);
end

%/* ************************************************** */
function list=leaves(this)
% function list=leaves(this)  
% collects all the leaves in a single array  
  list={};
  cnt=0;
  function visit(this,dummy)
    if this.isleaf
      cnt=cnt+1;
      list{cnt}=this;
    end
  end
  this.preorder(@visit,[],[]);
end
  

%/* ************************************************** */
function print(this)
% function print(this)  
% print information about a node  
  print_node = @(this,dummy)...
      fprintf('node  at level %d: lowerLeftCorner:[%1.4f %1.4f]\N',...
              this.level,this.lowerLeftCorner(1),this.lowerLeftCorner(2));
  this.preorder(print_node,[],[]);
end

%/* ************************************************** */
function print_mids(this,leaves_only)
% function print_mids(this,leaves_only)  
% print the Morton IDs of nodes (using the class morton_id)
% leaves_only : if true, it only prints the leaves.
  mid = morton_id;
  if nargin<2, leaves_only=false; end;
  function print_node(this,dummy)
    if leaves_only & ~this.isleaf, return; end;
    fprintf('mid:');
    id = mid.id(this.level,this.lowerLeftCorner);
    mid.print(id);
    fprintf(': %20u at level %2d: lowerLeftCorner:[%1.4f %1.4f]\N',...
            id, this.level,this.lowerLeftCorner(1),this.lowerLeftCorner(2));

  end
  this.preorder(@print_node,[],[]);    
end

%/* ************************************************** */
function depth=find_depth(this)
% function depth=find_depth(this)  
% function depth=find_depth(this)  
% depth : depth of the tree.  
   depth = 0;
   function v=fdt(this,dummy)
      depth=max(this.level,depth);
   end
   this.preorder(@fdt,[],[]);
end

%/* ************************************************** */
function  plotnode(this,width) 
% function  plotnode(this,width)   
% plot the rectangul corresponding to the node  
  if nargin<2,width=1;end;
  rectangle('Curvature',[0.00,0.00],...
     'Position', [this.lowerLeftCorner(1),this.lowerLeftCorner(2),1/2^this.level,1/2^this.level], 'LineWidth',width);
end

%/* ************************************************** */
function plottree(this,markersize)
% function plottree(this,markersize)  
%  function plottree(this,markersize)
%  plots the whole tree.  
  hold on;
  if nargin<2,markersize=2;end;
  plotnode=@(node,dummy)node.plotnode(markersize);
  this.preorder(plotnode,[],[]);
  hold off;
end

end % methods
end % classdef
