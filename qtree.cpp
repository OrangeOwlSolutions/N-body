#include "qtree.h"

#include <thrust\merge.h>
#include <thrust\execution_policy.h>

/**************/
/* DESTRUCTOR */
/**************/
qtree::~qtree()
{
	//// not leaf
	//if(!isleaf)
	//	delete[] children;

	return;
}

/*********************************/
/* CONSTRUCTOR WITH NO ARGUMENTS */
/*********************************/
qtree::qtree() {

	float2 lowerLeftCorner; lowerLeftCorner.x = 0.; lowerLeftCorner.y = 0.;
	
	this->parent = NULL;
	this->level = 0;
	this->lowerLeftCorner = lowerLeftCorner;
	this->children.resize(4);
	this->children[0] = NULL;
	this->children[1] = (qtree *)NULL;
	this->children[2] = (qtree *)NULL;
	this->children[3] = (qtree *)NULL;
	//this.data = [];
	this->isleaf = true;

}

/************************/
/* INSERT POINTS METHOD */
/************************/
void qtree::insertPoints(thrust::host_vector<int> &globalIDs) {

	// --- Function adding points to the quad tree.
	// --- particleCoordinates           : 2 x N array of the point coordinates (N is the number of particles)
	// --- globalIDs : ids of particleCoordinates(unsigned int)
	// --- maxNumPointsPerNode : used to decide wether to split the node or not
	//
	//     Note that this function will break if globalIDs is not a permuation of 1 : N

	if (globalIDs.size() == 0) return;

	printf("Sto qui\n");

	// --- If the node is a leaf...
	if (this->isleaf) {
		const int numParticlesToBeInserted = globalIDs.size();
		const int numAlreadyPresentParticles = this->globalIDs.size();
		printf("Number of Particles already present %d\n", numAlreadyPresentParticles);
		// ...the particle IDs are merged...
		this->globalIDs.reserve(this->globalIDs.size() + globalIDs.size());
		this->globalIDs.insert(this->globalIDs.end(), globalIDs.begin(), globalIDs.end());
		//this.globalIDs = unique([this.globalIDs(:); globalIDs(:)]);

		// ...if the node can host the already contained particles plus the new particles, or if the node level is the maximum possible,
		//then there is nothing to do and the routine returns...
		//	if (numAlreadyPresentParticles + numParticlesToBeInserted <= maxNumPointsPerNode || this.level == maxNumLevels)
		//		return;
		//end

		//	% ...if the node cannot host the already contained particles plus the new particles, or if the node level is less than maximum
		//	%    possible, then the node is split.
		//	splitNode(this);

		//% -- - The particle IDs are temporary saved and canceled from the current node that has become a parent.
		//	globalIDs = this.globalIDs;
		//this.globalIDs = [];
	}

	//	% -- - Now we have to insert the particle coordinates to the children.
	//	for k = 1 : 4
	//		% -- - First, we have to check which points belong to which child...
	//		idx = this.children{ k }.getPointIndicesInNode(particleCoordinates(:, globalIDs));
	//% ... and then the function recursively calls itself.
	//	this.children{ k }.insertPoints(globalIDs(idx), particleCoordinates, maxNumPointsPerNode, maxNumLevels);
	//end

	//get_centroid();
	//
	//if (isleaf){

	//	if((level==maxNumLevels) || (globalIDs.size()<maxNumPointsPerNode)){
	//		isleaf=1; // becomes leaf
	//		return 1;
	//	}
	//	
	//	else create_kids();

	//}

	//if (n_th > 1){
	//	// now insert points to children
	//	// for(int i=0; i<4; i++){
	//	#pragma omp parallel shared (n_th)
	//	{

	//	#pragma omp sections 
	//	{
	//		#pragma omp section 
	//		{
	//			#pragma omp atomic
	//			n_th--;
	//			// #pragma omp critical
	//			// cout<<omp_get_thread_num()<<" start"<<endl;
	//			children[0].points_in_node(globalIDs);
	//			
	//			if(children[0].insertPoints(n_threads/2))
	//				childrenExist[0] = 0;
	//			
	//			children[1].points_in_node(globalIDs);
	//			if(children[1].insertPoints(n_threads/2))
	//				childrenExist[1] = 0;
	//		
	//			// #pragma omp critical
	//			// cout<<omp_get_thread_num()<<" end"<<endl;
	//		
	//			#pragma omp atomic
	//			n_th++;
	//		}

	//		#pragma omp section
	//		{
	//			#pragma omp atomic
	//			n_th--;
	//			// #pragma omp critical
	//			// cout<<omp_get_thread_num()<<" start"<<endl;

	//			children[2].points_in_node(globalIDs);
	//			if(children[2].insertPoints(n_threads-n_threads/2))
	//				childrenExist[2] = 0;
	//			
	//			children[3].points_in_node(globalIDs);
	//			if(children[3].insertPoints(n_threads-n_threads/2))
	//				childrenExist[3] = 0;

	//			#pragma omp atomic
	//			n_th++;
	//			
	//			// #pragma omp critical
	//			// cout<<omp_get_thread_num()<<" end"<<endl;
	//			
	//		}

	//	} 
	//}
	//}
	//// serial
	//else{
	//	for(int i=0; i<4; i++){
	//		children[i].points_in_node(globalIDs);
	//		if(children[i].insertPoints(1))
	//			childrenExist[i] = 0;
	//	}
	//		
	//}
	//
	//
	//return 0;

}


//void qtree::qtree(qtree* ini_parent, int ini_level, point ini_anchor,
//	int ini_max_level, int ini_max_pts,
//	point ini_width, int ini_lid,
//	vector<particle>* points, int num_points,
//	long ini_gid, int st_point) {
//
//	localNodeID		= ini_lid;
//	globalNodeID	= ini_gid;
//	this.level			= ini_level;   // level of the node 
//	particleCoordinates = points;
//	// globalIDs
//	totalNumberOfParticles = num_points;
//	// children
//	// parent
//	isleaf = 1;
//	lowerLeftCorner = ini_anchor;  // [x;y] coordinates of lower left point of a node
//	nodeSize = ini_width; // nodeSize of a box
//	parent = ini_parent;  // parent of the node
//	maxNumLevels = ini_max_level;
//	maxNumPointsPerNode = ini_max_pts;
//	st_pt = st_point;
//
//	for (int i = 0; i<4; i++)
//		childrenExist[i] = 1;
//}

//void qtree::initialize(qtree* ini_parent, int ini_level, point ini_anchor,
//	int ini_max_level, int ini_max_pts,
//	point ini_width, int ini_lid,
//	vector<particle>* points, int num_points,
//	long ini_gid, int st_point) {
//
//	localNodeID = ini_lid;
//	globalNodeID = ini_gid;
//	level = ini_level;   // level of the node 
//	particleCoordinates = points;
//	// globalIDs
//	totalNumberOfParticles = num_points;
//	// children
//	// parent
//	isleaf = 1;
//	lowerLeftCorner = ini_anchor;  // [x;y] coordinates of lower left point of a node
//	nodeSize = ini_width; // nodeSize of a box
//	parent = ini_parent;  // parent of the node
//	maxNumLevels = ini_max_level;
//	maxNumPointsPerNode = ini_max_pts;
//	st_pt = st_point;
//
//	for (int i = 0; i<4; i++)
//		childrenExist[i] = 1;
//}

////	qt->initialize_root( NULL, 0, anch, maxNumLevels,
////						 maxNumPointsPerNode, nodeSize, 0,
////						 particleCoordinates, N, 0, st_pt, np_proc );
//
//void qtree::initialize_root( qtree* ini_parent, int ini_level, point ini_anchor,
//							 int ini_max_level, int ini_max_pts,
//							 point ini_width, int ini_lid,
//							 vector<particle>* points, int num_points,
//							 long ini_gid, int st_point,
//							 int np_proc )
//{
//	initialize( ini_parent, ini_level, ini_anchor,
//				ini_max_level, ini_max_pts,
//				ini_width, ini_lid, points, num_points, ini_gid, st_point );
//	globalIDs.resize(np_proc,0);
//	for(int i=0; i<np_proc; i++){
//		globalIDs[i]=st_point+i;
//	}
//	localNodeID = 0;
//	globalNodeID = 0;
//};
//
//
//void qtree::create_kids()
//{
//	int kid_level = level+1;
//	point kid_width;
//	kid_width.x = nodeSize.x/2.0;
//	kid_width.y = nodeSize.y/2.0;
//	point kid_anchor;
//
//	// create new children objects
//	children  = new qtree[4];
//	
//	for(int i=0; i<4; i++){
//		kid_anchor.x = lowerLeftCorner.x + kid_width.x*kids_pos[i][0];
//		kid_anchor.y = lowerLeftCorner.y + kid_width.y*kids_pos[i][1];
//
//		// kid_global_id = 4*parent_global_id + kid_local_id
//		children[i].initialize(this, kid_level, kid_anchor,
//						   maxNumLevels, maxNumPointsPerNode, kid_width, i, particleCoordinates, totalNumberOfParticles,
//						   i+globalNodeID*4, st_pt );
//	}
//
//	// not a leaf any more
//	isleaf = 0;
//			
//	// show_kids();
//	
//	return;
//}
//
//
//void qtree::show_kids()
//{
//	// for(int i=0; i<4; i++)
//	// 	cout<<"kid "<<i<<": "<<children[i].lowerLeftCorner.x<<" "<<children[i].lowerLeftCorner.y<<endl;
//
//}
//
//// check which pointes are in nodes
//void qtree::points_in_node(vector<int>& idx_parent)
//{
//	// number of points that a parent has 
//	int npts = idx_parent.size();
//		
//	for(int i=0; i<npts; i++){
//		// int mt_id = (*particleCoordinates)[real_idx].mt_id;
//		#ifndef LONG_LONG
//		unsigned int mt_id = (*particleCoordinates)[idx_parent[i]].mt_id;
//		#endif
//		#ifdef LONG_LONG
//		unsigned long long mt_id = (*particleCoordinates)[idx_parent[i]].mt_id;
//		#endif
//		
//		int offset = maxNumLevels-level;		
//		
//		if(( (mt_id & mt_checker[offset]) >> (offset*2)  ) == localNodeID){
//			globalIDs.push_back(idx_parent[i]);
//				// cout<<(*particleCoordinates)[i].mt_id<<" is in box "<<localNodeID
//				// <<" on level "<<level<<endl;
//		}
//	}
//	
//	
//}
//
//// show tree using preorder traversal
//void qtree::show_tree()
//{
//	// cout<<"level: "<<level<<endl
//	// 	// <<"localNodeID:   "<<localNodeID<<endl
//	// 	<<"globalNodeID:   "<<globalNodeID<<endl
//		// <<"n_pts: "<<globalIDs.size()<<endl<<endl;
//	for(int i=0; i<4; i++){
//		if(!children[i].isleaf){
//			// cout<<"not leaf"<<endl;
//			children[i].show_tree();
//		}
//		else{
//			cout<<"level: "<<children[i].level
//				<<" id: "<<children[i].globalNodeID<<endl;
//			for(int j=0; j<children[i].globalIDs.size(); j++)
//				cout<<children[i].globalIDs[j]<<endl;
//		}
//		// cout<<endl;
//	}
//
//	
//}
//
//void qtree::get_centroid( )
//{
//	centroid.x=0;
//	centroid.y=0;
//	totalMass=0;
//	
//	for(int i=0; i<globalIDs.size(); i++){
//		centroid.x += (*particleCoordinates)[globalIDs[i]].x * (*particleCoordinates)[globalIDs[i]].m;
//		centroid.y += (*particleCoordinates)[globalIDs[i]].y * (*particleCoordinates)[globalIDs[i]].m;
//		totalMass += (*particleCoordinates)[globalIDs[i]].m;
//	}
//
//	if(totalMass>0){
//		centroid.x /= totalMass;
//		centroid.y /= totalMass;
//	}
//	// need to think about the case where there is no points
//	else{
//		centroid.x = 0;
//		centroid.y = 0;
//	}
//	
//	// cout<<centroid.x<<" "<<centroid.y<<endl;
//	
//	return;
//}
