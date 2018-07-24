/******************************/
/* QUAD-TREE CLASS DEFINITION */
/******************************/
#ifndef QTREE_H
#define QTREE_H

//#include <iostream>
#include <thrust\host_vector.h>
//#include <cstdlib>
//#include <omp.h>
//#include "point.h"
//#include "particle.h"
//#include "bits_numbers.h"

using namespace std;

#define LONG_LONG


// keep track of number of threads at work
extern int n_th;

/*******************/
/* QUAD-TREE CLASS */
/*******************/
class qtree {

public:

	int localNodeID;											// --- Local node ID [0, 1, 2, 3]
	long globalNodeID;											// --- Global node ID
  	int level;													// --- Node level 
	
	thrust::host_vector<float2> *particleCoordinates;			// --- 2 x N array of the point coordinates (N is the number of particles)

	thrust::host_vector<int> globalIDs;							// --- Global particle IDs

	// starting point in the array for parallel
	int st_pt;
	
	// --- Pointers to relatives
	thrust::host_vector<qtree *> children;						// --- children[i], i = 1 : 4, are the children of the current node - It is a 4-size array
	qtree *parent;												// --- Node parent

	// --- Dimensional data
	float2 lowerLeftCorner;										// --- Coordinates of lower left point of a node
	float2 nodeSize;											// --- x and y dimensions of a box

	bool isleaf;												// --- True if node is a leaf or false otherwise
	bool childrenExist[4];										// --- True if corresponding child exists

	// --- Centroid data
	double totalMass;											// --- Node total mass
	float2 centroid;											// --- Centroid coordinates
	
	int totalNumberOfParticles;									// --- Total number of particles in node

	int maxNumLevels;											// --- Maximum tree depth
	int maxNumPointsPerNode;									// --- Maximum number of particles per node

	// --- Constructor with no arguments
	qtree();

	// --- Destructor
	~qtree();
	
	//void initialize ( qtree* ini_parent, int ini_level, point ini_anchor,
	//				  int ini_max_level, int ini_max_pts,
	//				  point ini_width, int ini_lid,
	//				  vector<particle>* points, int num_points,
	//				  long ini_gid, int st_point );
	//void initialize_root( qtree* ini_parent, int ini_level, point ini_anchor,
	//					  int ini_max_level, int ini_max_pts,
	//					  point ini_width, int ini_lid,
	//					  vector<particle>* points, int num_points,
	//					  long ini_gid, int st_point,
	//					  int np_local);
	//
	
	void insertPoints(thrust::host_vector<int> &globalIDs);
	//void create_kids();
	//void show_kids();
	//void points_in_node( vector<int>& idx_parent );
	//void show_tree();
	//void get_centroid( );

};

#endif 
