#include <iostream>

#include <thrust\host_vector.h>
#include <thrust\sequence.h>
#include <thrust\generate.h>
#include <thrust\for_each.h>
#include <thrust\execution_policy.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/constant_iterator.h>

#include "qtree.h"

/************************************/
/* RANDOM NUMBERS GENERATION STRUCT */
/************************************/
template<typename T>
struct rand_01 {

	T constant;

	rand_01(T _constant) : constant(_constant) {}
	__host__ void operator()(T& VecElem) const { VecElem = (T)rand() / RAND_MAX / constant; }
};

//#include <omp.h>
//#include <iostream>
//#include <fstream>
//#include <algorithm>
//#include <vector>
//#include <stack>
//#include <cstdlib>
//#include <cmath>
//#include <bitset>
//#include <iomanip>
//
//// #include "../msort/msort.h"
//#include "point.h"
//#include "particle.h"
//#include "qtree.h"
//
//
//using namespace std;
//
//// keep track of the number of threads at work
//int n_th;
//
//int compar_vector( const particle left, const particle right)
//{
//	return left.mt_id<right.mt_id;
//}
//
//// 2-norm
//double norm2( const double x1, const double y1,
//			  const double x2, const double y2)
//{
//	return sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
//}
//
//
//// kernel
//double g_kernel( const double x1, const double y1,
//				 const double x2, const double y2) 
//{
//	if ((x1!=x2) && (y1!=y2))
//		return -1/(2*pi) * log( norm2(x1, y1, x2, y2) );
//	
//	return 0;
//}
//
//void build_tree( const double xmin,
//				 const double ymin,
//				 const int maxNumPointsPerNode,
//				 const int maxNumLevels,
//				 point nodeSize,
//				 vector<particle>* particleCoordinates,
//				 int N,
//				 qtree* qt,
//				 const int st_pt,
//				 const int np_proc,
//				 const int nt )
//{
//	point anch(xmin, ymin);
//	qt->initialize_root( NULL, 0, anch, maxNumLevels,
//						 maxNumPointsPerNode, nodeSize, 0,
//						 particleCoordinates, N, 0, st_pt, np_proc );
//
//	qt->insertPoints(nt);	
//		
//	return;
//}
//	
//
//// read from the file
//int read_points( vector<particle>& particleCoordinates,
//				 const int N )
//{
//	ifstream file_in ("points.dat");
//	if (!file_in.is_open()) return 1;
//
//	particleCoordinates.resize(N);
//	
//	for(int i=0; i<N; i++){
//		file_in>>particleCoordinates[i].x;
//		file_in>>particleCoordinates[i].y;
//		file_in>>particleCoordinates[i].mt_id;
//	}
//	
//	file_in.close();
//	
//	return 0;
//}
//
//
//int write_points( const vector<particle>& particleCoordinates,
//				  const int N )
//{
//	ofstream ofile;
//	ofile.open ("points.dat");
//
//	for (int i=0; i<N; i++){
//
//		ofile<<particleCoordinates[i].x<<" "<<particleCoordinates[i].y<<" "<<setprecision(25)<<particleCoordinates[i].mt_id<<" "
//			 <<particleCoordinates[i].m<<" "<<particleCoordinates[i].u<<endl;
//	}
//	ofile.close();
//
//	return 0;
//}
//
//int write_boxes( const qtree* qt, ofstream& ofile, int parent )
//{
//	// only make ofstream instance for parent node
//	if(parent){
//		ofile.open ("boxes.dat");
//	}
//
//	// qt is not leaf
//	if (!qt->isleaf){
//		for(int i=0; i<4; i++){
//			write_boxes(&(qt->children[i]), ofile, 0);
//		}
//	}
//	// qt is leaf
//	else{
//		ofile<<qt->level<<" "
//			 <<qt->lowerLeftCorner.x<<" "<<qt->lowerLeftCorner.y<<" "
//			 <<qt->nodeSize.x<<" "<<qt->nodeSize.y<<endl;
//	}
//
//	// only close ofstream for parent node
//	if(parent){
//		ofile.close();
//	}
//
//	return 0;
//}
//
//// not necessary?
//void average_trees( qtree* qts, const int nt )
//{
//		
//	return;
//}
//
//// get which kid has the node with id
//qtree* get_local_node(const int id, qtree* qt, const int level )
//{
//	if ((level==0) || (qt->isleaf))
//		return qt;
//
//	// qtree qt;
//	// qt.initialize_root( NULL, 0, anch, maxNumLevels,
//	// 					maxNumPointsPerNode, nodeSize, 0,
//	// 					particleCoordinates, N, 0, st_pt, np_proc );
//
//	// for(int i=0; i<nt; i++){
//	// 	for(int j=0; j<4; j++){
//	// 		if
//	// 	}
//
//	// }
//	
//	for(int i=0; i<4; i++){
//		for(int j=0; j<qt->children[i].globalIDs.size(); j++){
//			if (id==qt->children[i].globalIDs[j]){
//				return &(qt->children[i]);
//				// return get_local_node(id, &(qt->children[i]), level-1);
//			}
//		}
//	}
//	
//	return NULL;
//}
//
//// get which kid has the node with id
//qtree* get_local_node(const int id, qtree* qt, const vector<particle>& particleCoordinates )
//{
//	if (qt->isleaf)
//		return qt;
//	
//	for(int i=0; i<4; i++){
//		if( particleCoordinates[id].x > qt->children[i].lowerLeftCorner.x &&
//			particleCoordinates[id].x < (qt->children[i].lowerLeftCorner.x+qt->children[i].nodeSize.x) &&
//			particleCoordinates[id].y > qt->children[i].lowerLeftCorner.y &&
//			particleCoordinates[id].y < (qt->children[i].lowerLeftCorner.y+qt->children[i].nodeSize.y) )
//
//			return &(qt->children[i]);
//			
//	}
//	
//	return NULL;
//}
//
//
//// check if the target and source are well-separated
//int well_separated( const double x_target,
//					const double y_target,
//					qtree* source )
//{
//	double xdiff = abs((source->lowerLeftCorner.x) - (x_target));
//	double ydiff = abs((source->lowerLeftCorner.y) - (y_target));
//
//	if ((xdiff>source->nodeSize.x) || (ydiff>source->nodeSize.y) )
//		return 1;
//	
//	return 0;
//}
//
//// evaluate the potential at target using individual particles
//double direct_evaluation( const int id, qtree* source,
//						const vector<particle>& particleCoordinates)
//{
//	double u=0.0;
//	
//	for(int i=0; i<source->globalIDs.size(); i++){
//		u += g_kernel( particleCoordinates[id].x, particleCoordinates[id].y,
//					   particleCoordinates[source->globalIDs[i]].x, particleCoordinates[source->globalIDs[i]].y )
//			* particleCoordinates[source->globalIDs[i]].m;
//	}
//		
//
//	return u;
//}
//
//// evaluate the potential at a target using a box
//double approximate_evaluation( const int id, qtree* source,
//							   const vector<particle>& particleCoordinates )
//{
//	double u = g_kernel( particleCoordinates[id].x , particleCoordinates[id].y,
//						 source->centroid.x, source->centroid.y )
//		* source->totalMass;
//
//	return u;
//}
//
//// evaluate the forces using quadtree
//void evaluate_trees( const int id,
//					 qtree* source, vector<particle>& particleCoordinates)
//{
//	// get qtree pointer of the kid that has the particleCoordinates[id]
//	// target = get_local_node(id, target, 1);
//	// target = get_local_node(id, target, particleCoordinates);
//
//	
//	// cout<<"x "<<particleCoordinates[id].x<<" y "<<particleCoordinates[id].y<<endl;
//	// cout<<"globalNodeID "<<target->globalNodeID<<" level "<<target->level<<endl;
//
//	// now compare the target with the children of the source
//	for(int i=0; i<4; i++){
//		// not separated enough
//		if(!well_separated(particleCoordinates[id].x, particleCoordinates[id].y, &(source->children[i]))){
//			// but the kid is a leaf
//			if(source->children[i].isleaf){
//				// cout<<target->globalNodeID<<" and "<<source->children[i].globalNodeID
//					// <<" need direct evaluation on level "
//					// <<level+1<<endl;	// cout<<"leaf"<<endl;
//				particleCoordinates[id].u += direct_evaluation(id, &(source->children[i]), particleCoordinates);
//			}
//			// not leaf
//			else{
//				// cout<<"not separated enough"<<endl;
//				// go one level down
//				evaluate_trees(id, &(source->children[i]), particleCoordinates);
//			}
//		}
//		else{
//			// cout<<target->globalNodeID<<" and "<<source->children[i].globalNodeID
//				// <<" are sufficiently SEPARATED on level "
//				// <<level+1<<" "<<target->level<<endl;
//			particleCoordinates[id].u += approximate_evaluation(id, &(source->children[i]), particleCoordinates );
//		}
//	}
//		
//	
//	return;
//}
//
/********/
/* MAIN */
/********/
int main() {

	const int N						= 2 << 5;			// --- Number of particles
	const int maxNumPointsPerNode	= 1;				// --- Maximum number of particles per node
	const int maxNumLevels			= 20;				// --- Maximum tree depth

	// --- Particle coordinates
	thrust::host_vector<double>		particleCoordinates(N * 2);
	thrust::for_each(particleCoordinates.begin(), particleCoordinates.end(), rand_01<double>(1.));

	thrust::host_vector<double>		particleMasses(N);
	thrust::for_each(particleMasses.begin(), particleMasses.end(), rand_01<double>((double)N));

	// --- Running the algorithm
	thrust::host_vector<int>		globalIDs(N);
	//thrust::transform(thrust::make_counting_iterator(0),
	//	thrust::make_counting_iterator(N),
	//	thrust::make_constant_iterator(1),
	//	globalIDs.begin(),
	//	thrust::plus<int>());

	thrust::sequence(thrust::host, globalIDs.begin(), globalIDs.end());
	for (int k = 0; k < N; k++) printf("globalIDs = %d\n", globalIDs[k]);

	qtree quadTreeObject;				// --- Initializing the quad tree object
	
	quadTreeObject.insertPoints(globalIDs);
	
	printf("level = %d\n", quadTreeObject.level);
	
	//quadTreeObject = qtree();  

	
	//globalIDs = 1 : N;  % -- - Global particle IDs



		//quadTreeObject = qtree;  % -- - Initializing the quad tree object


	//std::cout << "Values generated: " << std::endl;
	//for (int k = 0; k < N; k++)
	//	std::cout << particleMasses[k] << " : ";
	//std::cout << std::endl;

	//const point nodeSize(xrange, yrange);

	// thread nesting enabled
	//omp_set_nested(1);
	
	// number of threads
	//const int nt = 1;

	// domain
	const double xmin = 0;
	const double xmax = 1;
	const double ymin = 0;
	const double ymax = 1;

//	// mass
//	const double mmin = 1.0;
//	const double mrange = 10;
//	
//	// get domain range and grid size
//	const double xrange = xmax-xmin;
//	const double yrange = ymax-ymin;
//	const double x_grid_size = (xmax-xmin)/pow(2.,maxNumLevels);
//	
//	const double y_grid_size = (ymax-ymin)/pow(2.,maxNumLevels);
//
//	// for measuring wall clock time
//	double start, end;
//	
//	// generate points
//	vector<particle> particleCoordinates;
//	particleCoordinates.resize(N);
//	#pragma omp parallel for shared(particleCoordinates)
//	for(int i=0; i<N; i++){
//
//	if(i<N/10){
//			particleCoordinates[i].gen_coords_cluster(xmin, xrange, ymin, yrange,
//				mmin, mrange, xmin+xrange/4,
//				ymin+yrange/3, min(xrange,yrange)/5);
//		}
//		else if(i<4*N/10) {
//			particleCoordinates[i].gen_coords_cluster(xmin, xrange, ymin, yrange,
//				mmin, mrange, xmin+3*xrange/4,
//				ymin+2*yrange/3, min(xrange,yrange)/6);
//		}
//		else{
//			particleCoordinates[i].gen_coords_cluster(xmin, xrange, ymin, yrange,
//				mmin, mrange,
//				xmin+xrange/4, ymin+3*yrange/4,
//				min(xrange,yrange)/5);
//		}
//
//	// particleCoordinates[i].gen_coords(xmin, xrange, ymin, yrange, mmin, mrange);
//	}
//	
//	// compute morton ids
//	start=omp_get_wtime();	
//	#pragma omp parallel for default(none), shared(particleCoordinates, N, maxNumLevels)
//	for(int i=0; i<N; i++)
//		particleCoordinates[i].get_morton_id(xmin, ymin, x_grid_size, y_grid_size, maxNumLevels);
//	end=omp_get_wtime();	
//	cout<<"get_morton_id: "<<end-start<<endl;
//	
//	// parallel merge sort
//	// vector<particle> tmp;
//	// tmp.resize(N);
//	// mergesort<particle>(&particleCoordinates[0], nt, N, &tmp[0], compar_vector);
//
//	// read points data
//	// read_points(particleCoordinates, N);
//
//	cout<<endl;
//	for(int i=0; i<5; i++){
//		n_th=pow(2.0,i);
//		cout<<"proc: "<<n_th<<endl;
//			
//		qtree* qt1 = new qtree;
//		// parallel tree construction
//		start=omp_get_wtime();
//		for(int i=0; i<nt; i++){
//			// build tree
//			build_tree( xmin, ymin, maxNumPointsPerNode, maxNumLevels, nodeSize,
//						&particleCoordinates, N, qt1, 0, N, 1);
//		}
//		end=omp_get_wtime();
//		cout<<"build_tree: "<<end-start<<endl;
//
//		start=omp_get_wtime();	
//		#pragma omp parallel for shared(qt1, particleCoordinates) num_threads(n_th)
//		for(int i=0; i<N; i++)
//			evaluate_trees(i, qt1, particleCoordinates);
//		end=omp_get_wtime();
//		cout<<"evaluate_trees: "<<end-start<<endl;
//	
//		// ofstream ofs;
//		// write_boxes( qt1, ofs, 1 );
//		// write_points(particleCoordinates, N);
//	
//		delete qt1;
//		
//		cout<<endl;
//	}
//	
//	// n_th=2;
//	// qtree* qt2 = new qtree;
//	// // parallel tree construction
//	// start=omp_get_wtime();
//	// // #pragma omp parallel for shared(particleCoordinates, np_local)
//	// for(int i=0; i<nt; i++){
//	// 	// build tree
//	// 	build_tree( xmin, ymin, maxNumPointsPerNode, maxNumLevels, nodeSize,
//	// 	&particleCoordinates, N, qt2, 0, N, 2);
//	// }
//	// end=omp_get_wtime();
//	// cout<<"wall clock time = "<<end-start<<endl;
//
//	// n_th=4;
//	// qtree* qt4 = new qtree;
//	// // parallel tree construction
//	// start=omp_get_wtime();
//	// // #pragma omp parallel for shared(particleCoordinates, np_local)
//	// for(int i=0; i<nt; i++){
//	// 	// build tree
//	// 	build_tree( xmin, ymin, maxNumPointsPerNode, maxNumLevels, nodeSize,
//	// 	&particleCoordinates, N, qt4, 0, N, 4);
//	// }
//	// end=omp_get_wtime();
//	// cout<<"wall clock time = "<<end-start<<endl;
//
//	// n_th=8;
//	// qtree* qt8 = new qtree;
//	// // parallel tree construction
//	// start=omp_get_wtime();
//	// // #pragma omp parallel for shared(particleCoordinates, np_local)
//	// for(int i=0; i<nt; i++){
//	// 	// build tree
//	// 	build_tree( xmin, ymin, maxNumPointsPerNode, maxNumLevels, nodeSize,
//	// 	&particleCoordinates, N, qt8, 0, N, 8);
//	// }
//	// end=omp_get_wtime();
//	// cout<<"wall clock time = "<<end-start<<endl;
//
//	// delete qt1, qt2, qt4, qt8;
	
	return 0;
}
//
//
//// int main()
//// {
//// 		for(int i=1; i<67; i+=2)
//// 			cout<<setprecision(25)<<pow(2,(i-1)) + pow(2,i)<<", ";
//// 		// for(int i=1; i<33; i++)
//// 			// cout<<setprecision(25)<<2*pow(4,(i-1))<<", ";
//
//// }
