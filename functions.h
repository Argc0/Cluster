#ifndef __FUNCTIONS__
#define __FUNCTIONS__

#include <iostream>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <list>
#include <numeric> 
#include <map>
#include <algorithm> 
#include <ctime> 
#include <iterator>
#include <math.h> 
#include <random>

//a define variable that changes when the points are int or double
#define TP double
//a max number
#define INT_MAX 34456778
//a max number for ri of lsh euclidean hash function (g)
#define MAX 20
//a w number for h functions of lsh euclidean
#define W 2

using namespace std;
using std::string;

/* Get arguments */
void get_args( int argc, char** argv, string& input_name, string& configuration_name, string& output_name, string& metric);

/* Euclidean distance for points */
double euclidean_distance(vector<TP> &p, vector<TP> &v);

/* Cosine distance for points */
double cosine_distance(vector<TP> &x, vector<TP> &y);

/* Internal product for points */
double internal_product(vector<TP> &p, vector<double> &v);

/* Computes hamming distance of the neighbor probes of the Hypercube */
void hamming(int n, int i, int rest, vector<int>* num, int* counter, int max_probes);

/* Finds the nearest adjacent probes for the Hypercube */
void find_probes(int n, int probes, int d, vector<int> *num);

/* Computes the minimum distance that a point has from the others and returns the value and the index of the closest point of them */
double min_distance(vector<TP> &q, vector< vector<TP> > * points, int metric, int* index);

/* Calculates the minimum distance between the centroids */
double min_distance_between_centers(vector < vector<TP> > *centroids, int m, int num_clusters);

/* Read the configuration file */
void read_configuration_file(string& configuration_name, int& num_clusters, int& num_hash_functions, int& num_hash_tables);

/* Read the input file */
void read_input_file(string& input_name, vector<string> * in_ids, vector< vector<TP> > *points, int& n);

/* Simplest random selection of k points */
void initialization_random(vector< vector<TP> > *points, vector<string> * in_ids, vector < vector<TP> > *centroids, vector<string> * centroids_ids, int& num_clusters, int& n);

/* Initialization using the kmeans++ algorithm */
void initialization_kmeanspp(vector< vector<TP> > *points, vector<string> * in_ids, vector < vector<TP> > *centroids, vector<string> * centroids_ids, int& num_clusters, int& m, int& n);

/* Binary Search for the kmeans++ algorithm */
int binarySearch(double P[], int left, int right, double x);

#endif