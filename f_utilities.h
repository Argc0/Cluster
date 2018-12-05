#ifndef __F_UTILITIES__
#define __F_UTILITIES__

#include "cube.h"
#include "Dots.h"

/* Simple assignment with Lloyd's algorithm */
void Lloyds_assignment(Dots *dots, vector< vector<TP> > *points, vector < vector<TP> > *centroids, int& counter, int n, int metric);

/* Assignment by Range Search LSH or Hypercube */
void range_assignment(Dots *dots, vector < vector<TP> > *centroids, HashTable **table, GFunctions ** fs, Bcube * cube, int num_hash_tables, int num_hash_functions, int num_clusters, int& counter, int metric, int n, int choice);

/* Simple assignment for the rest points of the Range Search algorithm */
void rest_assignment(vector < vector<TP> > *centroids, Dots * dots, int& counter, int metric, int num_clusters);

/* A Range Search with the use of LSH */
void Range_Search(vector < vector<TP> > *q, Dots *dots, double r, int& counter, int& assignment, HashTable ** table, GFunctions ** g, int L, int num_clusters, int metric);

/* A Range Search with the use of Hypercube */
void Range_Search_cube(vector < vector<TP> > *q, Dots *dots, double r, int& counter, int& assignment, Bcube * cube, int probes, int num_clusters, int metric);

/* Create/make structs/classes and insert the points for utilization for the next steps */
void insert_in(Dots ** dots, HashTable *** table, GFunctions *** fs, Bcube ** cube, vector< vector<TP> > *points, vector<string> * in_ids, int n, int num_hash_tables, int num_hash_functions, int m);

/* Deletion of the structs/classes */
void deletion(double *sil, vector<string> * in_ids, vector< vector<TP> > *points, vector<string> *centroids_ids, vector < vector<TP> > *centroids, GFunctions ** fs, HashTable **table, Bcube *cube, Dots *dots, int n, int num_clusters, int num_hash_tables);

/* Update the centroids of a cluster with kmeans algorithm */
void update_kmeans(Dots * dots, vector < vector<TP> > *centroids, int num_clusters, int dimension);

/* Update the centroids of a cluster with Partitioning Around Medoids (PAM) improved like Lloyd’s (kmedoid) */
void update_kmedoid(Dots * dots, vector < vector<TP> > *centroids, vector<string> * centroids_ids, int num_clusters, int m);

/* Computation of the Silhouette function for each cluster and all */
void silhouette(vector < vector<TP> > *centroids, Dots * dots, double *sil, int metric, int num_clusters, int num);

/* Prints out the result of each iteration as requested by the exercise in the output file */
void print_out(Dots * dots, vector < vector<TP> > *centroids, vector<string> * centroids_ids, int metric, double secs, double * silhouette, int kmeans, int num_clusters, int i, int a, int u, ofstream& myfile_out);

#endif