#include "functions.h"
#include "f_utilities.h"

int main(int argc, char *argv[]) {
	string input_name, output_name, configuration_name;
	string metric;
  ofstream myfile_out;

  vector<string> * in_ids = new vector<string>();
  vector< vector<TP> > *points = new vector< vector<TP> >() ;
  vector < vector<TP> > *centroids = new vector< vector<TP> >();
  vector<string> * centroids_ids = new vector<string>();
 
	int num_clusters, num_hash_functions=4, num_hash_tables=5;
  int counter, n=0, m, kmeans, dimension;

  clock_t begin_, end_;
  double elapsed_secs;

  Dots * dots;
  HashTable **table;
  GFunctions ** fs;
  Bcube * cube;

	srand((unsigned)time(0));

  //get arguments
	get_args(argc, argv, input_name, configuration_name, output_name, metric);

	cout << "The input file is: " << input_name << endl << "The configuration file is: " << configuration_name << endl << "The output file is: " << output_name << endl;
	cout << "The metric is: " << metric << endl;
  
  if(metric == "euclidean") m=0; else m=1;

  //read configuration file
  read_configuration_file(configuration_name, num_clusters, num_hash_functions, num_hash_tables);

	cout << "Num of clusters: " << num_clusters << endl << "Num of hash functions: " << num_hash_functions << endl << "Num of hash tables: " << num_hash_tables << endl;

  //read input file
  read_input_file(input_name, in_ids, points, n);

  dimension = points[0][0].size();

  double *sil = new double[num_clusters+1];

  //clustering algorithms

  insert_in(&dots, &table, &fs, &cube, points, in_ids, n, num_hash_tables, num_hash_functions, m);

  myfile_out.open(output_name.c_str());

  for(int  i= 0; i < 2; i++){
    for(int a = 0; a < 3; a++){
      for(int u = 0; u < 2; u++){
        cout << "Initialization: " << i << ", Assignment: " << a << ", Update: " << u << endl;
        begin_ = clock();
        //initialization
        if(i == 0)
          initialization_random(points, in_ids, centroids, centroids_ids, num_clusters, n);
        else
          initialization_kmeanspp(points, in_ids, centroids, centroids_ids, num_clusters, m, n);
        //assignment
        counter = n;
        while( counter > 0){
          counter=0;
          if(a == 0){
            Lloyds_assignment(dots, points, centroids, counter, n, m);
          }else {
            range_assignment(dots, centroids, table, fs, cube, num_hash_tables, num_hash_functions, num_clusters, counter, m, n, a);
          }
          //update
          if(u == 0){
            update_kmeans(dots, centroids, num_clusters, dimension);
            kmeans = 1;
          }else{
            update_kmedoid(dots, centroids, centroids_ids, num_clusters, m);
            kmeans = 0;
          }
        }
        silhouette(centroids, dots, sil, m, num_clusters, n);

        end_= clock();
        elapsed_secs = double(end_ - begin_) / CLOCKS_PER_SEC;

        print_out(dots, centroids, centroids_ids, m, elapsed_secs, sil, kmeans, num_clusters, i, a, u, myfile_out);
        centroids->clear(); centroids_ids->clear(); dots->clear_clusters();
      }
    }
  }
  myfile_out.close();

  deletion(sil, in_ids, points, centroids_ids, centroids, fs, table, cube, dots, n, num_clusters, num_hash_tables);

  cout << "End\n";
}