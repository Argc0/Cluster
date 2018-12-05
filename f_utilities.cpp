#include "f_utilities.h"

//enimerwsi twn kentrwn me to meso tou kathe cluster (means)
void update_kmeans(Dots *dots, vector < vector<TP> > *centroids, int num_clusters, int dimension){
	vector <TP> new_centroid;
	vector<TP> *temp;
	int size_cluster;
	double mean = 0.0, sum = 0.0;

	for(int i=0; i<num_clusters; i++){
		for(int j=0; j<dimension; j++){
			sum = 0.0;
			size_cluster = 0;
			for(int k = 0; k < dots->getSize(); k++){
				if(i == dots->returnIndex(k)){
					temp = dots->returnPoint(k);
					sum += temp[0][j];
					size_cluster++;
				}
			}
			if(size_cluster == 0)
				mean=0;
			else
				mean = sum/size_cluster;
			new_centroid.push_back(mean);
		}
		centroids[0][i] = new_centroid;
		new_centroid.clear();
	}
}

//enimerwsi twn kentrwn me ti diameso(simeio pou to athroismma twn apostasewn apo ta alla simeia tou cluster einai to mikrotero) tou kathe cluster
void update_kmedoid(Dots * dots, vector < vector<TP> > *centroids, vector<string> * centroids_ids, int num_clusters, int m){
	vector <TP> *new_centroid;
	string new_id;
	vector<TP> *temp1, *temp2;
	vector<double> distances;
	double min_sum = double(INT_MAX);
	double dist = 0.0, sum = 0.0;

	for(int i=0; i < num_clusters; i++){
		min_sum = double(INT_MAX);
		for(int j=0; j < dots->getSize(); j++){
			sum=0.0;
			if(i == dots->returnIndex(j)){
				if(centroids_ids[0][i] != dots->returnIdPoint(j)){
					temp1 = dots->returnPoint(j);
					for(int k=0; k < dots->getSize(); k++){
						if(k != j){
							if(i == dots->returnIndex(k)){
								temp2 = dots->returnPoint(k);
								if(m)
									dist=cosine_distance(*temp1, *temp2);
								else
									dist=euclidean_distance(*temp1, *temp2);
								distances.push_back(dist);
							}
						}
					}
					sum = accumulate(distances.begin(), distances.end(), 0.0);
				}
				if(distances.size() == 0){
						new_centroid = dots->returnPoint(j);
						new_id = dots->returnIdPoint(j);
				}else{
					if(sum < min_sum){
						min_sum = sum;
						new_centroid = dots->returnPoint(j);
						new_id = dots->returnIdPoint(j);
					}
				}
				distances.clear();
			}
		}
		centroids[0][i] = *new_centroid;
		centroids_ids[0][i] = new_id;
	}
}

//sinartisi pou ilopoihei ton aplo algorithmo anathesis tou Lloyd's
void Lloyds_assignment(Dots *dots, vector< vector<TP> > *points, vector < vector<TP> > *centroids, int& counter, int n, int metric){
	int index;
	for(int i=0; i<n; i++){
		min_distance(points[0][i], centroids, metric, &index);
		//elenxos an to simeio anike sto idio cluster proigoumenws (simvalei sti sinthiki termatismou diladi ola ta simeia na vriskontai xwris allages se ena cluster)
		if(dots->returnIndex(i) != index)
			counter++;
		dots->setIndexforPoint(i, index);
	}
}

//sinartisi pou ilopoihei ton algoithmo range search me aktina to kentro tou cluster
void range_assignment(Dots *dots, vector < vector<TP> > *centroids, HashTable **table, GFunctions ** fs, Bcube * cube, int num_hash_tables, int num_hash_functions, int num_clusters, int& counter, int metric, int n, int choice){
	//korifes gia search
	int korifes = num_hash_functions / 2;
	double r = min_distance_between_centers(centroids, metric, num_clusters)/2;
    int assignment = n;

    dots->refreshMarked();
    //oi anatheseis na 3epernoun to 0.1% tou sinolikou plithous twn simeiwn
    while(assignment > 0.01*n ){
    	assignment = 0;
    	if(choice == 1)
    		Range_Search(centroids, dots, r, counter, assignment, table, fs, num_hash_tables, num_clusters, metric);
    	else if(choice == 2)
    		Range_Search_cube(centroids, dots, r, counter, assignment, cube, korifes, num_clusters, metric);
    	r = r*2;
	}
    rest_assignment(centroids, dots, counter, metric, num_clusters);
}


//anazitisi prosegistika geitonwn twn kentrwn twn clusters entos aktinas r
void Range_Search(vector < vector<TP> > *q, Dots *dots, double r, int& counter, int& assignment, HashTable ** table, GFunctions ** g, int L, int num_clusters, int metric){
	double dist;
	long long int num;
	//gia apofigi diplotipwn
	list<string> key_list;
	list<HashNode<TP>*> * temp;
	//domi pou voithaei stin prosorini apothikeusi tou deikti enos cluster pou anikei ena simeio
	map <string, int> map_index;

	for(int index=0; index < num_clusters; index++){
		for(int i=0; i < L; i++){
			num = g[i]->getNum(q[0][index]);
			temp = table[i]->return_bucket(g[i]->getBucket_Num(num));
			//if(temp->size() > 3*L) break; //trick stopping sooner the algorithm
			for (std::list< HashNode<TP>* >::iterator it = temp->begin(); it != temp->end(); ++it){
				string k = (*it)->getKey();
				long long int g = (*it)->getG();
				if(g != num) continue;
				vector <TP> *p = (*it)->getValue();
				if(!metric)
					dist = euclidean_distance(q[0][index], *p);
				else
					dist = cosine_distance(q[0][index], *p);	
				if( dist < r){
					list<string>::iterator result1 = find(key_list.begin(), key_list.end(), k);
    				if (result1 == key_list.end()) {
        				key_list.push_back(k);

        				if(dots->getMark(k)==0){
        					dots->setMark(k,1);
        					assignment++;
        					if(dots->returnIndexbyKey(k) == -1)
        						map_index[k] = index;
        					else
        						map_index[k] = dots->returnIndexbyKey(k);
        				}else{
        					vector < vector<TP> > n;
        					int new_index, final_index;
        					if(dots->returnIndexbyKey(k) == index) continue;
        					//se periptwsi conflict ipologizete to pragmatika kontinotero kentro ki anathetetai sto simeio
        					n.push_back(q[0][dots->returnIndexbyKey(k)]);
        					n.push_back(q[0][index]);
        					min_distance(*p, &n, metric, &new_index);
        					if(new_index == 0){
        						final_index = dots->returnIndexbyKey(k);
        					}else{
        						final_index = index;
        					}
        					map_index[k]=final_index;
        				}
    				}
				}	
			}
		}
	}
	//telikh anathesi twn deiktwn twn cluster sta simeia
	for ( map<string, int>::iterator iter = map_index.begin(); iter != map_index.end(); iter++ ){
		if(dots->returnIndexbyKey(iter->first) != iter->second)
        	counter++;
    	dots->setIndexbyKey(iter->first, iter->second);
	}
	
}

//evresi prosegistikwn geitonwn apo ta kentra twn clusters entos aktinas r me ton cube
void Range_Search_cube(vector < vector<TP> > *q, Dots *dots, double r, int& counter, int& assignment, Bcube * cube, int probes, int num_clusters, int metric){
	double dist;
	int counter_probes =0, n; 
	//gia apofigi diplotipwn
	list<string> key_list;
	list<HashNode<TP>*> * temp;
	//domi pou voithaei stin prosorini apothikeusi tou deikti enos cluster pou anikei ena simeio
	map <string, int> map_index;

	int d = cube->return_dimension();

	vector<int> *num = new vector<int>();

	for(int index=0; index < num_clusters; index++){
		//evresi korifwn tou cube gia anazitisi
		n = cube->insert_num(q[0][index]);
		num->push_back(n);
		find_probes(n, probes-1, d, num);
		counter_probes = 0;

		while( counter_probes < probes){
			temp = cube->return_bucket(num[0][counter_probes]);
			for (std::list< HashNode<TP>* >::iterator it = temp->begin(); it != temp->end(); ++it){
				string k = (*it)->getKey();
				vector <TP> *p = (*it)->getValue();
				if(!metric)
					dist = euclidean_distance(q[0][index], *p);
				else
					dist = cosine_distance(q[0][index], *p);

				if( dist < r){
					list<string>::iterator result1 = find(key_list.begin(), key_list.end(), k);
					if (result1 == key_list.end()) {
        				key_list.push_back(k);

        				if(dots->getMark(k) == 0){
        					dots->setMark(k, 1);
        					assignment++;
        					if(dots->returnIndexbyKey(k) == -1)
        						map_index[k] = index;
        					else{
        						map_index[k] = dots->returnIndexbyKey(k);
        					}
        				}else{
        					vector < vector<TP> > v;
        					int new_index, final_index;
        					if(dots->returnIndexbyKey(k) == index) continue;
        					//se periptwsi conflict ipologizete to pragmatika kontinotero kentro ki anathetetai sto simeio
        					v.push_back(q[0][dots->returnIndexbyKey(k)]);
        					v.push_back(q[0][index]);
        					min_distance(*p, &v, metric, &new_index);
        					if(new_index == 0){
        						final_index = dots->returnIndexbyKey(k);
        					}else{
        						final_index = index;
        					}
        					map_index[k] = final_index;
        				}	
    				}
				}
			}
			counter_probes++;
		}
		num->clear();
	}
	//telikh anathesi twn deiktwn twn cluster sta simeia
	for (map<string, int>::iterator iter = map_index.begin(); iter != map_index.end(); iter++ ){
		if(dots->returnIndexbyKey(iter->first) != iter->second)
        	counter++;
    	dots->setIndexbyKey(iter->first, iter->second);
	}
	delete num;
}

//sinartisi pou kanei apli anathesi twn clusters sta simeia meta to range search
void rest_assignment(vector < vector<TP> > *centroids, Dots * dots, int& counter, int metric, int num_clusters){
	int n = dots->getSize();
	int index;

	for(int i=0; i<n; i++){
		if(dots->get_Marked(i) == 0){
			min_distance((*dots->returnPoint(i)), centroids, metric, &index);
      		if(dots->returnIndex(i) != index)
        		counter++;
      		dots->setIndexforPoint(i, index);
		}
	}
}

//sinartisi pou epistrefei stin metavliti sil tis times ths silhouette gia kathe cluster kai gia ola ta cluster mazi
void silhouette(vector < vector<TP> > *centroids, Dots * dots, double* sil, int metric, int num_clusters, int num){
	vector<TP> *temp, *temp1, *temp2;
	int index, neigh_index;
	int size_a = 0, size_b=0;
	double dist, a = 0.0, b=0.0, s=0.0;

	//refresh the array of silhouettes
	for(int i=0; i<num_clusters +1; i++){
		sil[i]=0.0;
	}
	//computation of silhouette
	for(int i=0; i<num; i++){
		temp = dots->returnPoint(i);
		index = dots->returnIndex(i);
		vector<vector<TP> > elem(centroids[0]);
		elem.erase(elem.begin()+ index);
		min_distance(centroids[0][index], &elem, metric, &neigh_index);
		if(index <= neigh_index)
			neigh_index = neigh_index + 1;
		size_a=0; size_b=0; a=0.0; b=0.0;
		for(int j=0; j<num; j++){
			if(i != j){
				if(dots->returnIndex(j) == index){
					temp1 = dots->returnPoint(j);
					if(!metric)
						dist = euclidean_distance(*temp, *temp1);
					else
						dist = cosine_distance(*temp, *temp1);
					a += dist;
					size_a++;
				}
			}
			if(dots->returnIndex(j) == neigh_index){
				temp2 = dots->returnPoint(j);
				if(!metric)
					dist = euclidean_distance(*temp, *temp2);
				else
					dist = cosine_distance(*temp, *temp2);
				b += dist;
				size_b++;
			}
		}
		//gia size_a == 0 or size_b == 0 exei apotelesma NaN ki paei sto else diladi sto 0
		a = a/size_a;
		b = b/size_b;

		if(a < b)
			s = 1 - a/b;
		else if(a > b)
			s = b/a - 1;
		else
			s = 0.0;

		sil[index] += s;
		sil[num_clusters] += s;
	}
}

//ektipwsi twn apotelesmatwn simfwna me ta zitoumena tis ekfwnisis
void print_out(Dots * dots, vector < vector<TP> > *centroids, vector<string> * centroids_ids, int metric, double secs, double * silhouette, int kmeans, int num_clusters, int in, int a, int u, ofstream& myfile_out){
	vector <string> temp_cluster_ids[num_clusters];

	//prosorini apothikeusi twn ids twn simiwn gia kathe cluster gia kaliteri organwsi tis ektipwsis
	for(int index = 0; index < num_clusters; index++){
		for(int j=0; j<dots->getSize(); j++){
			if(dots->returnIndex(j) == index){
				temp_cluster_ids[index].push_back(dots->returnIdPoint(j));	
			}
		}
	}

	if (myfile_out.is_open()){
		//Algorithm
		myfile_out << "Algorithm: I" << in << "A" << a << "U" << u << endl;
		//Metric
		if(metric)
			myfile_out << "Metric: Cosine" << endl;
		else
			myfile_out << "Metric: Euclidean" << endl;
		//Centroids
		for(int index = 0; index < num_clusters; index++){
			myfile_out << "CLUSTER-" << index+1 << " {size: " << temp_cluster_ids[index].size() << ", centroid: ";
			if(kmeans){
    			for(vector<TP>::const_iterator j = centroids[0][index].begin(); j != centroids[0][index].end(); ++j) {
      				myfile_out << *j << " ";
    			}
    			myfile_out << "}" << endl;
			}else{
				myfile_out << centroids_ids[0][index] << " }" << endl;
			}
		}
		//Time
		myfile_out << "clustering_time: " << secs << endl;
		//Silhouette
		myfile_out << "Silhouette: [";
		for(int i=0; i < num_clusters + 1; i++){
			if( i == num_clusters)
				myfile_out << silhouette[i]/dots->getSize() << "]" << endl;
			else
				myfile_out << silhouette[i]/temp_cluster_ids[i].size() << ", ";

		}
		//Clusters
		for(int index = 0; index < num_clusters; index++){
			myfile_out << "CLUSTER-" << index+1 << " {";
			for(int j=0; j<(int)temp_cluster_ids[index].size(); j++){
				if(j == (int)temp_cluster_ids[index].size()-1)
					myfile_out << temp_cluster_ids[index][j] << "}" << endl;
				else
					myfile_out << temp_cluster_ids[index][j] << ", ";
			}
		}
	}else { 
      cout << "Unable to open output file";
      exit(-1);
    }

    myfile_out << endl;
}

//dimiourgia domwn voithitikwn domwn gia ton kwdika kai katanomi tous
void insert_in(Dots ** dots, HashTable *** table, GFunctions *** fs, Bcube ** cube, vector< vector<TP> > *points, vector<string> * in_ids, int n, int num_hash_tables, int num_hash_functions, int m){
  int i, j;
  int dimension = points[0][0].size();
  //for assignment
  *dots = new Dots(points, in_ids, n);

  //insert into hastable depending on metric
  *table = new HashTable*[num_hash_tables]();
  *fs = new GFunctions*[num_hash_tables]();
  int w = W;
  long long int g;
  int tablesz, in;

  if(m){
    tablesz = pow(2.0, double(num_hash_functions));
    for(i=0; i < num_hash_tables; i++)
      (*fs)[i] = new GCosine(num_hash_functions, dimension);

  }else{
    tablesz = n/8;
    for(i=0; i < num_hash_tables; i++)
      (*fs)[i] = new GFunction(num_hash_functions, dimension, w, tablesz);   
  }

  for(i=0; i < num_hash_tables; i++){
    (*table)[i] = new HashTable(tablesz);
    for (j = 0; j < n; j++){
      g = (*fs)[i]->getNum(points[0][j]);
      in = (*fs)[i]->getBucket_Num(g);
      
      (*table)[i]->add_item(in, in_ids[0][j], g, &points[0][j]);
    }
  }


  //insert into cube depending on metric
  *cube = new Bcube(num_hash_functions, dimension, w, m);

  for (i = 0; i < n; i++){
    g = (*cube)->insert_num(points[0][i]);
    (*cube)->add_item(g, in_ids[0][i], g, &points[0][i]); 
  }
}

//diagrafi twn domwn pou dimiourgithikan
void deletion(double *sil, vector<string> * in_ids, vector< vector<TP> > *points, vector<string> *centroids_ids, vector < vector<TP> > *centroids, GFunctions ** fs, HashTable **table, Bcube *cube, Dots *dots, int n, int num_clusters, int num_hash_tables){
  delete[] sil;

  in_ids->clear();
  delete in_ids;

  for(int i=0; i< n; i++)
    points[0][i].clear();
  points->clear();
  delete points;

  centroids_ids->clear();
  delete centroids_ids;

  for(int i=0; i< num_clusters; i++)
    centroids[0][i].clear();
  centroids->clear();
  delete centroids;

  for(int i=0; i<num_hash_tables; i++)
    delete fs[i];
  delete[] fs;

  for(int i=0; i<num_hash_tables; i++)
    delete table[i];
  delete[] table;

  delete cube;

  delete dots;
}