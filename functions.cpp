#include "functions.h"

//sinartisi pou epistrefei ta orismata apo tin grammi entolwn
void get_args( int argc, char** argv, string& input_name, string& configuration_name, string& output_name, string& metric){
	int i, j;

	if (argc > 1) {

		if(argc != 9){
			cout << "The arguments that given are not right!\n";
			exit(-1);
		}

    	for (i = 1; i < argc-1; i += 2) {

			// Check if argument is given again
			for( j = i+2; j < argc-1; j += 2 ) {

				if( !strcmp(argv[i], argv[j]) ) {
					cout << "Error: Argument given again.\nPlease check README for more info.\n";
					exit(-1);
				}
			}
			// Get arguments
			if( !strcmp(argv[i], "-i") ) {
				input_name=argv[i+1];
			} else if( !strcmp(argv[i], "-c") ) {
				configuration_name = argv[i+1];
			} else if( !strcmp(argv[i], "-o") ) {
				output_name = argv[i+1];
			} else if( !strcmp(argv[i], "-d") ) {
				metric = argv[i+1];
			} else {
				cout << "Error: Wrong argument given.\nPlease check README for more info.\n";
				exit(-1);
			}
		}
	} else {
		cout << "Error: Arguments are missing! Please check README for more info\n";
	}
	return;
}

//eswteriko ginomeno dianismatwn
double internal_product(vector<TP> &p, vector<double> &v){
	int size = p.size();
	double sum = 0.0;
	for(int i = 0; i < size; ++i) {
		sum += p[i]*v[i];
	}
	return sum;
}

//euclidean apostasi dianismatwn
double euclidean_distance(vector<TP> &p, vector<TP> &v){
	int size = p.size();
	double d = 0.0;

	for(int i = 0; i < size; ++i) {
		//(y1-x1)^2
		d += (v[i] - p[i])*(v[i]-p[i]);
	}
	return sqrt(d);
}

//apostasi sinimitonou dianismatwn
double cosine_distance(vector<TP> &x, vector<TP> &y){
	int size = x.size();
	double d = 0.0, in = 0, a = 0, b = 0;

	//internal_product
	for(int i = 0; i < size; ++i) {
		in += x[i]*y[i];
	}
	for(int i = 0; i < size; ++i) {
		a += x[i]*x[i];
		b += y[i] *y[i]; 
	}

	d = in/(sqrt(a)*sqrt(b));
	return 1 - d;

}

//sinartisi pou epistrefei tin apostasi tou pragmatika kontinoterou geitona kai ton deikti (metavliti index) autou
double min_distance(vector<TP> &q, vector< vector<TP> > * points, int metric, int* index){
	int n = points[0].size();
	double true_d = double(INT_MAX), d;
	if(metric){
		for(int x=0; x < n; x++){
            d = cosine_distance(q, points[0][x]);
            if( true_d > d){
            	true_d = d;
            	*index = x;
            }
        }
	}else{
		for(int x=0; x < n; x++){
            d = euclidean_distance(q, points[0][x]);
            if( true_d > d){
            	true_d = d;
            	*index = x;
            }
        }       
	}
	return true_d;
}

//https://www.geeksforgeeks.org/binary-search/
//diadiki anazitisi
int binarySearch(double P[], int left, int right, double x) { 
   if (right >= left) { 
        int mid = left + (right - left)/2; 
  
        if(mid == 0)
        	return mid;

        if ((P[mid-1] < x) && (x <= P[mid]))   
            return mid; 
  
        if (P[mid] > x)  
            return binarySearch(P, left, mid-1, x); 
        
        return binarySearch(P, mid+1, right, x); 
   }  
   return -1; 
}

//anagnwsi twn parametrwn apo to configuration file
void read_configuration_file(string& configuration_name, int& num_clusters, int& num_hash_functions, int& num_hash_tables){
	ifstream myfile_con;
	string line;

	myfile_con.open(configuration_name.c_str());

	if (myfile_con.is_open()){

		while ( getline (myfile_con,line) ){
          long unsigned int position;
          string num;
          if (line.find("number_of_clusters:") != std::string::npos){
            position = line.find(":");
            num = line.substr (position + 1);
           	num_clusters =  stoi(num);
           	continue;
          }
          if (line.find("number_of_hash_functions:") != std::string::npos){
            position = line.find(":");
            num = line.substr (position + 1);
           	num_hash_functions =  stoi(num);
           	continue;
          }
          if (line.find("number_of_hash_tables:") != std::string::npos){
            position = line.find(":");
            num = line.substr (position + 1);
           	num_hash_tables =  stoi(num);
           	continue;
          }
		}
		myfile_con.close();
	}else {
		cout << "Unable to open configuration file";
		exit(-1);
	}
}

//anagnwsi twn simeiwn apo to input file
void read_input_file(string& input_name, vector<string> * in_ids, vector< vector<TP> > *points, int& n){
	ifstream myfile_in;
	
	string line;
	istringstream iss;
  	vector<TP> inputs;
  	
  	double it;
	string id;
	//metavliti voithitiki gia tin anagnwsi ki apothikeusi tou id tou simeiou
	int flag = 1;

	myfile_in.open(input_name.c_str());

  	//read from input file the points
  	if (myfile_in.is_open()){
    	while ( getline (myfile_in,line) ){
      		iss.clear();
      		iss.str(line);

      		flag=1;
    		while (iss >> it){

    			if(flag){
    				ostringstream strs;
					strs << it;
					string str = strs.str();
    				in_ids->push_back(str);
    				flag = 0;
    			}else{
        			inputs.push_back(it);
            	}
            	if (iss.peek() == ',' || iss.peek() == ' ' || iss.peek() == '\t')
            			iss.ignore();
    		}

      		points->push_back (inputs);
 
      		inputs.clear();
      		n++;
    	}
    	myfile_in.close();
  	}else { 
    	cout << "Unable to open input file";
    	exit(-1);
  	}
}

//tixaia arxikopoihsh kentrwn apo simeia pou anikoun sto dataset
void initialization_random(vector< vector<TP> > *points, vector<string> * in_ids, vector < vector<TP> > *centroids, vector<string> * centroids_ids, int& num_clusters, int& n){
	vector <int> num_centroids;
    vector<int>::iterator iter;

    while((int) num_centroids.size() != num_clusters){
      
      //n == total number of points
      int rand_num = rand() % n;
      iter = find (num_centroids.begin(), num_centroids.end(), rand_num);
      //oles oi tixaies times pou tha paraxthoun na einai diaforetikes
      if (iter == num_centroids.end())
         num_centroids.push_back(rand_num);
    }

    for(int i=0; i < num_clusters; i++){
      centroids->push_back(points[0][num_centroids.at(i)]);
      centroids_ids->push_back(in_ids[0][num_centroids.at(i)]);
    }
}

//arxikopoihsh kentrwn apo simeia tou dataset pou einai apomakrismena meta3i tous
void initialization_kmeanspp(vector< vector<TP> > *points, vector<string> * in_ids, vector < vector<TP> > *centroids, vector<string> * centroids_ids, int& num_clusters, int& m, int& n){
	int t=1, rn, pi, i, j, temp;
    double rand_num;
    double *mx;

    //antigrafi twn simeiwn ki twn id tous se nees metavlitis gia kaliteri diaxeirisi tou algorithmou
    vector<string> elem_id(in_ids[0]);
    vector<vector<TP> > elem(points[0]);
    rn = rand() % n;

    centroids_ids->push_back(in_ids[0][rn]);
    elem_id.erase(elem_id.begin() + rn);
    //eisagwgi tixaiou simeiou ston vector twn centroid ki diagrafi apo ti thesi pou eixe sti nea metavliti gia na min ipologizetai 3ana
    centroids->push_back(points[0][rn]);
    elem.erase(elem.begin()+ rn);

    //algorithmos tou kmeans++
	while(t != num_clusters){
		double P[n-t + 1];
		double D[n-t + 1];
		
		D[0]=0.0;
		for(i=0; i < n-t ; i++){
			D[i+1]= min_distance(elem.at(i), centroids, m, &temp);
        }
		mx=max_element(D,D+(n-t));
		for(i=0; i< n-t+1; i++){
			P[i]=0.0;
			for(j=0; j < i; j++){
				P[i] += (D[j+1]/mx[0])*(D[j+1]/mx[0]);
			}
		}
		//[0,P[n-t]]
		rand_num = (static_cast <double> (rand()) / static_cast <double> (RAND_MAX + 1.0)) * (P[n-t] + 1.0);
		
		pi=binarySearch(P, 0, n-t, rand_num);
		//periptwsi akraiwn periptwsewn
		if(pi != 0){ 
			pi = pi-1;
			if(pi == -2)
				pi = n-t;
		}

		centroids->push_back(elem.at(pi));
		elem.erase(elem.begin()+ pi);

		centroids_ids->push_back(elem_id.at(pi));
		elem_id.erase(elem_id.begin() + pi);
		
		t = t + 1;
	}
}

//epistrefei tin mikroteri apostasi anamesa sta kentra twn clusters
double min_distance_between_centers(vector < vector<TP> > *centroids, int m, int num_clusters){
	vector<double> distances_between;
	double dist;

	for(int i=0; i<num_clusters; i++){
		for(int k=i+1; k < num_clusters; k++){
			if(m)
				dist=cosine_distance(centroids[0][i], centroids[0][k]);
			else
				dist=euclidean_distance(centroids[0][i], centroids[0][k]);
			distances_between.push_back(dist);
		}
	}
	if(distances_between.size() == 0 )
		return 0;
	else
		return *min_element(distances_between.begin(),distances_between.end());
}

//epistrefei olous tous geitones apo apostasi haming 1 mexri d
//https://stackoverflow.com/questions/40813022/generate-all-sequences-of-bits-within-hamming-distance-t
void hamming(int n, int i, int rest, vector<int>* num, int* counter, int max_probes) {
	if (rest == 0) {
		num->push_back(n);
		*counter = *counter + 1;
		return;
	}

	if (i < 0) return;

	if(*counter == max_probes) return;
	//change currebt bit
	n ^= 1U << i;
	hamming(n, i-1, rest-1, num, counter, max_probes);

	if(*counter == max_probes) return;
	//change it again(undo)
	n ^= 1U << i;
	hamming(n, i-1, rest, num, counter, max_probes);

}

//epistrefei sti domi num tous katallilous geitones apostasis hamming
void find_probes(int n, int probes, int d, vector<int> *num) {
	int *counter = new int;
	*counter=0;
	for (int i = 1 ; i <= d ; ++i) {
		hamming(n, d-1, i, num, counter, probes);
	}
	delete counter;
}