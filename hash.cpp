#include "hash.h"
//kedikas apo thn 1h ergasia

HCosine::HCosine(int dimension){
	random_device rd;
	mt19937 gen(rd());
	ri = new vector<double>();
	//kataskevazei dimension tixaious arithmous pou anikoun sth kanoniki katanomi
	normal_distribution<double> distribution(0.0,1.0);
	for(int i = 0; i < dimension; ++i) {
		double number = distribution(gen);
		ri->push_back(number);
	}
}

int HCosine::getResult(vector<TP> &p){
	double inner = internal_product(p, *ri);
    if( inner >= 0) 
		return 1;
	else
		return 0;
}

GCosine::GCosine(int k, int dimension){
	hc = new HCosine*[k]();
	for(int i = 0; i < k; ++i) {
		hc[i]= new HCosine(dimension);
	}
	GCosine::k = k;
}

long long int GCosine::getNum(vector<TP> &p){
	int res = 0;
	//ftiaxnei ton dekadiko arithmo gia to bucket apo ta 0 h 1 pou epistefei kathe hc sinartisi
	for(int i = 0; i < k; i++){
		res += hc[k-(i+1)]->getResult(p)*int(pow(2.0, double(i)));
	}
	return res;
}

HFunction::HFunction(int dimension, int w){
	random_device rd;
	mt19937 gen(rd());
	v = new vector<double>();
	//kataskevazei dimension tixaious arithmous pou anikoun stin kanoniki katanomi
	normal_distribution<double> distribution_n(0.0,1.0);
	for(int i = 0; i < dimension; ++i) {
		double number = distribution_n(gen);
		v->push_back(number);
	}
	HFunction::w = w;
	uniform_real_distribution<double> distribution_u(0.0,double(w));
	t = distribution_u(gen);
}

//lsh euclidean hi
int HFunction::getResult(vector<TP> &p){
	double product = internal_product(p, *v);
	double result = (product + t) / w ;
	return int(floor(result));
}

GFunction::GFunction(int k, int dimension, int w, int tablesize){
	GFunction::k=k;
	GFunction::tablesz=tablesize;
	r = new vector<int>();
	h = new HFunction*[k]();
	//kataskevazei k omoiomorfa tixaious arithmous gia ti sinartisi g
	for(int i = 0; i < k; ++i) {
		int num = (-MAX) + (rand() / (RAND_MAX + 1.0))*(MAX - (-MAX) +1);
		r->push_back(num);
	}
	for(int j = 0; j < k; ++j) {
		h[j] = new HFunction(dimension, w);
	}
}

//metatropi tou apotelesmatos tis g wste na min iparxei iperxilisi
int GFunction::getBucket_Num(long long int g){
	int in=0, in1=0;
	M = int(pow(2.0, 32.0)) - 5;
	in1 = (g % M + M) % M;
	in = (in1 % tablesz + tablesz) % tablesz;
	return in;
}

//epistrefei to apotelesma tis g se long long gia na min ipar3ei iperxilisi
long long int GFunction::getNum(vector<TP> &p){
	long long int sum = 0;
	for(int i = 0; i < k; ++i) {
		sum += (long long int)h[i]->getResult(p)*r[0][i];
	}
	return sum;
}

//voithitiki sinartisi ektiposis
void HashBucket::print_list(){
	cout << "Print list\n";
	for (std::list< HashNode<TP>* >::iterator it = hashlist->begin(); it != hashlist->end(); ++it){
		vector <TP> *vit = (*it)->getValue();
		for(vector<TP>::const_iterator i = vit->begin(); i != vit->end(); ++i) {
			cout << *i << " ";
		}
		cout << endl;
	}
}

//kataskevh hashtable
HashTable::HashTable(int &size) : size(size){
	array_of_buckets = new HashBucket*[size]();
	for (int i = 0; i < size; ++i){
		array_of_buckets[i] = new HashBucket();
	}
}

//prosthiki simeiou sti lista tou bucket
void HashTable::add_item(int key, string &id, long long int g, vector<TP>* v){
	HashNode<TP> *node = new HashNode<TP>(id, g, v);
    array_of_buckets[key]->add_element(node);
}

//voithitiki sinartisi ektiposis
void HashTable::print_table(){
	cout << "Print Hash Table\n";
	for(int i = 0; i < size; ++i){
		cout << "Bucket " << i << endl;
		array_of_buckets[i]->print_list();
	}
}