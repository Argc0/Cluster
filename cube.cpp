#include "cube.h"
//kwdikas apo tin prwti ergasia

//kataskevi tou yperkivou simfwna me ti theoria
Bcube::Bcube(int &d, int dimension, int w, int &choice) : d(d),choice(choice), dimension(dimension){
  int c = int(pow(2.0, double(d)));
  buckets = new HashBucket*[c]();
  for (int i = 0; i < c; ++i){
    buckets[i] = new HashBucket();
  }
  hfs = new HFunctions*[d]();
  //1 for cosine 0 for lsh euclidean
  if(choice == 1){
    for (int i = 0; i < d; ++i){
      hfs[i] = new HCosine(dimension);
    }
  }else{
    for (int i = 0; i < d; ++i){
      hfs[i] = new HFunction(dimension, w);
    }
  }
}

//katastrofi tou iperkivou
Bcube::~Bcube(){
  int c = int(pow(2.0, double(d)));
  for (int i = 0; i < c; ++i){
    delete buckets[i];
  }
  delete[] buckets;
  for (int i = 0; i < d; ++i){
    delete hfs[i];
  }
  delete[] hfs;
}

//epistrefei tin katallili timi gia na isaxthei to stoixeio sto katallillo bucket
int Bcube::insert_num(vector<TP>& v){
  int key = 0;
  if(choice){
    for(int i = 0; i < d; i++){
      key += hfs[d-(i+1)]->getResult(v)*int(pow(2.0, double(i)));
    }
  }else{
    int res = 0,coin;
    map<int,int>::iterator it;
    for(int i = 0; i < d; i++){
      res = hfs[d-(i+1)]->getResult(v); 
      coin = rand() % 2;
      //ripsi nomismatos ki apothikeusi timis se periptwsi pou o arithmos epanalipthei na paei to idio apotelesma ripsis tou nomismatos
      it = pam.find(res);
      if (it == pam.end())
        pam[res]=coin;
      else
        coin = it->second;
      key += coin*int(pow(2.0, double(i)));
    }
  }
  return key;
}

//prosthiki stoixeiou me tis domis hashnode sto bucket
void Bcube::add_item(int key, string &id, long long int g, vector<TP>* v){
  HashNode<TP> *node = new HashNode<TP>(id, g, v);
  buckets[key]->add_element(node);
}

//voithitiki sinartisi ektiposis
void Bcube::print_cube(){
  cout << "Print Binary Cube\n";
  int c = int(pow(2.0, double(d)));
  for(int i = 0; i < c; ++i){
    cout << "Bucket " << i << endl;
    buckets[i]->print_list();
  }
}
