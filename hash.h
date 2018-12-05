#ifndef __HASH__
#define __HASH__

#include "functions.h"

//kedikas apo thn ergasia 1

//abstract classes
class GFunctions {
    public:
        virtual ~GFunctions() {}
        virtual long long int getNum(vector<TP> &p) = 0;
        virtual int getBucket_Num(long long int g) = 0;
};

class HFunctions {
    public:
        virtual ~HFunctions() {}
        virtual int getResult(vector<TP> &p) = 0;
};

//classes for hashtable and lsh hash functions

class HCosine : public HFunctions{
    private:
        vector<double> *ri;
    public:
        HCosine(int dimension);
        ~HCosine(){
            ri->clear();
            delete ri;
        }
        vector<double> * getRi() {
            return ri;
        }
        void setRi(vector<double>* ri) {
            HCosine::ri = ri;
        }
        int getResult(vector<TP> &p);
};

class GCosine : public GFunctions{
    private:
        HCosine ** hc;
        int k;
    public:
        GCosine(int k, int dimension);
        ~GCosine(){
            for(int i = 0; i < k; ++i) {
                delete hc[i];
            }
            delete[] hc;
        }
        HCosine ** getHc() {
            return hc;
        }
        void setHc(HCosine** hc) {
            GCosine::hc = hc;
        }
        int getBucket_Num(long long int g){
            return g;
        }
        long long int getNum(vector<TP> &p);
};

class HFunction : public HFunctions{
	private:
		double t;
		vector<double> *v;
        int w;
	public:
        HFunction(int dimension, int w);
        ~HFunction(){
            v->clear();
            delete v;
        }
        double getT() {
            return t;
        }
        vector<double> * getV() {
            return v;
        }
        void setT(int t) {
            HFunction::t = t;
        }
        void setV(vector<double>* v) {
            HFunction::v = v;
        }
        int getResult(vector<TP> &p);
};

class GFunction : public GFunctions {
    private:
        vector<int> *r;
        HFunction ** h;
        int M;
        //gia na ginei to mod gia to katallilo bucket
        int tablesz;
        int k;
    public:
        GFunction(int k, int dimension, int w, int tablesize);
        ~GFunction(){
            r->clear();
            delete r;
            for(int j = 0; j < k; ++j) {
                delete h[j];
            }
            delete[] h;
        }
        vector<int> * getR() {
            return r;
        }
        HFunction**  getH() {
            return h;
        }
        void setV(vector<int>* r) {
            GFunction::r = r;
        }
        void setV(HFunction** h) {
            GFunction::h = h;
        }
        int getBucket_Num(long long int g);
        long long int getNum(vector<TP> &p);
};


template <typename T>
class HashNode {
	private:
        //unique key
    	string key;
    	vector<T> *p;
        long long int g;

	public:
    	HashNode(string &key, long long int &g, vector<T> *p) : key(key), g(g) {
            HashNode::p = p;
        }
    	string getKey() {
        	return key;
    	}
        long long int getG() {
            return g;
        }
    	vector<T> * getValue() {
        	return p;
    	}
    	void setKey(string key) {
    		HashNode::key = key;
    	}
        void setG(long long int g) {
            HashNode::g = g;
        }
    	void setValue(vector<T>* p) {
        	HashNode::p = p;
    	}
};


class HashBucket {
	private:
		list< HashNode<TP>* > *hashlist;
	public:
        HashBucket(){
            hashlist = new list< HashNode<TP>* >();
        }
        ~HashBucket(){
            for(list<HashNode<TP>*> ::const_iterator j = hashlist->begin(); j != hashlist->end(); ++j) {
                delete *j;
            }
            hashlist->clear();
            delete hashlist;
        }
		list< HashNode<TP>* >* getList() {
        	return hashlist;
    	}
    	void setList(list< HashNode<TP>* >* list) {
    		HashBucket::hashlist = list;
    	}
        void add_element(HashNode<TP>* node){
            hashlist->push_back(node);
        }
        void print_list();
};

class HashTable{
    private:
        HashBucket **array_of_buckets;
        int size;
    public:
        HashTable(int &size);
        ~HashTable(){
            for (int i = 0; i < size; ++i){
                delete array_of_buckets[i];
            }
            delete[] array_of_buckets;
        }
        void add_item(int key, string &id, long long int g, vector<TP>* v);

        list<HashNode<TP> *> * return_bucket(int key){
            return array_of_buckets[key]->getList();
        }
        int return_size(){
            return size;
        }
        void print_table();
};

#endif