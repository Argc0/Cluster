#ifndef __CUBE__
#define __CUBE__

#include "hash.h"
//kwdikas ths 1hs ergasias

//klasi gia ton yperkuvo

class Bcube{
    private:
        HashBucket **buckets;
        int d;
        HFunctions **hfs;
        int choice;
        int dimension;
        map <int, int> pam;        // empty map container for "throwing coins"
    public:
        Bcube(int &d, int dimension, int w, int &choice);
        ~Bcube();
        int insert_num(vector<TP>& v);
        void add_item(int key, string &id, long long int g, vector<TP>* v);
        list<HashNode<TP> *> * return_bucket(int key){
            return buckets[key]->getList();
        }
        int return_size_bucket(){
            return int(pow(2.0, double(d)));
        }
        int return_dimension(){
            return d;
        }
        void print_cube();
};

#endif