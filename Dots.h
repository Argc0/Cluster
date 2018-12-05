#ifndef __DOTS__
#define __DOTS__

//domi gia apothikeusi pliroforiwn twn simeiwn kai se poio cluster anikoun

class Dot{

  private:
    string id;
    //to simio se poio (deikti) cluster anikei
    int index_cluster;
    //an to simeio einai markarismeno
    int marked;
    vector<TP> *point;
  
  public:
    Dot(string& id, vector<TP> *v){
      Dot::id = id;
      Dot::index_cluster = -1;
      Dot::point = v;
    }
    void setIndex_Cluster(int index){
      Dot::index_cluster = index;
    }
    void setMarked(int marked){
      Dot::marked = marked;
    }
    int getIndex_Cluster(){
      return index_cluster;
    }
    int getMarked(){
      return marked;
    }
    string getId(){
      return id;
    }
    vector<TP> * getPoint(){
      return point;
    }
};

class Dots{
  
  private:
    Dot ** dots;
    //voithitiki domi gia tin evresi enos simeiou me key
    map<string, int> hmap;
    int size;
  
  public:
    Dots(vector< vector<TP> > *points, vector<string> * in_ids, int& n){
      dots = new Dot*[n]();
      for(int i=0; i<n; i++){
        dots[i] = new Dot(in_ids[0][i], &points[0][i]);
        hmap[in_ids[0][i]]=i;
      }
      size=n;
    }
    ~Dots(){
      hmap.clear();
      for(int i=0; i<size; i++){
        delete dots[i];
      }
      delete[] dots;
    }
    void clear_clusters(){
      for(int i=0; i<size; i++){
        dots[i]->setIndex_Cluster(-1);
      }
    }
    void refreshMarked(){
      for(int i=0; i<size; i++){
        dots[i]->setMarked(0);
      }
    }
    void setIndexforPoint(int i, int index){
      dots[i]->setIndex_Cluster(index);
    }
    void setIndexbyKey(string key, int index){
      dots[hmap[key]]->setIndex_Cluster(index);
    }
    void setMark(string& key, int m){
      dots[hmap[key]]->setMarked(m);
    }
    int getMark(string& k){
      return dots[hmap[k]]->getMarked();
    }
    int get_Marked(int i){
      return dots[i]->getMarked();
    }
    int getSize(){
      return size;
    }
    vector<TP> * returnPoint(int i){
      return dots[i]->getPoint();
    }
    string returnIdPoint(int i){
      return dots[i]->getId();
    }
    int returnIndex(int j){
      return dots[j]->getIndex_Cluster();
    }
    int returnIndexbyKey(string key){
      return dots[hmap[key]]->getIndex_Cluster();
    }
    void print_clusters(int& n, int num_clusters){
      for(int i=0; i<num_clusters; i++){
        cout << "Cluster " << i+1 << endl;
        for(int j=0; j<n; j++){
          if(i == dots[j]->getIndex_Cluster()){
            cout << dots[j]->getId() << endl;
          }
        }
      }
    }

};

#endif