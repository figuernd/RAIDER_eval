#ifndef LMER_VECTOR_H
#define LMER_VECTOR_H
#include <vector>
#include <stdexcept>
#include <iostream>

using namespace std;

#ifndef uint
typedef unsigned int uint;
#endif

class Family;

class LmerVector {
public:
  
  LmerVector(uint position) {
    family = nullptr;
    push_back(position);
  }
  
  ~LmerVector() {
  }
  
  void push_back(uint val) {
    if(lmers.size() > 0){
      previous = back();
    }
    lmers.push_back(val);
  }
  
  uint operator [](uint i) const { return lmers[i]; }
  uint & operator [](uint i)     { return lmers[i]; }
  
  uint prev()  const { return previous; }
  uint back()  const { return lmers.back(); }
  uint front() const { return lmers.front(); }
  
  uint getOff() const { return off; }
  void setOff(uint i) { off = i;    }
  
  uint size() const { return lmers.size(); }
  
  Family* getFamily() const   { return family; }
  void setFamily(Family *fam) {
    if (family) {
      prevFamily = family;
    }
    family = fam;
  }
  
  Family* getPrevFamily() const   { return prevFamily; }
  void setPrevFamily(Family *fam) { prevFamily = fam;  }
  
  friend ostream &operator<<( ostream &output, const LmerVector &L) {
    output  << "(\tFront:\t" << L.front()
    << "\tBack:\t" << L.back()
    << "\tSize:\t" << L.size()
    << "\tOffset In Family\t" << L.getOff()
    << "\tFamily:\t" << L.getFamily() << "\t)";
    return output;
  }
  
  vector<uint> lmers;
  uint previous;
  uint off;
  Family* family;
  Family* prevFamily;
};

#endif //LMER_VECTOR_H
