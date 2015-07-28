#ifndef FAMILY_H
#define FAMILY_H
#include <vector>
#include "LmerVector.h"
#include <iostream>
#include <assert.h>
#include <algorithm>

class Family {
public:
  void adopt(LmerVector *v) {
    v->setFamily(this);
    vectors.push_back(v);
    setLast(v);
    setExpectedEnd(v->back() + 1);
  }
  
  void adopt(LmerVector *v, uint L){
    v->setFamily(this);
    if (vectors.size() > 0){
      uint off = v->front() - vectors.back()->front();
      v->setOff(off - 1);
    }
    else{
      v->setOff(0);
    }
    vectors.push_back(v);
    setLast(v);
    setExpectedEnd(v->back() + L);
  }
  
  uint size() const {
    if (vectors.size() == 0) {
      return 0;
    }
    return vectors.front()->size();
  }
  
  uint repeatLength(uint L) const {
    return vectors.back()->front() - vectors.front()->front() + L;
  }
  
  vector<LmerVector*>* getLmers() { return &vectors; }
  
  LmerVector* at(uint index) { return vectors.at(index); }
  LmerVector* see(uint index) const { return vectors.at(index); }
  void push_back(LmerVector* v) { vectors.push_back(v); }
  
  LmerVector* findByPos(uint pos){
    uint sumOff = 0;
    if (pos == 0){
      return vectors.front();
    }
    for (LmerVector* v : vectors){
      sumOff += v->getOff();
      if (pos == sumOff){
        return v;
      }
      sumOff++;
    }
    return vectors.front();
  }
    
  LmerVector* getPrefix() const { return vectors.front(); }
  LmerVector* getSuffix() const { return vectors.back();  }
  
  LmerVector* getLast() { return last; }
  uint getLastIndex() const { return last_index; }
  bool lastRepeatComplete() const { return last == getSuffix(); }

  void setLast(LmerVector* v) {
    assert(v->getFamily() == this);
    last = v;
    last_index = v->back();
  }
  
  LmerVector* getOneAfterLast(){
    std::vector<LmerVector*>::iterator start = std::find(vectors.begin(), vectors.end(), last);
    return *(start+1);
  }
  
  uint getExpectedEnd() const { return expected_end; }
  void setExpectedEnd(uint expected) { expected_end = expected; }
  
  std::vector<LmerVector*>* getSkipped() { return &skipped; }

  std::vector<LmerVector*> popSkipped() {
    std::vector<LmerVector*> ret(skipped);
    skipped.clear();
    return ret;
  }
  
  void setRemainingSkipped(){
    std::vector<LmerVector*>::iterator start = std::find(vectors.begin(), vectors.end(), last);
    std::move(start+1, vectors.end(), std::back_inserter(skipped));
    vectors.erase(start+1, vectors.end()); // no longer part of backbone of this family
    vectors.shrink_to_fit();
  }
  
  
  void moveSkippedRange(LmerVector* lastSeen){
    std::vector<LmerVector*>::iterator start = std::find(vectors.begin(), vectors.end(), last);
    std::vector<LmerVector*>::iterator end = std::find(vectors.begin(), vectors.end(), lastSeen);
    std::move(start+1, end, std::back_inserter(skipped));
    vectors.erase(start+1, end); // no longer part of backbone of this family
    vectors.shrink_to_fit();
    lastSeen->setOff(lastSeen->front() - last->front());  // adjust offset
  }
  
  void removeOne(LmerVector* v){
    std::vector<LmerVector*>::iterator start = std::find(vectors.begin(), vectors.end(), v);
    (*(start+1))->setOff((*(start+1))->front() - (*(start-1))->front());
    vectors.erase(start, start+1, vectors.end()); // no longer part of backbone of this family
    vectors.shrink_to_fit();
  }

  
  friend std::ostream &operator<< (std::ostream &output, const Family &F){
    output  << "Size: " << F.vectors.size()
    << " Expected End: " << F.expected_end
    << "\t\tLast:\t" << *(F.last) << std::endl
    << "\t\tPrefix:\t" << *(F.getPrefix()) << std::endl
    << "\t\tSuffix:\t" << *(F.getSuffix()) << std::endl;
    for(uint i = 0; i < F.vectors.size(); i++){
      output << "\t\t\t"<< i << ":\t" << *(F.see(i)) << std::endl;
    }
    return output;
  }
  
  LmerVector* last;
  uint last_index;
  uint expected_end;
  vector<LmerVector*> vectors;
  std::vector<LmerVector*> skipped;
};

#endif //FAMILY_H