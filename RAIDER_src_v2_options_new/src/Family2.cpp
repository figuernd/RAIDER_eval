#include "Family.h"

void Family::adopt(LmerVector *v) {
  v->setFamily(this);
  vectors.push_back(v);
  setLast(v);
  setExpectedEnd(v->back() + 1);
}

void Family::adopt(LmerVector *v, uint L){
  v->setFamily(this);
  //cout << L << endl;
  if (vectors.size() > 0){
    uint off = v->front() - vectors.back()->front();
    addToOffset(off - 1);
    setRepeatLength(getRepeatLength() + off);
    v->setOff(off - 1);
  }
  else{
    setRepeatLength(L);
    setOffset(0);
    v->setOff(0);
  }
  vectors.push_back(v);
  setLast(v);
  setExpectedEnd(v->back() + L);
}

uint Family::size() const {
  if (vectors.size() == 0) {
    return 0;
  }
  return vectors.front()->size();
}

void Family::addToOffset(uint off){
  offset += off;
}

uint Family::repeatLength(uint L) const {
  return repeatlength + L - L;
}

LmerVector* Family::findByPos(uint pos){
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

void Family::setLast(LmerVector* v) {
  assert(v->getFamily() == this);
  last = v;
  last_index = v->back();
}

bool Family::lastRepeatComplete() const {
  return last == getSuffix();
}

std::ostream &operator<<( std::ostream &output, const Family &F) {
  output  << "Size: " << F.vectors.size()
          << " Repeat Length: " << F.repeatlength
          << " Expected End: " << F.expected_end
          << " Offset: " << F.offset << std::endl
          << "\t\tLast:\t" << *(F.last) << std::endl
          << "\t\tPrefix:\t" << *(F.getPrefix()) << std::endl
          << "\t\tSuffix:\t" << *(F.getSuffix()) << std::endl;
  for(uint i = 0; i < F.vectors.size(); i++){
    output << "\t\t\t"<< i << ":\t" << *(F.see(i)) << std::endl;
  }
  return output;
}
