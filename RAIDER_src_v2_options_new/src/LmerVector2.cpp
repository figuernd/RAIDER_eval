#include "LmerVector.h"

LmerVector::LmerVector(uint position) {
  family = nullptr;
  push_back(position);
}

void LmerVector::push_back(uint val) {
  if(lmers.size() > 0){
    previous = back();
  }
  lmers.push_back(val);
}

std::ostream &operator<<( std::ostream &output, const LmerVector &L) {
  output  << "(\tFront:\t" << L.front()
          << "\tBack:\t" << L.back()
          << "\tSize:\t" << L.size()
          << "\tOffset In Family\t" << L.getOff()
          << "\tFamily:\t" << L.getFamily() << "\t)";
  return output;
}
