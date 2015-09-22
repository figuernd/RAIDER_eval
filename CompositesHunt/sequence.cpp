#include "sequence.hh"

using namespace std;

Sequence::Sequence(const Sequence& seq) {
  this->s = seq_copy(s);
}

string Sequence::subsequence(int start, int length) {
  int begin = MAX(start, 0);
  int end = MIN(start+length, size());
  char c = s->seq[end];
  s->seq[end] = '\0';
  string result = (char*)(s->seq + begin);
  s->seq[end] = c;
  return result;
}

Sequence& Sequence::operator=(const Sequence& seq) {
  if (this != &seq) {
    this->s = seq_copy(s);
  } 
  return *this;
}


