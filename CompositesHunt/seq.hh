#ifndef JK_SEQ_HH
#define JK_SEQ_HH

#include <string>
extern "C" {
  #include "seq.h"
}


class Sequence {
public:
  Sequence(std::string fileName) {s = seq_open(fileName.c_str());}
  Sequence(char* fileName) {s = seq_open(fileName);}
  Sequence(const Sequence& seq);
  ~Sequence() {seq_close(s); read();}


  bool read() {return seq_read(s);}
  
  std::string name() {return SEQ_NAME(s);}
  unsigned size() {return SEQ_LEN(s);}
  std::string header() {return SEQ_HEAD(s);}
  
  uchar* seqPointer() {return s->seq;}

  std::string subsequence(int start, int length);
  void rev_comp() {seq_revcomp_inplace(s);}
  
  
  uchar operator[](int i) {return s->seq[i];}
  Sequence& operator=(const Sequence& seq);
  

private: 
  SEQ* s;
};

#endif
