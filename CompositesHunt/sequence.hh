#ifndef JK_SEQ_HH
#define JK_SEQ_HH

#include <iostream>
#include <string>
extern "C" {
  #include "seq.h"
}


class Sequence {
public:
  Sequence(std::string fileName) {s = seq_open(fileName.c_str()); read();}
  Sequence(const char* fileName) {s = seq_open(fileName); read();}
  Sequence(const Sequence& seq);
  ~Sequence() {seq_close(s);}


  bool read() {return seq_read(s);}
  
  std::string name() {return SEQ_NAME(s);}
  unsigned size() {return SEQ_LEN(s);}
  std::string header() {return SEQ_HEAD(s);}
  
  uchar* seqPointer() {return s->seq;}

  std::string subsequence(int start, int length);
  void revComplement() {seq_revcomp_inplace(s);}
  
  
  char operator[](int i) {return (char)s->seq[i];}
  Sequence& operator=(const Sequence& seq);
  
private: 
  SEQ* s;
};

  

#endif
