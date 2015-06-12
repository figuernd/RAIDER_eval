// interval_list.h
#include <list>
#include <iostream>

typedef std::pair<int,int> interval;

class interval_list {
 public:
  interval_list() {;};
  void add(const interval& I);
  void add(int s, int f) {add(std::make_pair(s,f));}
  int num_covered();
  void print(std::ostream& out);

 private:
  std::list<interval> ilist;
  
};

std::ostream operator>>(std::ostream& out, const interval_list& IL);
