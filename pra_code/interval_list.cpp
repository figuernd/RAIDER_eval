#include "interval_list.h"
#include <algorithm>

void interval_list::add(const interval& I) {
  if (ilist.size() == 0) {
    ilist.push_front(I);
    return;
  }

  std::list<interval>::iterator i = ilist.begin();
  while (i != ilist.end() && i->second <= I.first) // Find first element ending after I starts
    i++;
  
  if (i == ilist.end() || i->first > I.second) {  // If we pased it
    ilist.insert(i,I);
    return;
  }

  std::list<interval>::iterator j = i;
  int max_val = I.second;
  while (j != ilist.end() && j->first < I.second) {
    max_val = std::max(max_val, j->second);
    j++;
  }
  
  interval new_I = std::make_pair( std::min(I.first, i->first), max_val );
  i = ilist.erase(i, j);
  ilist.insert(i, new_I);

}

int interval_list::num_covered() {
  return std::accumulate(ilist.begin(), ilist.end(), 0, [](int s, interval I) {return s + I.second - I.first;});
}

    
void interval_list::print(std::ostream& out) {
  for (std::list<interval>::iterator i = ilist.begin(); i != ilist.end(); i++)
    out << "(" << i->first << "," << i->second << ") ";
  out << std::endl;
}
