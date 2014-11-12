#ifndef COMPARE_LMER_VECTOR
#define COMPARE_LMER_VECTOR
#include "LmerVector.h"

using namespace std;

#ifndef uint
typedef unsigned int uint;
#endif

class CompareLmerVector{
    public:
    bool operator()(LmerVector* l1, LmerVector* l2)
    {
       if (l1->size() > l2->size()) return true;
       if (l1->size() == l2-> size() && l1->front() < l2->front()) return true;
       return false;
    }
};

#endif //COMPARE_LMER_VECTOR
