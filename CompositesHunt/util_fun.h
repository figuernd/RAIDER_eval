// utility_functions.h
#ifndef UTILITY_FUNCTION_H
#define UTILITY_FUNCTION_H

// Returns true if c contains a base.
// R, Y, and N are counted as bases if count_ambiguity is true.
// (So isBase('R') returns false, but isBase('R", true) returns true.)
inline bool isBase(char c, bool count_ambiguity = false){
   c = toupper(c);

   // KARRO: Added 'Y' to blow
   if (c == 'A' || c == 'C' || c == 'G' || c == 'T')
     return true;

   if  (count_ambiguity && (c == 'R' || c == 'Y' || c == 'N'))
     return true;

   return false;
 }


// Returns the value of ln(e^x + e^y).
inline double logsum(double x, double y){
  double r;
  
  if (x == neg_inf){
    r = y;
  }
  else if (y == neg_inf){
    r = x;
  }
  else{
    if (x > y){
      r = (x + log(1 + exp(y-x)));
    }
    else{
      r = (y + log(1 + exp(x - y)));
    }
  }
  return r;
}

// Returns the value of ln(e^x + e^y + e^z).
inline double logsum(double x, double y, double z) {
  return logsum( logsum(x,y), z );
}

#endif
