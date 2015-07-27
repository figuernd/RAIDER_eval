#ifndef RAIDER2_SEEDCHAIN_H
#define RAIDER2_SEEDCHAIN_H

#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <stdexcept>
#include <vector>
#include "Family.h"
#include <unordered_map>
#include <ctime>
#include <string>

using namespace std;

typedef unordered_map<size_t, LmerVector*> LmerMap;

int p = getpid();

void prettyPrintTabbing(uint t){
  for(uint i = 0; i < t; i++){
    cout << '\t';
  }
}

void prettyPrintLmer(int t, LmerVector* v){
  prettyPrintTabbing(t);
  cout << "Lmer:\t" << *v << endl << endl;
}


void prettyPrintFamily(uint t, Family* fam, bool isNew, bool isMod){
  prettyPrintTabbing(t);
  if (isNew) {
    cout << "New family:\t" << *fam << endl << endl;
  }
  else if (isMod) {
    cout << "Modified family:\t" << *fam << endl;
  }
  else{
    cout << "Family:\t" << *fam << endl;
  }
}

void prettyPrintMethodState(uint t, string methodname, LmerVector* v, Family* fam, bool isNew, bool isMod){
  prettyPrintTabbing(t);
  cout << "--- " << methodname << " ---" << endl;
  prettyPrintLmer(t + 1, v);
  prettyPrintFamily(t + 1, fam, isNew, isMod);
}


// Determines numeric representation of base. Only works for meaningful bases
// ('A', 'C', 'G', 'T'). Helper for hashing function.
size_t baseToInt(char b) {
  switch (toupper(b)) {
    case 'A':
      return 0;
    case 'C':
      return 1;
    case 'G':
      return 2;
    case 'T':
      return 3;
    default:
      throw std::out_of_range("Unknown base type");
      return -1;
  }
}


// Hashes seed (char array), considering only matching ('1') positions in mask.
size_t seedToInt(const char* seed, const seqan::CharString &mask) {
  //TODO use bitset or vector<bool> to allower for potentially large bitsets
  size_t result = 0;
  uint L = seqan::length(mask);
  for (uint i = 0; i < L; i++) {
    if (mask[i] == '1') {
      result <<= 2;
      result |= baseToInt(seed[i]);
    }
  }
  return result;
}


// Determines whether there is another Lmer remaining in the sequence at this index.
// Modifies seed (char array) to contain the next L meaningful characters of the
// sequence, filtering out 'N' characters.
bool getNextSeed(seqan::Dna5String &sequence, uint &index, uint seqLength, uint L, char* seed) {
  seed[L - 1] = 0;
  for (uint i = 0; i < L && index + i < seqLength; i++) {
    if (toupper((char) sequence[index + i]) == 'N') {
      index = index + i + 1;
      i = -1;
    } else {
      seed[i] = sequence[index + i];
    }
  }
  return seed[L - 1];
}


// Returns lmervector object that the seed maps to in lmers (hash map).
// If seed has never been seen before, creates new lmervector object.
// Adds most recent index to the lmervector's list of locations.
LmerVector* getAndInsert(LmerMap &lmers, size_t seed, uint index) {
  LmerVector*& v = lmers[seed];
  if (v == 0) {
    v = new LmerVector(index);
  } else {
    v->push_back(index);
  }
  return v;
}


// Determines whether lmervector v belongs to family fam.
bool isNewFamily(Family* fam, LmerVector* v, bool overlaps) {
  // Make sure that fam is a viable family of size exactly 2
  assert(!v->getFamily());
  if (!fam || fam->size() > 2) {
    return true;
  }
  assert(fam->size() == 2);
  
  LmerVector* suffix = fam->getSuffix();
  uint s1 = suffix->front();
  uint s2 = suffix->back();
  uint prev = v->front();
  uint index = v->back();
  
  // Ensure both lmers' location pairs are the same distance apart
  if (index - prev != s2 - s1) {
    return true;
  }
  // Ensure v's most recent location is less than L from the last location added to fam
  // (or equal to L if overlaps not required)
  if (index > fam->getExpectedEnd() || (overlaps && index == fam->getExpectedEnd())) {
    return true;
  }
  return false;
}


// Returns the index of the family in curr_fams to which lmervector v belongs.
// If there is no such family, returns L to signal that none exists.
uint getFamilyIndex(LmerVector* v, Family* curr_fams[], int curr_index, int L, bool overlaps) {
  assert(!v->getFamily());
  // Go through families in array from most recent to least (left from curr_index)
  for(int i = 0; i < L; i++){
    // Determine family index (allows for array wrap around)
    int fam_index = (curr_index - i < 0)? L + (curr_index-i) : curr_index - i;
    Family* fam = curr_fams[fam_index];
    
    if (fam) {
      if (isNewFamily(fam, v, overlaps))
        continue;
      else
        return fam_index;
    }
    else {
      return L;
    }
  }
  return L;
}


Family* exciseRepeatsByLmer(Family* fam, LmerVector *v, uint L, uint verbosity) {
  if(verbosity > 2) prettyPrintMethodState(1, "Excise Repeats by Lmer", v, fam, false, false);
  
  assert(fam == v->getFamily());
  vector<LmerVector*> *lmers = fam->getLmers();
  LmerVector* currLast = fam->getLast();
  Family *newFam = new Family;
  uint numRemoved = 0;
  uint oldLength = lmers->size();
  uint offset = oldLength;
  vector<LmerVector*> *toKeep = new vector<LmerVector*>;
  
  for (uint i = 0; i < oldLength; i++) {
    LmerVector *u = lmers->at(i);
    if (u == v) {
      assert(i != 0);
      offset = i;
      newFam->adopt(u, L);
      numRemoved++;
    }
    if (i > offset){
      if (u->size() == v->size() - 1 && u->back() <= newFam->getLast()->back() + L ) {
        newFam->adopt(u, L);
        numRemoved++;
      }
      else {
        toKeep->push_back(u);
      }
    }
  }
  
  lmers->resize(offset);
  lmers->shrink_to_fit();
  for(uint i = 0; i < toKeep->size(); i++){
    fam->adopt(toKeep->at(i), L);
  }
  fam->setLast(currLast);
  fam->setExpectedEnd(fam->getPrefix()->back() + fam->repeatLength(L));
  fam->setOffset(fam->getSuffix()->front() - fam->getPrefix()->front());
  fam->setRepeatLength(L + fam->getOffset());
  
  newFam->setLast(v);
  newFam->setExpectedEnd(v->back() + newFam->repeatLength(L));
  
  if (verbosity > 2){
    prettyPrintFamily(tabbing+1, fam, false, true);
    prettyPrintFamily(tabbing+1, newFam, true, false);
  }
  return newFam;
  
  
  
}

Family* splitRepeatsByLmer(Family* fam, LmerVector *v, bool keepV, uint L, uint verbosity, bool excising) {
  if(verbosity > 2) prettyPrintMethodState(1, "Split Repeats by Lmer", v, fam, false, false);
  
  assert(fam == v->getFamily());
  vector<LmerVector*> *lmers = fam->getLmers();
  Family *newFam = new Family;
  
  uint oldLength = lmers->size();
  uint offset = oldLength;
  
  // Go through the lmers of the family until find v.
  // Move all lmers following v (and v if keepV is false) into a new family
  for (uint i = 0; i < oldLength; i++) {
    LmerVector *u = lmers->at(i);
    
    if (u == v) {
      // assert(i != 0); // CARLY: commented this out because it seems incorrect...
      offset = i;
    }
    
    if ((!keepV && i >= offset) || i > offset) {
      newFam->adopt(u, L);
    }
  }
  
  // Resize lmers in original family to new size. Take back memory.
  lmers->resize(keepV ? offset + 1 : offset);
  lmers->shrink_to_fit();
  
  // Reset last and expected end of original family
  fam->setLast(fam->getSuffix());
  fam->setExpectedEnd(fam->getSuffix()->back() + L);
  
  //TODO: shouldn't this be everywhere???
  if (excising) {
    fam->setOffset(fam->getSuffix()->front() - fam->getPrefix()->front());
    fam->setRepeatLength(L + fam->getOffset());
  }
  
  // Reset last and expected end of new family if it will contain v
  if (!keepV) {
    newFam->setLast(v);
    newFam->setExpectedEnd(v->back() + newFam->repeatLength(L));
  }
  
  if (verbosity > 2){
    prettyPrintFamily(tabbing+1, fam, false, true);
    prettyPrintFamily(tabbing+1, newFam, true, false);
  }
  
  return newFam;
}

bool repeatExpects(Family *fam, LmerVector* v, uint L, uint verbosity, bool excising) {
  if(verbosity > 2) prettyPrintMethodState(1, "Repeat Expects", v, fam, false, false);

  const uint index = v->back();
  const uint length = fam->repeatLength(1);
  uint dist_to_end = fam->getExpectedEnd() - index;
  
  // The last location of v cannot be beyond the expected end of the family and
  // the distance between the two cannot exceed the repeat length
  if (index > fam->getExpectedEnd() || dist_to_end > length) {
    return false;
  }

  // If excising and the previous location of v is before the last index of the family
  //else if (excising && v->prev() < fam->getLastIndex()){
  //  if (verbosity > 2){
  //    prettyPrintTabbing(tabbing + 1);
  //    cout << "The previous location of lmer is before the last index of the family." << endl;
  //  }
  //  return false;
  //}

  // The lmer's alleged position in the family can be calculated. It must match
  // the lmer found in the family at this position
  const int pos = length - dist_to_end;
  if (pos >= 0 && pos < (int) length && fam->findByPos(pos) == v) {
    return true;
  }
  return false;
}


bool fragmentSplit(LmerVector* v, uint L, vector<Family*> &families, uint verbosity, bool excising) {
  Family* fam = v->getFamily();
  assert(fam);
  if(verbosity > 2) prettyPrintMethodState(0, "Fragment Split", v, fam, false, false);
  
  // If expecting to start a new repeat instance for this family, this lmer must
  // be the family prefix. If not, split family into prefix...(v-1) and v...suffix
  if (fam->lastRepeatComplete()) {
    if (v != fam->getPrefix()) {
      Family* newFam = splitRepeatsByLmer(fam, v, false, L, verbosity, excising);
      families.push_back(newFam);
      return true;
    }
  }
  
  // If v occurs after it is supposed to (over L past the last index of the family),
  // split the family into prefix---last, (last+1)---suffix
  else if (v->back() > fam->getLastIndex() + L) {
    // **************FROM NEW:
    //// if v was placed in newFamily but should have stayed where it was
    //if (v->back() <= prevFam->getLastIndex() + L && v->size() == prevFam->size()) {
    //    if (verbosity > 2){
    //        prettyPrintTabbing(tabbing + 1);
    //        cout << "The last location of the lmer is over L (" << L << ") past the last index of new family, but not prev. Revert." << endl;
    //    }
    //    v->setFamily(prevFam);
    //    return false;
    //}
    Family* newFam = splitRepeatsByLmer(fam, fam->getLast(), true, L, verbosity, excising);
    families.push_back(newFam);
    //NOT IN OLD OR OLDEST
    //Family* newFam2 = splitRepeatsByLmer(newFam, newFam->getLast(), false, L, verbosity, excising);
    //families.push_back(newFam2);
    return true;
  }
  
  // if it was previously skipped in forming a family instance
  else if (excising && v->prev() < fam->getLast()->prev() && v->size() < fam->size() ) {
    if (verbosity > 2){
      prettyPrintTabbing(tabbing + 1);
      cout << "The lmer was skipped in a previous iteration" << endl;
    }
    // OLDEST (1):
    // Family* newFam = splitRepeatsByLmer(fam, v, false, L, verbosity, excising);
    Family* newFam = exciseRepeatsByLmer(fam, v, L, verbosity);
    families.push_back(newFam);
    return true;
  }
  
  // OLDEST
  // if it occurs after it is supposed to
  //else if (!excising && v->prev() < fam->getLast()->prev()) {
  //    if (verbosity > 2){
  //        prettyPrintTabbing(tabbing + 1);
  //        cout << "The previous location of the lmer is before the last index of the family." << endl;
  //    }
  
  //    Family* newFam = splitRepeatsByLmer(fam, v, false, L, verbosity, excising);
  //    families.push_back(newFam);
  //    return true;
  //}
  return false;
}

uint maxLength(const vector<seqan::CharString> &masks) {
  uint max = 0;
  for (uint i = 0; i < masks.size(); i++) {
    uint len = seqan::length(masks[i]);
    if (len > max) {
      max = len;
    }
  }
  return max;
}

bool isPrefix(LmerVector *v) {
  return v->getFamily()->getPrefix() == v;
}

void tieLooseEnds(vector<Family*> &families, uint L, uint verbosity, bool excising) {
  if (verbosity > 2){
    cout << "--- Tying Loose Ends ---" << endl;
  }
  for (auto fam : families) {
    if (!fam->lastRepeatComplete()) {
      Family* newFam = splitRepeatsByLmer(fam, fam->getLast(), true, L, verbosity, excising);
      families.push_back(newFam);
    }
  }
}


void getElementaryFamilies(seqan::Dna5String &sequence, vector<seqan::CharString> &masks, vector<Family*> &families, int verbosity, bool family_array, bool excising, bool overlaps, bool tieup) {
  
  const uint seqLength = seqan::length(sequence);
  const seqan::CharString mask = masks.front();
  const uint L = seqan::length(mask);
  const uint MAX_L = maxLength(masks);
  
  uint index = -1;
  char unmasked[MAX_L];
  LmerMap lmers;
  Family* fam = nullptr;
  
  // Only used if family_array switch set
  Family* curr_families[L];
  for(uint i = 0; i < L; i++){
    curr_families[i] = NULL;
  }
  int curr_fam = 0;
  int signed_L = seqan::length(mask);
  
  while (getNextSeed(sequence, ++index, seqLength, MAX_L, unmasked)) {
    size_t seed = seedToInt(unmasked, mask);
    LmerVector* v = getAndInsert(lmers, seed, index);
    
    // If lmer has exactly 2 occurrences, it is a potential elementary repeat.
    // Needs to be placed with a family (either existing or new).
    if (v->size() == 2) {
      // If family_array enabled, search for corresponding family within L of this instance.
      if(family_array){
        int fam_index = getFamilyIndex(v, curr_families, curr_fam, signed_L, overlaps);
        if (fam_index == signed_L) {
          fam = new Family();
          curr_fam = (curr_fam + 1) % L;
          curr_families[curr_fam] = fam;
          families.push_back(fam);
        }
        else {
          fam = curr_families[fam_index];
        }
      }
      // Else, check to see if this lmer is a part of most recent family
      else{
        if (isNewFamily(fam, v, overlaps)) {
          fam = new Family();
          families.push_back(fam);
        }
      }
      fam->adopt(v, L);
    }
    
    // If lmer has more than 2 occurrences, it is already associated with a family.
    // Check that this latest occurrence does not violate defn of elementary repeat.
    else if (v->size() > 2 && !fragmentSplit(v, L, families, verbosity, excising)) {
      Family* fam = v->getFamily();
      
      if (isPrefix(v)) {
        fam->setLast(v);
        fam->setExpectedEnd(v->back() + fam->repeatLength(L));
      }
      
      else if (repeatExpects(fam, v, L, verbosity, excising)) {
        fam->setLast(v);
      }
    }
  }
  
  tieLooseEnds(families, L, verbosity, excising);
}

#endif //RAIDER2_SEEDCHAIN_H
