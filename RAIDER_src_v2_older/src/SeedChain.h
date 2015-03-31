#ifndef SCANER_SEEDCHAIN
#define SCANER_SEEDCHAIN

#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <stdexcept>
#include <vector>
#include "Family.h"
#include <unordered_map>
#include <ctime>

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

// NOT IN OLDEST (func)
uint getFamilyIndex(LmerVector* v, Family* curr_fams[], int curr_index, int L) {
	assert(!v->getFamily());
    for(int i = 0; i < L; i++){ 
        int fam_index = (curr_index - i < 0)? L + (curr_index-i) : curr_index - i;    
        //fam_index = fam_index < 0? fam_index * -1 : fam_index;
        Family* fam = curr_fams[fam_index];
        // if array has no family here, went through all viable families and found none.
        if(!fam || fam->size() > 2){
            cout << "No suitable family" << endl;
            return L;
        }
        
        assert(fam->size() == 2);
	    LmerVector* suffix = fam->getSuffix();
	    uint s1 = suffix->front();
	    uint s2 = suffix->back();
	    uint prev = v->front();
	    uint index = v->back();
	    // Ensure both Lmer pairs are the same distance apart
	    if (index - prev != s2 - s1) {
		    continue;
	    }
        // Ensure second Lmer pair occurs within L of first Lmer pair
        if (index > fam->getExpectedEnd()) {
            continue;
        }
        cout << "Found suitable family" << endl;
        return fam_index;
    }


    cout << "No suitable family" << endl;
	return L;
}

bool isNewFamily(Family* fam, LmerVector* v) {
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
	// Ensure both Lmer pairs are the same distance apart
	if (index - prev != s2 - s1) {
		return true;
	}
	// Ensure first Lmer pair overlaps second Lmer pair
	// CES: took this out because spaced seed version should not require 
    //if (index != s2 + 1) {
	//	return true;
	//}
    // Ensure second Lmer pair occurs within L of first Lmer pair
    if (index > fam->getExpectedEnd()) {
        return true;
    }
	return false;
}

LmerVector* getAndInsert(LmerMap &lmers, size_t seed, uint index) {
	LmerVector*& v = lmers[seed];
	if (v == 0) {
		v = new LmerVector(index);
	} else {
		v->push_back(index);
	}
	return v;
}

// NOT IN OLDEST (func)
Family* exciseRepeatsByLmer(Family* fam, LmerVector *v, uint L, uint verbosity) {
	uint tabbing = 1;
    if(verbosity > 2){
        prettyPrintTabbing(tabbing);
        cout << "--- Excise Repeats by Lmer ---" << endl;
        prettyPrintLmer(tabbing + 1, v);
        prettyPrintFamily(tabbing + 1, fam, false, false);
    }
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
    //lmers->erase(lmers->begin() + offset, lmers->begin() + offset + newFam->getLmers()->size());
    //cout << "\tResize to size\t" << offset;
    //for (uint i = oldLength - 1; i >= offset; i--){
    //    cout << i <<"\t\t" << *fam << endl;
    //    lmers->pop_back();
    //}
    
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
	uint tabbing = 1;
    if(verbosity > 2){
        prettyPrintTabbing(tabbing);
        cout << "--- Split Repeats by Lmer ---" << endl;
        prettyPrintLmer(tabbing + 1, v);
        prettyPrintFamily(tabbing + 1, fam, false, false);
    }
    assert(fam == v->getFamily());
	vector<LmerVector*> *lmers = fam->getLmers();
	Family *newFam = new Family;

	uint oldLength = lmers->size();
	uint offset = oldLength;

	for (uint i = 0; i < oldLength; i++) {
		LmerVector *u = lmers->at(i);
		if (u == v) {
			// CS: commented out because seems wrong...
            assert(i != 0); // so prefix never split?
			offset = i;
		}
		if ((!keepV && i >= offset) || i > offset) {
            newFam->adopt(u, L);
		}
	}
	lmers->resize(keepV ? offset + 1 : offset);
	lmers->shrink_to_fit();

	fam->setLast(fam->getSuffix());
	fam->setExpectedEnd(fam->getSuffix()->back() + L);
    // NOT IN OLDEST (2)
    if (excising) {
        fam->setOffset(fam->getSuffix()->front() - fam->getPrefix()->front());
        fam->setRepeatLength(L + fam->getOffset());
    }
	if (!keepV) {
		newFam->setLast(v);
        // set end to the last location of v + length of repeats
		newFam->setExpectedEnd(v->back() + newFam->repeatLength(L));
	}
    if (verbosity > 2){
        prettyPrintFamily(tabbing+1, fam, false, true);
        prettyPrintFamily(tabbing+1, newFam, true, false);
	}
    return newFam;
}

bool repeatExpects(Family *fam, LmerVector* v, uint L, uint verbosity, bool excising) {
    uint tabbing = 1;
    if(verbosity > 2){
        prettyPrintTabbing(tabbing);
        cout << "--- Repeat Expects ---" << endl;
        prettyPrintLmer(tabbing + 1, v);
        prettyPrintFamily(tabbing + 1, fam, false, false);
    }
    const uint index = v->back(); 
	const uint length = fam->repeatLength(1);
	uint dist_to_end = fam->getExpectedEnd() - index;
    if (index > fam->getExpectedEnd() || dist_to_end > length) {
        if (verbosity > 2){
		    prettyPrintTabbing(tabbing + 1);
            cout << "Last location of lmer is too far past the expected end of family OR the distance to end is greater than repeat length" << endl;
        }
        return false;
	}
    // IN OLDEST (else)
    else if (excising && v->prev() < fam->getLastIndex()){
        if (verbosity > 2){
            prettyPrintTabbing(tabbing + 1);
            cout << "The previous location of lmer is before the last index of the family." << endl;
        }
        return false;
    }
	// pos based on whether excising (OLD) or not excising (OLDEST)
    const int pos = excising? v->getPosition() : length - dist_to_end - L;
    if (pos >= 0 && pos < (int) length && fam->at(pos) == v) {
        if (verbosity > 2){
            prettyPrintTabbing(tabbing + 1);
            cout << "The lmer matches family lmer at position " << pos << endl;
		}
        return true;
	}
    return false;
}


bool fragmentSplit(LmerVector* v, uint L, vector<Family*> &families, uint verbosity, bool excising) {
    Family* fam = v->getFamily();
    assert(fam);
    uint tabbing = 0;
    if(verbosity > 2){
        prettyPrintTabbing(tabbing);
        cout << "--- Fragment Split ---" << endl;
        prettyPrintLmer(tabbing + 1, v);
        prettyPrintFamily(tabbing + 1, fam, false, false);
    }
    if (fam->lastRepeatComplete()) {
        // if starting new repeat instance
        if (v != fam->getPrefix()) {
            // if starting new repeat instance but lmer is not family prefix
            if (verbosity > 2){
                prettyPrintTabbing(tabbing + 1);
                cout << "Last Repeat Complete but Lmer is not Family Prefix." << endl;
			}
            Family* newFam = splitRepeatsByLmer(fam, v, false, L, verbosity, excising);
			families.push_back(newFam);
            return true;
		}
        // NOT IN OLD or OLDEST(else)
        //else if (v->prev() < fam->getLast()->prev() && v->size() != fam->getLast()->size() + 1) {
        //    if ( verbosity > 2){
        //        prettyPrintTabbing(tabbing + 1);
        //        cout << "Last Repeat Complete but prefix is not correct frequency." << endl;
        //        Family* newFam =  exciseRepeatsByLmer(fam, v, L, verbosity);
        //        families.push_back(newFam);
        //        return true;
        //    }
        //}
	}  
    // if it occurs after it is supposed to
    else if (v->back() > fam->getLastIndex() + L) {
        if (verbosity > 2){
            prettyPrintTabbing(tabbing + 1);
            cout << "The last location of the lmer is over L (" << L << ") past the last index of family." << endl;
	    }	
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
    // if it occurs before it is supposed to
    else if (!excising && v->prev() < fam->getLastIndex()) {
        if (verbosity > 2){
            prettyPrintTabbing(tabbing + 1);
            cout << "The previous location of the lmer is before the last index of the family." << endl;
        }
        Family* newFam = splitRepeatsByLmer(fam, v, false, L, verbosity, excising);
        families.push_back(newFam);
        return true;
    }
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

void getElementaryFamilies(seqan::Dna5String &sequence, vector<seqan::CharString> &masks, vector<Family*> &families, int verbosity, int age) {
    const bool family_array = age < 2 ? true : false;
    const bool excising = age < 2 ? true : false;
	const uint seqLength = seqan::length(sequence);
	const seqan::CharString mask = masks.front();
	const uint L = seqan::length(mask);
	const uint MAX_L = maxLength(masks);

	uint index = -1;
	char unmasked[MAX_L];
	LmerMap lmers;
	Family* fam = nullptr;
    // OLD
    //if (family_array) {
        Family* curr_families[L];
        for(uint i = 0; i < L; i++){
            curr_families[i] = NULL;
        }
        int curr_fam = 0;
	    int signed_L = seqan::length(mask);
    //}
    while (getNextSeed(sequence, ++index, seqLength, MAX_L, unmasked)) {
		size_t seed = seedToInt(unmasked, mask);
		LmerVector* v = getAndInsert(lmers, seed, index);
		if (v->size() == 2) {
            // OLD
            if(family_array){
                int fam_index = getFamilyIndex(v, curr_families, curr_fam, signed_L);
                if (fam_index == signed_L){
                    fam = new Family();
                    curr_fam = (curr_fam + 1) % L;
                    curr_families[curr_fam] = fam;
                    families.push_back(fam);
                }
                else{
                    fam = curr_families[fam_index];
                }
			}
            // OLDEST
            else{
                if (isNewFamily(fam, v)) {
                    fam = new Family();
                    families.push_back(fam);
                }
            }
            fam->adopt(v, L);
		} else if (v->size() > 2 && !fragmentSplit(v, L, families, verbosity, excising)) {
            Family* fam = v->getFamily();
			if (isPrefix(v)) {
                fam->setLast(v);
                // CS: modify expected end to include L
				fam->setExpectedEnd(v->back() + fam->repeatLength(L));
			} else if (repeatExpects(fam, v, L, verbosity, excising)) {
				fam->setLast(v);
			}
		}
	}
	tieLooseEnds(families, L, verbosity, excising);
}

#endif //SCANER_SEEDCHAIN
