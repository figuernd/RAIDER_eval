#ifndef LIU_RAIDER_FAMILY
#define LIU_RAIDER_FAMILY

#include <iostream>
#include <vector>
#include <iterator>
#include <algorithm>
#include <map>
#include <unordered_map>
//#include <unordered_set>
#include <tuple>

class Family;
typedef std::tuple<Family*, Family*, Family*> tpl_Fam;
typedef Family* famPr;

struct key_hash : public std::unary_function<tpl_Fam, std::size_t> {
	std::size_t operator()(const tpl_Fam & tf) const {
		return (size_t)std::get<0>(tf) * (size_t)std::get<1>(tf) * (size_t)std::get<2>(tf);
	}
};

struct key_equal : public std::binary_function<tpl_Fam, tpl_Fam, bool> {
	bool operator()(const tpl_Fam& v0, const tpl_Fam& v1) const {
		return (
			std::get<0>(v0) == std::get<0>(v1) &&
			std::get<1>(v0) == std::get<1>(v1) &&
			std::get<2>(v0) == std::get<2>(v1)
			);
	}
};

class Family
{
public:
	Family();
	Family(Family* otherFam);
	void increaseCount();
	int getCount() const;
	int getTotalLen() const;
	void setMinL(int L);
	void addStartPosition(int startPos);
	void addPosition(int startPos, int endPos);
	int getID() const;
	void setID(int ID);
	void dfsHunt(Family* preNode); 
	void dfsHuntRePath(Family* nxtNode);
	void clear();

	bool operator<(const Family& otherFam) const;
	//bool operator=(const Family& otherFam) const;
	
	int cnt;
	int id;
	int minL;
	float maxWgt;
	int reWgt;
	int nxtNoneInterLen; 
	int status;
//	Family *parent;
	Family *nextNode;
	bool hasParent;
	std::vector<int> startPositions;
	std::vector<int> endPositions;
	std::vector<Family*> parents;
	std::unordered_map<Family*, float> descendants;
	std::unordered_map<Family*, int> deNoneInterLen;
};

#endif

