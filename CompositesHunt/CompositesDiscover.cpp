#include "CompositeFind.h"

int main(int argc, char** argv)
{
	CompositeFind compositeFind = CompositeFind();
	compositeFind.readInFamilies(argv[1]);
	compositeFind.buildGraph(argv[2], argv[3]);
	return 0;
// 	CompositeFind compositeFind = CompositeFind();
// 	compositeFind.readInFamilies("C:\\Users\\alpha\\Desktop\\Research\\Com_Finder\\hg18\\ele\\elements");
// 	compositeFind.buildGraph("", "C:\\Users\\alpha\\Desktop\\Research\\Com_Finder\\composites_cpp.txt");
// 	return 0;
}

