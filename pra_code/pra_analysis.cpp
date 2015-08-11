// Analyze a pra file
#include <cassert>
#include <string>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <unordered_map>
#include "interval_list.h"
#include <ctime>

using namespace std;

typedef unordered_map<string,interval_list> si_map;


int main(int argc, char** argv) {
  assert(argc == 5);
  string consensus_file = argv[1];
  string rm_fa_file = argv[2];
  string database_file = argv[3];
  string output_file = argv[4];
  string blast_file = consensus_file.substr(0, consensus_file.rfind(".")) + ".blast.6.txt";

  clock_t start_time = clock();
  string blast_format = "\"6 qseqid sseqid qstart qend qlen sstart send slen\"";
  string cmd = (string)"blastn -out " + blast_file + " -outfmt " + blast_format + " -query " + consensus_file + " -db " + database_file + 
																					" -evalue 0.001  -task blastn-short ";  // -word_size 7 -reward 4 -ungapped
  int v = system(cmd.c_str());
  clock_t blast_time = clock();

  ifstream fin(blast_file);
  ofstream fout(output_file);

  string qid;
  string sid;
  int qstart;
  int qend;
  int qlen;
  int sstart;
  int send;
  int slen;

  si_map qmap;
  si_map tmap;
  unordered_map<string,int> lenmap;
  
  fout << "## Blast time: " << (blast_time - start_time)/CLOCKS_PER_SEC << endl;
  while (fin >> qid >> sid >> qstart >> qend >> qlen >> sstart >> send >> slen) {
    if (qend < qstart)
      std::swap(qstart, qend);
    if (send < sstart)
      std::swap(send, sstart);

    
    si_map::iterator i;
    i = qmap.find(qid);
    if (i == qmap.end()) {
      i = qmap.insert(make_pair(qid, interval_list())).first;
      lenmap[qid] = qlen;
    }
    i->second.add(qstart, qend);

    i = tmap.find(sid);
    if (i == tmap.end()) {
      i = tmap.insert(make_pair(sid, interval_list())).first;
      lenmap[sid] = slen;
    }
    i->second.add(sstart, send);
  }
		      

  // Get stats on query sequences
  fout << "# Consesnsus converage (id, num covered bases, length, ratio)\n";
  double sum = 0;
  int count = 0;
  
  for (si_map::iterator i = qmap.begin(); i != qmap.end(); i++) {
    string id = i->first;
    int len = lenmap[id];
    int covered = i->second.num_covered();
    double p = ((double)covered)/len;
    fout << id << " " << covered << " " << len << " " << ((double)covered)/len << endl;

    sum += p;
    count++;
  }
  fout << "# Average consensus coverage: " << sum / count << endl;

  fout << "# Query coverage (id, num covered bases, length, ratio)\n";
  
  sum = 0.0;
  count = 0;

  unordered_map<string,pair<double,int>> D;
  for (si_map::iterator i = tmap.begin(); i != tmap.end(); i++) {
    string id = i->first;
    int len = lenmap[id];
    int covered = i->second.num_covered();
    double p = ((double)covered)/len;
    
    fout << id << " " << covered << " " << len << " " << ((double)covered)/len << endl;

    string s = id.substr(0, id.find("|"));
    unordered_map<string,pair<double,int>>::iterator j = D.find(s);
    if (j == D.end())
      j = D.insert(make_pair(s,make_pair(0.0,0))).first;
    j->second.first += p;
    j->second.second += 1;
  }

  fout << ("# Average family coverage (family, average coverage) -- not sure if this statistic is meaningful\n");
  
  for (auto i=D.begin(); i != D.end(); i++)
    fout << i->first << "\t" << i->second.first << "\t" << i->second.second << endl;

  clock_t finish_time = clock();
  fout << "## Process time: " << ((finish_time - blast_time) / CLOCKS_PER_SEC) << endl;

  return 0;
}
