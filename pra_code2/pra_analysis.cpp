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

bool debug;

// Get the number of sequences in an fa file
int fa_size(string file) {
  ifstream fin(file);

  if (not fin) {
    cerr << "pra_analysis: File not found (" << file << ")" << endl;
    exit(0);
  }

  int count = 0;
  string line;
  while (getline(fin, line))
    if (line[0] == '>')
      count++;
  return count;
}

int main(int argc, char** argv) {
  int marker = 0;
  if ((string)argv[1] == (string)"-d") {
    debug = true;
    marker++;
  }
  else {
    debug = false;
  }

  string output_file = argv[1];

 
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
  
  while (cin >> qid >> sid >> qstart >> qend >> qlen >> sstart >> send >> slen) {
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


  return 0;
}
