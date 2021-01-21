#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <stdlib.h>

using namespace std;

const string complexfile="/expanse/lustre/scratch/ewbell/temp_project/SPRINGDB/70CDHITstruct.txt";
const int maxmono=5000;

struct namedScore{
  string name;
  double score;

  bool operator<(namedScore other) const {
    return score > other.score;
  }
};

unordered_map<string,vector<string>> constructComplexes(string complexfilename){
  ifstream complexf(complexfilename.c_str());
  unordered_map<string, vector<string>> complexlist;
  for (string line; getline(complexf,line);){
    int dashind=line.find("-");
    string firstid=line.substr(0,dashind);
    string secondid=line.substr(dashind+1);
    unordered_map<string,vector<string>>::const_iterator fit=complexlist.find(firstid);
    if (fit==complexlist.end()){
      vector<string> toadd;
      toadd.push_back(secondid);
      complexlist[firstid]=toadd;
    } else {
      complexlist[firstid].push_back(secondid);
    }

    unordered_map<string,vector<string>>::const_iterator sit=complexlist.find(secondid);
    if (sit==complexlist.end()){
      vector<string> toadd;
      toadd.push_back(firstid);
      complexlist[secondid]=toadd;
    } else {
      complexlist[secondid].push_back(firstid);
    }
  }
  return complexlist;
}

vector<namedScore> fetchHits(string hhrfilename){
  vector<namedScore> hits;

  ifstream hhrfile(hhrfilename.c_str());
  vector<string> templates;
  vector<double> scores;

  double total=0.0;
  for (string line; getline(hhrfile,line);){
    if (line.substr(0,1).compare(">") == 0){
      string target=line.substr(1);
      if (find(templates.begin(),templates.end(),target) != templates.end()){
        continue;
      }
      templates.push_back(target);
      getline(hhrfile,line);
      double score=atof(line.substr(line.find("Sum_probs=")+10).c_str());
      //double score=atof(line.substr(line.find("Score=")+6,6).c_str());
      //cout<<line<<endl<<score<<endl;
      total+=score;
      scores.push_back(score);
    }
  }

  if (scores.size() != templates.size()){
    cout << "Array sizes are not equal" << endl;
    return hits;
  }

  //cout << templates.size() << endl;
  double mean=total/scores.size();
  //cout << mean << endl;
  double stdev=0.0;
  for (vector<double>::iterator it=scores.begin(); it!=scores.end(); ++it){
    stdev+=pow(*it-mean,2.0)/scores.size();
  }
  stdev=pow(stdev,0.5);

  vector<string>::iterator tit=templates.begin();
  vector<double>::iterator sit=scores.begin();
  while (tit!=templates.end() || sit!=scores.end()){
    double z=((*sit-mean)/stdev);
    namedScore hit;
    hit.name=*tit;
    hit.score=z;
    hits.push_back(hit);
    ++tit;
    ++sit;
  }
  
  sort(hits.begin(),hits.end());
  return hits;
}

vector<namedScore> fetchDimers(vector<namedScore> hhr1hits, vector<namedScore> hhr2hits, unordered_map<string,vector<string>> complexlist){
  vector<namedScore> dimerlist;
  for (vector<namedScore>::iterator iti=hhr1hits.begin(); iti!=hhr1hits.end(); ++iti){
    string name1=iti->name;
    //string pdbid1=name1.substr(0,4);
    vector<string> partners=complexlist[name1];
    for (vector<string>::iterator itj=partners.begin(); itj!=partners.end(); ++itj){
      string pname=*itj;
      for (vector<namedScore>::iterator itk=hhr2hits.begin(); itk!=hhr2hits.end(); ++itk){
	string name2=itk->name;
	if (name2.compare(pname) == 0){
	  double score;
	  if (iti->score < itk->score){
	    score=iti->score;
	  } else {
	    score=itk->score;
	  }
	  namedScore dimerpair;
	  dimerpair.name=name1+" "+name2;
	  dimerpair.score=score;
	  dimerlist.push_back(dimerpair);
	  break;
	}
      }
    }
  }
  sort(dimerlist.begin(),dimerlist.end());
  return dimerlist;
}

int main(int argc, char** argv){
  
  unordered_map<string,vector<string>> complexlist=constructComplexes(complexfile);
  string hhr1=argv[1];
  string hhr2=argv[2];

  vector<namedScore> hhr1hits=fetchHits(hhr1);
  vector<namedScore> hhr2hits=fetchHits(hhr2);

  cout<<hhr1hits.begin()->name<<endl<<hhr2hits.begin()->name<<endl;
  vector<namedScore> tophhr1(hhr1hits.begin(),hhr1hits.begin()+maxmono);
  vector<namedScore> tophhr2(hhr2hits.begin(),hhr2hits.begin()+maxmono);

  vector<namedScore> dimers=fetchDimers(tophhr1,tophhr2,complexlist);

  for (vector<namedScore>::iterator it=dimers.begin(); it!=dimers.end(); ++it){
    cout << it->name << " " << it->score << endl;
  }
  return 0;
}
