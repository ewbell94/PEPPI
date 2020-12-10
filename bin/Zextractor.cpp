#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <stdlib.h>

using namespace std;

int main(int argc, char** argv){
  ifstream hhrfile(argv[1]);
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
      total+=score;
      scores.push_back(score);
    }
  }

  if (scores.size() != templates.size()){
    cout << "Array sizes are not equal" << endl;
    return 1;
  }

  double mean=total/scores.size();
  double stdev=0.0;
  for (vector<double>::iterator it=scores.begin(); it!=scores.end(); ++it){
    stdev+=pow(*it-mean,2.0)/scores.size();
  }
  stdev=pow(stdev,0.5);
  
  vector<string>::iterator tit=templates.begin();
  vector<double>::iterator sit=scores.begin();
  while (tit!=templates.end() || sit!=scores.end()){
    double z=((*sit-mean)/stdev);
    cout << *tit << " " << z << endl;
    ++tit;
    ++sit;
  }
}
