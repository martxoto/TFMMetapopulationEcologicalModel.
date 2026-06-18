#ifndef MODEL_H
#define MODEL_H

#include <iostream>
#include <cmath>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include <tuple>
#include <algorithm>
#include <iomanip>
#include <random>

using namespace std;

extern double rp;
extern double rv;
extern double alphap;
extern double alphav;
extern double D;
constexpr double ha = 1.0;
constexpr double ap = 0.1;
constexpr double av = 0.1;
constexpr double viability = 0.01;




bool plantExistsInPatch(int plantID, int site, int insectCount, const vector<vector<vector<double>>>& gamma);
bool insectExistsInPatch(int insectID, int site, int plantCount, const vector<vector<vector<double>>>& gamma);

void evaluaFp(double p, const vector<vector<double>>& v, double &fp, const vector<vector<vector<double>>>& gamma, int pindex, int site, int insectCount);

void evaluaFv(const vector<vector<double>>& p, const vector<vector<double>>& v, double &fv, const vector<vector<vector<double>>>& gamma, int vindex, int site, int plantCount, int numpatch, double D);

void rungekutta(vector<vector<double>>& p, vector<vector<double>>& v, const vector<vector<vector<double>>>& gamma, double h, int plantCount, int insectCount, int numpatch, double D);

void findSteadyState(double& t, vector<vector<double>>& p, vector<vector<double>>& v, ofstream& fichp, ofstream& fichv, int plantCount, int insectCount, int numpatch,  const vector<vector<vector<double>>>& gamma, double h, double D);

void loadGamma(const string &filename, map<string, int>& plantIndex, map<string, int>& insectIndex, int& plantCount, int& insectCount, int& numpatch, vector<vector<vector<double>>>& gamma);

void runExtinctionExperiment(const vector<vector<double>>& p, const vector<vector<double>>& v, const vector<vector<vector<double>>>& gamma, double h, int plantCount, int insectCount, int numpatch, double D);


void runRandomExtinctionExperiment(const vector<vector<double>>& p, const vector<vector<double>>& v, const vector<vector<vector<double>>>& gamma, double h, int plantCount, int insectCount, int numpatch, double D);

void runTargetExtinctionExperiment(const vector<vector<double>>& p, const vector<vector<double>>& v, const vector<vector<vector<double>>>& gamma, double h, int plantCount, int insectCount, int numPatch, double D);

void runLocalTargetExtinctionExperiment(const vector<vector<double>>& p, const vector<vector<double>>& v, const vector<vector<vector<double>>>& gamma, double h, int plantCount, int insectCount, int numPatch, double D);

#endif
