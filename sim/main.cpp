//This program simulates a metapopulation model of pollinators-plants ecosystem

#include <iostream>
#include <cmath>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include <tuple>

#define r 2.0
#define Kp 100
#define Kv 1000
#define ha 1
#define d 0.3
#define D 2.5

using namespace std;

void evaluaFp(double p, const vector<vector<double>>& v, double &fp, double** gamma, int pindex, int site, int insectCount)
{
    fp = 0.0;
    double sum = 0.0;
    for( int j=0 ; j<insectCount ; j++)
        sum += (gamma[pindex][j]*v[j][site]) / (1.0 + (ha*gamma[pindex][j]*v[j][site])); 
    fp = (r*p*(1.0-(p/Kp))) + (p*sum); 
    return;
}

void evaluaFv(const vector<vector<double>>& p, const vector<vector<double>>& v, double &fv, double** gamma, int vindex, int site, int plantCount, int numpatch)
{
    fv = 0.0;
    double sum = 0.0;
    double sumD = 0.0;
    for( int i=0 ; i<plantCount ; i++)
        sum += (gamma[i][vindex]*p[i][site]*v[vindex][site]) / (1.0 + (ha*gamma[i][vindex]*p[i][site]));
    for( int i=0 ; i<numpatch ; i++)
        if ( i!=site)
            sumD += (v[vindex][i]-v[vindex][site]); 
    fv = (-1.0*d*v[vindex][site]*(1.0+(v[vindex][site]/Kv))) + sum + (D*sumD); 
    return;
}

void rungekutta(vector<vector<double>>& p, vector<vector<double>>& v, double** gamma, double h, int plantCount, int insectCount, int numpatch)
{   
    double sumav, sumap;
    
    vector<vector<double>> k1v(insectCount, vector<double>(numpatch));
    vector<vector<double>> k2v(insectCount, vector<double>(numpatch));
    vector<vector<double>> k3v(insectCount, vector<double>(numpatch));
    vector<vector<double>> k4v(insectCount, vector<double>(numpatch));
    
    vector<vector<double>> k1p(plantCount, vector<double>(numpatch));
    vector<vector<double>> k2p(plantCount, vector<double>(numpatch));
    vector<vector<double>> k3p(plantCount, vector<double>(numpatch));
    vector<vector<double>> k4p(plantCount, vector<double>(numpatch)); 
    
    vector<vector<double>> fv(insectCount, vector<double>(numpatch));
    vector<vector<double>> fp(plantCount, vector<double>(numpatch));
    
    vector<vector<double>> auxp(plantCount, vector<double>(numpatch));
    vector<vector<double>> auxv(insectCount, vector<double>(numpatch));
    
    //k1p k1v
    for(int i=0 ; i<plantCount ; i++)
    {
        for(int site=0 ; site<numpatch ; site++)
        {
            k1p[i][site] = 0.0;
            evaluaFp(p[i][site],v,fp[i][site],gamma,i,site,insectCount);
            k1p[i][site] = h*fp[i][site];
        }
    }
    for(int i=0 ; i<insectCount ; i++)
    {
        for(int site=0 ; site<numpatch ; site++)
        {
            k1v[i][site] = 0.0;
            evaluaFv(p,v,fv[i][site],gamma,i,site,plantCount,numpatch);
            k1v[i][site] = h*fv[i][site];
        }
    }
    //k2p k2v
    for(int i=0 ; i<plantCount ; i++)
    {
        for(int site=0 ; site<numpatch ; site++)
        {
            auxp[i][site] = 0.0;
            auxp[i][site] = p[i][site]+(0.5*k1p[i][site]);
        }
    }
    
    for(int i=0 ; i<insectCount ; i++)
    {
        for(int site=0 ; site<numpatch ; site++)
        {
            auxv[i][site] = 0.0; 
            auxv[i][site] = v[i][site]+(0.5*k1v[i][site]);
        }
    }
    
    for(int i=0 ; i<plantCount ; i++)
    {
        for(int site=0 ; site<numpatch ; site++)
        {
            k2p[i][site] = 0.0;
            evaluaFp(auxp[i][site],auxv,fp[i][site],gamma,i,site,insectCount);
            k2p[i][site] = h*fp[i][site];
        }
    }
    
    for(int i=0 ; i<insectCount ; i++)
    {
        for(int site=0 ; site<numpatch ; site++)
        {
            k2v[i][site] = 0.0;
            evaluaFv(auxp,auxv,fv[i][site],gamma,i,site,plantCount,numpatch);
            k2v[i][site] = h*fv[i][site];
        }
    }
    //k3p k3v
    for(int i=0 ; i<plantCount ; i++)
    {
        for(int site=0 ; site<numpatch ; site++)
        {
            auxp[i][site] = 0.0;
            auxp[i][site] = p[i][site]+(0.5*k2p[i][site]);
        }
    }
    
    for(int i=0 ; i<insectCount ; i++)
    {
        for(int site=0 ; site<numpatch ; site++)
        {
            auxv[i][site] = 0.0;
            auxv[i][site] = v[i][site]+(0.5*k2v[i][site]);
        }
    }
    
    for(int i=0 ; i<plantCount ; i++)
    {
        for(int site=0 ; site<numpatch ; site++)
        {
            k3p[i][site] = 0.0;
            evaluaFp(auxp[i][site],auxv,fp[i][site],gamma,i,site,insectCount);
            k3p[i][site] = h*fp[i][site];
        }
    }
    
    for(int i=0 ; i<insectCount ; i++)
    {
        for(int site=0 ; site<numpatch ; site++)
        {
            k3v[i][site] = 0.0;
            evaluaFv(auxp,auxv,fv[i][site],gamma,i,site,plantCount,numpatch);
            k3v[i][site] = h*fv[i][site];
        }
    }
    //k4p k4v
    for(int i=0 ; i<plantCount ; i++)
    {
        for(int site=0 ; site<numpatch ; site++)
        {
            auxp[i][site] = 0.0;
            auxp[i][site] = p[i][site] + k3p[i][site];
        }
    }
    
    for(int i=0 ; i<insectCount ; i++)
    {
        for(int site=0 ; site<numpatch ; site++)
        {
            auxv[i][site] = v[i][site] + k3v[i][site];
        }
    }
    
    for(int i=0 ; i<plantCount ; i++)
    {
        for(int site=0 ; site<numpatch ; site++)
        {
            k4p[i][site] = 0.0;
            evaluaFp(auxp[i][site],auxv,fp[i][site],gamma,i,site,insectCount);
            k4p[i][site] = h*fp[i][site];
            sumap=0.0;
            sumap=k1p[i][site]+(2*k2p[i][site])+(2*k3p[i][site])+k4p[i][site];
            p[i][site]+=(sumap/6.0);
        }
    }
    
    for(int i=0 ; i<insectCount ; i++)
    {
        for(int site=0 ; site<numpatch ; site++)
        {
            k4v[i][site] = 0.0;
            evaluaFv(auxp,auxv,fv[i][site],gamma,i,site,plantCount,numpatch);
            k4v[i][site] = h*fv[i][site];
            sumav=0.0;
            sumav=k1v[i][site]+(2*k2v[i][site])+(2*k3v[i][site])+k4v[i][site];
            v[i][site]+=(sumav/6.0);
        }
    }
    return;
}

int main(void)
{
    int n, numpatch;
    double h, t;
    ofstream fichp, fichv;
    cout << "Introduce the step: " << endl;
    cin >> h;
    
    cout << "Introduce the number of iterations: ";
    cin >> n;
    
    cout << "Introduce the number of patches: ";
    cin >> numpatch;
    
    //Gamma
    
    ifstream intfich("interactionsBurnham-on-Sea.txt");
    
    map<string, int> plantIndex;
    map<string, int> insectIndex;
    string plant, insect;
    double weight;
    
    int plantCount = 0, insectCount = 0;    
    
    string line;
    
    vector<tuple<string, string, double>> data;
    
    while (getline(intfich, line)) 
    {
        istringstream iss(line);
        if (!(iss >> plant >> insect >> weight)) continue;

        if (plantIndex.find(plant) == plantIndex.end())
            plantIndex[plant] = plantCount++;
        if (insectIndex.find(insect) == insectIndex.end())
            insectIndex[insect] = insectCount++;

        data.push_back({plant, insect, weight});
    }
    
    intfich.close();
    
    double** gamma = new double*[plantCount];
    
    for (int i = 0; i < plantCount; ++i) 
    {
        gamma[i] = new double[insectCount];
        for (int j = 0; j < insectCount; ++j)
            gamma[i][j] = 0.0;
    }
    
    for (const auto& item : data) 
    {
        string pl = get<0>(item);
        string ins = get<1>(item);
        double w = get<2>(item);
        int i = plantIndex[pl];
        int j = insectIndex[ins];
        gamma[i][j] = w;
    }
    
    vector <vector<double>> p(plantCount,vector<double>(numpatch));
    vector <vector<double>> v(insectCount,vector<double>(numpatch));
    //Initialization
    
    for(int i=0 ; i<plantCount ; i++)
    {
        for(int site=0 ; site<numpatch ; site++)
            p[i][site] = 100;   
    }
    
    for(int i=0 ; i<insectCount ; i++)
    {
        for(int site=0 ; site<numpatch ; site++)
            v[i][site] = 500;    
    }
    
    cout << "Number of plants: " << plantCount << endl;
    cout << "Number of insects: " << insectCount << endl;
    //Run
    fichp.open("evolutionp.txt");
    fichv.open("evolutionv.txt");

    t = 0.0;
    for(int it=0 ; it<n ; it++)
    {   
        fichp << t << " ";
        fichv << t << " ";
        
        for(int i=0 ; i<plantCount ; i++)
        {
            for(int site=0 ; site<numpatch ; site++)
            {
                fichp << p[i][site] << " ";
            }
        }
        
        for(int i=0 ; i<insectCount ; i++)
        {
            for(int site=0 ; site<numpatch ; site++)
            {
                fichv << v[i][site] << " ";
            }
        }
        fichp << endl;
        fichv << endl;


        rungekutta(p,v,gamma,h,plantCount,insectCount,numpatch); 
        t+=h;   
    }
    fichp.close();
    fichv.close();
    
    for (int i = 0; i < plantCount; ++i)
        delete[] gamma[i];
    delete[] gamma;
    return 0;
}
