//First experiment

#include "model.h"
#include <clocale>

double rp = -0.1;
double rv = -0.1;
double alphap = 8.0;
double alphav = 8.0;
double h = 0.01;
double D = 0.01;
int main(int argc, char* argv[])
{
    setlocale(LC_ALL, "C");
    if (argc == 5) 
    {
        std::istringstream iss_rp(argv[1]);
        iss_rp.imbue(std::locale("C"));
        iss_rp >> rp;
         
        std::istringstream iss_rv(argv[2]); 
        iss_rv.imbue(std::locale("C"));
        iss_rv >> rv;
        
        std::istringstream iss_ap(argv[3]); 
        iss_ap.imbue(std::locale("C"));
        iss_ap >> alphap;
        
        std::istringstream iss_av(argv[4]); 
        iss_av.imbue(std::locale("C"));
        iss_av >> alphav;
        
    }
    
    else if (argc == 1)
    {
        cout << "Manual mode: executing with default values." << endl;
        cout << "rp=" << rp << ", rv=" << rv << ", ap=" << alphap << ", av=" << alphav << endl;
    }
    
    else
    {
        cout << "ERROR: Incorrect number of parameters." << endl;
        cout << "Manual use: ./exe" << endl;
        cout << "Script use: ./exe <rp> <rv> <alphap> <alphav>" << endl;
        return 1;
    }
    //Gamma
    
    map<string, int> plantIndex;
    map<string, int> insectIndex;
    
    int plantCount = 0, insectCount = 0, numPatch = 0;    
    
    vector<vector<vector<double>>> gamma;
    loadGamma("interactions/3interactions_Oxwich_Bay_patches.txt", plantIndex,insectIndex,plantCount,insectCount,numPatch,gamma);
    
    cout << "Number of plants: " << plantCount << endl;
    cout << "Number of insects: " << insectCount << endl;
    
    //Initialization
    
    vector <vector<double>> p(plantCount,vector<double>(numPatch));
    vector <vector<double>> v(insectCount,vector<double>(numPatch));
    
    for(int i=0 ; i<plantCount ; i++)
    {
        for(int site=0 ; site<numPatch ; site++)
        {
            if (plantExistsInPatch(i, site, insectCount, gamma)) 
            {
                p[i][site] = 10.0;
            } 
            else 
            {
                p[i][site] = 0.0; 
            } 
        }
    }
    
    for(int i=0 ; i<insectCount ; i++)
    {
        for(int site=0 ; site<numPatch ; site++)
        {
            if (insectExistsInPatch(i, site, plantCount, gamma)) 
            {
                v[i][site] = 10.0;
            } 
            else 
            {
                v[i][site] = 0.0;
            }
        }
    }   
    
    //Run
    
    runExtinctionExperiment(p, v, gamma, h, plantCount, insectCount, numPatch, D);
    
    return 0;
}
