//First experiment

#include "model.h"
#include <clocale>

double rp = -0.5;
double rv = -0.8;
double alphap = 5.0;
double alphav = 9.0;
double h = 0.01;
double D = 5.0;

int main(int argc, char* argv[])
{
    setlocale(LC_ALL, "C");
    if (argc == 6) 
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

        std::istringstream iss_D(argv[5]); 
        iss_D.imbue(std::locale("C"));
        iss_D >> D;
        
    }
    
    else if (argc == 1)
    {
        cout << "Manual mode: executing with default values." << endl;
        cout << "rp=" << rp << ", rv=" << rv << ", ap=" << alphap << ", av=" << alphav << ", D=" << D << endl;
    }
    
    else
    {
        cout << "ERROR: Incorrect number of parameters." << endl;
        cout << "Manual use: ./exe" << endl;
        cout << "Script use: ./exe <rp> <rv> <alphap> <alphav> <D>" << endl;
        return 1;
    }
    //Gamma
    
    map<string, int> plantIndex;
    map<string, int> insectIndex;
    
    int plantCount = 0, insectCount = 0, numPatch = 0;    
    
    vector<vector<vector<double>>> gamma;
    
    loadGamma("interactions/3interactions_Woolacombe_patches.txt", plantIndex,insectIndex,plantCount,insectCount,numPatch,gamma);
    
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
                v[i][site] = 4.0;
            } 
            else 
            {
                v[i][site] = 0.0;
            }
        }
    }   
    
    //Run
    
    runTargetExtinctionExperiment(p, v, gamma, h, plantCount, insectCount, numPatch, D);
    
    return 0;
}
