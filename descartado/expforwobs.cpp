//First experiment

#include "model.h"

int main(void)
{
    double h, D;
    
    /*
    cout << "Introduce the step: " << endl;
    cin >> h;

    cout << "Introduce Dispersal (D): " << endl;
    cin >> D;
    */
    
    h = 0.01;
    D = 0.0;
    //Gamma
    
    map<string, int> plantIndex;
    map<string, int> insectIndex;
    
    int plantCount = 0, insectCount = 0, numPatch = 0;    
    
    vector<vector<vector<double>>> gamma;
    
    loadGamma("interactions/interactions_Penhale_Sands_patches.txt", plantIndex,insectIndex,plantCount,insectCount,numPatch,gamma);
    
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
    
    runExtinctionExperiment(p, v, gamma, h, plantCount, insectCount, numPatch, D, rp, rv);
    
    return 0;
}
