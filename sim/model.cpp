//This program simulates a metapopulation model of pollinators-plants ecosystem

#include "model.h"

bool plantExistsInPatch(int plantID, int site, int insectCount, const vector<vector<vector<double>>>& gamma)
{
    for (int j = 0; j < insectCount; ++j) {
        if (gamma[site][plantID][j] > 0.0) {
            return true;
        }
    }
    return false;
}

bool insectExistsInPatch(int insectID, int site, int plantCount, const vector<vector<vector<double>>>& gamma)
{
    for (int i = 0; i < plantCount; ++i) {
        if (gamma[site][i][insectID] > 0.0) {
            return true;
        }
    }
    return false; 
}

void evaluaFp(double p, const vector<vector<double>>& v, double &fp, const vector<vector<vector<double>>>& gamma, int pindex, int site, int insectCount)
{
    fp = 0.0;
    double sum = 0.0;
    
    //Compute degree
    double degree = 0.0;
    for(int i=0 ; i<insectCount ; i++)
    {
        if (gamma[site][pindex][i] > 0.0)
            degree += 1.0;
    }   
    
    if (degree == 0.0) degree = 1.0;
    
    for(int i=0 ; i<insectCount ; i++)
    {
        if (gamma[site][pindex][i] > 0.0)
        {
            double gammaNorm = gamma[site][pindex][i] / degree;
            sum += (gammaNorm * v[i][site])/(1.0+( ha * gammaNorm * v[i][site]));
        }     
    }    
    fp = p*(rp - (ap*p) + (alphap*sum));
    return;
}

void evaluaFv(const vector<vector<double>>& p, const vector<vector<double>>& v, double &fv, const vector<vector<vector<double>>>& gamma, int vindex, int site, int plantCount, int numPatch, double D)
{
    fv = 0.0;
    double sum = 0.0, sumD = 0.0;
    
    //Compute degree
    double degree = 0.0;
    for(int i=0 ; i<plantCount ; i++)
        if (gamma[site][i][vindex] > 0.0) degree += 1.0;
    
    if (degree == 0.0) degree = 1.0;
    
            
    for(int i=0 ; i<plantCount ; i++)
    {
        if (gamma[site][i][vindex] > 0.0)
        {
            double gammaNorm = gamma[site][i][vindex] / degree;
            sum += (gammaNorm * p[i][site])/(1.0+(ha*gammaNorm * p[i][site]));
        }
    }     
    for( int i=0 ; i<numPatch ; i++)
        if ( i!=site)
            sumD += (v[vindex][i]-v[vindex][site]); 
            
    fv = v[vindex][site]*(rv-(av*v[vindex][site]) + (alphav*sum)) + (D*sumD);
    return;
}


void rungekutta(vector<vector<double>>& p, vector<vector<double>>& v, const vector<vector<vector<double>>>& gamma, double h, int plantCount, int insectCount, int numPatch, double D)
{   
    double sumav, sumap;
    
    vector<vector<double>> k1v(insectCount, vector<double>(numPatch));
    vector<vector<double>> k2v(insectCount, vector<double>(numPatch));
    vector<vector<double>> k3v(insectCount, vector<double>(numPatch));
    vector<vector<double>> k4v(insectCount, vector<double>(numPatch));
    
    vector<vector<double>> k1p(plantCount, vector<double>(numPatch));
    vector<vector<double>> k2p(plantCount, vector<double>(numPatch));
    vector<vector<double>> k3p(plantCount, vector<double>(numPatch));
    vector<vector<double>> k4p(plantCount, vector<double>(numPatch)); 
    
    vector<vector<double>> fv(insectCount, vector<double>(numPatch));
    vector<vector<double>> fp(plantCount, vector<double>(numPatch));
    
    vector<vector<double>> auxp(plantCount, vector<double>(numPatch));
    vector<vector<double>> auxv(insectCount, vector<double>(numPatch));
    
    //k1p k1v
    for(int i=0 ; i<plantCount ; i++)
    {
        for(int site=0 ; site<numPatch ; site++)
        {
            k1p[i][site] = 0.0;
            evaluaFp(p[i][site],v,fp[i][site],gamma,i,site,insectCount);
            k1p[i][site] = h*fp[i][site];
        }
    }
    for(int i=0 ; i<insectCount ; i++)
    {
        for(int site=0 ; site<numPatch ; site++)
        {
            k1v[i][site] = 0.0;
            evaluaFv(p,v,fv[i][site],gamma,i,site,plantCount,numPatch,D);
            k1v[i][site] = h*fv[i][site];
        }
    }
    //k2p k2v
    for(int i=0 ; i<plantCount ; i++)
    {
        for(int site=0 ; site<numPatch ; site++)
        {
            auxp[i][site] = 0.0;
            auxp[i][site] = p[i][site]+(0.5*k1p[i][site]);
        }
    }
    
    for(int i=0 ; i<insectCount ; i++)
    {
        for(int site=0 ; site<numPatch ; site++)
        {
            auxv[i][site] = 0.0; 
            auxv[i][site] = v[i][site]+(0.5*k1v[i][site]);
        }
    }
    
    for(int i=0 ; i<plantCount ; i++)
    {
        for(int site=0 ; site<numPatch ; site++)
        {
            k2p[i][site] = 0.0;
            evaluaFp(auxp[i][site],auxv,fp[i][site],gamma,i,site,insectCount);
            k2p[i][site] = h*fp[i][site];
        }
    }
    
    for(int i=0 ; i<insectCount ; i++)
    {
        for(int site=0 ; site<numPatch ; site++)
        {
            k2v[i][site] = 0.0;
            evaluaFv(auxp,auxv,fv[i][site],gamma,i,site,plantCount,numPatch,D);
            k2v[i][site] = h*fv[i][site];
        }
    }
    //k3p k3v
    for(int i=0 ; i<plantCount ; i++)
    {
        for(int site=0 ; site<numPatch ; site++)
        {
            auxp[i][site] = 0.0;
            auxp[i][site] = p[i][site]+(0.5*k2p[i][site]);
        }
    }
    
    for(int i=0 ; i<insectCount ; i++)
    {
        for(int site=0 ; site<numPatch ; site++)
        {
            auxv[i][site] = 0.0;
            auxv[i][site] = v[i][site]+(0.5*k2v[i][site]);
        }
    }
    
    for(int i=0 ; i<plantCount ; i++)
    {
        for(int site=0 ; site<numPatch ; site++)
        {
            k3p[i][site] = 0.0;
            evaluaFp(auxp[i][site],auxv,fp[i][site],gamma,i,site,insectCount);
            k3p[i][site] = h*fp[i][site];
        }
    }
    
    for(int i=0 ; i<insectCount ; i++)
    {
        for(int site=0 ; site<numPatch ; site++)
        {
            k3v[i][site] = 0.0;
            evaluaFv(auxp,auxv,fv[i][site],gamma,i,site,plantCount,numPatch,D);
            k3v[i][site] = h*fv[i][site];
        }
    }
    //k4p k4v
    for(int i=0 ; i<plantCount ; i++)
    {
        for(int site=0 ; site<numPatch ; site++)
        {
            auxp[i][site] = 0.0;
            auxp[i][site] = p[i][site] + k3p[i][site];
        }
    }
    
    for(int i=0 ; i<insectCount ; i++)
    {
        for(int site=0 ; site<numPatch ; site++)
        {
            auxv[i][site] = v[i][site] + k3v[i][site];
        }
    }
    
    for(int i=0 ; i<plantCount ; i++)
    {
        for(int site=0 ; site<numPatch ; site++)
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
        for(int site=0 ; site<numPatch ; site++)
        {
            k4v[i][site] = 0.0;
            evaluaFv(auxp,auxv,fv[i][site],gamma,i,site,plantCount,numPatch,D);
            k4v[i][site] = h*fv[i][site];
            sumav=0.0;
            sumav=k1v[i][site]+(2*k2v[i][site])+(2*k3v[i][site])+k4v[i][site];
            v[i][site]+=(sumav/6.0);
        }
    }
    return;
}

void findSteadyState(double& t, vector<vector<double>>& p, vector<vector<double>>& v, ofstream& fichp, ofstream& fichv, int plantCount, int insectCount, int numPatch,  const vector<vector<vector<double>>>& gamma, double h, double D)
{
    //Stationary state detection
    double max_delta = 1.0;
    const double TOLERANCE = 1e-6;
    int iter_count = 0;
    int max_iter = 1000000;
    
    vector<vector<double>> p_prev(plantCount, vector<double>(numPatch));
    vector<vector<double>> v_prev(insectCount, vector<double>(numPatch));
    
    while(iter_count < max_iter)
    {   
        fichp << t << " ";
        fichv << t << " ";
        
        for(int i=0 ; i<plantCount ; i++)
            for(int site=0 ; site<numPatch ; site++)
                fichp << p[i][site] << " ";
        
        for(int i=0 ; i<insectCount ; i++)
            for(int site=0 ; site<numPatch ; site++)
                fichv << v[i][site] << " ";

        fichp << endl;
        fichv << endl;

        p_prev = p;
        v_prev = v;

        rungekutta(p,v,gamma,h,plantCount,insectCount,numPatch,D); 
        t+=h;   
        iter_count++;
        
        max_delta = 0.0;
        for(int i=0 ; i<plantCount ; i++)
        {
            for(int site=0 ; site<numPatch ; site++)
            {
                double delta = abs(p[i][site] - p_prev[i][site]);
                if(delta > max_delta)
                    max_delta = delta;
            }
        }
        
        for(int i=0 ; i<insectCount ; i++)
        {
            for(int site=0 ; site<numPatch ; site++)
            {
                double delta = abs(v[i][site] - v_prev[i][site]);
                if(delta > max_delta)
                    max_delta = delta;
            }
        }
        
        if (max_delta < TOLERANCE && iter_count > 1000) 
        {
            cout << "Stationary state at t = " << t << " (iter " << iter_count << ")" << endl;
            break;
        }
    }
    
    if (iter_count == max_iter)
        cout << "  No convergence." << endl;
    
    fichp << t << " ";
    fichv << t << " ";
    for(int i=0 ; i<plantCount ; i++)
        for(int site=0 ; site<numPatch ; site++)
            fichp << p[i][site] << " ";
    for(int i=0 ; i<insectCount ; i++)
        for(int site=0 ; site<numPatch ; site++)
            fichv << v[i][site] << " ";
    fichp << endl;
    fichv << endl;
    
    return;
}

void loadGamma(const string &filename, map<string, int>& plantIndex, map<string, int>& insectIndex, int& plantCount, int& insectCount, int& numPatch, vector<vector<vector<double>>>& gamma)
{
    ifstream intfich(filename);

    string plant, insect;
    double weight;
    int site;

    string line;
    
    vector<tuple<int, string, string, double>> data;
    
    while (getline(intfich, line)) 
    {
        istringstream iss(line);
        if (!(iss >> site >> plant >> insect >> weight)) continue;

        if (plantIndex.find(plant) == plantIndex.end())
            plantIndex[plant] = plantCount++;
        if (insectIndex.find(insect) == insectIndex.end())
            insectIndex[insect] = insectCount++;
        
        if (site + 1 > numPatch)
            numPatch = site + 1;

        data.push_back({site, plant, insect, weight});
    }
    
    intfich.close();
    
    cout << endl << numPatch << " patches." << endl << endl;
    
    gamma.resize(numPatch, vector<vector<double>>(plantCount, vector<double>(insectCount, 0.0)));
    
    for (const auto& item : data) 
    {
        int p_id = get<0>(item);
        string pl = get<1>(item);
        string ins = get<2>(item);
        double w = get<3>(item);
        
        int i = plantIndex[pl];
        int j = insectIndex[ins];
        gamma[p_id][i][j] = w;
    }
    
    return;
}

void runExtinctionExperiment(const vector<vector<double>>& p, const vector<vector<double>>& v, const vector<vector<vector<double>>>& gamma, double h, int plantCount, int insectCount, int numPatch, double D)
{
    cout << "\n---Initializing extinction experiment ---" << endl;  
    
    //Preparation
    ofstream experimentFile("results.txt");
    ofstream fichp("evolutionp.txt");
    ofstream fichv("evolutionv.txt");
    
    experimentFile << "# Num_Extinctions Robustness_Ratio Surv_Plants Surv_Insects Pollination_Service Gini_Plants Gini_Insects Norm_Biomass" << endl;
    
    //ofstream dummy_p("/dev/null");
    //ofstream dummy_v("/dev/null");
    
    vector<vector<double>> pCurrent = p;
    vector<vector<double>> vCurrent = v;
    
    double tDummy = 0.0;
    
    findSteadyState(tDummy,pCurrent,vCurrent,fichp,fichv,plantCount,insectCount,numPatch,gamma,h,D);
    
    int initialSurvPlants = 0;
    int initialSurvInsects = 0;
    double initialBiomass = 0.0;

    vector<bool> isDead(plantCount, false);
    
    for(int i=0; i<plantCount; ++i) 
    {
        double total = 0.0;
        for(int s=0; s<numPatch; ++s) 
            total += pCurrent[i][s];
            
        if(total > viability)
        {
            initialSurvPlants++;
            initialBiomass += total;
        }
        
        else isDead[i] = true;
    }
    
    for(int i=0; i<insectCount; ++i) 
    {
        double total = 0.0;
        for(int s=0; s<numPatch; ++s) 
            total += vCurrent[i][s];
            
        if(total > viability)
        {
            initialSurvInsects++;
            initialBiomass += total;
        }
    }
    
    int totalInitialSpecies = initialSurvPlants + initialSurvInsects;
    
    cout << "Especies Iniciales Vivas: " << totalInitialSpecies << " (Plantas: " << initialSurvPlants << ", Insectos: " << initialSurvInsects << ")" << endl;
    
    double sumRobustness = 0.0;
    int kEffective = 0, steps = 0;
    //Experiment
    for(int k=0 ; k<=plantCount ; k++)
    {
        if (k>0)
        {
            int plantToRemove = -1;
            double minAbundance = 1.0e20;
            
            for(int i=0 ; i<plantCount ; i++)
            {
                if (isDead[i]) continue;
                
                double currentAbundance = 0.0;
                for(int s=0 ; s<numPatch ; s++) currentAbundance += pCurrent[i][s];
                
                if (currentAbundance <= viability)
                {
                    isDead[i] = true;
                    continue;
                }
                
                if (currentAbundance < minAbundance)
                {
                    minAbundance = currentAbundance;
                    plantToRemove = i;
                }
            }
            
            if (plantToRemove == -1) break;
            
            for(int s=0 ; s<numPatch ; s++) pCurrent[plantToRemove][s] = 0.0;
            
            isDead[plantToRemove] = true;
            
            kEffective++;
            
            findSteadyState(tDummy, pCurrent, vCurrent, fichp, fichv, plantCount, insectCount, numPatch, gamma, h,D);
        }
        
        //Metrics
        int survPlants = 0;
        int survInsects = 0;
        double plantBiomass = 0.0, insectBiomass = 0.0;
        double sump = 0.0, sumv = 0.0;
        
        for (int i=0 ; i<plantCount ; i++)
        {
            double total = 0.0;
            for(int site=0 ; site<numPatch ; site++)
                total += pCurrent[i][site];
            
            if(total > viability) 
            {
                survPlants++;
                plantBiomass += total;
            }            
        }
        
        if (plantBiomass > 0)
        {
            for(int i=0 ; i<plantCount ; i++)
            {
                double total = 0.0;
                for(int site=0 ; site<numPatch ; site++)
                    total += pCurrent[i][site];
                if(total > viability)
                {
                    double prob = total / plantBiomass;
                    if (prob > 0) sump -= (prob * log(prob));
                }
            }
        }
        
        for (int i=0 ; i<insectCount ; i++)
        {
            double total = 0.0;
            for(int site=0 ; site<numPatch ; site++)
                total += vCurrent[i][site];
            
            if(total > viability) 
            {
                survInsects++;
                insectBiomass += total;  
            }          
        }
        
        if (insectBiomass > 0)
        {
            for(int i=0 ; i<insectCount ; i++)
            {
                double total = 0.0;
                for(int site=0 ; site<numPatch ; site++)
                    total += vCurrent[i][site];
                if(total > viability)
                {
                    double prob = total / insectBiomass;
                    if (prob > 0) sumv -= (prob * log(prob));
                }
            }
        }
        
        double robustness = (double)(survPlants + survInsects) / (1.0*totalInitialSpecies);
        double pollinationService = insectBiomass;
        double shannonP = sump;
        double shannonV = sumv;
        double totalBiomass = plantBiomass + insectBiomass;
        if (initialBiomass > 0)
            totalBiomass /= initialBiomass;
        
        sumRobustness += robustness;
        steps++;
        experimentFile << kEffective << " " << fixed << setprecision(6) << robustness << " " << survPlants << " " << survInsects << " " << pollinationService << " " << shannonP << " " << shannonV << " " << totalBiomass << endl;
    }
    
    double Rint = 0.0;
    if(steps > 0) Rint = sumRobustness / steps;
    cout << " R (robustness) = " << Rint << endl;
    
    ofstream rFile("robustnessR.txt");
    rFile << Rint << endl;
    rFile.close();
    
    experimentFile.close();
    fichp.close();
    fichv.close();
    
    cout << "---Extinction experiment complete ---" << endl;
    
    return;
}
/*
void runRandomExtinctionExperiment(const vector<vector<double>>& p, const vector<vector<double>>& v, const vector<vector<vector<double>>>& gamma, double h, int plantCount, int insectCount, int numPatch, double D)
{
    int numRealizations = 50;
    cout << "\n---Initializing Random extinction experiment ---" << endl;
    
    //Initial stationary state
    cout << "Computing initial stationary state..." << flush;
    
    vector<vector<double>> pSteadyInitial = p;
    vector<vector<double>> vSteadyInitial = v;
    double tDummy = 0.0;
    
    ofstream fichp("dummy_p.txt");
    ofstream fichv("dummy_v.txt");
    
    findSteadyState(tDummy,pSteadyInitial,vSteadyInitial,fichp,fichv,plantCount,insectCount,numPatch,gamma,h,D);
    
    int initialSurvPlants = 0;
    int initialSurvInsects = 0;
    
    for(int i=0; i<plantCount; ++i) 
    {
        double total = 0.0;
        for(int s=0; s<numPatch; ++s) 
            total += pSteadyInitial[i][s];
        if(total > viability) 
            initialSurvPlants++;
    }
    
    for(int i=0; i<insectCount; ++i) 
    {
        double total = 0.0;
        for(int s=0; s<numPatch; ++s) 
            total += vSteadyInitial[i][s];
        if(total > viability) 
            initialSurvInsects++;
    }
    
    int totalInitialSpecies = initialSurvPlants + initialSurvInsects;
    
    cout << "Especies Iniciales Vivas: " << totalInitialSpecies << " (Plantas: " << initialSurvPlants << ", Insectos: " << initialSurvInsects << ")" << endl;
    
    vector<double> avgRobustness(plantCount + 1, 0.0);
    vector<double> avgSurvPlants(plantCount + 1, 0.0);
    vector<double> avgSurvInsects(plantCount + 1, 0.0);
    vector<double> avgPollination(plantCount + 1, 0.0);
    vector<double> avgGiniP(plantCount + 1, 0.0);
    vector<double> avgGiniV(plantCount + 1, 0.0);
    
    for(int r=0 ; r<numRealizations ; r++)
    {
        cout << "Realization " << r + 1 << " / " << numRealizations << "...\r" << flush;
        
        //Plant ranking
        vector<int> plantIndices(plantCount);
        
        for(int i=0; i<plantCount ; i++)
            plantIndices[i] = i;
    
        random_device rd;
        auto rng = default_random_engine(rd());
        shuffle(plantIndices.begin(), plantIndices.end(), rng);
    
        //Preparation
        ofstream experimentFile("resultsRandom.txt");
    
        experimentFile << "# Num_Extinctions Robustness_Ratio Surv_Plants Surv_Insects Pollination_Service Gini_Plants Gini_Insects" << endl;
    
        vector<vector<double>> pCurrent = pSteadyInitial;
        vector<vector<double>> vCurrent = vSteadyInitial;
    
        double tLocal = tDummy;
    
        double sumRobustness = 0.0;
        int kEffective = 0, steps = 0;
    
        //Experiment
        for(int k=0; k<=plantCount ; k++)
        {
            if (k>0)
            {
                int plantToRemove = plantIndices[k-1];
                double currentAbundance = 0.0;
                
                for(int site=0 ; site<numPatch ; site++)
                    currentAbundance += pCurrent[plantToRemove][site];
            
                if (currentAbundance > viability)
                {                
                    for(int site=0 ; site<numPatch ; site++)
                        pCurrent[plantToRemove][site] = 0.0;
                
                    findSteadyState(tLocal, pCurrent, vCurrent, fichp, fichv, plantCount, insectCount, numPatch, gamma, h, D);
                }    
            }
        
            //Metrics
            int survPlants = 0;
            int survInsects = 0;
            double plantBiomass = 0.0, insectBiomass = 0.0;
            double sum2p = 0.0, sum2v = 0.0;
        
            for (int i=0 ; i<plantCount ; i++)
            {
                double total = 0.0;
                for(int site=0 ; site<numPatch ; site++)
                    total += pCurrent[i][site];
            
                if(total > viability) 
                {
                    survPlants++;
                    plantBiomass += total;
                }            
            }
        
            if (plantBiomass > 0)
            {
                for(int i=0 ; i<plantCount ; i++)
                {
                    double total = 0.0;
                    for(int site=0 ; site<numPatch ; site++)
                        total += pCurrent[i][site];
                    if(total > viability)
                    {
                        double prob = total / plantBiomass;
                        sum2p += (prob * prob);
                    }
                }
            }
            
            for (int i=0 ; i<insectCount ; i++)
            {
                double total = 0.0;
                for(int site=0 ; site<numPatch ; site++)
                    total += vCurrent[i][site];
                
                if(total > viability) 
                {
                    survInsects++;
                    insectBiomass += total;  
                }          
            }
        
            if (insectBiomass > 0)
            {
                for(int i=0 ; i<insectCount ; i++)
                {
                    double total = 0.0;
                    for(int site=0 ; site<numPatch ; site++)
                        total += vCurrent[i][site];
                    if(total > viability)
                    {
                        double prob = total / insectBiomass;
                        sum2v += (prob * prob);
                    }
                }
            }
        
            double robustness = (double)(survPlants + survInsects) / (totalInitialSpecies*1.0);
            double pollinationService = insectBiomass;
            double giniP = 1.0 - sum2p;
            double giniV = 1.0 - sum2v;
        
        sumRobustness += robustness;
        steps++;
        
        experimentFile << kEffective << " " << fixed << setprecision(6) << robustness << " " << survPlants << " " << survInsects << " " << pollinationService << " " << giniP << " " << giniV << endl;
    }
    
    double Rint = sumRobustness / steps;
    cout << " R (robustness) = " << Rint << endl;
    }
    ofstream rFile("robustnessD.txt");
    rFile << Rint << endl;
    rFile.close();
    
    experimentFile.close();
    fichp.close();
    fichv.close();
    
    cout << "---Extinction experiment complete ---" << endl;
    
    return;
}
*/
void runTargetExtinctionExperiment(const vector<vector<double>>& p, const vector<vector<double>>& v, const vector<vector<vector<double>>>& gamma, double h, int plantCount, int insectCount, int numPatch, double D)
{
    cout << "\n---Initializing Target extinction experiment ---" << endl;  
    
    //Preparation
    ofstream experimentFile("resultsTarget.txt");
    ofstream fichp("evolutionp.txt");
    ofstream fichv("evolutionv.txt");
    
    experimentFile << "# Num_Extinctions Robustness_Ratio Surv_Plants Surv_Insects Pollination_Service Gini_Plants Gini_Insects Norm_Biomass" << endl;

    //ofstream dummy_p("/dev/null");
    //ofstream dummy_v("/dev/null");
    
    vector<vector<double>> pCurrent = p;
    vector<vector<double>> vCurrent = v;
    
    double tDummy = 0.0;
    
    findSteadyState(tDummy,pCurrent,vCurrent,fichp,fichv,plantCount,insectCount,numPatch,gamma,h,D);
    
    int initialSurvPlants = 0;
    int initialSurvInsects = 0;
    double initialBiomass = 0.0;
    
    vector<bool> isDead(plantCount, false);
    
    for(int i=0; i<plantCount; ++i) 
    {
        double total = 0.0;
        for(int s=0; s<numPatch; ++s) 
            total += pCurrent[i][s];
            
        if(total > viability)
        {
            initialSurvPlants++;
            initialBiomass += total;
        }
        
        else isDead[i] = true;
    }
    
    for(int i=0; i<insectCount; ++i) 
    {
        double total = 0.0;
        for(int s=0; s<numPatch; ++s) 
            total += vCurrent[i][s];
            
        if(total > viability)
        {
            initialSurvInsects++;
            initialBiomass += total;
        }
    }
    
    int totalInitialSpecies = initialSurvPlants + initialSurvInsects;
    
    cout << "Especies Iniciales Vivas: " << totalInitialSpecies << " (Plantas: " << initialSurvPlants << ", Insectos: " << initialSurvInsects << ")" << endl;
    
    double sumRobustness = 0.0;
    int kEffective = 0, steps = 0;
    //Experiment
    for(int k=0 ; k<=plantCount ; k++)
    {
        if (k>0)
        {
            int plantToRemove = -1;
            double maxAbundance = -1.0;
            
            for(int i=0 ; i<plantCount ; i++)
            {
                if (isDead[i]) continue;
                
                double currentAbundance = 0.0;
                for(int s=0 ; s<numPatch ; s++) currentAbundance += pCurrent[i][s];
                
                if (currentAbundance <= viability)
                {
                    isDead[i] = true;
                    continue;
                }
                
                if (currentAbundance > maxAbundance)
                {
                    maxAbundance = currentAbundance;
                    plantToRemove = i;
                }
            }
            
            if (plantToRemove == -1) break;
            
            for(int s=0 ; s<numPatch ; s++) pCurrent[plantToRemove][s] = 0.0;
            
            isDead[plantToRemove] = true;
            
            kEffective++;
            
            findSteadyState(tDummy, pCurrent, vCurrent, fichp, fichv, plantCount, insectCount, numPatch, gamma, h,D);
        }
        
        //Metrics
        int survPlants = 0;
        int survInsects = 0;
        double plantBiomass = 0.0, insectBiomass = 0.0;
        double sum2p = 0.0, sum2v = 0.0;
        
        for (int i=0 ; i<plantCount ; i++)
        {
            double total = 0.0;
            for(int site=0 ; site<numPatch ; site++)
                total += pCurrent[i][site];
            
            if(total > viability) 
            {
                survPlants++;
                plantBiomass += total;
            }            
        }
        
        if (plantBiomass > 0)
        {
            for(int i=0 ; i<plantCount ; i++)
            {
                double total = 0.0;
                for(int site=0 ; site<numPatch ; site++)
                    total += pCurrent[i][site];
                if(total > viability)
                {
                    double prob = total / plantBiomass;
                    sum2p += (prob * prob);
                }
            }
        }
        
        for (int i=0 ; i<insectCount ; i++)
        {
            double total = 0.0;
            for(int site=0 ; site<numPatch ; site++)
                total += vCurrent[i][site];
            
            if(total > viability) 
            {
                survInsects++;
                insectBiomass += total;  
            }          
        }
        
        if (insectBiomass > 0)
        {
            for(int i=0 ; i<insectCount ; i++)
            {
                double total = 0.0;
                for(int site=0 ; site<numPatch ; site++)
                    total += vCurrent[i][site];
                if(total > viability)
                {
                    double prob = total / insectBiomass;
                    sum2v += (prob * prob);
                }
            }
        }
        
        double robustness = (double)(survPlants + survInsects) / (1.0*totalInitialSpecies);
        double pollinationService = insectBiomass;
        double giniP = 1.0 - sum2p;
        double giniV = 1.0 - sum2v;
        double totalBiomass = plantBiomass + insectBiomass;
        if (initialBiomass > 0)
            totalBiomass /= initialBiomass;
        
        sumRobustness += robustness;
        steps++;
        experimentFile << kEffective << " " << fixed << setprecision(6) << robustness << " " << survPlants << " " << survInsects << " " << pollinationService << " " << giniP << " " << giniV << " " << totalBiomass << endl;
    }
    
    double Rint = 0.0;
    if(steps > 0) Rint = sumRobustness / steps;
    cout << " R (robustness) = " << Rint << endl;
    
    ofstream rFile("robustnessR.txt");
    rFile << Rint << endl;
    rFile.close();
    
    experimentFile.close();
    fichp.close();
    fichv.close();
    
    cout << "---Extinction experiment complete ---" << endl;
    
    return;
}

void runLocalTargetExtinctionExperiment(const vector<vector<double>>& p, const vector<vector<double>>& v, const vector<vector<vector<double>>>& gamma, double h, int plantCount, int insectCount, int numPatch, double D)
{
    cout << "\n---Initializing Local Target extinction experiment ---" << endl;  
    
    //Preparation
    ofstream experimentFile("resultsLocalTarget.txt");
    ofstream fichp("evolutionp.txt");
    ofstream fichv("evolutionv.txt");
    
    experimentFile << "# Num_Extinctions Robustness_Ratio Surv_Plants Surv_Insects Pollination_Service Shannon_Plants Shannon_Insects Norm_Biomass" << endl;

    vector<vector<double>> pCurrent = p;
    vector<vector<double>> vCurrent = v;
    
    double tDummy = 0.0;
    
    findSteadyState(tDummy,pCurrent,vCurrent,fichp,fichv,plantCount,insectCount,numPatch,gamma,h,D);
    
    int initialSurvPlants = 0;
    int initialSurvInsects = 0;
    double initialBiomass = 0.0;
    
    // Ya no usamos isDead global, comprobaremos viabilidad localmente
    
    for(int i=0; i<plantCount; ++i) 
    {
        double total = 0.0;
        for(int s=0; s<numPatch; ++s) 
            total += pCurrent[i][s];
            
        if(total > viability)
        {
            initialSurvPlants++;
            initialBiomass += total;
        }
    }
    
    for(int i=0; i<insectCount; ++i) 
    {
        double total = 0.0;
        for(int s=0; s<numPatch; ++s) 
            total += vCurrent[i][s];
            
        if(total > viability)
        {
            initialSurvInsects++;
            initialBiomass += total;
        }
    }
    
    int totalInitialSpecies = initialSurvPlants + initialSurvInsects;
    
    cout << "Especies Iniciales Vivas: " << totalInitialSpecies << " (Plantas: " << initialSurvPlants << ", Insectos: " << initialSurvInsects << ")" << endl;
    
    double sumRobustness = 0.0;
    int kEffective = 0, steps = 0;
    
    // Experiment: Limite ampliado para permitir extinciones por población local
    int maxExtinctions = plantCount * numPatch;
    
    for(int k=0 ; k<=maxExtinctions ; k++)
    {
        if (k>0)
        {
            int plantToRemove = -1;
            int patchToRemove = -1;
            
            // TARGET ATTACK: Buscamos la MAXIMA abundancia local
            double maxAbundance = -1.0;
            
            for(int i=0 ; i<plantCount ; i++)
            {
                for(int s=0 ; s<numPatch ; s++)
                {
                    double localAbundance = pCurrent[i][s];
                    
                    // Solo consideramos poblaciones vivas
                    if (localAbundance > viability)
                    {
                        if (localAbundance > maxAbundance)
                        {
                            maxAbundance = localAbundance;
                            plantToRemove = i;
                            patchToRemove = s;
                        }
                    }
                }
            }
            
            // Si no quedan plantas vivas en ningún parche, paramos
            if (plantToRemove == -1) break;
            
            // ATAQUE LOCAL QUIRÚRGICO:
            // Eliminamos la planta MÁS abundante ÚNICAMENTE en su parche correspondiente
            pCurrent[plantToRemove][patchToRemove] = 0.0;
            
            kEffective++;
            
            findSteadyState(tDummy, pCurrent, vCurrent, fichp, fichv, plantCount, insectCount, numPatch, gamma, h, D);
        }
        
        //Metrics
        int survPlants = 0;
        int survInsects = 0;
        double plantBiomass = 0.0, insectBiomass = 0.0;
        
        // Cambiamos sum2p y sum2v (Gini) por shannonP y shannonV
        double shannonP = 0.0, shannonV = 0.0;
        
        for (int i=0 ; i<plantCount ; i++)
        {
            double total = 0.0;
            for(int site=0 ; site<numPatch ; site++)
                total += pCurrent[i][site];
            
            if(total > viability) 
            {
                survPlants++;
                plantBiomass += total;
            }            
        }
        
        if (plantBiomass > 0)
        {
            for(int i=0 ; i<plantCount ; i++)
            {
                double total = 0.0;
                for(int site=0 ; site<numPatch ; site++)
                    total += pCurrent[i][site];
                if(total > viability)
                {
                    double prob = total / plantBiomass;
                    if (prob > 0.0) shannonP -= (prob * log(prob)); // Índice de Shannon
                }
            }
        }
        
        for (int i=0 ; i<insectCount ; i++)
        {
            double total = 0.0;
            for(int site=0 ; site<numPatch ; site++)
                total += vCurrent[i][site];
            
            if(total > viability) 
            {
                survInsects++;
                insectBiomass += total;  
            }          
        }
        
        if (insectBiomass > 0)
        {
            for(int i=0 ; i<insectCount ; i++)
            {
                double total = 0.0;
                for(int site=0 ; site<numPatch ; site++)
                    total += vCurrent[i][site];
                if(total > viability)
                {
                    double prob = total / insectBiomass;
                    if (prob > 0.0) shannonV -= (prob * log(prob)); // Índice de Shannon
                }
            }
        }
        
        double robustness = (double)(survPlants + survInsects) / (1.0*totalInitialSpecies);
        double pollinationService = insectBiomass;
        double shannon_P = shannonP; // Renombrado para mantener consistencia
        double shannon_V = shannonV;
        double totalBiomass = plantBiomass + insectBiomass;
        if (initialBiomass > 0)
            totalBiomass /= initialBiomass;
        
        sumRobustness += robustness;
        steps++;
        experimentFile << kEffective << " " << fixed << setprecision(6) << robustness << " " << survPlants << " " << survInsects << " " << pollinationService << " " << shannon_P << " " << shannon_V << " " << totalBiomass << endl;
    }
    
    double Rint = 0.0;
    if(steps > 0) Rint = sumRobustness / steps;
    cout << " R (robustness) = " << Rint << endl;
    
    ofstream rFile("robustnessR.txt");
    rFile << Rint << endl;
    rFile.close();
    
    experimentFile.close();
    fichp.close();
    fichv.close();
    
    cout << "---Extinction experiment complete ---" << endl;
    
    return;
}
