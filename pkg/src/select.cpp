#include "select.hpp"

Select::Select(){}

Select::Select(Vect v, CritClust b, SelectReg sReg, int packSize)
{
    this->v = v;
    this->b = b;
    this->sReg = sReg;
    this->packSize = packSize;
}

Select::Select(Vect v, SelectReg sReg, int packSize)
{
    this->v = v;
    this->sReg = sReg;
    this->packSize = packSize;
}

List Select::selectS(vector<int> Order)
{
    int  InitialProjectsNb = this->v.experiments.size();
    //string m;
    double CritValue;
    const int numeromodeleaux = 1;
    int firstIndex, lastIndex, idx;
    firstIndex = 1;
    lastIndex = firstIndex + packSize;
    int ClustVar = 1000;
    vector<int> aux, varSelectReg, varSelectClust_aux, varSelectClust;
    double critClustaux = 0.0, critDiffClust = 0.0;
    List Mylist, Mylistaux;
    
    varSelectClust.clear();
    varSelectClust.push_back(Order[0]);
    Mylist = b.ClustBestModel(varSelectClust);
    CritValue = as<double>(Mylist["criterionValue"]);
    
    while((ClustVar > 0) && (firstIndex < (int)Order.size()))
    {
        ClustVar = 0;
        for(idx = firstIndex; idx < lastIndex; ++idx)
        {
            if(idx < (int)Order.size())
            {
                aux.clear(); varSelectReg.clear(); varSelectClust_aux.clear();
                aux.push_back(Order[idx]);
                varSelectReg = (this->sReg).selectReg(varSelectClust,aux,InitialProjectsNb);
                varSelectClust_aux = (this->v).ajouter_var(varSelectClust,aux);
                Mylistaux = b.ClustBestModel(varSelectClust_aux);
                critClustaux = as<double>(Mylistaux["criterionValue"]);
                critDiffClust = critClustaux - CritValue - v.bicReggen(aux, varSelectReg, numeromodeleaux);
                if(critDiffClust > 0)
                {
                    varSelectClust =  varSelectClust_aux;
                    Mylist = Mylistaux;
                    CritValue = as<double>(Mylist["criterionValue"]);
                    ClustVar++;
                    
                }
            }
        }
        firstIndex = lastIndex;
        lastIndex += packSize;
        
    }
    
    //cout << " .... S is selected .... " << endl;
    return List::create(Named("S") = wrap(varSelectClust),
                        Named("model") = Mylist["model"],
                        Named("criterionValue") = Mylist["criterionValue"],
                        Named("criterion") = Mylist["criterion"],
                        Named("nbCluster") = Mylist["nbCluster"],
                        Named("proba") = Mylist["proba"],
                        Named("partition") = Mylist["partition"]);
}

vector<int> Select::selectW(vector<int> Order, vector<int>OtherVar)
{
    vector<int> varIndep;
    int InitialProjectsNb = this->v.experiments.size();
    int firstIndex, lastIndex, idx;
    lastIndex = Order.size();
    firstIndex = lastIndex - packSize;
    int ClustVar = 1000;
    
    //vector<int> OtherVar, varIndep_Role(Order.size()), varRegAux, aux;
    vector<int> varIndep_Role(Order.size()), varRegAux, aux;
    for(int l = 0; l < (int)Order.size(); ++l)
        varIndep_Role[l] = 0;
    
    varIndep.clear();
    while((ClustVar > 0) && (firstIndex >= 0))
    {
        ClustVar =0;
        for(idx = (lastIndex - 1); idx >= firstIndex; idx--)
        {
            if(idx >= 0)
            {
                varRegAux.clear(); aux.clear();
                //OtherVar.clear();
                //for(int l = 0; l < idx; ++l)
                //    if((varIndep_Role[l] == 0) && (l != idx))
                //        OtherVar.push_back(Order[l]);
                
                aux.push_back(Order[idx]);
                varRegAux=sReg.selectReg(OtherVar,aux,InitialProjectsNb);              	    
                if(varRegAux.empty())
                {
                    varIndep.push_back(Order[idx]);
                    ClustVar++;
                    varIndep_Role[idx] = 1; 
                } 
                
            } 
        } 
        lastIndex = firstIndex; 
        firstIndex -= packSize; 
    }
    //cout << " .... W is selected .... " << endl;
    return(varIndep);
};




