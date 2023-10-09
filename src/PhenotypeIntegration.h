#ifndef READ_H
#define READ_H

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <math.h>
#include <algorithm>

using namespace std;

class PhenotypeIntegration
{
private:
  string inpObo;
  string inpVarTab;
  string inpUserHPOs;
  string inpGenToPhen;
  string outGenNamScoreDisID;
  double mapGenToOnt;
  double mapGenToHPOkidz;
  double maxRes;
  double maxEric;
  double meanRes;
  double meanEric;
  double BMARes;
  double BMAEric;
  double BMSRes;
  double BMSEric;
  double FunSimAvRes;
  double FunSimAvEric;
  double FunSimMaxRes;
  double FunSimMaxEric;
  vector<string> vPathUp;
  vector<string> vPathDown;
  vector<string> vPathDownSortUniq;
  vector<string> vHPOs;
  vector<string> vInpUserHPOs;
  vector<string> vUTInpUserHPOs;
  vector<string> vInpVarTabGenName;
  vector<double> vInpVarTabScore;
  vector<string> vGenToPhenGenID;
  vector<string> vGenToPhenDisHPOs;
  vector<string> vGenToPhenGenName;
  vector<string> vGenToPhenDisID;
  vector<string> vGenToPhenHPOterm;
  vector<string> vGenToPhenHPOtermName;
  vector<string> vUnderontGenTerm;
  vector<vector<string>> vvGenToPhenDisHPOs;
  vector<vector<string>> vvGenNamHPOsDiseaseID;
  map<string, vector<string>> mHPOandPars;
  map<string, vector<string>> mHPOandKids;
  map<string, vector<string>> mHPOsGen;
  map<string, vector<vector<string>>> mGeneNamDiseHPOs;
  map<string, vector<double>> mHPOsSign;

  void readObo();
  void readVarTab();
  void readUserHPOs();
  void frameObo();
  void readGenToPhen();
  void frameGenToPhen();
  void treeWalkUp(string qHPO);
  void treeWalkDown(string qHPO);
  void uniqNodAllDesc();
  void findHPOassGen(string qHPO);
  void calcSemSimSet(vector<string> HPOSet1, vector<string> HPOSet_2);
  void getAssHPOs(string genNam);
  void preCompICs();
  void buildOut(string UTgenNam); 
  string findMICA(string qHPO1, string qHPO2);
  double calcIC(string qHPO);

public:
  PhenotypeIntegration()
  {
    init();
  }
  ~PhenotypeIntegration() {}

  void setupHPO(const string &pathHPO)
  {
    string fnObo = pathHPO + "/hp.obo";
    inpObo = fnObo;
    readObo();

    string fnGenPhen = pathHPO + "/genes_to_phenotype.txt";
    inpGenToPhen = fnGenPhen;
    readGenToPhen();
  }

  void read_variant_table(const string &pathVarTab)
  {
    string fnVarTab = pathVarTab;
    inpVarTab = fnVarTab;
    readVarTab();
  }

  void read_user_HPOs(const string &pathUserHPOs)
  {
    string fnUserHPOs = pathUserHPOs;
    inpUserHPOs = fnUserHPOs;
    readUserHPOs();
  }

  vector<string> getInpContObo()
  {
    return vHPOs;
  }
  vector<string> getInpContGenPhen()
  {
    return vGenToPhenHPOterm;
  }

  void init()
  {
    outGenNamScoreDisID.clear();
    inpObo.clear();
    vInpUserHPOs.clear();
    vUTInpUserHPOs.clear();
    vHPOs.clear();
    vGenToPhenGenID.clear();
    vGenToPhenGenName.clear();
    vGenToPhenDisID.clear();
    vGenToPhenHPOterm.clear();
    vGenToPhenHPOtermName.clear();
    vGenToPhenDisHPOs.clear();
    vInpVarTabGenName.clear();
    vInpVarTabScore.clear();
    vvGenToPhenDisHPOs.clear();
    vvGenNamHPOsDiseaseID.clear();
    vUnderontGenTerm.clear();
    mHPOandPars.clear();
    mHPOandKids.clear();
    mHPOsGen.clear();
    mGeneNamDiseHPOs.clear();
    mHPOsSign.clear();
    vPathUp.clear();
    vPathDown.clear();
    vPathDownSortUniq.clear();
    mapGenToOnt = -1;
    mapGenToHPOkidz = -1;
    maxRes = -1;
    maxEric = -1;
    meanRes = -1;
    meanEric = -1;
    BMARes = -1;
    BMAEric = -1;
    BMSRes = -1;
    BMSEric = -1;
    FunSimAvRes = -1;
    FunSimAvEric = -1;
    FunSimMaxRes = -1;
    FunSimMaxEric = -1;
  }

  void checkQueryInHPOdata(vector<string> HPOSet);
  void printTermPar();
  void printTermKids();
  void printInpNetForCytoscape();
  void queryTermForPar(string qHPO);
  void queryTermForKids(string qHPO);

  void set_treeWalkUp(string qHPO)
  {
    vPathUp.clear();
    vPathDown.clear();
    vPathUp.push_back(qHPO);
    treeWalkUp(qHPO);
  }
  vector<string> get_treeWalkUp()
  {
    return vPathUp;
  }
  void showPath_treeWalkUp()
  {
    for (unsigned int i = 0; i < vPathUp.size(); i++)
    {
      cout << " " << vPathUp[i] << endl;
    }
    cout << "" << endl;
  }

  void set_treeWalkDown(string qHPO)
  {
    vPathDown.clear();
    vPathUp.clear();
    vPathDown.push_back(qHPO);
    treeWalkDown(qHPO);
    uniqNodAllDesc(); // calc number of descendants here
  }
  vector<string> get_treeWalkDown()
  {
    return vPathDown;
  }
  void showPath_treeWalkDown()
  {
    for (unsigned int i = 0; i < vPathDown.size(); i++)
    {
      cout << " " << vPathDown[i] << endl;
    }
    cout << "" << endl;
  }

  string get_findMICA(string qHPO1, string qHPO2)
  {
    string MICA = findMICA(qHPO1, qHPO2);
    return MICA;
  }

  void set_uniqNodAllDesc() // to get access from unit tests
  {
    uniqNodAllDesc();
  }
  vector<string> get_uniqNodAllDesc()
  {
    return vPathDownSortUniq;
  }
  int get_uniqNodAllDesc_numberOf()
  {
    return vPathDownSortUniq.size();
  }

  void set_findHPOassGen(string qHPO)
  {
    findHPOassGen(qHPO);
  }
  int get_findHPOassGen()
  {
    return mapGenToHPOkidz;
  }


  void set_getAssHPOs(string genNam){
    getAssHPOs(genNam);
  }
  vector<vector<string>> get_getAssHPOs()
  {
    return vvGenNamHPOsDiseaseID;
  }

  void set_buildOut_unit_test_inpUserHPOs(vector<string>  UnitTest_UserInput){
    vUTInpUserHPOs = UnitTest_UserInput;
  }
  void set_buildOut(string UTgenNam){
    buildOut(UTgenNam);
  }
  string get_buildOut(){
    return outGenNamScoreDisID;
  }

  double get_InformationContent(string qHPO)
  {
    map<string, vector<double>>::iterator it;
    it = mHPOsSign.find(qHPO);
    if (mHPOsSign.count(qHPO) == 0)
    {
      // base case: reached the bottom of the tree
      return 0;
    }
    else
    {
      // cout<<"n: "<<(*it).second[0]<<endl;
      // cout<<"N: "<<(*it).second[1]<<endl;
      // cout<<"IC: "<<(*it).second[2]<<endl;
      return (*it).second[2];
    }
  }

  vector<double> get_calcSemSimSet(vector<string> HPOSet1, vector<string> HPOSet_2)
  {
    calcSemSimSet(HPOSet1, HPOSet_2);
    vector<double> SemSim;
    SemSim = {
        maxRes,
        maxEric,
        meanRes,
        meanEric,
        BMARes,
        BMAEric,
        BMSRes,
        BMSEric,
        FunSimAvRes,
        FunSimAvEric,
        FunSimMaxRes,
        FunSimMaxEric};
    return SemSim;
  }

  void set_preCompICs()
  {
    preCompICs();
  }
};

#endif
