#include <gtest/gtest.h>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>
#include "../src/PhenotypeIntegration.h"
using namespace std;

class semSimAndICTest : public ::testing::Test
{
protected:
    void SetUp() override {}
    void TearDown() override {}
};

TEST_F(semSimAndICTest, TestSemSimAndICTest)
{
    PhenotypeIntegration phen;

    phen.init();
    phen.setupHPO("./HPO/HPO_release20220414/");

    // Make comparisons where two HPO terms are one the ancestor of the other, so that the ancestor is automatically the MICA. Consequently, the following equalities have to hold: IC(MICA)=IC(ancestor)=ResnikSim=ERICSim.

    // "HP:0002450", "HP:0012639" (Ancestor/MICA)
    // "HP:0032794", "HP:0001250" (Ancestor/MICA)
    // "HP:0009099", "HP:0000175" (Ancestor/MICA)
    // "HP:0009099", "HP:0000163" (Ancestor/MICA)

    vector<string> HPOSet1;
    vector<string> HPOSet_2;
    double IC;
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

    HPOSet1 = {"HP:0002450"};
    HPOSet_2 = {"HP:0012639"};
    IC = phen.get_InformationContent(HPOSet_2[0]);
    maxRes = phen.get_calcSemSimSet(HPOSet1, HPOSet_2)[0];
    maxEric = phen.get_calcSemSimSet(HPOSet1, HPOSet_2)[1];
    meanRes = phen.get_calcSemSimSet(HPOSet1, HPOSet_2)[2];
    meanEric = phen.get_calcSemSimSet(HPOSet1, HPOSet_2)[3];
    BMARes = phen.get_calcSemSimSet(HPOSet1, HPOSet_2)[4];
    BMAEric = phen.get_calcSemSimSet(HPOSet1, HPOSet_2)[5];
    BMSRes = phen.get_calcSemSimSet(HPOSet1, HPOSet_2)[6];
    BMSEric = phen.get_calcSemSimSet(HPOSet1, HPOSet_2)[7];
    FunSimAvRes = phen.get_calcSemSimSet(HPOSet1, HPOSet_2)[8];
    FunSimAvEric = phen.get_calcSemSimSet(HPOSet1, HPOSet_2)[9];
    FunSimMaxRes = phen.get_calcSemSimSet(HPOSet1, HPOSet_2)[10];
    FunSimMaxEric = phen.get_calcSemSimSet(HPOSet1, HPOSet_2)[11];
    // cout<<"IC = "<<IC<<endl;
    // cout<<"Resnik = "<<Resnik<<endl;
    // cout<<"ERIC = "<<ERIC<<endl;
    ASSERT_EQ(IC, maxRes);
    ASSERT_EQ(maxRes, maxEric);
    ASSERT_EQ(maxEric, meanRes);
    ASSERT_EQ(meanRes, meanEric);
    ASSERT_EQ(meanEric, BMARes);
    ASSERT_EQ(BMARes, BMAEric);
    ASSERT_EQ(BMAEric, BMSRes);
    ASSERT_EQ(BMSRes, BMSEric);
    ASSERT_EQ(BMSEric, FunSimAvRes);
    ASSERT_EQ(FunSimAvRes, FunSimAvEric);
    ASSERT_EQ(FunSimAvEric, FunSimMaxRes);
    ASSERT_EQ(FunSimMaxRes, FunSimMaxEric);

    HPOSet1 = {"HP:0032794"};
    HPOSet_2 = {"HP:0001250"};
    IC = phen.get_InformationContent(HPOSet_2[0]);
    maxRes = phen.get_calcSemSimSet(HPOSet1, HPOSet_2)[0];
    maxEric = phen.get_calcSemSimSet(HPOSet1, HPOSet_2)[1];
    meanRes = phen.get_calcSemSimSet(HPOSet1, HPOSet_2)[2];
    meanEric = phen.get_calcSemSimSet(HPOSet1, HPOSet_2)[3];
    BMARes = phen.get_calcSemSimSet(HPOSet1, HPOSet_2)[4];
    BMAEric = phen.get_calcSemSimSet(HPOSet1, HPOSet_2)[5];
    BMSRes = phen.get_calcSemSimSet(HPOSet1, HPOSet_2)[6];
    BMSEric = phen.get_calcSemSimSet(HPOSet1, HPOSet_2)[7];
    FunSimAvRes = phen.get_calcSemSimSet(HPOSet1, HPOSet_2)[8];
    FunSimAvEric = phen.get_calcSemSimSet(HPOSet1, HPOSet_2)[9];
    FunSimMaxRes = phen.get_calcSemSimSet(HPOSet1, HPOSet_2)[10];
    FunSimMaxEric = phen.get_calcSemSimSet(HPOSet1, HPOSet_2)[11];
    ASSERT_EQ(IC, maxRes);
    ASSERT_EQ(maxRes, maxEric);
    ASSERT_EQ(maxEric, meanRes);
    ASSERT_EQ(meanRes, meanEric);
    ASSERT_EQ(meanEric, BMARes);
    ASSERT_EQ(BMARes, BMAEric);
    ASSERT_EQ(BMAEric, BMSRes);
    ASSERT_EQ(BMSRes, BMSEric);
    ASSERT_EQ(BMSEric, FunSimAvRes);
    ASSERT_EQ(FunSimAvRes, FunSimAvEric);
    ASSERT_EQ(FunSimAvEric, FunSimMaxRes);
    ASSERT_EQ(FunSimMaxRes, FunSimMaxEric);

    HPOSet1 = {"HP:0009099"};
    HPOSet_2 = {"HP:0000175"};
    IC = phen.get_InformationContent(HPOSet_2[0]);
    maxRes = phen.get_calcSemSimSet(HPOSet1, HPOSet_2)[0];
    maxEric = phen.get_calcSemSimSet(HPOSet1, HPOSet_2)[1];
    meanRes = phen.get_calcSemSimSet(HPOSet1, HPOSet_2)[2];
    meanEric = phen.get_calcSemSimSet(HPOSet1, HPOSet_2)[3];
    BMARes = phen.get_calcSemSimSet(HPOSet1, HPOSet_2)[4];
    BMAEric = phen.get_calcSemSimSet(HPOSet1, HPOSet_2)[5];
    BMSRes = phen.get_calcSemSimSet(HPOSet1, HPOSet_2)[6];
    BMSEric = phen.get_calcSemSimSet(HPOSet1, HPOSet_2)[7];
    FunSimAvRes = phen.get_calcSemSimSet(HPOSet1, HPOSet_2)[8];
    FunSimAvEric = phen.get_calcSemSimSet(HPOSet1, HPOSet_2)[9];
    FunSimMaxRes = phen.get_calcSemSimSet(HPOSet1, HPOSet_2)[10];
    FunSimMaxEric = phen.get_calcSemSimSet(HPOSet1, HPOSet_2)[11];
    ASSERT_EQ(IC, maxRes);
    ASSERT_EQ(maxRes, maxEric);
    ASSERT_EQ(maxEric, meanRes);
    ASSERT_EQ(meanRes, meanEric);
    ASSERT_EQ(meanEric, BMARes);
    ASSERT_EQ(BMARes, BMAEric);
    ASSERT_EQ(BMAEric, BMSRes);
    ASSERT_EQ(BMSRes, BMSEric);
    ASSERT_EQ(BMSEric, FunSimAvRes);
    ASSERT_EQ(FunSimAvRes, FunSimAvEric);
    ASSERT_EQ(FunSimAvEric, FunSimMaxRes);
    ASSERT_EQ(FunSimMaxRes, FunSimMaxEric);

    HPOSet1 = {"HP:0009099"};
    HPOSet_2 = {"HP:0000163"};
    IC = phen.get_InformationContent(HPOSet_2[0]);
    maxRes = phen.get_calcSemSimSet(HPOSet1, HPOSet_2)[0];
    maxEric = phen.get_calcSemSimSet(HPOSet1, HPOSet_2)[1];
    meanRes = phen.get_calcSemSimSet(HPOSet1, HPOSet_2)[2];
    meanEric = phen.get_calcSemSimSet(HPOSet1, HPOSet_2)[3];
    BMARes = phen.get_calcSemSimSet(HPOSet1, HPOSet_2)[4];
    BMAEric = phen.get_calcSemSimSet(HPOSet1, HPOSet_2)[5];
    BMSRes = phen.get_calcSemSimSet(HPOSet1, HPOSet_2)[6];
    BMSEric = phen.get_calcSemSimSet(HPOSet1, HPOSet_2)[7];
    FunSimAvRes = phen.get_calcSemSimSet(HPOSet1, HPOSet_2)[8];
    FunSimAvEric = phen.get_calcSemSimSet(HPOSet1, HPOSet_2)[9];
    FunSimMaxRes = phen.get_calcSemSimSet(HPOSet1, HPOSet_2)[10];
    FunSimMaxEric = phen.get_calcSemSimSet(HPOSet1, HPOSet_2)[11];
    ASSERT_EQ(IC, maxRes);
    ASSERT_EQ(maxRes, maxEric);
    ASSERT_EQ(maxEric, meanRes);
    ASSERT_EQ(meanRes, meanEric);
    ASSERT_EQ(meanEric, BMARes);
    ASSERT_EQ(BMARes, BMAEric);
    ASSERT_EQ(BMAEric, BMSRes);
    ASSERT_EQ(BMSRes, BMSEric);
    ASSERT_EQ(BMSEric, FunSimAvRes);
    ASSERT_EQ(FunSimAvRes, FunSimAvEric);
    ASSERT_EQ(FunSimAvEric, FunSimMaxRes);
    ASSERT_EQ(FunSimMaxRes, FunSimMaxEric);
}
