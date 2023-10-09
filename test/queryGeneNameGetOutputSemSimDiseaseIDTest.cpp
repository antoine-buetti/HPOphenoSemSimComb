#include <gtest/gtest.h>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>
#include "../src/PhenotypeIntegration.h"
using namespace std;

// Define a test fixture class for the Read class
class queryGeneNameGetOutputSemSimDiseaseIDTest : public ::testing::Test
{
protected:
    void SetUp() override {}
    void TearDown() override {}
};

TEST_F(queryGeneNameGetOutputSemSimDiseaseIDTest, TestQueryGeneNameGetOutputSemSimDiseaseIDTest)
{

    string actual_output;
    string expected_output;
    string str;
    vector<string> UserInput;
    string GeneNameInput;
    PhenotypeIntegration phen;

    phen.init();
    phen.setupHPO("./HPO/HPO_release20220414/");


    GeneNameInput = "MKKS";
    UserInput = {"HP:0001249", "HP:0000003", "HP:0000135", "HP:0006101", "HP:0000426", "HP:0000822", "HP:0008736", "HP:0001513", "HP:0000470", "HP:0000368", "HP:0000494", "HP:0004322", "HP:0000639", "HP:0002167", "HP:0001162", "HP:0001395", "HP:0000100", "HP:0000365"};

    phen.set_buildOut_unit_test_inpUserHPOs(UserInput);
    phen.set_buildOut(GeneNameInput);
    str = phen.get_buildOut();

    actual_output = str;
    expected_output = "MKKS { 2.810634 [ORPHA:110 (Generalized hirsutism)] | 1.426470 [OMIM:605231 (Polydactyly)] | 1.112169 [ORPHA:2473 (High palate)] | 0.807745 [OMIM:236700 (Syndactyly)] }";
    EXPECT_EQ(expected_output, actual_output);





    GeneNameInput = "MKKS";
    UserInput = {"HP:0000047", "HP:0000510", "HP:0001249", "HP:0003241", "HP:0000819", "HP:0000007", "HP:0000107", "HP:0001513", "HP:0001159", "HP:0010442"};

    phen.set_buildOut_unit_test_inpUserHPOs(UserInput);
    phen.set_buildOut(GeneNameInput);
    str = phen.get_buildOut();

    actual_output = str;
    expected_output = "MKKS { 2.919124 [OMIM:605231 (Polydactyly)] | 1.262325 [ORPHA:2473 (High palate)] | 1.245842 [ORPHA:110 (Generalized hirsutism)] | 1.203254 [OMIM:236700 (Syndactyly)] }";
    EXPECT_EQ(expected_output, actual_output);




    GeneNameInput = "TSR2";
    UserInput = {"HP:0001419", "HP:0011800", "HP:0000653", "HP:0001972", "HP:0000494", "HP:0000347", "HP:0000405", "HP:0008551"};

    phen.set_buildOut_unit_test_inpUserHPOs(UserInput);
    phen.set_buildOut(GeneNameInput);
    str = phen.get_buildOut();

    actual_output = str;
    expected_output = "TSR2 { 3.686166 [OMIM:300946 (Microtia)] | 0.565993 [ORPHA:124 (High palate)] }";
    EXPECT_EQ(expected_output, actual_output);

    


    phen.init();
    phen.setupHPO("./HPO/HPO_release20230127/");


    GeneNameInput = "MKKS";
    UserInput = {"HP:0000047", "HP:0000510", "HP:0001249", "HP:0003241", "HP:0000819", "HP:0000007", "HP:0000107", "HP:0001513", "HP:0001159", "HP:0010442"};

    phen.set_buildOut_unit_test_inpUserHPOs(UserInput);
    phen.set_buildOut(GeneNameInput);
    str = phen.get_buildOut();

    actual_output = str;
    expected_output = "MKKS { 2.583563 [OMIM:605231 (Pigmentary retinopathy)] | 1.247046 [ORPHA:2473 (Urethral stricture)] | 1.234232 [ORPHA:110 (Skeletal muscle atrophy)] | 1.190181 [OMIM:236700 (Hydrometrocolpos)] }";
    EXPECT_EQ(expected_output, actual_output);




    GeneNameInput = "TSR2";
    UserInput = {"HP:0001419", "HP:0011800", "HP:0000653", "HP:0001972", "HP:0000494", "HP:0000347", "HP:0000405", "HP:0008551"};

    phen.set_buildOut_unit_test_inpUserHPOs(UserInput);
    phen.set_buildOut(GeneNameInput);
    str = phen.get_buildOut();

    actual_output = str;
    expected_output = "TSR2 { 2.575602 [OMIM:300946 (Cleft palate)] | 0.556915 [ORPHA:124 (Absent thumb)] }";
    EXPECT_EQ(expected_output, actual_output);

    


}
