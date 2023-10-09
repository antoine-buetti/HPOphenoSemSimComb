#include <gtest/gtest.h>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>
#include "../src/PhenotypeIntegration.h"
using namespace std;

// Define a test fixture class for the Read class
class getAssHPOsTest : public ::testing::Test
{
protected:
    void SetUp() override {}
    void TearDown() override {}
};

TEST_F(getAssHPOsTest, TestgetAssHPOsTest)
{

    vector<string> actual_output;
    vector<string> expected_output;
    vector<vector<string>> _vvGenNamHPOsDiseaseID;

    PhenotypeIntegration phen;

    phen.init();
    phen.setupHPO("./HPO/HPO_release20220414/");

    phen.set_getAssHPOs("MKKS");
    _vvGenNamHPOsDiseaseID = phen.get_getAssHPOs();

    actual_output = _vvGenNamHPOsDiseaseID[0];
    expected_output = {"HP:0000047", "HP:0000510", "HP:0001249", "HP:0003241", "HP:0000819", "HP:0000007", "HP:0000107", "HP:0001513", "HP:0001159", "HP:0010442", "MKKS", "OMIM:605231", "Polydactyly"};
    EXPECT_EQ(expected_output, actual_output);

    actual_output = _vvGenNamHPOsDiseaseID[1];
    expected_output = {"HP:0001162", "HP:0008678", "HP:0000028", "HP:0006101", "HP:0000003", "HP:0004322", "HP:0008368", "HP:0001830", "HP:0001636", "HP:0000126", "HP:0002251", "HP:0002023", "HP:0100779", "HP:0030010", "HP:0001643", "HP:0004383", "HP:0000807", "HP:0001631", "HP:0005916", "HP:0001263", "HP:0012227", "HP:0001629", "HP:0004397", "HP:0001508", "HP:0000175", "HP:0001249", "HP:0001156", "HP:0000218", "MKKS", "ORPHA:2473", "High_palate"};
    EXPECT_EQ(expected_output, actual_output);

    actual_output = _vvGenNamHPOsDiseaseID[2];
    expected_output = {"HP:0001249", "HP:0000003", "HP:0000135", "HP:0006101", "HP:0000426", "HP:0000822", "HP:0008736", "HP:0001513", "HP:0000470", "HP:0000368", "HP:0000494", "HP:0004322", "HP:0000639", "HP:0002167", "HP:0001162", "HP:0001395", "HP:0000100", "HP:0000365", "HP:0000028", "HP:0003202", "HP:0008724", "HP:0000512", "HP:0010747", "HP:0000580", "HP:0002230", "MKKS", "ORPHA:110", "Generalized_hirsutism"};
    EXPECT_EQ(expected_output, actual_output);

    actual_output = _vvGenNamHPOsDiseaseID[3];
    expected_output = {"HP:0000145", "HP:0001374", "HP:0000113", "HP:0000969", "HP:0000007", "HP:0002089", "HP:0002023", "HP:0010741", "HP:0000143", "HP:0001586", "HP:0001162", "HP:0000028", "HP:0000148", "HP:0006159", "HP:0030010", "HP:0000072", "HP:0002251", "HP:0000126", "HP:0030680", "HP:0001159", "MKKS", "OMIM:236700", "Syndactyly"};
    EXPECT_EQ(expected_output, actual_output);

    
    phen.set_getAssHPOs("TSR2");
    _vvGenNamHPOsDiseaseID = phen.get_getAssHPOs();


    actual_output = _vvGenNamHPOsDiseaseID[0];
    expected_output = {"HP:0002863", "HP:0001199", "HP:0000104", "HP:0001631", "HP:0001518", "HP:0020118", "HP:0001790", "HP:0000431", "HP:0012410", "HP:0001873", "HP:0000185", "HP:0000085", "HP:0000912", "HP:0012758", "HP:0000508", "HP:0002669", "HP:0001254", "HP:0000252", "HP:0001894", "HP:0001087", "HP:0001895", "HP:0000486", "HP:0005518", "HP:0005532", "HP:0001896", "HP:0040276", "HP:0000047", "HP:0009944", "HP:0000980", "HP:0000286", "HP:0000294", "HP:0000347", "HP:0006758", "HP:0410030", "HP:0001227", "HP:0000470", "HP:0030270", "HP:0004808", "HP:0001680", "HP:0011904", "HP:0000316", "HP:0001875", "HP:0005280", "HP:0012133", "HP:0000465", "HP:0001629", "HP:0000519", "HP:0008551", "HP:0000369", "HP:0001882", "HP:0004322", "HP:0009778", "HP:0009777", "HP:0000218", "TSR2", "ORPHA:124", "High_palate"};
    EXPECT_EQ(expected_output, actual_output);

    actual_output = _vvGenNamHPOsDiseaseID[1];
    expected_output = {"HP:0001419", "HP:0011800", "HP:0000653", "HP:0001972", "HP:0000494", "HP:0000347", "HP:0000405", "HP:0008551", "TSR2", "OMIM:300946", "Microtia"};
    EXPECT_EQ(expected_output, actual_output);

    phen.init();
    phen.setupHPO("./HPO/HPO_release20230127/");

    phen.set_getAssHPOs("MKKS");
    _vvGenNamHPOsDiseaseID = phen.get_getAssHPOs();

    actual_output = _vvGenNamHPOsDiseaseID[0];
    expected_output = {"HP:0000218", "HP:0001156", "HP:0001830", "HP:0001508", "HP:0008678", "HP:0006101", "HP:0001643", "HP:0004383", "HP:0000807", "HP:0001629", "HP:0000126", "HP:0008368", "HP:0001263", "HP:0000003", "HP:0001636", "HP:0004397", "HP:0001631", "HP:0005916", "HP:0002023", "HP:0000175", "HP:0002251", "HP:0001249", "HP:0100779", "HP:0001162", "HP:0004322", "HP:0000028", "HP:0030010", "HP:0012227", "MKKS", "ORPHA:2473", "Urethral_stricture"};
    EXPECT_EQ(expected_output, actual_output);

    actual_output = _vvGenNamHPOsDiseaseID[1];
    expected_output = {"HP:0000126", "HP:0000143", "HP:0001586", "HP:0006159", "HP:0002023", "HP:0000072", "HP:0000007", "HP:0000113", "HP:0001159", "HP:0010741", "HP:0000145", "HP:0000969", "HP:0030680", "HP:0002089", "HP:0000028", "HP:0001374", "HP:0000148", "HP:0001162", "HP:0002251", "HP:0030010", "MKKS", "OMIM:236700", "Hydrometrocolpos"};
    EXPECT_EQ(expected_output, actual_output);

    actual_output = _vvGenNamHPOsDiseaseID[2];
    expected_output = {"HP:0000007", "HP:0000107", "HP:0001249", "HP:0001159", "HP:0001513", "HP:0003241", "HP:0100259", "HP:0000819", "HP:0003577", "HP:0000047", "HP:0000148", "HP:0000510", "HP:0000580", "MKKS", "OMIM:605231", "Pigmentary_retinopathy"};
    EXPECT_EQ(expected_output, actual_output);

    actual_output = _vvGenNamHPOsDiseaseID[3];
    expected_output = {"HP:0000470", "HP:0002167", "HP:0000822", "HP:0008724", "HP:0000580", "HP:0000368", "HP:0000100", "HP:0001513", "HP:0000365", "HP:0001162", "HP:0004322", "HP:0000512", "HP:0001395", "HP:0001249", "HP:0000426", "HP:0008736", "HP:0000028", "HP:0006101", "HP:0010747", "HP:0000135", "HP:0000494", "HP:0002230", "HP:0000639", "HP:0000003", "HP:0003202", "MKKS", "ORPHA:110", "Skeletal_muscle_atrophy"};
    EXPECT_EQ(expected_output, actual_output);

    phen.set_getAssHPOs("TSR2");
    _vvGenNamHPOsDiseaseID = phen.get_getAssHPOs();

    actual_output = _vvGenNamHPOsDiseaseID[0];
    expected_output = {"HP:0001631", "HP:0001199", "HP:0000470", "HP:0005280", "HP:0006758", "HP:0009778", "HP:0000294", "HP:0000465", "HP:0004808", "HP:0009944", "HP:0002669", "HP:0001882", "HP:0000286", "HP:0000980", "HP:0001629", "HP:0001894", "HP:0002863", "HP:0000369", "HP:0001254", "HP:0011904", "HP:0000347", "HP:0000218", "HP:0020118", "HP:0012133", "HP:0410030", "HP:0040276", "HP:0004322", "HP:0000486", "HP:0001873", "HP:0005518", "HP:0001087", "HP:0030270", "HP:0001227", "HP:0001518", "HP:0000252", "HP:0001875", "HP:0012410", "HP:0000431", "HP:0001896", "HP:0005532", "HP:0001680", "HP:0000912", "HP:0000316", "HP:0000047", "HP:0008551", "HP:0001790", "HP:0000085", "HP:0012758", "HP:0000508", "HP:0000519", "HP:0001895", "HP:0000104", "HP:0000185", "HP:0009777", "TSR2", "ORPHA:124", "Absent_thumb"};
    EXPECT_EQ(expected_output, actual_output);

    actual_output = _vvGenNamHPOsDiseaseID[1];
    expected_output = {"HP:0030270", "HP:0001419", "HP:0011800", "HP:0011904", "HP:0000413", "HP:0003577", "HP:0005518", "HP:0012741", "HP:0000653", "HP:0008551", "HP:0001972", "HP:0000494", "HP:0000405", "HP:0000347", "HP:0000175", "TSR2", "OMIM:300946", "Cleft_palate"};
    EXPECT_EQ(expected_output, actual_output);
}
