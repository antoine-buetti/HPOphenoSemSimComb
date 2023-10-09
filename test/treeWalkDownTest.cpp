#include <gtest/gtest.h>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>
#include "../src/PhenotypeIntegration.h"
using namespace std;

class treeWalkDownTest : public ::testing::Test
{
protected:
    void SetUp() override {}
    void TearDown() override {}
};

TEST_F(treeWalkDownTest, TestTreeWalkDownTest)
{

    vector<string> actual_output;
    vector<string> expected_output;

    PhenotypeIntegration phen;

    phen.init();
    // phen.setupHPO("./data/hp_small.obo");
    phen.setupHPO("./HPO/HPO_release20220414/");

    // // case 1: HPO root node for "Mendelian inheritance" branch
    // phen.set_treeWalkDown("HP:0034345");
    // actual_output = phen.get_treeWalkDown();
    // expected_output = {
    //     "HP:0034345",
    //     "HP:0000006",
    //     "HP:0000007",
    //     "HP:0001417",
    //     "HP:0001419",
    //     "HP:0001423",
    //     "HP:0001427",
    //     "HP:0001450",
    //     "HP:0032113",
    //     "HP:0034339",
    //     "HP:0034340",
    //     "HP:0034341"};
    // EXPECT_EQ(expected_output, actual_output);
    // // -> this whole branch of "Mendelian inheritance" was only added later than 2022, see below tests with newer release

    // case 2: HPO root node for "Mortality/Aging" branch
    phen.set_treeWalkDown("HP:0040006");
    actual_output = phen.get_treeWalkDown();
    expected_output = {
        "HP:0040006",
        "HP:0011420",
        "HP:0001522",
        "HP:0003811",
        "HP:0003819",
        "HP:0011421",
        "HP:0033763",
        "HP:0033764",
        "HP:0033765",
        "HP:0100613",
        "HP:0034241",
        "HP:0003826",
        "HP:0005268"
        };
    EXPECT_EQ(expected_output, actual_output);
    // case 3: with branching point -> warning that this make the branching convergence node appear twice
    phen.set_treeWalkDown("HP:0002450");
    actual_output = phen.get_treeWalkDown();
    expected_output = {
        "HP:0002450",
        "HP:0002127",
        "HP:0002366",
        "HP:0006802",
        "HP:0002398",
        "HP:0007277",
        "HP:0007373",
        "HP:0002398",
        "HP:0007354"
        };
    EXPECT_EQ(expected_output, actual_output);

    // case 4: another branching point (with more complex branching example in TER-760 slide deck)
    phen.set_treeWalkDown("HP:0002197");
    actual_output = phen.get_treeWalkDown();
    expected_output = {
        "HP:0002197",
        "HP:0002121",
        "HP:0007270",
        "HP:0011147",
        "HP:0011152",
        "HP:0011149",
        "HP:0011150",
        "HP:0032677",
        "HP:0002123",
        "HP:0001327",
        "HP:0032678",
        "HP:0011149",
        "HP:0010818",
        "HP:0011169",
        "HP:0011170",
        "HP:0025190",
        "HP:0007207",
        "HP:0032661",
        "HP:0032795",
        "HP:0032855",
        "HP:0032842",
        "HP:0032887"
        };
    EXPECT_EQ(expected_output, actual_output);


    // case 5: select an end leaf
    phen.set_treeWalkDown("HP:0032887");
    actual_output = phen.get_treeWalkDown();
    expected_output = {
        "HP:0032887"
        };
    EXPECT_EQ(expected_output, actual_output);

    // case 6: select an other end leaf
    phen.set_treeWalkDown("HP:0007354");
    actual_output = phen.get_treeWalkDown();
    expected_output = {
        "HP:0007354"
        };
    EXPECT_EQ(expected_output, actual_output);



// // Repeat tests with another release of the HPO data



    phen.init();
    phen.setupHPO("./HPO/HPO_release20230127/");

    // case 1: HPO root node for "Mendelian inheritance" branch
    phen.set_treeWalkDown("HP:0034345");
    actual_output = phen.get_treeWalkDown();
    expected_output = {
        "HP:0034345",
        "HP:0000006",
        "HP:0000007",
        "HP:0001417",
        "HP:0001419",
        "HP:0001423",
        "HP:0001427",
        "HP:0001450",
        "HP:0032113",
        "HP:0034339",
        "HP:0034340",
        "HP:0034341"
        };
    EXPECT_EQ(expected_output, actual_output);

    // case 2: HPO root node for "Mortality/Aging" branch
    phen.set_treeWalkDown("HP:0040006");
    actual_output = phen.get_treeWalkDown();
    expected_output = {
        "HP:0040006",
        "HP:0011420",
        "HP:0001522",
        "HP:0003811",
        "HP:0003819",
        "HP:0011421",
        "HP:0033763",
        "HP:0033764",
        "HP:0033765",
        "HP:0100613",
        "HP:0034241",
        "HP:0003826",
        "HP:0005268"
        };
    EXPECT_EQ(expected_output, actual_output);
    // case 3: with branching point -> warning that this make the branching convergence node appear twice
    phen.set_treeWalkDown("HP:0002450");
    actual_output = phen.get_treeWalkDown();
    expected_output = {
        "HP:0002450",
        "HP:0002127",
        "HP:0002366",
        "HP:0006802",
        "HP:0002398",
        "HP:0007277",
        "HP:0007373",
        "HP:0002398",
        "HP:0007354"
        };
    EXPECT_EQ(expected_output, actual_output);

    // case 4: another branching point (with more complex branching example in TER-760 slide deck)
    phen.set_treeWalkDown("HP:0002197");
    actual_output = phen.get_treeWalkDown();
    expected_output = {
        "HP:0002197",
        "HP:0002121",
        "HP:0007270",
        "HP:0011147",
        "HP:0011152",
        "HP:0011149",
        "HP:0011150",
        "HP:0032677",
        "HP:0002123",
        "HP:0001327",
        "HP:0032678",
        "HP:0011149",
        "HP:0010818",
        "HP:0011169",
        "HP:0011170",
        "HP:0025190",
        "HP:0007207",
        "HP:0032661",
        "HP:0032795",
        "HP:0032855",
        "HP:0032842",
        "HP:0032887"
        };
    EXPECT_EQ(expected_output, actual_output);


    // case 5: select an end leaf
    phen.set_treeWalkDown("HP:0032887");
    actual_output = phen.get_treeWalkDown();
    expected_output = {
        "HP:0032887"
        };
    EXPECT_EQ(expected_output, actual_output);

    // case 6: select an other end leaf
    phen.set_treeWalkDown("HP:0007354");
    actual_output = phen.get_treeWalkDown();
    expected_output = {
        "HP:0007354"
        };
    EXPECT_EQ(expected_output, actual_output);



}
