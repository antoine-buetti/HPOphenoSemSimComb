#include <gtest/gtest.h>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>
#include "../src/PhenotypeIntegration.h"
using namespace std;

class treeWalkUpTest : public ::testing::Test
{
protected:
    void SetUp() override {}
    void TearDown() override {}
};

TEST_F(treeWalkUpTest, TestTreeWalkUpTest)
{
    vector<string> actual_output;
    vector<string> expected_output;


    PhenotypeIntegration phen;

    phen.init();
    // phen.setupHPO("./data/hp_small.obo");
    phen.setupHPO("./HPO/HPO_release20220414/");


    // case 1: single up path
    phen.set_treeWalkUp("HP:0020219");
    actual_output = phen.get_treeWalkUp();
    expected_output = {
        "HP:0020219",
        "HP:0001250",
        "HP:0012638",
        "HP:0000707",
        "HP:0000118",
        "HP:0000001"};
    EXPECT_EQ(expected_output, actual_output);

    // case 2: same single up path as case 1
    phen.set_treeWalkUp("HP:0002197");
    actual_output = phen.get_treeWalkUp();
    expected_output = {
        "HP:0002197",
        "HP:0001250",
        "HP:0012638",
        "HP:0000707",
        "HP:0000118",
        "HP:0000001"};
    EXPECT_EQ(expected_output, actual_output);

    // case 3: common branching point just one level below case 1 and 2
    phen.set_treeWalkUp("HP:0032677");
    actual_output = phen.get_treeWalkUp();
    expected_output = {
        "HP:0032677",
        "HP:0002197",
        "HP:0001250",
        "HP:0012638",
        "HP:0000707",
        "HP:0000118",
        "HP:0000001",
        "HP:0020219",
        "HP:0001250",
        "HP:0012638",
        "HP:0000707",
        "HP:0000118",
        "HP:0000001"};
    EXPECT_EQ(expected_output, actual_output);

    // case 4: another branching point one more level below
    phen.set_treeWalkUp("HP:0002123");
    actual_output = phen.get_treeWalkUp();
    expected_output = {
        "HP:0002123",
        "HP:0032677",
        "HP:0002197",
        "HP:0001250",
        "HP:0012638",
        "HP:0000707",
        "HP:0000118",
        "HP:0000001",
        "HP:0020219",
        "HP:0001250",
        "HP:0012638",
        "HP:0000707",
        "HP:0000118",
        "HP:0000001",
        "HP:0032794",
        "HP:0020219",
        "HP:0001250",
        "HP:0012638",
        "HP:0000707",
        "HP:0000118",
        "HP:0000001"};
    EXPECT_EQ(expected_output, actual_output);
}
