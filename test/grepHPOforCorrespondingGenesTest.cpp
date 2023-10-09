#include <gtest/gtest.h>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>
#include "../src/PhenotypeIntegration.h"
using namespace std;

// Define a test fixture class for the Read class
class findHPOassGenTest : public ::testing::Test
{
protected:
    void SetUp() override {}
    void TearDown() override {}
};

TEST_F(findHPOassGenTest, TestfindHPOassGenTest)
{

    int actual_output;
    int expected_output;
    PhenotypeIntegration phen;

    phen.init();
    phen.setupHPO("./HPO/HPO_release20220414/");

    phen.set_findHPOassGen("HP:0000041");
    actual_output = phen.get_findHPOassGen();
    expected_output = 16;
    EXPECT_EQ(expected_output, actual_output);

    phen.set_findHPOassGen("HP:0000042");
    actual_output = phen.get_findHPOassGen();
    expected_output = 3;
    EXPECT_EQ(expected_output, actual_output);

    // there was the term in the tree, but no associated to genes, while in later releases there was a gene association
    phen.set_findHPOassGen("HP:0000811"); // MICA of the 2 above
    actual_output = phen.get_findHPOassGen();
    expected_output = 1039;
    EXPECT_EQ(expected_output, actual_output);

    phen.set_findHPOassGen("HP:0002450");
    actual_output = phen.get_findHPOassGen();
    expected_output = 65;
    EXPECT_EQ(expected_output, actual_output);

    phen.set_findHPOassGen("HP:0011147");
    actual_output = phen.get_findHPOassGen();
    expected_output = 16;
    EXPECT_EQ(expected_output, actual_output);

    phen.set_findHPOassGen("HP:0000707"); // MICA of the 2 above
    actual_output = phen.get_findHPOassGen();
    expected_output = 5409;
    EXPECT_EQ(expected_output, actual_output);



    phen.init();
    phen.setupHPO("./HPO/HPO_release20230127/");

    phen.set_findHPOassGen("HP:0000041");
    actual_output = phen.get_findHPOassGen();
    expected_output = 17;
    EXPECT_EQ(expected_output, actual_output);

    phen.set_findHPOassGen("HP:0000042");
    actual_output = phen.get_findHPOassGen();
    expected_output = 3;
    EXPECT_EQ(expected_output, actual_output);

    phen.set_findHPOassGen("HP:0000811");
    actual_output = phen.get_findHPOassGen();
    expected_output = 1096;
    EXPECT_EQ(expected_output, actual_output);

    phen.set_findHPOassGen("HP:0002450");
    actual_output = phen.get_findHPOassGen();
    expected_output = 65;
    EXPECT_EQ(expected_output, actual_output);

    phen.set_findHPOassGen("HP:0011147");
    actual_output = phen.get_findHPOassGen();
    expected_output = 19;
    EXPECT_EQ(expected_output, actual_output);

    phen.set_findHPOassGen("HP:0000707"); // MICA of the 2 above
    actual_output = phen.get_findHPOassGen();
    expected_output = 5614;
    EXPECT_EQ(expected_output, actual_output);
}
