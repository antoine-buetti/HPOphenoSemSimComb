#include <gtest/gtest.h>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>
#include "../src/PhenotypeIntegration.h"
using namespace std;

// Define a test fixture class for the Read class
class findMICATest : public ::testing::Test
{
protected:
    void SetUp() override {}
    void TearDown() override {}
};

TEST_F(findMICATest, TestfindMICATest)
{
    string actual_output;
    string expected_output;

    PhenotypeIntegration phen;

    phen.init();
    // phen.setupHPO("./data/hp_small.obo");
    phen.setupHPO("./HPO/HPO_release20220414/");

    // Different cases where a pair of HPO terms is injected into the function, and with the test HPO dataset ("hp_release20220414.obo") the expected MICA is checked using as ground truth the verification with Cytoscape UFO app and literature cases.

    // Simple cases where the path to the root follows HPO terms with descending indexes

    // case 0a (perfect match user-disease term)
    actual_output = phen.get_findMICA("HP:0000100", "HP:0000100");
    expected_output = "HP:0000100";
    EXPECT_EQ(expected_output, actual_output);

    // case 0b (perfect match user-disease term)
    actual_output = phen.get_findMICA("HP:0000200", "HP:0000200");
    expected_output = "HP:0000200";
    EXPECT_EQ(expected_output, actual_output);

    // case 0c (perfect match user-disease term)
    actual_output = phen.get_findMICA("HP:0000001", "HP:0000001");
    expected_output = "HP:0000001";
    EXPECT_EQ(expected_output, actual_output);

    // case 1
    actual_output = phen.get_findMICA("HP:0000100", "HP:0000200");
    expected_output = "HP:0000118";
    EXPECT_EQ(expected_output, actual_output);

    // case 2
    actual_output = phen.get_findMICA("HP:0000300", "HP:0000200");
    expected_output = "HP:0000271";
    EXPECT_EQ(expected_output, actual_output);

    // case 3
    actual_output = phen.get_findMICA("HP:0000300", "HP:0000400");
    expected_output = "HP:0000118";
    EXPECT_EQ(expected_output, actual_output);

    // case 4
    actual_output = phen.get_findMICA("HP:0000501", "HP:0000400");
    expected_output = "HP:0000118";
    EXPECT_EQ(expected_output, actual_output);

    // case 5
    actual_output = phen.get_findMICA("HP:0002450", "HP:0020219");
    expected_output = "HP:0000707";
    EXPECT_EQ(expected_output, actual_output);

    // case 5b (when one term is a direct ancestor)
    actual_output = phen.get_findMICA("HP:0002450", "HP:0012639");
    expected_output = "HP:0012639";
    EXPECT_EQ(expected_output, actual_output);

    // case 5c (when one term is a direct ancestor)
    actual_output = phen.get_findMICA("HP:0032794", "HP:0001250");
    expected_output = "HP:0001250";
    EXPECT_EQ(expected_output, actual_output);

    // case 6
    actual_output = phen.get_findMICA("HP:0002450", "HP:0020219");
    expected_output = "HP:0000707";
    EXPECT_EQ(expected_output, actual_output);

    // // Following cases are tricky because the path to the root follows HPO terms with non-descending indexes

    // case 7
    actual_output = phen.get_findMICA("HP:0020219", "HP:0002197");
    expected_output = "HP:0001250";
    EXPECT_EQ(expected_output, actual_output);

    // case 8
    actual_output = phen.get_findMICA("HP:0020219", "HP:0033259");
    expected_output = "HP:0001250";
    EXPECT_EQ(expected_output, actual_output);

    // case 9
    actual_output = phen.get_findMICA("HP:0032677", "HP:0002121");
    expected_output = "HP:0002197";
    EXPECT_EQ(expected_output, actual_output);

    // // Following below are examples from the scheme in Fig.2 of
    // // Carmody et al., Orphanet Journal of Rare Diseases (2020) https://doi.org/10.1186/s13023-020-1313-0
    // // (yes the scheme match with the HPO release used for testing)

    // case 10
    actual_output = phen.get_findMICA("HP:0009099", "HP:0000161");
    expected_output = "HP:0000202";
    EXPECT_EQ(expected_output, actual_output);

    // case 11 (when one term is a direct ancestor)
    actual_output = phen.get_findMICA("HP:0009099", "HP:0000175");
    expected_output = "HP:0000175";
    EXPECT_EQ(expected_output, actual_output);

    // case 12 (when one term is a direct ancestor)
    actual_output = phen.get_findMICA("HP:0009099", "HP:0000163");
    expected_output = "HP:0000163";
    EXPECT_EQ(expected_output, actual_output);

    // case 13
    actual_output = phen.get_findMICA("HP:0009099", "HP:0000159");
    expected_output = "HP:0000163";
    EXPECT_EQ(expected_output, actual_output);

    // case 14
    actual_output = phen.get_findMICA("HP:0009099", "HP:0000177");
    expected_output = "HP:0000163";
    EXPECT_EQ(expected_output, actual_output);

    // case 15
    actual_output = phen.get_findMICA("HP:0000175", "HP:0000177");
    expected_output = "HP:0000163";
    EXPECT_EQ(expected_output, actual_output);

    // case 16: tricky case with confusing loop
    actual_output = phen.get_findMICA("HP:0000041", "HP:0000042");
    expected_output = "HP:0000811"; // and not HP:0012243
    EXPECT_EQ(expected_output, actual_output);

    // PhenotypeIntegration phen;

    phen.init();
    phen.setupHPO("./HPO/HPO_release20230127/");

    // // Same as above but with another release

    // case 0a (perfect match user-disease term)
    actual_output = phen.get_findMICA("HP:0000100", "HP:0000100");
    expected_output = "HP:0000100";
    EXPECT_EQ(expected_output, actual_output);

    // case 0b (perfect match user-disease term)
    actual_output = phen.get_findMICA("HP:0000200", "HP:0000200");
    expected_output = "HP:0000200";
    EXPECT_EQ(expected_output, actual_output);

    // case 0c (perfect match user-disease term)
    actual_output = phen.get_findMICA("HP:0000001", "HP:0000001");
    expected_output = "HP:0000001";
    EXPECT_EQ(expected_output, actual_output);

    // case 1
    actual_output = phen.get_findMICA("HP:0000100", "HP:0000200");
    expected_output = "HP:0000118";
    EXPECT_EQ(expected_output, actual_output);

    // case 2
    actual_output = phen.get_findMICA("HP:0000300", "HP:0000200");
    expected_output = "HP:0000271";
    EXPECT_EQ(expected_output, actual_output);

    // case 3
    actual_output = phen.get_findMICA("HP:0000300", "HP:0000400");
    expected_output = "HP:0000118";
    EXPECT_EQ(expected_output, actual_output);

    // case 4
    actual_output = phen.get_findMICA("HP:0000501", "HP:0000400");
    expected_output = "HP:0000118";
    EXPECT_EQ(expected_output, actual_output);

    // case 5
    actual_output = phen.get_findMICA("HP:0002450", "HP:0020219");
    expected_output = "HP:0000707";
    EXPECT_EQ(expected_output, actual_output);

    // case 5b (when one term is a direct ancestor)
    actual_output = phen.get_findMICA("HP:0002450", "HP:0012639");
    expected_output = "HP:0012639";
    EXPECT_EQ(expected_output, actual_output);

    // case 5c (when one term is a direct ancestor)
    actual_output = phen.get_findMICA("HP:0032794", "HP:0001250");
    expected_output = "HP:0001250";
    EXPECT_EQ(expected_output, actual_output);

    // case 6
    actual_output = phen.get_findMICA("HP:0002450", "HP:0020219");
    expected_output = "HP:0000707";
    EXPECT_EQ(expected_output, actual_output);

    // // Following cases are tricky because the path to the root follows HPO terms with non-descending indexes

    // case 7
    actual_output = phen.get_findMICA("HP:0020219", "HP:0002197");
    expected_output = "HP:0001250";
    EXPECT_EQ(expected_output, actual_output);

    // case 8
    actual_output = phen.get_findMICA("HP:0020219", "HP:0033259");
    expected_output = "HP:0001250";
    EXPECT_EQ(expected_output, actual_output);

    // case 9
    actual_output = phen.get_findMICA("HP:0032677", "HP:0002121");
    expected_output = "HP:0002197";
    EXPECT_EQ(expected_output, actual_output);

    // // Following below are examples from the scheme in Fig.2 of
    // // Carmody et al., Orphanet Journal of Rare Diseases (2020) https://doi.org/10.1186/s13023-020-1313-0
    // // (yes the scheme match with the HPO release used for testing)

    // case 10: tricky casa that with this one release was misidentified with HP:0000163 due to an additional link
    actual_output = phen.get_findMICA("HP:0009099", "HP:0000161");
    expected_output = "HP:0000202";
    EXPECT_EQ(expected_output, actual_output);

    // case 11 (when one term is a direct ancestor)
    actual_output = phen.get_findMICA("HP:0009099", "HP:0000175");
    expected_output = "HP:0000175";
    EXPECT_EQ(expected_output, actual_output);

    // case 12 (when one term is a direct ancestor)
    actual_output = phen.get_findMICA("HP:0009099", "HP:0000163");
    expected_output = "HP:0000163";
    EXPECT_EQ(expected_output, actual_output);

    // case 13
    actual_output = phen.get_findMICA("HP:0009099", "HP:0000159");
    expected_output = "HP:0000163";
    EXPECT_EQ(expected_output, actual_output);

    // case 14
    actual_output = phen.get_findMICA("HP:0009099", "HP:0000177");
    expected_output = "HP:0000163";
    EXPECT_EQ(expected_output, actual_output);

    // case 15
    actual_output = phen.get_findMICA("HP:0000175", "HP:0000177");
    expected_output = "HP:0000163";
    EXPECT_EQ(expected_output, actual_output);

    // case 16: tricky case with confusing loop
    actual_output = phen.get_findMICA("HP:0000041", "HP:0000042");
    expected_output = "HP:0000811"; // and not HP:0012243
    EXPECT_EQ(expected_output, actual_output);
}
