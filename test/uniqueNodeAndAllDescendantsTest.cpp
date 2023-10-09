#include <gtest/gtest.h>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>
#include "../src/PhenotypeIntegration.h"
using namespace std;

// Define a test fixture class for the Read class
class uniqNodAllDescTest : public ::testing::Test
{
protected:
    void SetUp() override {}
    void TearDown() override {}
    PhenotypeIntegration phen;
};

TEST_F(uniqNodAllDescTest, TestuniqNodAllDescTest)
{
    phen.init();
    phen.setupHPO("./HPO/HPO_release20220414/");

    // these number of descendants can be validated also with the jax HPO browser (https://hpo.jax.org/app/)
    // (they report a total number of HPO terms inferior by 1 because they do not include the search term itself)


    // check the sub-ontology sizes

    //  Mode of inheritance 
    phen.set_treeWalkDown("HP:0000005");
    phen.set_uniqNodAllDesc();
    ASSERT_EQ(phen.get_uniqNodAllDesc_numberOf(), 32);

    // Phenotypic abnormality 
    phen.set_treeWalkDown("HP:0000118");
    phen.set_uniqNodAllDesc();
    ASSERT_EQ(phen.get_uniqNodAllDesc_numberOf(),  16321);

    // Clinical modifier 
    phen.set_treeWalkDown("HP:0012823"); 
    phen.set_uniqNodAllDesc();
    ASSERT_EQ(phen.get_uniqNodAllDesc_numberOf(), 192);

    // Blood group 
    phen.set_treeWalkDown("HP:0032223"); 
    phen.set_uniqNodAllDesc();
    ASSERT_EQ(phen.get_uniqNodAllDesc_numberOf(), 9);

    // Past medical history
    phen.set_treeWalkDown("HP:0032443"); 
    phen.set_uniqNodAllDesc();
    ASSERT_EQ(phen.get_uniqNodAllDesc_numberOf(), 19);

    // Frequency
    phen.set_treeWalkDown("HP:0040279"); 
    phen.set_uniqNodAllDesc();
    ASSERT_EQ(phen.get_uniqNodAllDesc_numberOf(), 7);



    phen.init();
    phen.setupHPO("./HPO/HPO_release20230127/");

    //  Mode of inheritance 
    phen.set_treeWalkDown("HP:0000005");
    phen.set_uniqNodAllDesc();
    ASSERT_EQ(phen.get_uniqNodAllDesc_numberOf(), 41);

    // Phenotypic abnormality 
    phen.set_treeWalkDown("HP:0000118");
    phen.set_uniqNodAllDesc();
    ASSERT_EQ(phen.get_uniqNodAllDesc_numberOf(),  16562);

    // Clinical modifier 
    phen.set_treeWalkDown("HP:0012823"); 
    phen.set_uniqNodAllDesc();
    ASSERT_EQ(phen.get_uniqNodAllDesc_numberOf(), 210);

    // Blood group 
    phen.set_treeWalkDown("HP:0032223"); 
    phen.set_uniqNodAllDesc();
    ASSERT_EQ(phen.get_uniqNodAllDesc_numberOf(), 9);

    // Past medical history
    phen.set_treeWalkDown("HP:0032443"); 
    phen.set_uniqNodAllDesc();
    ASSERT_EQ(phen.get_uniqNodAllDesc_numberOf(), 44);

    // Frequency
    phen.set_treeWalkDown("HP:0040279"); 
    phen.set_uniqNodAllDesc();
    ASSERT_EQ(phen.get_uniqNodAllDesc_numberOf(), 7);




    // Other cases

    phen.init();
    phen.setupHPO("./HPO/HPO_release20220414/");

    // Inheritance modifier (only present in release20230127)
    phen.set_treeWalkDown("HP:0034335"); 
    phen.set_uniqNodAllDesc();
    ASSERT_EQ(phen.get_uniqNodAllDesc_numberOf(), 1); // 1 because the query is added to the set, but during runtime would not happen as the queries are checked for presence in the HPO data very initially.

    // A term reported twice by treeWalkDown function (9->8)
    phen.set_treeWalkDown("HP:0002450");
    phen.set_uniqNodAllDesc();
    ASSERT_EQ(phen.get_uniqNodAllDesc_numberOf(), 8);

    // Also with a term reported twice by treeWalkDown
    phen.set_treeWalkDown("HP:0002197");
    phen.set_uniqNodAllDesc();
    ASSERT_EQ(phen.get_uniqNodAllDesc_numberOf(), 21);



    phen.init();
    phen.setupHPO("./HPO/HPO_release20230127/");

    // Inheritance modifier (only present in release20230127)
    phen.set_treeWalkDown("HP:0034335"); 
    phen.set_uniqNodAllDesc();
    ASSERT_EQ(phen.get_uniqNodAllDesc_numberOf(), 22);



}
