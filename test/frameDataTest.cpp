#include <gtest/gtest.h>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>
#include "../src/PhenotypeIntegration.h"
using namespace std;

class frameDataTest : public ::testing::Test
{
protected:
    void SetUp() override {}
    void TearDown() override {}
};

TEST_F(frameDataTest, TestframeDataTest)
{
    PhenotypeIntegration phen; 

    phen.init();
    phen.setupHPO("./HPO/HPO_release20220414/");
    ASSERT_EQ(phen.getInpContObo().size(), 983450);
    ASSERT_EQ(phen.getInpContGenPhen().size(), 259829);

    phen.init();
    phen.setupHPO("./HPO/HPO_release20230127/");
    ASSERT_EQ(phen.getInpContObo().size(), 1004349);
    ASSERT_EQ(phen.getInpContGenPhen().size(), 268983);

}

