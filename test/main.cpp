#include <cstdlib>
#include <gtest/gtest.h>
using namespace std;

int main(int argc, char **argv)
{
    char cwd[256];
    getcwd(cwd, 255);

    cout << "Executing tests from: " << cwd << endl;

    ::testing::InitGoogleTest(&argc, const_cast<char **>(argv));
    return RUN_ALL_TESTS();
}
