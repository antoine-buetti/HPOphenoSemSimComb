sudo apt-get install cmake

mkdir lib
cd lib
git rm --cached ./googletest
rm -rf ./googletest
git submodule add https://github.com/google/googletest.git
git submodule update --init --recursive

cd ..
