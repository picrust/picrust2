sh ro_submodules.sh 
git submodule init
git submodule update
mkdir build
cd build/
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j2
cp papara ..
cd ..

