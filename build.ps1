mkdir build
cd build
cmake .. -G "MinGW Makefiles" -DCMAKE_MAKE_PROGRAM=mingw32-make
mingw32-make