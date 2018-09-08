#!/bin/bash
g++ -std=c++11 -o out main.cpp interface.cpp image.cpp atom.cpp bv.cpp bvv.cpp lj12.cpp ewald.cpp
chmod +x out
./out
