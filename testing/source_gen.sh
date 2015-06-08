#!/bin/bash

cur_dir=$(pwd -P)

cd ..

find src \( -iname "*.cpp" -o -iname "*.h" -o -iname "*.ipp" \) -exec bash -c "sha512sum '{}' >> $cur_dir/source.sha512" ';'
find src -maxdepth 1 -iname "Makefile" -exec bash -c "sha512sum '{}' >> $cur_dir/source.sha512" ';'
find .make \( -iname "src4obj.h" -o -iname "dependency.h" \) -exec bash -c "sha512sum '{}' >> $cur_dir/source.sha512" ';'
chmod 600 "$cur_dir/source.sha512"
