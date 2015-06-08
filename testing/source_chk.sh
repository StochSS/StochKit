#!/bin/bash

cur_dir=$(pwd -P)

cd ..

sha512sum -w -c "$cur_dir/source.sha512" | grep -v "OK"
