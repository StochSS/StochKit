#!/bin/bash

cur_dir=$(pwd -P)

cd ..

sha512sum -w -c "$cur_dir/result.sha512" | grep -v "OK"
sha512sum -w -c "$cur_dir/result2.sha512" | grep -v "OK"
