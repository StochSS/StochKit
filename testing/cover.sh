#!/bin/bash

pushd ".." > /dev/null
par_dir=$(pwd -P)
popd > /dev/null

lcov --checksum --directory "$par_dir/obj" --directory "$par_dir/src" --base-directory "$par_dir/src" --no-external -c -o coverage_temp.lcov

genhtml --legend --demangle-cpp -o html/ coverage_temp.lcov
