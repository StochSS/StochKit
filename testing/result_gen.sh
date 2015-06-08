#!/bin/bash

cur_dir=$(pwd -P)

cd ..

find models/examples/ \( -ipath "*output_case2*" -o -ipath "*output_case3*" -o -ipath "*output_case5*" -o -ipath "*output_case6*" \) \( -iname "*variance*txt" -o -iname "*mean*txt" -o -iname "*hist*dat" -o -iname "*trajectory*txt" \) -exec bash -c "sha512sum '{}' >> $cur_dir/result.sha512" ';'
find models/examples/ \( -ipath "*output_case2*" -o -ipath "*output_case3*" -o -ipath "*output_case5*" -o -ipath "*output_case6*" \) \( \! -ipath "*schlogl*" \! -ipath "*events*" \! -ipath "*mixed*" \) \( -iname "*variance*txt" -o -iname "*mean*txt" -o -iname "*hist*dat" -o -iname "*trajectory*txt" \) -exec bash -c "sha512sum '{}' >> $cur_dir/result2.sha512" ';'
chmod 600 "$cur_dir/result.sha512" "$cur_dir/result2.sha512"
