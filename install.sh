#!/bin/bash
#perform installation steps

#set -x

for shell_argv in "$@"
do
  if [ "$shell_argv" == "--omit_dir_diff" ]; then
    omit_dir_diff="YES"
  fi
done

if [ "$STOCHKIT_HOME" == "" ]; then
  pushd "$(dirname "$0")" > /dev/null
  export STOCHKIT_HOME="`pwd -P`"
  popd > /dev/null
elif [ "$STOCHKIT_HOME" != "`pwd -P`" ]; then
  if [ "$omit_dir_diff" != "YES" ]; then
    echo "\$STOCHKIT_HOME=$STOCHKIT_HOME is not current directory, do you wish to continue? Y/N:"
    while read user_choice
    do
      case $user_choice in
        [yY]* ) break;;
        [nN]* ) exit;;
        * )     echo "Please enter Y or N:";;
      esac
    done
  fi
fi

case $STOCHKIT_HOME in
  *\ * )
    stochkittmpdir=$(mktemp -d /tmp/tmp.stochkit.XXXXX)
    trap 'rm -rf "$stochkittmpdir"' EXIT INT TERM HUP
    ln -s "$STOCHKIT_HOME" $stochkittmpdir/stochkit_home
    export STOCHKIT_ORGININAL_HOME="$STOCHKIT_HOME"
    export STOCHKIT_HOME=$stochkittmpdir/stochkit_home;;
esac

cd "$STOCHKIT_HOME/libs/boost_1_53_0/"

./bootstrap.sh

if [ "$?" -ne 0 ]; then
   echo "bootstrap.sh failed to run" 
   exit 1
fi

./b2 --with-thread

if [ "$?" -ne 0 ]; then
   echo "b2 failed to run" 
   exit 1
fi

./b2 --with-program_options

if [ "$?" -ne 0 ]; then
   echo "b2 failed to run" 
   exit 1
fi

./b2 --with-filesystem

if [ "$?" -ne 0 ]; then
   echo "b2 failed to run" 
   exit 1
fi

cd "$STOCHKIT_HOME"

if [ "$#" -ge "1" -a "$1" == "MPI" ]
then
   make -j 4 COMPILE_MPI="1" EXTRA_OPTIONS="-DCOMPILE_MPI"
else
   make -j 4
fi

if [ "$?" -eq 2 ]; then
   echo "make failed" 
   exit 1
fi

cd "$STOCHKIT_HOME/custom_drivers/single_trajectory"
if [ "$#" -ge "1" -a "$1" == "MPI" ]
then
   make -j 4 COMPILE_MPI="1" EXTRA_OPTIONS="-DCOMPILE_MPI"
else
   make -j 4
fi

if [ "$?" -eq 2 ]; then
   echo "make failed" 
   exit 1
fi

cd "$STOCHKIT_HOME/custom_drivers/custom_user_output"
if [ "$#" -ge "1" -a "$1" == "MPI" ]
then
   make -j 4 COMPILE_MPI="1" EXTRA_OPTIONS="-DCOMPILE_MPI"
else
   make -j 4
fi

if [ "$?" -eq 2 ]; then
   echo "make failed" 
   exit 1
fi

#echo "**********************************************"
#echo "Please remember to add $STOCHKIT_HOME/libs/boost_1_53_0/stage/lib/ to the LD_LIBRARY_PATH environment variable on Linux machines, or the DYLD_LIBRARY_PATH environment variable on Mac machines."
echo "You may add the home directory of StochKit to the PATH environment variable to run ssa or tau_leaping from any directories."
#echo "**********************************************"

