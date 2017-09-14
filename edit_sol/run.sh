#!/bin/sh

run()
{
   soldir=$1
   datadir=$2

   dir=$3
   dirname=${3##*/}

   for subdir in `ls -d ${dir}/*`
   do
      subrun $soldir $datadir $dirname $subdir
   done
}

soldir=$1
datadir=$2

for dir in `ls -d ${soldir}/spiked*`
do
   run $soldir $datadir $dir
done
