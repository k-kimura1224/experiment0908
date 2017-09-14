#!/bin/sh

edit()
{
   solfile=$1
   buf=${1##*/}
   datadir=$2
   dirname=$3
   subdirname=$4

   cp $solfile cp.sol
   numline=$(cat $solfile | wc -l)
   for i in `grep -e objective -n cp.sol | sed -e 's/:.*//g'`
   do
      x=`expr $numline - $i + 1`
      tail -n $x $solfile > buf.sol
   done
   rm -f cp.sol
   mv buf.sol $buf
   mv $buf ${datadir}/${dirname}/${subdirname}
}

subrun()
{
   soldir=$1
   datadir=$2
   dirname=$3

   subdir=$4
   subdirname=${4##*/}

   for solfile in `ls -d $subdir/*.sol`
   do
      edit $solfile $datadir $dirname $subdirname
   done
}

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
