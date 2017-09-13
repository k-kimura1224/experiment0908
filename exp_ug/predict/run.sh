#!/bin/sh

ug()
{
   file=$1
   sth=$2
   dirname=$3
   subdirname=$4
   buf=${1##*/}
   filename=${buf%.*}

   ../bin/fscip ../settings/exp10.set $file -q -sth $sth
   #../bin/fscip ../settings/exp20.set $file -q -sth $sth
   #../bin/fscip ../settings/exp30.set $file -q -sth $sth
   mv ${filename}* ${dirname}/${subdirname}
}

subrun()
{
   subdir=$1
   sth=$2
   dirname=$3
   subdirname=${1##*/}

   mkdir ${dirname}/${subdirname}

   for file in `ls ${subdir}/*_sample*.linereg`
   do
      ug $file $sth $dirname $subdirname
   done
}

run()
{
   dir=$1
   sth=$2
   dirname=${1##*/}

   rm -rf ${dirname}
   mkdir ${dirname}

   for subdir in `ls -d ${dir}/*`
   do
      subrun ${subdir} ${sth} ${dirname}
   done
}

datadir=$1
sth=$2

for dir in `ls -d ${datadir}/spiked*`
do
   run $dir $sth
done

