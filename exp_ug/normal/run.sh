#/bin/sh
#chmod u+x gen.sh
run_ug()
{
   file=$1
   buf=${1##*/}
   name=${buf%.*}
   sth=$2
   ../bin/fscip ../settings/exp10.set $file -q -sth $sth
   mv ${name}* exp10
   ../bin/fscip ../settings/exp20.set $file -q -sth $sth
   mv ${name}* exp20
   ../bin/fscip ../settings/exp30.set $file -q -sth $sth
   mv ${name}* exp30
}

datadir=$1
sth=$2

rm -rf exp10
rm -rf exp20
rm -rf exp30

mkdir exp10
mkdir exp20
mkdir exp30

for file in `ls ${datadir}/*.linereg`
do
   run_ug $file $sth
done

