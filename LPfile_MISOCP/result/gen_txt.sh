#/bin/sh
#chmod u+x gen.sh
generate(){
   file=$1
   buf=${1##*/}
   name=${buf%.*}
   txtname=$2
   echo "read $file"
   echo "read $file" >> $txtname
   echo "opt"
   echo "opt" >> $txtname
   echo "write ${name}.sol"
   echo "write ${name}.sol" >> $txtname
}

time=5000
txtname=run.txt

rm -f $txtname
echo "generate $txtname"
echo "set timelimit $time"
echo "set timelimit $time" >> $txtname

for file in `ls ../output/*.lp`
do
   generate $file $txtname
done

echo "q"
echo "q" >> $txtname
