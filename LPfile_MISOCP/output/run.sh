#/bin/sh
#chmod u+x gen.sh
generate(){
   file=$1
   buf=${1##*/}
   name=${buf%.*}
   k=$2
   ../bin/genelpfile $file $k > ${name}_${k}.lp
}

k=30

for file in `ls ../data/*.linereg`
do
   generate $file $k
   echo "generate lp-file from ${file}"
done

k=20

for file in `ls ../data/*.linereg`
do
   generate $file $k
   echo "generate lp-file from ${file}"
done

k=10

for file in `ls ../data/*.linereg`
do
   generate $file $k
   echo "generate lp-file from ${file}"
done
