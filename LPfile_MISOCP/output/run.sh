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

#for file in `ls ../data/*.data`
for file in `ls ../../../data_scip/*.linereg`
do
   generate $file $k
   echo "generate lp-file from ${file}"
done


#generate housing
#generate servo
#generate auto-mpg
#generate automobile
#generate breastcancer
#generate crime
#generate forestfires
#generate solarflareC
#generate solarflareX
#generate solarflareM

