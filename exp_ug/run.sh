#/bin/sh
#chmod u+x gen.sh
generate(){
   for i in `seq 0 $2`
   do
      echo "$1"
   done
}

num=$1

for file in `ls ../data/*.linereg`
do
   generate $file $num
done

