#!/bin/sh

run()
{
   dir=$1
   echo "** dir = $1 **********************"

   ct=0
   sum_sens=0
   sum_spec=0
   for file in `ls -d ${dir}/*/*.pred`
   do
      sens=$(sed -n 10p $file)
      sum_sens=$(echo "$sum_sens + $sens" | bc -l)

      spec=$(sed -n 12p $file)
      sum_spec=$(echo "$sum_spec + $spec" | bc -l)

      ct=`expr $ct + 1`
   done

   ave_sens=$(echo "$sum_sens / $ct" | bc -l)
   ave_spec=$(echo "$sum_spec / $ct" | bc -l)

   echo "Sensitivity: $ave_sens"
   echo "Specificity: $ave_spec"
}

for dir in `ls -d ./spiked*`
do
   run $dir
done
