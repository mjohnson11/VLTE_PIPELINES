MT=1
for ((i=$((${MT}*10));i<$((${MT}*10+10));i++));
do
  SAMP=$(sed -n ${i}'{p;q}' ../accessory_files/Samples.txt)
  echo $i
  echo $SAMP
done
