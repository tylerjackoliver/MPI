for i in 1 2 4 8 16 24
do
    echo "${i}" >> ring_data/data.txt
    aprun -n $i ./ring >> ring_data/data.txt
    echo "Done ${i}"

done
