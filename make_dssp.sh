fpid=()
counter=0
cores=100
# for f in pdb_data/PDBs/*.cif.gz
# do
# 	f_out=${f##*/}
# 	f_out=pdb_data/dssp_temp/${f_out%.cif.gz}.dssp
# 	# echo $f_out
# 	rm -rf $f_out
# 	mkdssp -i $f -o $f_out &
# 	fpid+=($!)
# 	((counter=counter+1))
# 
# 	if (( $counter % $cores == 0 ));then
# 	  for f in ${fpid[*]}; do
#       # echo $f
#       wait $f
#     done
#   fi
# done
# for f in ${fpid[*]}; do
#       echo $f
#       wait $f
# done

fpid=()
counter=0
cores=100
for f in pdb_data/PDBs/*.cif.gz
do
	bash clean_dssp.sh $f &
	fpid+=($!)
	((counter=counter+1))

	if (( $counter % $cores == 0 ));then
	  for f in ${fpid[*]}; do
      # echo $f
      wait $f
    done
    fpid=()
  fi
done