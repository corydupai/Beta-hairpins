fpid=()
counter=0
cores=1000
f=$1
f_out=${f##*/}
f_out2=pdb_data/dssp_temp/${f_out%.cif.gz}.dssp
f_out=pdb_data/dssp/${f_out%.cif.gz}.dssp
# echo $f
cp $f_out2 $f_out
sed -i '1,27d' $f_out
sed -i 's/.*\!.*//' $f_out
sed -i -r 's/(.{15}).{19}/\1/' $f_out
sed -i  -r 's/-([[:digit:]])/ -\1/g' $f_out
# sed -i -r 's/(.{19}).{19}/\1/g' $f_out #115 19

sed -i 's/RESIDUE/Position  Chain1/' $f_out
sed -i 's/AUTHCHAIN/Chain/' $f_out
sed -i 's/ AA / aa /' $f_out
sed -i 's/ BP1/ bp1/' $f_out
sed -i 's/ BP2/ bp2/' $f_out
sed -i 's/ ACC/ acc/' $f_out
sed -i 's/\#/No./' $f_out
sed -i '/^\s*$/d' $f_out
sed -i 's/, /,/g' $f_out

sed -i -e 's/\([0-9][A-Z][A-Z]\)/\1\#/g' $f_out # These two handle alternative amino acids
sed -i -e 's/\([A-Z]\#\)/ \1/g' $f_out # These two handle alternative amino acids

sed -i -e 's/[[:space:]]\+/\t/g' $f_out

# sed -i -e 's/\t/,/g' $f_out
# sed -i 's/\#//' $f_out

sed -i -e 's/^\t//g' $f_out