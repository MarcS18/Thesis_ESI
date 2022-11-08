#! /bin/bash

cd ML_structures_pyr || exit

find . \( -name "*.cpcm" -o -name "*.smd.out" -o -name "*.xtbtopo.mol" -o -name "*.engrad" -o -name "*.charges" -o -name "*.gbw" -o -name "*.gradient" -o -name "*.lastxtb" -o -name "*.lastxtberr" -o -name "*.proc*" -o -name "*.prop" -o -name "xtb.err" -o -name "*.sh.o*" -o -name "*.sh.e*" -o -name "*.scfp" -o -name "*.molinp" -o -name "*.report" -o -name "*.wbo" -o -name "*.xtbrestart" -o -name "jobscript" -o -name "*.log" -o -name "*.tmp" -o -name "*_property.txt" -o -name "core*" -o -name ".*.sccnotconverged" -o -name ".*.xtboptok" -o -name "coord" -o -name "coord.original" -o -name "cre_members" -o -name "gfnff_topo" -o -name "struc.xyz" -o -name "wbo" -o -name ".history*" -o -name "*.ges" -o -name "*.smd.grd" -o -name "*_atom53.out" -o -name "*.opt" -o -name "*.densities" -o -name "*.gbw1" -o -name "*.hostnames" -o -name "*.lastgrad" -o -name "*.lastint" -o -name "*.lastscf" -o -name "*.Dipoles" -o -name "*.Gradients" -o -name "*.scfgrad.inp" -o -name "*.loc" \) -type f -delete -print

echo "The script has completed"

