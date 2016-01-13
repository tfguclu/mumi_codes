package require autopsf

set filelist [glob protein_res*_mutated_wb.pdb]
set sortedfilelist [lsort -dictionary $filelist]
foreach file $sortedfilelist { 
mol new $file
set strmolid [molinfo top]
autopsf -mol $strmolid -top top_all27_prot_lipid.inp
resetpsf
mol delete $strmolid
}
