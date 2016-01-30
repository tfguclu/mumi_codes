package require autopsf

set filelist [glob 1lit.pdb]
set sortedfilelist [lsort -dictionary $filelist]
foreach file $sortedfilelist {
set filewhext [file rootname $file]
mol new $file
set strmolid [molinfo top]
autopsf -mol $strmolid -top top_all27_prot_lipid.inp
package require solvate
solvate ${filewhext}_autopsf.psf ${filewhext}_autopsf.pdb -t 10 -o ${filewhext}_autopsf_wb
package require autoionize
autoionize -psf ${filewhext}_autopsf_wb.psf -pdb ${filewhext}_autopsf_wb.pdb  -o ${filewhext}_autopsf_wb_ionized -seg ION -sc 0.15 -cation POT -anion CLA -from 5 -between 5
resetpsf
mol delete all
}
