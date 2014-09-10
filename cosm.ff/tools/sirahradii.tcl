#
# TCL script to set VdW radii of SIRAH atom types in VMD
# Definitions will be applyed to existing molecules and
# to new loaded onces. Use in combination with an X-PLOR
# PSF or AMBER topology file.
#
# AUTHOR:  Matias Machado
# E-MAIL:  mmachado@pasteur.edu.uy
# 
# SIRAH version [Aug 2013]
#
# USAGE
# VMD Tk console: source sirahradii.tcl
# VMD start up:   vmd -e sirahradii.tcl
#

proc def_sirah_radii {args} {

     lassign $args fname molid 
     
                         # type  sigma(nm)
     set sirah_atomtypes {
                           GNaz  0.65000
                           GObz  0.65000
                           GNan  0.40000
                           GObn  0.40000
                           GN    0.40000
                           GO    0.40000
                           GC    0.40000
                           Y1C   0.40000
                           Y2Ca  0.41000
                           Y3Sm  0.47000
                           Y4Cv  0.40000
                           Y5Sx  0.50000
                           Y6Cp  0.43000
                           A1C   0.35000
                           A1Cw  0.35000
                           A2C   0.45000
                           A3P   0.35000
                           A4O   0.35000
                           A5No  0.40000
                           A6No  0.40000
                           A5Nu  0.40000
                           A6Nu  0.40000
                           A7N   0.40000
                           A8P   0.35000
                           P1O   0.50000
                           P1S   0.50000
                           P2P   0.35000
                           P3Cn  0.40000
                           P3Cq  0.40000
                           P4O   0.35000
                           P5N   0.35000
                           C1Ck  0.47000
                           C2Cr  0.47000
                           C3Cr  0.40000
                           C4Cd  0.40000
                           C4Ce  0.40000
                           C5N   0.45000
                           C6O   0.45000
                           C7Nk  0.55000
                           C8C   0.40000
                           PX    0.46327
                           KX    0.42906
                           KN    0.33997
                           CX    0.26698
                           NF    0.32500
                           NL    0.32500
                           NR    0.32500
                           NS    0.32500
                           NU    0.32500
                           NW    0.32500
                           NX    0.32500
                           OV    0.29599
                           OX    0.29599
                           OY    0.29599
                           OZ    0.29599
                           WT    0.42000
                           NaW   0.58000
                           KW    0.64500
                           ClW   0.68000
                         }

     foreach {type sigma} $sirah_atomtypes {
             
             set my_sel [atomselect $molid "type $type"]
             
             #       vdw radius = rm/2 = 2^(1/6) * sigma / 2
             $my_sel set radius [ expr pow(2.0,(1.0/6))*$sigma*10.0/2.0 ]
            
             $my_sel delete
     }

}

# Set SIRAH radii for existing molecules
foreach mol [array name vmd_initialize_structure] {

        if ($vmd_initialize_structure($mol)) { def_sirah_radii vmd_initialize_structure $mol w }
}

# Set SIRAH radii for new molecules
trace variable vmd_initialize_structure w def_sirah_radii

puts "SIRAH VdW radii have been set! No need to source the script again during the current VMD session" 
