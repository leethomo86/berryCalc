# berryCalc

A program that can read a series of Gaussian matrix files containing a determinant expansion (molecular orbitals and coefficients) to determine the geometric phase of a closed loop in nuclear configuration space. Perl scripts to generate closed loops using several programs can be found in this repository. Currently, the following methods are implemented:
* SCF (using Gaussian 16)
* Wheeler-Hill projection (https://github.com/thompsonresearchgroup/Wheeler-HIll-PUHF)
* NOCI (https://github.com/thompsonresearchgroup/ResHF)

To compile the program, edit the makefile included as required and type make.

For further information on running the program type ./berryCalc.exe --help
