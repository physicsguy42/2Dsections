		FILES IN THIS DIRECTORY:
                                            
AAreadme.txt - this file                                            

>>>>>>> Two fortran source codes (for special versions, see under amatrix.dat):

random_sectioning.f 
	- Program generates nr random cuts (apparent diameters) per initial
	  particle size given a 3d true diameter distribution (Nv) 
	  with a given number of bins (nb).
	- Produces a histogram of particle # (Na) per bin of apparent diameter.
	- Applies a realistic particle size sampling bias factor.
	- Can be used for forward modeling comparison with 2d data .

inverse.f - 
	- Program that takes a 2d apparent diameter histogram (Na) and 
	  "inverts" it to recover a 3d true diameter distribution (Nv).
	- User provides input file of particle number vs apparent diameter.
	- Arithmetic or geometric binning can be used, arbitrary # bins (nb).
	
>>>>>>> Two sample output files from random sectioning code:

random_arith_160_160_500k_25.dat 
	- Sample output from random sectioning program with arithmetic	
	  binning. The parameters used were delta = 25 microns, nb = 80 bins,
 	  nr=500000 sections, mean = width = 160 microns.
	- The values of Nv are normalized by the total number of counts.
	- This file can be read directly by the inversion codes as a check

random_geom_160_160_500k_1_15.dat 
	- Sample output from random sectioning program with geometric
	  binning. The parameters used were c = 1.15 (ratio between bins), 
	  nb = 80 bins, nr=500000 sections, mean = width = 160 microns. 
	- The values of Nv are normalized by the total number of counts.
	- This file can be read directly by the inversion codes as a check
                                      
	Input and output files for inverse code, for arith or geom binning

inverse_lognorm_arith_160_160_25_in.dat 
	- Sample input file for inverse.f with arithmetic binning.
	- columns: apparent diam Dj at upper edge of bin, and # counts per bin (Na).
	- The input file provided is simply a trimmed version of the output
          from the random sectioning code, so the output Nv (below) can be 
          compared with the original Nv in random_arith_160_160_500k_25.dat.

inverse_lognorm_arith_160_160_25_out.dat 
	- Sample output file for inverse.f with arithmetic binning
	  using inverse_lognorm_arith_160_160_25_in.dat as the input file.
	- columns: diam Dj at lower edge, midpoint, & upper edge of bin; 
		then # counts per appar. diam bin (Na), 
		then recovered # counts per true diameter bin (Nv).
	- The values of Nv are normalized by the total number of counts.

inverse_lognorm_geom_160_160_1_15_in.dat 
	- Sample input file for inverse_export with geometric binning.
	- columns: apparent diam Dj at upper edge of bin, and # counts per bin (Na).
	- The input file provided is simply a trimmed version of the output
          from the random sectioning code, so the output Nv (below) can be 
          compared with the original Nv in random_geom_160_160_500k_1_15.dat.

inverse_lognorm_geom_160_160_1_15_out.dat 
	- Sample output file for inverse_export with geometric binning
	  using inverse_lognorm_geom_160_160_1_15_in.dat as input file.
	- columns: diam Dj at lower edge, midpoint, & upper edge of bin; 
		then # counts per appar. diam bin (Na), 
		then recovered # counts per true diameter bin (Nv).
	- The values of Nv are normalized by the total number of counts.

>>>>>>> F_ij matrix (eqns 1, 6, 9) for testing or comparison purposes

amatrix.dat - 
	- written by 'inverse.f'
	- Contains the values of F_ij matrix elements for an 80 bin case. 
	- Used for checking that program is working properly.
	- A specific version 'inverse_test_15.f' is provided to illustrate this;
	  this version outputs k_ij = F_ij/delta as 'kmatrix_15.dat', which agrees 
	  with published tables. See "The Inversion Technique" section of paper).
	- A comparable 15 bin version of 'random_sectioning.f' is also included.

			GENERAL COMMENTS: 

Best to read the paper first to clarify notation.
    
It is necessary to have a Fortran compiler installed on your system. 
We recommend the GNU Fortran package 
(http://directory.fsf.org/wiki/Gfortran),
which runs on unix (including Mac OSX) and linux. One test user with limited
programming experience also recommended "Simply Fortran":
(http://simplyfortran.com), 
which is a commercial product inclusing a user interface, contains basically 
the GNU compiler, and runs on Windows machines. For users unfamiliar with 
running code from the command line this might be a good alternate. The directions 
below assume the code is being run from a command line in a unix (or Mac OSX 
terminal window) environment.

Running RANDOM_SECTIONING.f requires no input file, but you need to define a few
parameters within the source code, using any text editor. The user defines the 
mean radius and width of a lognormal size distribution, the number of bins in 
the histogram of number vs (true) diameter, the upper diameter of the smallest 
bin, the bin width (for arithmetic binning) OR the ratio between bin sizes (for 
geometric binning), and the number of cuts per bin. The output histogram is 
normalized by the total number of counts.  See comments within source
code for details.

For the INVERSE.f code, which recovers the true diameter distribution from some input
histogram of number vs. apparent diameter, the first step is to provide the input 
data file containing your own histogram of counts vs apparent diameter. The structure 
of the input files and meaning of the columns is described above, and in the comment 
lines (lines preceded by "c") of the enclosed codes. These comment lines also 
provide directions as to the commands you need to enter (on the command line) 
to run the codes. The name of an output file must also be chosen.

In the program 'inverse.f' for instance, the input and output files are specified 
within single quotes in the lines:
      PARAMETER (inputFile =
     &'inverse_lognorm_arith_160_160_25_in.dat')
and
      PARAMETER(outputFile=
     & 'inverse_lognorm_arith_160_160_25_out.dat')

you can replace the example file names with whatever filename you want. 
In Fortran it is required to start actual commands on or after column 7 
(a continuation character is added in column 6, as '&' here); these requirements 
should be maintained to keep you out of trouble. It's an old language, but rugged. 

In INVERSE.f, you also need to edit the source code to tell it the number of diameter 
bins in your histogram (nb), but that's it! The code automatically determines if it is 
geometric/arithmetic binning and everything else. 

In running the codes it is assumed here that the source code and input/output files 
are all in the same directory. 

Typically in Fortran one compiles the (human-readable) source code, eg test.f, 
by typing a compile command after the system prompt:
prompt> gfortran test.f -o test   (or, eg., prompt> gfortran inverse.f -o inverse1)
where the "gfortran" command might be called something different on your system
(like f77, f90, f95, ftn, etc..).
The piece "-o test" defines an executable program called "test" which you then 
execute simply by typing its name preceded by characters "./" to tell the computer 
where to look (in this case the current directory):
prompt> ./test                    (or, eg., prompt> ./inverse1)
The output text file appears in the directory and can be accessed in your editor. 


For more details see source code or:

"Obtaining 3D particle size distributions from 2D section counts”", 
Cuzzi and Olson, MAPS 2016

