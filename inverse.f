      PROGRAM INVERSE
c  Note - "c" in column 1 denotes a comment line
c  For details see  "Obtaining 3D particle size distributions from 2D 
c                   section counts", J. Cuzzi and D. Olson, MAPS 2016
c
c This code inverts the binned apparent diameter histogram to get the 
c true diameter distribution, for an arbitrary number of bins.
c User specifies the input and output files as "inputFile" and "outputFile".
c User supplies the number of size bins (nb) in the input data file 
c (not including the two line header).
c
c  We recommend having the latest version of gfortran, although it 
c  isn't critical. This code has been tested on multiple versions 
c  of gfortran.  
c
c How to compile and run the code from a command line in unix or linux:
c relocate into the directory where this file resides and compile by typing,
c 	gfortran inverse.f -o inverse  
c Your compile command may differ from "gfortran" depending on your system.
c The executeable file "inverse" can be called anything.
c then execute the program by typing at the prompt,
c 	./inverse 

c >>> here you need to define your number of diameter bins, nb >>>>>>>
      INTEGER nb
      PARAMETER (nb = 80)
      CHARACTER(LEN=70)inputFile
      CHARACTER(LEN=70)outputFile
c>>>> then define the input and output filenames below >>>
c>>>> filenames are inside single quotes; can be anything (< 70 characters)>>>
c>>>> here filenames are mnemonic for size distribution parameters >>>
      PARAMETER (inputFile =
     &           'inverse_lognorm_arith_160_160_25_in.dat')     
      PARAMETER(outputFile=
     &           'inverse_lognorm_arith_160_160_25_out.dat')
c>>>>    the code does the rest    >>>>

C23456789012345678901234567890123456789012345678901234567890123456789012
      REAL*8 a(nb,nb),y(nb,nb),yy(nb,nb),prod(nb,nb),adup(nb,nb),nvT
      REAL*8 dj(nb),ns(nb),nv(nb),na(nb),nt(nb),nq(nb),da, c, f1,f2      
      INTEGER indx(nb)           
      CHARACTER(LEN=90) dummy1, dummy2
      CHARACTER(LEN=21) matrixFile
      PARAMETER (matrixFile = 'amatrix.dat')      
      CHARACTER(LEN=20) Binning

           
      open(8,file=inputFile,status='unknown')
      open(9,file=outputFile,status='unknown')
      open(10,file=matrixFile,status='unknown')      

      n = nb
c read in histogram data
      read(8,6001)     dummy1 ! Skip one line
      read(8,6001)     dummy2 ! Skip one line
      do j = 1,nb
      read(8,*)     dj(j), na(j)
      end do
      
c Code automatically determines if binning is geometric or arithmetic
      if ( dj(3)-dj(2) .gt. dj(2) - dj(1) ) then       
          Binning = 'GEOMETRIC'
      else      
          Binning = 'ARITHMETIC'
      end if        
c use upper bin (second column of file) since these are apparent diameters     
      da = dj(2) - dj(1) ! assume uniform sampling for arithmetic Binning      
      c = dj(2)/dj(1)    ! geometric binning ratio   
      
c the unfolding will replace apparent diameter histograms with recovered nv,
c so we copy the original values into new array.
      do j = 1,n
      
         ns(j) = na(j)
         
      end do      
      
c here we calculate the F matrix 
c using eqns 6 (for arithmetic binning) or 8 (geometric) of the paper, 
c In the paper, roles of j and i differ from standard array notation

      if(Binning .eq. 'ARITHMETIC') then

      do j=1,n   ! using equation 6
            do i=1,j-1
               f1 = sqrt((j-0.5)**2 -(i-1)**2)  
               f2 = sqrt((j-0.5)**2 - i**2)
               a(i,j)= da*(f1-f2)
            end do
            
            do i=j+1,n
               a(i,j)=0.
            end do
      a(j,j)= da*sqrt(j-0.75)
      end do
      elseif (Binning .eq. 'GEOMETRIC') then
      
      do j=1,n  ! Using equation 8
          do i=1,j-1
              f1 = (c**(j-1))*sqrt(1-c**(2*(i-j-1)))
              f2 = (c**(j-1))*sqrt(1-c**(2*(i-j)))
              a(i,j) = dj(1)*(f1-f2)
              
          end do
          do i=j+1,n
               a(i,j)=0.0
          end do
      a(j,j) = dj(1)*c**(j-1)*sqrt(1-c**(-2))    
      end do
      end if
      
c The matrix A will be destroyed; thus duplicate original above
c for testing purposes
      do i=1,n
      do j=1,n
      adup(i,j)=a(i,j)
      end do
      end do
c A is upper triangular so we use back substitution to solve the system.
c Routine is based on lubksb from Numerical recipes.
c (section 2.3 page 39 of Numerical Recipes in Fortran 77 Press etal 1992).
c The system is upper trangular so we don't need to do 
c a LU decomposition and just use the A matrix itself.

      call mybacksub(a,n,nb,ns)
            
c ns now has the components of the recovered number density nv.

       do i = 1,n
       write (10,1000), (adup(i,j),j=1,n)

       end do

c recall that ns was replaced by the recovered values by the 
c backsubstition method. We now copy those into the nv array.

       
       nvT=0.
       do j = 1, n
           nv(j) = ns(j)
           nvT = nvT + nv(j)      
       end do
       do j=1,n
       nv(j) = nv(j)/nvT
       end do
c nv(j) is normalized; retaining units in the transform left to "the reader"       
       write(9,6001) dummy1
       write(9,6003) 
       DO j = 1, nb

c print bin boundaries and number densities       
       if (Binning .eq. 'GEOMETRIC' ) then       
           write(9,6000) dj(j)/c, dj(j)/dsqrt(c), dj(j),na(j),nv(j)
       else if(Binning .eq. 'ARITHMETIC') then    
           write(9,6000) dj(j)-DA,dj(j)-DA/2.0d0,dj(j),na(j),nv(j)
       end if
       end do  
       
 1000 format(15(1pe8.1,1x))
 6000 format(15(1pe11.4,1x))
 6001 format(A80)
 6003 format(3x,'D_low',7x,'D_mid',8x,'D_up',7x,'Na_obs',
     &3x,'Nv_recovered')
      close(8)
      close(9)
      close(10)
      end

c >>>>>>>>>>>>>> Other notes:
c To obtain the recovered nv vector we solve the inverse of the system 
c A*nv = na, where A is an nb by nb matrix, and na and nv are vectors of size nb. 
c The inverse is determined using back substitution in the matrix A (called 'amatrix'). 
c
c "amatrix.dat" contains the values of the matrix A specified in the paper as F_ij.
c This is generated for testing purposes. 80 arithmetic bin case is enclosed.
c This instance of amatrix contains 80 rows of 80 values each. 
c For the 15 bin case mentioned in the paper as tabulated by Weibel in section
c "The Inversion Technique", amatrix is 15 rows by 15 columns. 
c Geometric binning needs a different A than for arithmetic bin spacing. 
c
c
c Input file has two columns; looks like:
c
c mean =160.000, width =160.000,lognormal, with 500000 sections - arithmetic spacing
c      Dj        Na      
c 2.5000E+01 2.5468E+00
c 5.0000E+01 6.0868E+00
c 7.5000E+01 7.3422E+00
c 1.0000E+02 7.3188E+00
c .
c .
c 1st column is the upper bin diameter 
c 2nd column is the number of samples in that (apparent) diameter bin.
c
c Output file will have 5 columns; looks like:
c
c mean =    160.0, width =    160.0, lognormal, with 500000 sections - arithmetic 
c    D_low       D_mid        D_up       Na_obs   Nv_recovered
c  0.0000E+00  1.2500E+01  2.5000E+01  2.5468E+00  9.3902E-02
c  2.5000E+01  3.7500E+01  5.0000E+01  6.0868E+00  1.4752E-01
c  5.0000E+01  6.2500E+01  7.5000E+01  7.3422E+00  1.4282E-01
c .
c .
c 1st column is the lower bin boundary,
c 2nd column is the midpoint (arithmetic or geometric mean)
c 3rd column is the upper bin boundary, 
c 4th column is original number per apparent diameter bin (na), 
c 5th column is the recovered number per true diameter bin (nv)

c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      SUBROUTINE  mybacksub(a,n,nb,b)
     
      integer  n, nb
     
      real*8   a(n,nb),b(n),sum
     
c Solve system AX = B with the asumption that A is upper triangular 

c the components of B will be replaced with those of the solution x
c For our particular case A is the matrix F in the paper and B is vector 
c containing na values. The vector x is the unfolded nv values.

c addapted from Numerical recipies routine lubksub 
c (in our case LU decompision is not needed).

      do i=n,1,-1
         sum = b(i)
         do j=i+1,n
            sum=sum-a(i,j)*b(j)
         end do
         b(i) = sum/a(i,i)
      end do
     
      end
     
          
