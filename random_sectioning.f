      PROGRAM RANDOM_SECTIONING
      PARAMETER (nb = 80)
      PARAMETER (nr = 500000)
c  Note - "c" in column 1 denotes a comment line
c  For details see  "Obtaining 3D particle size distributions from 2D 
c                   section counts", J. Cuzzi and D. Olson, MAPS 2016
c
c  This code calculates counts in nb bins for a lognormal size distribution,
c  generates nr random cuts or apparent diameters per initial particle 
c  size, and rebins samples into histograms also with nb bins.
c
c  We recommend having the latest version of gfortran, although it 
c  isn't critical. This code has been tested on multiple versions 
c  of gfortran.  
c
c How to compile and run the code from a command line in unix or linux:
c relocate into the directory where this file resides and compile by typing,
c     gfortran random_sectioning.f -o random_sectioning  
c Your compile command may differ from "gfortran" depending on your system.
c The executeable file "random_sectioning" can be called anything.
c then execute the program by typing at the prompt,
c    ./random_sectioning
c Lines below need to stay there
      CHARACTER(LEN=120) dummy
      character*50 fileout
      character*20 binning
      REAL*8 x, dj, da(nb), mean, width, mu, sigma, d(nb)
      REAL*8 nv(nb), na(nb), nvT
      real*8 ratio,amp,expon,nv1(nb), delta

c >>>>>>>>>> start parameters you need to change >>>>>>>>>> 
c First set the # of size bins (nb) and # random cuts per bin (nr).
c These are set in the PARAMETER statements at the top of the code.
c Below you set the binning type to 'geometric' or 'arithmetic':
      parameter (binning ='arithmetic')                 
c 
c Below you name the output file you want to create;
c the filename below (in single quotes) hints at input parameters 
c but you can call it anything. 
      parameter (fileout=
     &          'random_arith_160_160_500k_25.dat')
      open(8,file=fileout,status='unknown')
c
c
c Adjustable parameters below determine the initial size distribution.     
c Sampled output apparent diameters will be binned on the same grid.      
      d1    = 25.0     ! = upper boundary of smallest bin     
      ratio = 1.15     ! = geometric spacing (dj+1/dj)
      delta = 25.0     ! = arithmetic spacing (dj+1 - dj)     
c
c Lognormal initial size distribution is assumed; arbitrary units, here microns      
      mean  = 160.	! lognormal mean size
      width = 160.	! lognormal effective width parameter
c >>>>>>>>>> end of parameters you need to change >>>>>>>>>> 


c mu and sigma are parameters called for by lognormal distribution
      mu =log(   mean**2/sqrt(width**2+mean**2)  ) 
      sigma = sqrt(log(1.+ width**2/mean**2))
c writes two-line  header inside file
      write(8,6002) 'mean =',mean,'width =',width,
     &', lognormal, with', nr,'sections -', binning, 'spacing'
      write(8,6003)

c starts loop to determine size distribution
      do ib = 1,nb
      na(ib)=0.
c Bin labeled by upper boundary      
      if (binning == 'geometric') then
          d(ib) = d1*ratio**(float(ib)-1.0d0) ! Geometric binning
      else
          d(ib) = delta*float(ib)
      endif    
      amp = 1.0d3/(sigma*d(ib)*dsqrt(6.28d0))
      expon = dexp(-((dlog(d(ib))-mu)**2)/(2.0d0*sigma**2))
      nv(ib) = amp*expon
c nv(ib) is the initial lognormal size distribution
     
c in Fortran LOG is natural log, LOG10 is base 10
      end do
      
      do ib = 1,nb    
         if (binning == 'geometric') then
           if (ib .eq. 1 ) then
              da(ib) = d(ib) - d(1)/(ratio)
           else   
              da(ib) = d(ib) - d(ib-1)
           end if      
         else  
           da(ib) = delta
         endif  
         
c here we do the random sectioning
         do ir = 1,nr
c rand is the built in fortran random number generator. 
c Returns a uniform random number between 0 and 1
         x=rand(0)  
         dj = d(ib)*sqrt(1.-x**2)
         if(binning == 'geometric') then
             if (dj .ge. d1/ratio) then
              iib = int(log(dj/d1)/log(ratio) + 2.0D0) !geom binning
             endif   
         else 
             iib = int(dj/da(ib) +1.)                  !arith binning
c            write(*,*) ib, iib, dj, da(ib)
         endif 
         
c Here we accumulate the number density of the 2d sections. 
c The factor d(ib)/d(nb) accounts for the sampling bias,  
c or the larger probability of hitting larger spheres.
c For instance if there is one each of a small and a large sphere
c per unit volume, the probability of getting a sample of the 
c large sphere is larger by the ratio of their diameters. 
c The nv factor applies the underlying per-unit-volume distribution
         
         na(iib) = na(iib)+da(ib)*nv(ib)*(d(ib)/d(nb))/float(nr)
           
      end do
c output Nv multiplied by bin width because lognormal above is per-unit-radius 
      nv1(ib) = da(ib)*nv(ib)  
      end do
c normalizes nv by the total number of counts
      nvT=0.
      do j = 1,nb
      nvT = nvT + nv1(j)
      end do

c write out results
      do iib = 1,nb
      write(8,1000), d(iib), na(iib), nv1(iib)/nvT
      enddo

c An example of the output file: 
c 1st column is upper diam (dj) of the bin
c 2nd column is the number of objects with apparent diameter in the bin (Na)
c 3rd column is normalized number density of actual particles in the bin (Nv) 
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c mean =    160.0, width =    160.0, lognormal, with 500000 sections - arithmetic spacing
c      Dj        Na      Nv_orig
c 2.5000E+01 2.5468E+00 9.2382E-02
c 5.0000E+01 6.0868E+00 1.4781E-01
c 7.5000E+01 7.3422E+00 1.4111E-01
c 1.0000E+02 7.3188E+00 1.1825E-01
c .
c .
c .

 1000 format(15(1pe10.4,1x))
 6000 format(15(1pe11.4,1x))
 6001 format(A80)
 6002 format(a6,1x,F8.1,',',1x,a7,1x,F8.1,a17,1x,I6,1x,a10,1x,
     & a10,1x,a7)
 6003 format(5x,'Dj',8x,'Na',6x'Nv_orig')
      close(8)
      end
