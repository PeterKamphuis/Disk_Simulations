	program veldust
C-----------------------------------------------------------------------------
C
C    Program to project exponential disks inclined under various
C    position angles. Includes uniform dust layer. This program calculates
C    the projected velocity line profiles. For the moment only the major
C    axis line profiles are calculated. Output will be an IRAF map.
C    Soon an upgrade will follow to a datacube.
C    Compile with g77 -o emissionmodel EmissionModel.f -L/opt/local/lib/ -lcfitsio
C-----------------------------------------------------------------------------

	Parameter (MAXNPIX=2048*2048*100)
	Parameter (MAXNPIXD2=2048*2048)
	Implicit Real*8 (A-H,O-Z)
	real*4 Inc,R0,z0,Rdi,velscale,vellow
	real*8 Incd,R0d,z0d,Rdid,lag
	real*8 dveldisp,cenvd,Distance,kpcscale
	real*8 bdratd,b0d,dbratio,dbuldisp
	real*4 Xlo,Xhi,Zlo,Zhi,t,tau,rdust,veldisp,cenvel
	real*4 veldispchange,veldispouter
	real*8 dveldispchange,dveldispouter,trRdi,bounem(2)
	real*4 rad(1000),rotvel(1000),scalefactor
	real*4 boun(4),bdrat,B0,bratio,buldisp
	real*8 Xlod,Xhid,Zlod,Zhid,vlow,vscale
	real*4 Outdisk(MAXNPIX)
	real*4 Outdiskdusttwodim(MAXNPIXD2)

	integer nxn,nzn,f,i,p
	integer nx,nz,nv,N,nvor
	integer blocksize,status,exists,unit,fpixel
	integer bitpix,naxis,naxes(3),pcount,gcount,group
	integer blocksize2,status2,exists2
	integer bitpix2,naxis2,naxes2(2),group2,unit2,fpixel2

	logical simple2,extend2
	logical Log,Fold,Unproj, Create_Cube
	logical simple,extend

	Character*80 outname, infile, outname2
	character*1 c
	Character*80 rotcur
	blocksize2=1
	blocksize=1
	Pi=3.1415926535897D0
	call get_command_argument(1,infile)
	write (6,*) ' file ',infile, ' detected as input parameter file. '
c	write (6,*) ' Give file with input parameters '
c	read (5,99999) infile
c	infile='Input_Parameters.txt'
99999	format (A)
c	write (6,*) infile
	open (unit=1,file=infile,status='old')
c	write (6,*) ' file ',infile, ' opened '
	read (1,99999) c
	read (1,99999) outname2
	read (1,99999) c
	read (1,*) inc
	read (1,99999) c
	read (1,*) R0
	read (1,99999) c
	read (1,*) z0
	read (1,99999) c
	read (1,*) Rdi
	read (1,99999) c
	read (1,*) t
	read (1,99999) c
	read (1,*) tau
	read (1,99999) c
	read (1,*) bdrat
	read (1,99999) c
	read (1,*) b0
	read (1,99999) c
	read (1,*) bratio
	read (1,99999) c
	read (1,99999) c
	if ((c.eq.'N').or.(c.eq.'n')) then
	   Log=.FALSE.
	else
	   Log=.TRUE.
	end if
	read (1,99999) c
	read (1,99999) c
	if ((c.eq.'N').or.(c.eq.'n')) then
	   Fold=.FALSE.
	else
	   Fold=.TRUE.
	end if
	read (1,99999) c
	read (1,*) trRdi
	read (1,99999) c
	read (1,*) Xlo
	read (1,99999) c
	read (1,*) Xhi
	read (1,99999) c
	read (1,*) Zlo
	read (1,99999) c
	read (1,*) Zhi
	read (1,99999) c
	read (1,*) nx
	read (1,99999) c
	read (1,*) nz
	read (1,99999) c
	read (1,*) N
	read (1,99999) c
	read (1,*) bounem(1)
	read (1,99999) c
	read (1,*) bounem(2)
	read (1,99999) c
	read (1,*) boun(1)
	read (1,99999) c
	read (1,*) boun(2)
	read (1,99999) c
	read (1,*) boun(3)
	read (1,99999) c
	read (1,*) boun(4)
	read (1,99999) c
	read (1,*) Distance
	read (1,99999) c
	read (1,99999) c
	if ((c.eq.'N').or.(c.eq.'n')) then
		 Create_Cube=.FALSE.
	else
		 Create_Cube=.TRUE.
	end if
	If (Create_Cube) then
		read (1,99999) c
		read (1,99999) outname
		read (1,99999) c
		read (1,99999) rotcur
		read (1,99999) c
		read (1,*) lag
		read (1,99999) c
		read (1,*) veldisp
		read (1,99999) c
		read (1,*) veldispouter
		read (1,99999) c
		read (1,*) buldisp
		read (1,99999) c
		read (1,*) vellow
		read (1,99999) c
		read (1,*) velscale
		read (1,99999) c
		read (1,*) cenvel
		write (6,*) cenvel
		read (1,99999) c
		write (6,*) c
		read (1,*) nv
		read (1,99999) c
		write (6,*) nv
	else
		write (6,*) 'Not requiring velocity information'
		outname = 'no_cube.fits'
		rotcur = 'no_rc.txt'
		lag = 0.
		veldisp = 0.
		veldispouter = 0.
		buldisp = 0.
		vellow = 0.
		velscale = 0.
		cenvel = 0.
		nv = 0.
	End if
	close(1)

C The following parameters were in the original input but not used anywhere
	rdust = 0.
C	#  Spatial extent of dust
	Unproj = .FALSE.
C # Test Parameter
C
C	Now adjust the logarithmic boundaries
C
	Distance=Distance*1E3
	kpcscale=Distance*((1./(360.*60.))*2*Pi)
c	write (6,*),Distance,kpcscale,(60./360.)*2*Pi
	write (6,*) 'You have given a Distance of', distance, ' kpc'

	vlow=vellow
	bdratd=bdrat
	b0d=B0
	dbratio=bratio
	dbuldisp=buldisp
	taud=tau
	Incd=Inc
	rdustd=rdust
	R0d=R0
	Rdid=Rdi
	z0d=z0
	td=t

	If (Log) then
	   If (Xlo.le.0.) then
	      Xlod=1.d0
	   else
	      Xlod=Xlo
	   end if
	   If (Zlo.le.0.) then
	      Zlod=1.d0
	   else
	      Zlod=Zlo
	   end if
	   Xhid=Xhi
	   Zhid=Zhi
	   Xlod=Log10(Xlod)
	   Xhid=Log10(Xhid)
	   Zlod=Log10(Zlod)
	   Zhid=Log10(Zhid)
	else
	   Xlod=Xlo
	   Xhid=Xhi
	   Zlod=Zlo
	   Zhid=Zhi
	End if

	Vlow=vellow
	vscale=velscale
	cenvd=cenvel
	do i=1,1000
	   rad(i)=0.
	end do

	If (Create_Cube) then
C	read the file with the rotation curve
C
	   dveldisp=veldisp
	   dveldispchange=veldispchange
	   dveldispouter=veldispouter
	   open(unit=65,file=rotcur,status='old')
	   nr=1
11	   read (65,*,end=21) rad(nr),rotvel(nr)
	   nr=nr+1
	   goto 11
 21	   nr=nr-1

C	Now the main integration loop, using Gauss-Legendre
C	integration.
C
c	write(6,*) rad,'bla',rotvelor
	   nvp2=nv+1
	   nvor=nv
	   write (6,*) 'The dimensions of x=',nx,'z=',nz,' v=',nvor
c	write (6,*) 'bla',nvor

	else
	   dveldisp=0
	   dveldispchange=0
	   dveldispouter=0
	   do i=1,1000
	      rotvel(i)=0.
	   end do


C	Now the main integration loop, using Gauss-Legendre
C	integration.
C

	   nvp2=nv+1
	   nvor=nv
	   write (6,*) 'The dimensions of x=',nx,'z=',nz
	end if

	call calccube(nx,nz,nvp2,outdisk,outdiskdusttwodim,Xlod,Xhid,
     1                Zlod,Zhid,vlow,vscale,cenvd,Incd,R0d,z0d,Rdid,td,
     2                taud,rdustd,Log,Fold,N,Unproj,dveldisp,
     3	              rad,rotvel,nr,boun,bdratd,b0d,dbratio,
     4                dbuldisp,outname,outname2,scalefactor,
     5                dveldispouter,dveldispchange,trRdi,kpcscale,nvor
     6	              ,lag,bounem,rotcur,Distance,Create_Cube)







	return
	end

	subroutine calccube(nx,nz,nvp2,Res,Diskdusttwodim,Xlo,Xhi
     1                ,Zlo,Zhi,vlow,vscale,cenvel,Inc,R0,z0,
     2                Rdi,td,tau,rdust,Log,Fold,N,Unproj,dveldisp,
     3		      rad,rotvelor,nr,boun, bdrat,b0,bratio,buldisp,
     4                outname,outname2,scalefactor,dveldispouter,
     5                dveldispchange,trRdi,kpcscale,nvor,lag,bounem,
     6                rotcur,distance,c_cube)

	implicit REAL*8 (A-H,O-Z)
	real*8 W(2048),XX(2048)
	real*4 rad(nr),rotvelor(nr),rotvel(nr),xxx,scalefactor
	real*8 Inc,R0,z0,csttau,Rdi,dveldisp,cenvel,lag,incout
	real*8 bdrat,b0,bratio,buldisp,dustint,kpcscale
	real*4 Xint(2048),Yint(2048),y2(2048)
	real*4 Res(nx,nz,nvor),Diskdusttwodim(nx,nz)
c	real*4 Res3(nx,nz,nvp2)
	real*8 vprof(2048),vprof2(2048),vadd(2048),vlow,vscale
	real*4 vv(2048),vprn(2048),vpro(2048),result
	real*8 dveldispchange,dveldispouter,trRdi,bounem(2)
	real*4 boun(4)
	real*8 distance

	integer nxn,nzn,N,nv,f,i,p
	integer nx,nz,nr,nvp2,nvor
	integer blocksize,status ,exists ,unit,fpixel
	integer bitpix,naxis,naxes(3),pcount,gcount,group
	integer blocksize2,status2,exists2
	integer bitpix2,naxis2,naxes2(2),group2,unit2,fpixel2

	logical simple2,extend2
	logical Log,Fold,Unproj,c_cube
	logical simple,extend

	character*128 string
	character*80 outname
	character*80 outname2
	character*80 rotcur
	write (6,*) 'You have provided the following input'
	if (c_cube) then
	   write (6,*) 'The dimensions of x=',nx,'z=',nz,' v=',nvor
	else
	   write (6,*) 'The dimensions of x=',nx,'z=',nz
	End if
	write (6,*) 'Physically ranging in X from ', Xlo,'to'
     1	   ,Xhi, 'kpc in X'
 	write (6,*) 'Physically ranging in Z from ', Zlo,'to',
     1	Zhi, 'kpc in Z'
	if (c_cube) then
	   write (6,*) 'Physically ranging in vel from ',
     1	   vlow,'to',vscale*nvor+vlow , ' km/s in velocity'
	End if
c	write (6,*) 'Because nvor=',nvor,'vscale=',vscale
 	write (6,*) 'A Hole from ',bounem(1),' to ',bounem(2),
     1 	' kpc'
	write (6,*) 'The emission disk has a scale length ',R0 ,
     1	' and height',z0
	write (6,*) 'The emission disk has a inclination of ',inc
	write (6,*) 'The disk is truncated at',trRdi ,'kpc'
	write (6,*) 'You have given a Distance of', distance, ' kpc'
	write (6,*) 'The dust disk has a scale length ',Rdi ,
     1	'and height', td
	write (6,*) 'With an optical depth of', tau,
     1	' and a spatial extend from', boun(1), ' to ', boun(2)
	write (6,*)  ' and a second part from', boun(3), ' to ', boun(4)
	write (6,*) 'The bulge to disk ratio is',bdrat,
     1	' with an ellipticity of', bratio
	write (6,*) 'With a scale length of', b0
	write (6,*) 'You have logarithmic binning?', Log
	write (6,*) 'You fold the disk?',Fold
	write (6,*) 'The Number of steps in Gauss-legendre integration'
     1	,N
	if (c_cube) then
	   write (6,*) 'You have a rotation curve'
	   do i=1,nr
	      write (6,*) rad(i),rotvelor(i)
	   end do
	   write (6,*) 'The disk has a central velocity of',cenvel
	   write (6,*) 'With central dispersion ',dveldisp,
     1	   'km/s declining to',dveldispouter
	   write (6,*) 'And a bulge dispersion of', buldisp
	   write (6,*) 'And is lagging with ', lag,'km/s/kpc'
	End if

c	open (unit=77,name='N2048/checkdimake5vel.txt',status='new')
	incout=Inc
	nv=nvp2-1
c	real*4 Res3(nx,nz,nv)
	blocksize=1
	Pi=3.1415926535897D0
	Inc=Abs(Inc)
	If (Inc.gt.90.) Inc=90.
	Inc=Inc*Pi/180.D0
	cosi=cos(Inc)
	sini=sin(Inc)
	if (Fold) then
	   nxn=(nx)/2
	   nzn=(nz)/2
	else
	   nxn=nx+1
	   nzn=nz+1
	end if
	xs=(Xhi-Xlo)/(nxn-1.D0)
	zs=(Zhi-Zlo)/(nzn-1.D0)
	NN=128
	csttau=1.D0
	y=0.D0
	z=0.D0
	if (c_cube) then
	   do i=1,nv
	      vv(i)=vlow+(i-1)*vscale
	   end do
	end if

	call dimake(NN,td,csttau,Rdi,y,z,sini,cosi,pi,
     1		xint,yint,trRdi)
	NA=NN

	csttau=tau/yint(1)

	IF (Fold) then
	   dveldispchange=(dveldisp-dveldispouter)/nxn
	   dveldisp=dveldisp+dveldispchange
	else
	   dveldispchange=((dveldisp-dveldispouter)/((nxn-1)/2.))
	   dveldisp=dveldispouter-dveldispchange
	end if

c	do 2000 p=1,nzn
c	write (6,*) tau,csttau,'Here'
c	open (unit=88,name='N64/check3.txt',status='new')
c	open (unit=77,name='N2048/checkxi20yi20vel.txt',status='new')
c	write (6,*) dveldispchange

	do 1000 i=1,nxn
	if ((5*(i/5)).eq.i) then
	      write (6,*) 'We are at x pixel',i
	end if
	if (c_cube) then
	   IF (Fold) then
	      dveldisp=dveldisp-dveldispchange
	   else
	      IF (i.lt.nxn/2.+1) then

		 dveldisp=dveldisp+dveldispchange

	      else
		 dveldisp=dveldisp-dveldispchange
	      end if
	   end if
	end if
	do 2000 p=1,nzn
	   y=Xlo+(i-1.D0)*xs
	   z=Zlo+(p-1.D0)*zs
	   do f=0,nr
	      if (rotvelor(f).lt.0) then
		 rotvel(f)=rotvelor(f)+(lag*(abs(z)))
	      else
		 rotvel(f)=rotvelor(f)-(lag*(abs(z)))
	      end if
	      if (rotvelor(f).eq.0) then
		 rotvel(f)=0.
	      endif
	   end do

	   if (N.eq.0) then
	         frac=(Zhi-z)/(Zhi-Zlo)
	      NN=int(64. *frac)
	      If (NN.lt.16) NN=16
	   else
	      NN=N
	   end if
	   call dimake(NN,td,csttau,Rdi,y,z,sini,cosi,pi,
     1	        xint,yint,trRdi)

	   NA=NN
	   dum=Pi/2.D0

	   call Gauleg(-dum,dum,XX,W,NN)

	   do ii=1,nv
	      vprof(ii)=0.D0
	      vprof2(ii)=0.D0
	   end do
	   Sum3=0.
	   Do 1600 k=1,NN
	      f2=Func3(XX(k),y,z,R0,z0,sini,cosi,td,tau,csttau,
     1	            xint,yint,y2,NN,rad,rotvel,nr,dveldisp,vadd,
     2		    vlow,vscale,nv,boun,bdrat,b0,bratio,buldisp,
     3              p,k,i,trRdi,bounem,c_cube)


	      Sum3=Sum3+W(k)*f2
	      if (c_cube) then
		 do ii=1,nv
		    vprof2(ii)=vprof2(ii)+W(k)*vadd(ii)
		 end do
		 do ii=1,nv
		    vprof(ii)=vprof(ii)+W(k)*vadd(ii)
		 end do
	      end if

1600	   continue
	   test=0

	   Val=Sum2
	   Val2=Sum3

	   If (Fold) then
	         Diskdusttwodim(nx-nxn+i,nz-nzn+p)=Sum3
	         Diskdusttwodim(nx-nxn+i,nzn-p+2)=Sum3
	         Diskdusttwodim(nxn-i+2,nz-nzn+p)=Sum3
	         Diskdusttwodim(nxn-i+2,nzn-p+2)=Sum3
	   else
	         Diskdusttwodim(i,p)=Sum3
	   end if
	   if (c_cube) then
	      do ii=1,nv
		 vpro(ii)=vprof(ii)
	      end do

	      If (Fold) then
c	      write (6,*) nxn,nzn,p
		 do j=1,nv
		    Res(nx-nxn+i,nz-nzn+p,j)=vprof(j)
		    Res(nx-nxn+i,nzn-p+2,j)=vprof(j)
		 end do
c	   	   	   write (6,*) cenvel,nzn
		 do ii=1,nv
		    xxx=0-vv(ii)
		    call tplint(vv,vpro,nv,xxx,result)
		    vprn(ii)=result

		 end do

		 Val2=0
		 do ii=1,nv
		    Val2=Val2+vprn(ii)
		 end do
		 do j=1,nv
		    Res(nxn-i+2,nz-nzn+p,j)=vprn(j)
		    Res(nxn-i+2,nzn-p+2,j)=vprn(j)

		 end do
	      else

		 do j=1,nv
		    Res(i,p,j)=vprof(j)

		 end do
	      end if
	 end if
2000	 continue
1000	continue
c2000	continue
c	close (unit=88)
	close (unit=77)

c       Check if the file already exists

	naxes(1)=nx
	naxes(2)=nz
	halfax1=(nx/2.)+1
	halfax2=(nz/2.)+1
	cdeltscale1=(xs/kpcscale)*60
	xcrval=0
	cdeltscale2=(zs/kpcscale)*60
	ycrval=0

	if (c_cube) then

	   status=0
	   exist=0
	   blocksize=1
	   unit=0

	   call ftexist(outname, exists, status);
	   write (6,*) 'Does ',outname, ' exists y=1 n=0',exists
	   if (exists.eq.1) then
	      call ftopen(unit,outname,rwmode,
     1	           blocksize,status)
	      call ftdelt(unit, status)
	      call ftexist(outname, exists, status);
	      write (6,*) 'Do we remove it? y=0 n=1',exists
	   end if

C
C	create a new empty fitsfile
C

	   call ftinit(unit,outname,blocksize,status)
	   write (6,*) 'Writing the file ',
     1	              outname,status
C     initialize parameters about the FITS image
	   simple=.true.

	   naxes(3)=nvor
	   halfax3=(nvor/2.)+1

	   intnum=abs(N)
	   naxis=3
	   bitpix=-32
	   extend=.true.
	   call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
	   call ftdkey(unit,'COMMENT',status)
	   call ftdkey(unit,'COMMENT',status)
	   call ftpkyd(unit,'CRPIX1',halfax1,3,'/ ',status)
	   call ftpkyd(unit,'CRVAL1',xcrval,3,'/ ',status)
	   call ftpkyd(unit,'CDELT1',cdeltscale1,3,'/ ',status)
	   call ftpkys(unit,'CTYPE1','ARCSEC1',' ',status)
	   call ftpkys(unit,'CUNIT1','ARCSEC',' ',status)
	   call ftpkyd(unit,'CRPIX2',halfax2,3,'/ ',status)
	   call ftpkyd(unit,'CRVAL2',ycrval,3,'/ ',status)
	   call ftpkyd(unit,'CDELT2',cdeltscale2,3,'/ ',status)
	   call ftpkys(unit,'CTYPE2','ARCSEC2',' ',status)
	   call ftpkys(unit,'CUNIT2','ARCSEC',' ',status)
	   call ftpkyd(unit,'CRPIX3',halfax3,3,'/ ',status)
	   call ftpkyd(unit,'CRVAL3',cenvel,3,'/ ',status)
	   call ftpkyd(unit,'CDELT3',vscale,3,
	1	'PRIMARY PIXEL SEPARATION',status)
c	call ftpkys(unit,'CTYPE3','VELO-HEL',' ',status)
	   call ftpkys(unit,'CTYPE3','VELOCITY',' ',status)
	   call ftpkys(unit,'CUNIT3','KM/S',' ',status)
	   call ftpkyd(unit,'COMMENT',tau,3,'optical depth ',status)
	   call ftpkyd(unit,'COMMENT',incout,3,'inclination ',status)
	   call ftpkyj(unit,'COMMENT',N,'Number of intergrations',status)
	   call ftpkyd(unit,'COMMENT',tau,3,'optical depth ',status)
	   call ftpkyd(unit,'COMMENT',Incout,3,'inclination ',status)
	   call ftpkyd(unit,'COMMENT',R0,3,'Disk scale length',status)
	   call ftpkyd(unit,'COMMENT',z0,3,
	1	'Disk scale height (kpc) ',status)
	   call ftpkyd(unit,'COMMENT',Rdi,3,'Dust scale length',status)
	   call ftpkyd(unit,'COMMENT',td,3,'Dust scale height',status)
	   call ftpkyd(unit,'COMMENT',rdust,3,
	1	'Spatial extent of dust ',status)
	   call ftpkyd(unit,'COMMENT',brat,3,'Bulge to disk ratio '
	1	,status)
	   call ftpkyd(unit,'COMMENT',b0,3,'B0',status)
	   call ftpkyd(unit,'COMMENT',bratio,3,
	1	'Ellipticity of bulge ',status)
	   If (Log) then
	      call ftpkys(unit,'COMMENT','Y',
	1	   ' Logarithmic binning of disk',status)
	   else
	      call ftpkys(unit,'COMMENT','N',
	1	   ' Logarithmic binning of disk',status)
	   endif
	   If (Fold) then
	      call ftpkys(unit,'COMMENT','Y',
	1	   ' Fold disk to make disk symmetric',status)
	   else
	      call ftpkys(unit,'COMMENT','N',
	1	   'Fold disk to make disk symmetric',status)
	   endif
	   If (Unproj) then
	      call ftpkys(unit,'COMMENT','Y',
	1	   'Test parameter ',status)
	   else
	      call ftpkys(unit,'COMMENT','N',
	1	   'Test parameter',status)
	   endif
	   call ftpkys(unit,'COMMENT',rotcur,
	1	' File with rotation curve',status)
	   call ftpkyd(unit,'COMMENT',lag,3,'lag in km/s/kpc  ',status)
	   call ftpkyd(unit,'COMMENT',dveldisp,3,
	1	'Velocity dispersion (isotropic)(central)',status)
	   call ftpkyd(unit,'COMMENT',dveldispouter,3,
	1	'Velocity dispersion (isotropic)(outer) ',status)
	   call ftpkyd(unit,'COMMENT',trRdi,3,
	1	'Truncation Radius (kpc)',status)
	   call ftpkyd(unit,'COMMENT',buldisp,3,
	1	'Velocity dispersion of bulge',status)
	   call ftpkyd(unit,'COMMENT',cenvel,3,
	1	'  Central velocity of galaxy (km/s)',status)
	   call ftpkyd(unit,'COMMENT',Xlo,3,'X lower limit',status)
	   call ftpkyd(unit,'COMMENT',Xhi,3,'X upper limit ',status)
	   call ftpkyd(unit,'COMMENT',Zlo,3,'Z lower limit',status)
	   call ftpkyd(unit,'COMMENT',Zhi,3,'Z upper limit',status)
	   call ftpkyd(unit,'COMMENT',bounem(1),3,
	1	' Inner Hole radius (kpc)',status)
	   call ftpkyd(unit,'COMMENT',bounem(2),3,
	1	' Outer Hole radius (kpc)',status)
	   call ftpkyf(unit,'COMMENT',boun(1),3,'Boundary 1 ',status)
	   call ftpkyf(unit,'COMMENT',boun(2),3,'Boundary 2',status)
	   call ftpkyf(unit,'COMMENT',boun(3),3,'Boundary 3',status)
	   call ftpkyf(unit,'COMMENT',boun(4),3,'Boundary 4 ',status)
	   call ftpkyd(unit,'COMMENT',distance,3,'Distance(kpc) ',status)




c	write (6,*) status,' 3'
	   group=1
	   fpixel=1
	   nelements=naxes(1)*naxes(2)*naxes(3)
c	write (6,*) Res(20,20,25),'check this'
	   call ftppre(unit,group,fpixel,nelements,Res,status)




	   call ftclos(unit, status)

	   call ftfiou(unit, status)
	endif
c	write (6,*) status,' 6'
c	call ftpkyj(unit,'EXPOSURE',1500,'Totalexposure time',status)

c       Check if the second file already exists

	status2=0
	exist2=0
	blocksize2=1
	unit2=0

	call FTEXIST(outname2, exists2, status2);
	write (6,*) 'Does ',outname2, ' exists? y=1 n=0',exists2


	if (exists2.eq.1) then
	   call FTOPEN(unit2,outname2,rwmode2, blocksize2,status2)
	   call FTDELT(unit2, status2)
           call FTEXIST(outname2, exists2, status2);
	   write (6,*) 'Do we remove it? y=0 n=1', exists2
	end if

	call ftinit(unit2,outname2,blocksize2,status2)
	write (6,*) 'Writing the file ',
     1	outname2,status2
C     initialize parameters about the FITS image
	simple2=.true.
	naxes2(1)=nx
	naxes2(2)=nz
	naxis2=2
	bitpix2=-32
	extend2=.true.
	call ftphpr(unit2,simple2,bitpix2,
     1	naxis2,naxes2,0,1,extend2,status2)
	call ftdkey(unit2,'COMMENT',status2)
	call ftdkey(unit2,'COMMENT',status2)
	call ftpkyd(unit2,'CRPIX1',halfax1,3,'/ ',status2)
	call ftpkyd(unit2,'CRVAL1',xcrval,3,'/ ',status2)
	call ftpkyd(unit2,'CDELT1',cdeltscale1,3,'/ ',status2)
	call ftpkys(unit2,'CTYPE1','ARCSEC1',' ',status2)
	call ftpkys(unit2,'CUNIT1','ARCSEC',' ',status2)
	call ftpkyd(unit2,'CRPIX2',halfax2,3,'/ ',status2)
	call ftpkyd(unit2,'CRVAL2',ycrval,3,'/ ',status2)
	call ftpkyd(unit2,'CDELT2',cdeltscale2,3,'/ ',status2)
	call ftpkys(unit2,'CTYPE2','ARCSEC2',' ',status2)
	call ftpkys(unit2,'CUNIT2','ARCSEC',' ',status2)
	call ftpkys(unit2,'COMMENT',outname,'3D companion cube ',status2)


c       	write (6,*) status2,' 3'
	group2=1
	fpixel2=1
	nelements2=naxes2(1)*naxes2(2)
	write (6,*) unit2,group2,fpixel2,nelements2,Diskdusttwodim(50,50)
	call ftppre(unit2,group2,fpixel2,nelements2,
     1	Diskdusttwodim,status2)
	call ftclos(unit2, status2)
	call ftfiou(unit2, status2)
c	write (6,*) status2,' 6'

	return
	end

	function func3(t,y,z,R0,z0,sini,cosi,td,tau,csttau,
     1		xint,yint,y2,NA,rad,rotvel,nr,dveldisp,vadd,
     2	 	vlow,vscale,nv,boun,bdrat,b0,bratio,buldisp,
     3          p,k,l,trRdi,bounem,c_cube)
	implicit REAL*8 (A-H,O-Z)
	character*128 string
	integer NA,nr,nv,p
	logical c_cube
	real*8 t,func3,y,z,R0,z0,sini,cosi,td,tau,tt,csttau
	real*8 dveldisp,vadd(2048),vlow,vscale,noemer,dustint,d1
	real*8 bdrat,b0,bratio,buldisp,noemer1,noemer2,test
	real*8 hulp(2048),trRdi,bounem(2),radialpos
	real*4 boun(4)
	real*4 rad(nr),rotvel(nr),result,y1,res1
	real*4 xint(NA),yint(NA),y2(NA)

        bz0=b0*bratio
	pi=3.1415926535897D0
	radialpos=sqrt(y*y+(z/cosi+tan(t)*sini)*(z/cosi+tan(t)*sini))
	if ((t.le.(-pi/2.D0+1.D-4)).or.(t.ge.(pi/2.D0-1.D-4))
     1	.or.(radialpos.gt.trRdi).or.((radialpos.ge.bounem(1))
     2  .and.(radialpos.le.bounem(2)))) then
c	   write (6,*) 'Radius=',radialpos,'hole',bounem(1),bounem(2)
c	   write (6,*) 'trRadius=',trRdi
	   f2b=0.D0
	   f2d=0.D0
	else
	   tt=tan(t)
	   if ((abs(tau)).lt.1.e-4) then
	      z1=z
c	      z=0
	     fd=exp(-(sqrt(y*y+(z/cosi+tt*sini) *
     1          (z/cosi+tt*sini)))/R0 -
     2           abs(tt*cosi/z0))
c	     z=z1
	     fb=exp(-(sqrt(y*y+(z/cosi+tt*sini) *
     1           (z/cosi+tt*sini)))/b0 -
     2           abs(tt*cosi/bz0))
	     f2d=fd*(1+tt*tt)
	     f2b=bdrat*fb*(1+tt*tt)

	   else

	     d1=dustint(t,NA,xint,yint,y2,boun,y,z,cosi,sini)



	     z1=z

c	     z=0
	     dd2=abs(tt*cosi/z0)
	     dd3=sqrt(y*y+(z/cosi+tt*sini)*(z/cosi+tt*sini))/R0
c	     z=z1
	     db2=abs(tt*cosi/bz0)
	     db3=sqrt(y*y+(z/cosi+tt*sini)*(z/cosi+tt*sini))/b0
	     arg1 =  d1
	     arg2 = d1 + dd2 + dd3
	     arg3 = db2 + db3
C	   write (string,'('' d1 = '',2F10.3)') d1,tt
C	   call gwrits(string,2)
	     f2d=exp(-arg2)*(1+tt*tt)
	     f2b=exp(-arg1)*(bdrat*exp(-arg3)) *
     4		(1+tt*tt)

	   end if
	end if
C	write (string,'('' in func3 '',2F10.3)') f2b,f2d
C	call gwrits(string,2)
C
C	now the velocity profile
C
C	first get the velocity at this radius (y)
C
	if (f2d.lt.1.D-50) then
	   f2d=0.
	end if
	if (c_cube) then
	   y1=y
c	f2d=20
	   call tplint(rad,rotvel,nr,y1,res1)

	   noemer=y*y+(tt*sini+z/cosi)**2

	   if (noemer.eq.0.) then
	      result=0.
	   else
	      result=res1*sini*abs(y)/sqrt(noemer)
	   end if

	   test=0.

	   do i=1,nv
	      vel=vlow+(i-1)*vscale
	      vadd(i)=dexp(-5.D-1*(((vel-result)/dveldisp)**2.))

	   end do
c	test=f2d
	   test=0
	   do i=1,nv
	      test=test+vadd(i)
	   end do
c	test=0
	   do i=1,nv
	      vadd(i)=(vadd(i)*f2d)/test
	   end do



C
C now the bulge, with central velocity 0
C
C	result=0.
	   do i=1,nv
	      vel=vlow+(i-1)*vscale
	      hulp(i)=dexp(-5.D-1*(((vel-result)/buldisp)**2.))
	   end do
	   test=0
	   do i=1,nv
	      test=test+hulp(i)
	   end do
c	write (6,*) f2d,f2b,vadd(10),test

	   do i=1,nv
	      vadd(i)=vadd(i)+((hulp(i)*f2b)/test)
	   end do
	end if
c	write (6,*) f2d,f2b,vadd(10),test
	func3=f2d




	return
	end

	subroutine dimake(NN,td,csttau,Rdi,y,z,sini,cosi,pi,
     1	     xint,yint,trRdi)
	implicit REAL*8 (A-H,O-Z)
	real*8 tt,td,tau,Rdi,y,z,sini,cosi,pi
	real*8 XX(2048),W(2048),csttau,f2,trRdi
	character string*128
	integer NN
	real*4 xint(NN),yint(NN)

c	write (6,*) NN,xint(10),yint(10)

	do i=1,NN
	   Sum2=0.
	   xst=-pi/2.+(i-1.)*pi/NN
	   xend=xst+pi/NN
	   call Gauleg(xst,xend,XX,W,16)
	   Do 160 k=1,NN
c	      if ((XX(k).le.(-pi/2.D0+1.D-4)).or.
c     1           (XX(k).ge.(pi/2.D0-1.D-4))) then
c	         Sum2=Sum2+0.D0
	      if ((XX(k).le.(-pi/2.D0+1.D-4)).or.
     1           (XX(k).ge.(pi/2.D0-1.D-4)).or.
     2	(sqrt(y*y+(z/cosi+dtan(XX(k))*sini)
     3	*(z/cosi+dtan(XX(k))*sini)).gt.trRdi)) then

		 Sum2=Sum2+0.D0
	      else
	         tt=dtan(XX(k))
	         f2=dexp(-(dsqrt(y*y+(z/cosi+tt*sini) *
     1              (z/cosi+tt*sini)))/Rdi -
     2              abs(tt*cosi/td)) *
     3		    (1.+tt*tt)
c    		 write (6,*) k,f2,xx(k),tt,rdi
	         Sum2=Sum2+W(k)*f2
	      end if
160        continue
c	write (6,*) ' after gauleg ',NN,yint(100),csttau,sum2
	   xint(i)=xst
	   yint(i)=Sum2*csttau
	enddo
c	write (6,*) ' end loop of dimake ',NN,yint(10)
	do i=NN-1,1,-1
	   yint(i)=yint(i)+yint(i+1)
	   xint(i)=xint(i)
	enddo
	do i=1,NN-1
	   xint(i)=(xint(i)+xint(i+1))/2.
	enddo
c	write (6,*) ' end of dimake ',NN,yint(10)

C	write (string,'(A,I)') ' In dimake: NN = ',NN
C	call gwrits(string,2)
	return
	end


	function dustint(t,NN,xint,yint,y2,boun,y,z,cosi,sini)
	implicit REAL*8 (A-H,O-Z)
	character*128 string
	real*8 t,dustint,y,z
	real*4 boun(4)
	real*4 xint(NN),yint(NN),y2(NN),tr,result

c	rad=dsqrt(t*t+y*y+z*z)
	rad=sqrt(y*y+(z/cosi+tan(t)*sini)*(z/cosi+tan(t)*sini))
	if (((rad.ge.boun(1)).and.(rad.le.boun(2))).or.
     1      ((rad.ge.boun(3)).and.(rad.le.boun(4)))) then
	   tr=t
	   call tplint(xint,yint,NN,tr,result)
	   if (result.lt.0.) then
	      result=0.
	   end if
	   dustint=result
c	   write (6,*) dustint
	else
	   dustint=0.D0
	end if
C	write (string,'(A,4F)') ' In dustint: NN = ',rad,dustint
C	call gwrits(string,2)
	return
	end

	SUBROUTINE GAULEG(X1,X2,X,W,N)
	IMPLICIT REAL*8 (A-H,O-Z)
	REAL*8 X1,X2,X(N),W(N)
	Character*128 String
	PARAMETER (EPS=3.D-14)
	M=(N+1)/2
	XM=0.5D0*(X2+X1)
	XL=0.5D0*(X2-X1)
	DO I=1,M
		Z=COS(3.141592654D0*(I-.25D0)/(N+.5D0))
c		write (6,*) z,i,n
1		CONTINUE
			P1=1.D0
			P2=0.D0
			DO J=1,N
				P3=P2
				P2=P1
				P1=((2.D0*J-1.D0)*Z*P2-(J-1.D0)*P3)/J
c				write (6,*) j,z,p1,p2,p3
			ENDDO
			PP=N*(Z*P1-P2)/(Z*Z-1.D0)
			Z1=Z
			Z=Z1-P1/PP
c			write (6,*) z,z1,z-z1
		IF(ABS(Z-Z1).GT.EPS)GO TO 1
		X(I)=XM-XL*Z
		X(N+1-I)=XM+XL*Z
c		write (6,*) M,N,I
		W(I)=2.D0*XL/((1.D0-Z*Z)*PP*PP)
		W(N+1-I)=W(I)
	ENDDO
C	DO I=1,N
C  	   write (string,'(I,2F,A)') i,W(I),X(I),' in Gauleg '
C    	   call gwrits(string,2)
C	ENDDO
	RETURN
	END


	SUBROUTINE SPLINE(X,Y,N,YP1,YPN,Y2)
	PARAMETER (NMAX=100)
	DIMENSION X(N),Y(N),Y2(N),U(NMAX)
	IF (YP1.GT..99E30) THEN
		Y2(1)=0.
		U(1)=0.
	ELSE
		Y2(1)=-0.5
		U(1)=(3./(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
	ENDIF
	DO I=2,N-1
		SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
		P=SIG*Y2(I-1)+2.
		Y2(I)=(SIG-1.)/P
		U(I)=(6.*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))
     *			/(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
	ENDDO
	IF (YPN.GT..99E30) THEN
		QN=0.
		UN=0.
	ELSE
		QN=0.5
		UN=(3./(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
	ENDIF
	Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1.)
	DO K=N-1,1,-1
		Y2(K)=Y2(K)*Y2(K+1)+U(K)
	ENDDO
	RETURN
	END

	SUBROUTINE SPLINT(XA,YA,Y2A,N,X,Y)
	DIMENSION XA(N),YA(N),Y2A(N)
	character*128 string
	KLO=1
	KHI=N
1	IF (KHI-KLO.GT.1) THEN
		K=(KHI+KLO)/2
		IF(XA(K).GT.X)THEN
			KHI=K
		ELSE
			KLO=K
		ENDIF
	GOTO 1
	ENDIF
C	H=XA(KHI)-XA(KLO)
C	IF (H.EQ.0.) PAUSE 'Bad XA input.'
C	A=(XA(KHI)-X)/H
C	B=(X-XA(KLO))/H
C	Y=A*YA(KLO)+B*YA(KHI)+
C     *		((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6.
C	if ((y.le.0.).and.(x.le.0.)) then
C	write (string,'(2I,2F10.3,2G10.3,F10.3)') khi,klo,ya(khi),
C     1	      ya(klo),y2a(khi),y2a(klo),y
C	call gwrits(string,2)
C	end if
	Y=YA(KLO)+(X-XA(KLO))*(YA(KHI)-YA(KLO))/(XA(KHI)-XA(KLO))
	RETURN
	END




	SUBROUTINE TPLINT(XA,YA,N,X,Y)
	DIMENSION XA(N),YA(N)
	character*128 string
	KLO=1
	KHI=N
1	IF (KHI-KLO.GT.1) THEN
		K=(KHI+KLO)/2
		IF(XA(K).GT.X)THEN
			KHI=K
		ELSE
			KLO=K
		ENDIF
	GOTO 1
	ENDIF
C	H=XA(KHI)-XA(KLO)
C	IF (H.EQ.0.) PAUSE 'Bad XA input.'
C	A=(XA(KHI)-X)/H
C	B=(X-XA(KLO))/H
C	Y=A*YA(KLO)+B*YA(KHI)+
C     *		((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6.
C	if ((y.le.0.).and.(x.le.0.)) then
C	end if
	Y=YA(KLO)+(X-XA(KLO))*(YA(KHI)-YA(KLO))/(XA(KHI)-XA(KLO))
	RETURN
	END
