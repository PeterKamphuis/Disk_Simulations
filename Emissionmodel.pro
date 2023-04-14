Pro EmissionModel


disp=''
WHILE (disp EQ '') DO BEGIN
   print, 'Give file with input parameters '
   disp='../Exampleparameter5.txt'
;   read,disp
ENDWHILE
filename=disp
filelength=FILE_LINES(filename)

h=''
counter=0.
close,1
openr,1,filename
for i=0,filelength-1 do begin
   readf,1,h
   tmp = STRMID( h,0,1 )
   IF tmp NE '#' then begin
      case counter of
         0: outname=strtrim(strcompress(h),2)
         1: outname2=strtrim(strcompress(h),2)
         2: inc=double(h)
         3: R0=double(h)
         4: z0=double(h)
         5: Rdi=double(h)
         6: t=double(h)
         7: tau=double(h)
         8: rdust=double(h)
         9: bdrat=double(h)
         10: b0=double(h)
         11: bratio=double(h)
         12: Log=STRUPCASE(strtrim(strcompress(h),2))
         13: Fold=STRUPCASE(strtrim(strcompress(h),2))
         14: rotcur=strtrim(strcompress(h),2)
         15: lag=double(h)
         16: veldisp=double(h)
         17: veldispouter=double(h)
         18: trRdi=double(h)
         19: buldisp=double(h)
         20: vellow=double(h)
         21: velscale=double(h)
         22: cenvel=double(h)
         23: Unproj=STRUPCASE(strtrim(strcompress(h),2))
         24: Xlo=double(h)
         25: Xhi=double(h)
         26: Zlo=double(h)
         27: Zhi=double(h)
         28: nx=double(h)
         29: nz=double(h)
         30: nv=double(h)
         31: N=double(h)
         32: begin
            bounem=dblarr(2)
            bounem[0]=double(h)
         end
         33: bounem[1]=double(h)
         34: begin
            boun=dblarr(4)
            boun[0]=double(h)
         end
         35: boun[1]=double(h)
         36: boun[2]=double(h)
         37: boun[3]=double(h)
         38: Distance=double(h)
         else: break
      endcase
      counter++
   endif
endfor
close,1
Distance=Distance*1E3
kpcscale=Distance*((1.d0/(360.d0*60.d0))*2.d0*!dpi)
print, 'You have given a Distance of', distance, ' kpc'
vlow=vellow




td=t
dveldisp=veldisp
dveldispouter=veldispouter
If Log NE 'N' then begin
   If Xlo le 0. then begin
      Xlod=1.d0
   endif else begin
      Xlod=double(Xlo)
   endelse
   If Zlo le 0.  then begin
      Zlod=1.d0
   endif else begin
      Zlod=double(Zlo)
   endelse
   Xhid=double(Xhi)
   Zhid=double(Zhi)
   Xlo=ALOG10(Xlod)
   Xhi=ALOG10(Xhid)
   Zlo=ALOG10(Zlod)
   Zhi=ALOG10(Zhid)
endif
Vlow=vellow
vscale=velscale
readcol,rotcur,rad,rotvelor,FORMAT=('D,D')
nr=n_elements(rad)
nvp2=nv+1
nvor=nv
Res =dblarr(nx,nz,nvor)
Diskdusttwodim=dblarr(nx,nz)
print,'You have provided the following input'
print,'The dimensions of x=',nx,'z=',nz,' v=',nvor
print,'Physically ranging in X from ', Xlo,'to',Xhi, 'kpc in X'
 print,'Physically ranging in Z from ', Zlo,'to', Zhi, 'kpc in Z'
print,'Physically ranging in vel from ', vlow,'to',vscale*nvor+vlow , ' km/s in velocity'
print,'Because nvor=',nvor,'vscale=',vscale
 print,'A Hole from ',bounem[0],' to ',bounem[1], ' kpc'
print,'The emission disk has a scale length ',R0 , ' and height',z0
print,'The emission disk has a inclination of ',inc
print,'The disk is truncated at',trRdi ,'kpc'
print,'You have given a Distance of', distance, ' kpc'
print,'The dust disk has a scale length ',Rdi , 'and height', td
print,'With an optical depth of', tau, 'and a spatial extend of', rdust
print,'The bulge to disk ratio is',bdrat,' with an ellipticity of', bratio
print,'With a scale length of', b0
print,'You have logarithmic binning?', Log
print,'You fold the disk?',Fold
print,'You have a rotation curve'
for i=0,nr-1 do begin
   print,rad[i],rotvelor[i]
endfor
print,'The disk has a central velocity of',cenvel
print,'With central dispersion ',dveldisp, 'km/s declining to',dveldispouter
print,'And a bulge dispersion of', buldisp
print,'And is lagging with ', lag,'km/s/kpc'
print,'The Boundaries are',boun
print,'The Number of steps in Gauss-legendre integration' ,N

incout=Inc
nv=nvp2-1
Inc=Abs(Inc)
If Inc gt 90. then Inc=90.d0
Inc=(Inc*!Pi)/180.D0
cosi=cos(Inc)
sini=sin(Inc)
if Fold NE 'N' then begin
   nxn=(nx)/2.d0
   nzn=(nz)/2.d0
endif else begin
   nxn=nx+1.d0
   nzn=nz+1.d0
endelse
xs=(Xhi-Xlo)/(nxn-1.D0)
zs=(Zhi-Zlo)/(nzn-1.D0)
NN=128.D0
csttau=1.D0
y=0.D0
z=0.D0
vv=dblarr(nv)
for i=0.d0,nv-1 do begin
   vv[i]=vlow+(i-1)*vscale
endfor

dimake,NN,td,csttau,Rdi,y,z,sini,cosi,xint,yint,trRdi

NA=NN
csttau=tau/yint[0]

IF Fold NE 'N' then begin
   dveldispchange=double((dveldisp-dveldispouter)/nxn)
   dveldisp=double(dveldisp+dveldispchange)
endif else begin
   dveldispchange=double(((dveldisp-dveldispouter)/((nxn-1)/2.)))
   dveldisp=double(dveldispouter-dveldispchange)
endelse
for i=0.d0,nxn-2 do begin
   if ((5*fix(i/5)) eq i) then print,'We are at x pixel',i
   IF Fold NE 'N' then begin
      dveldisp=dveldisp-dveldispchange
   endif else begin
      IF  i lt nxn/2.+1 then begin
         dveldisp=dveldisp+dveldispchange
      endif else begin
         dveldisp=dveldisp-dveldispchange
      endelse
   endelse
   for p=0.d0,nzn-2 do begin
      y=double(Xlo+(i)*xs)
      z=double(Zlo+(p)*zs)
      rotvel=dblarr(n_elements(rotvelor))
      for f=0.d,nr-1 do begin
         if (rotvelor[f] lt 0) then $
            rotvel[f]=rotvelor[f]+(lag*(abs(z))) else $
               rotvel[f]=rotvelor[f]-(lag*(abs(z)))
         if (rotvelor[f] eq 0) then rotvel[f]=0.
      endfor
      if (N eq 0) then begin
         frac=(Zhi-z)/(Zhi-Zlo)
         NN=int(64. *frac)
         If (NN lt 16) then NN=16
      endif else NN=N


      dimake,NN,td,csttau,Rdi,y,z,sini,cosi,xint,yint,trRdi
      NA=NN
      dum=!dpi/2.d0
      Gauleg,-dum,dum,XX,W,NN

      vprof=dblarr(nv)
      vprof2=dblarr(nv)
      Sum3=0.d0
      for k=0,NN-1 do begin
         f2=0.D0
         Func3,f2,XX(k),y,z,R0,z0,sini,cosi,td,tau,csttau,$
               xint,yint,y2,NN,rad,rotvel,nr,dveldisp,vadd,$
               vlow,vscale,nv,boun,bdrat,b0,bratio,buldisp,$
               p,k,i,trRdi,bounem
         Sum3=Sum3+W[k]*f2
         vprof[*]=vprof[*]+W[k]*vadd[*]
      ;   for ii=0,nv-1 do begin
      ;      vprof[ii]=vprof[ii]+ W[k]*vadd[ii]
      ;      print, W[k] , vadd[ii],vprof[ii] , k,ii,i,p ,'here'
;         print,W[k],vadd[2],k,i,p,'here'
      ;   endfor

      endfor
      if i eq 0 then stop
      Val=0.D0
      Val2=double(Sum3)
      If Fold NE 'N' then begin
         Diskdusttwodim[nx-nxn+i,nz-nzn+p]=Sum3
         Diskdusttwodim[nx-nxn+i,nzn-p+2]=Sum3
         Diskdusttwodim[nxn-i+2,nz-nzn+p]=Sum3
         Diskdusttwodim[nxn-i+2,nzn-p+2]=Sum3
      endif else begin
         Diskdusttwodim[i,p]=Sum3
      endelse
      vpro=dblarr(n_elements(vprof))
      vpro[*]=vprof[*]
      If Fold NE 'N' then begin
         for j=0,nv-1 do begin
            Res[nx-nxn+i,nz-nzn+p,j]=vprof[j]
            Res[nx-nxn+i,nzn-p+2,j]=vprof[j]
         endfor
         for ii=0,nv-1 do begin
            xxx=0-vv[ii]
            tplint,vv,vpro,nv,xxx,result
            vprn[ii]=result
         endfor
         Val2=0.D0
         for ii=0,nv-1 do begin
            Val2=Val2+vprn[ii]
         endfor
         for j=0,nv-1 do begin
            Res[nxn-i+2,nz-nzn+p,j]=vprn[j]
            Res[nxn-i+2,nzn-p+2,j]=vprn[j]
         endfor
      endif else begin
         for j=0,nv-1 do begin
            Res[i,p,j]=vprof[j]
           ; print,vprof[j],j,p,i,'here'
         endfor
      endelse
   endfor
endfor
print,XX[NN-1]
MKHDR, header,Res
halfax1=(nx/2.)+1
halfax2=(nz/2.)+1
halfax3=(ABS(vlow)/vscale)+1
cdeltscale1=(xs/kpcscale)*60
xcrval=0
cdeltscale2=(zs/kpcscale)*60
ycrval=0
intnum=abs(N)
sxaddpar,header,'CRPIX1',halfax1
sxaddpar,header,'CRVAL1',xcrval
sxaddpar,header,'CDELT1',cdeltscale1
sxaddpar,header,'CTYPE1','ARCSEC1'
sxaddpar,header,'CUNIT1','ARCSEC'
sxaddpar,header,'CRPIX2',halfax2
sxaddpar,header,'CRVAL2',ycrval
sxaddpar,header,'CDELT2',cdeltscale2
sxaddpar,header,'CTYPE2','ARCSEC2'
sxaddpar,header,'CUNIT2','ARCSEC'
sxaddpar,header,'CRPIX3',halfax3
sxaddpar,header,'CRVAL3',cenvel
sxaddpar,header,'CDELT3',vscale,'PRIMARY PIXEL SEPARATION'
sxaddpar,header,'CTYPE3','VELOCITY'
sxaddpar,header,'CUNIT3','KM/S'
sxaddpar,header,'COMMENT',tau,'optical depth '
sxaddpar,header,'COMMENT',incout,'inclination '
sxaddpar,header,'COMMENT',N,'Number of intergrations'
sxaddpar,header,'COMMENT',tau,'optical depth '
sxaddpar,header,'COMMENT',Incout,'inclination '
sxaddpar,header,'COMMENT',R0,'Disk scale length'
sxaddpar,header,'COMMENT',z0,'Disk scale height (kpc) '
sxaddpar,header,'COMMENT',Rdi,'Dust scale length'
sxaddpar,header,'COMMENT',td,'Dust scale height'
sxaddpar,header,'COMMENT',rdust,'Spatial extent of dust '
sxaddpar,header,'COMMENT',bdrat,'Bulge to disk ratio '
sxaddpar,header,'COMMENT',b0,'B0'
sxaddpar,header,'COMMENT',bratio,'Ellipticity of bulge '
sxaddpar,header,'COMMENT',Log,' Logarithmic binning of disk'
sxaddpar,header,'COMMENT',Fold,' Fold disk to make disk symmetric'
sxaddpar,header,'COMMENT',Unproj,'Test parameter'
sxaddpar,header,'COMMENT',rotcur,' File with rotation curve'
sxaddpar,header,'COMMENT',lag,'lag in km/s/kpc  '
sxaddpar,header,'COMMENT',dveldisp,'Velocity dispersion (isotropic)(central)'
sxaddpar,header,'COMMENT',dveldispouter,'Velocity dispersion (isotropic)(outer) '
sxaddpar,header,'COMMENT',trRdi,'Truncation Radius (kpc)'
sxaddpar,header,'COMMENT',buldisp,'Velocity dispersion of bulge'
sxaddpar,header,'COMMENT',cenvel,'  Central velocity of galaxy (km/s)'
sxaddpar,header,'COMMENT',Xlo,'X lower limit'
sxaddpar,header,'COMMENT',Xhi,'X upper limit '
sxaddpar,header,'COMMENT',Zlo,'Z lower limit'
sxaddpar,header,'COMMENT',Zhi,'Z upper limit'
sxaddpar,header,'COMMENT',bounem[0],' Inner Hole radius (kpc)'
sxaddpar,header,'COMMENT',bounem[1], ' Outer Hole radius (kpc)'
sxaddpar,header,'COMMENT',boun[0],'Boundary 1 '
sxaddpar,header,'COMMENT',boun[1],'Boundary 2'
sxaddpar,header,'COMMENT',boun[2],'Boundary 3'
sxaddpar,header,'COMMENT',boun[3],'Boundary 4 '
sxaddpar,header,'COMMENT',distance,'Distance(kpc) '
writefits,outname,res,header



MKHDR, header,Diskdusttwodim
sxaddpar,header,'CRPIX1',halfax1
sxaddpar,header,'CRVAL1',xcrval
sxaddpar,header,'CDELT1',cdeltscale1
sxaddpar,header,'CTYPE1','ARCSEC1'
sxaddpar,header,'CUNIT1','ARCSEC'
sxaddpar,header,'CRPIX2',halfax2
sxaddpar,header,'CRVAL2',ycrval
sxaddpar,header,'CDELT2',cdeltscale2
sxaddpar,header,'CTYPE2','ARCSEC2'
sxaddpar,header,'CUNIT2','ARCSEC'
sxaddpar,header,'COMMENT',outname2,'3D companion cube '


writefits,outname2,Diskdusttwodim,header

end


Pro Func3,f2,t,y,z,R0,z0,sini,cosi,td,tau,csttau,$
               xint,yint,y2,NA,rad,rotvel,nr,dveldisp,vadd,$
               vlow,vscale,nv,boun,bdrat,b0,bratio,buldisp,$
               p,k,l,trRdi,bounem

  bz0=double(b0*bratio)
  radialpos=double(SQRT(y*y+(z/cosi+TAN(t)*sini)*(z/cosi+TAN(t)*sini)))

  if ((t le (-!dpi/2.D0+1.D-4)) or (t ge (!dpi/2.D0-1.D-4)) $
      or (radialpos gt trRdi) or $
      ((radialpos ge bounem[0]) and (radialpos le bounem[1]))) then begin
     f2b=0.D0
     f2d=0.D0
     tt=0.D0
  endif else begin
     tt=double(TAN(t))
     if ((ABS(tau)) lt 1.e-4) then begin
        z1=double(z)

        fd=double(EXP(-(SQRT(y*y+(z/cosi+tt*sini)* $
                      (z/cosi+tt*sini)))/R0 - $
               abs(tt*cosi/z0)))
        fb=double(EXP(-(SQRT(y*y+(z/cosi+tt*sini)* $
                      (z/cosi+tt*sini)))/b0 - $
               abs(tt*cosi/bz0)))
        f2d=double(fd*(1+tt*tt))
        f2b=double(bdrat*fb*(1+tt*tt))

     endif else begin

        dustint,d1,t,NA,xint,yint,y2,boun,y,z,cosi,sini

        z1=double(z)
        dd2=double(ABS(tt*cosi/z0))
        dd3=double(sqrt(y*y+(z/cosi+tt*sini)*(z/cosi+tt*sini))/R0)

        db2=double(abs(tt*cosi/bz0))
        db3=double(sqrt(y*y+(z/cosi+tt*sini)*(z/cosi+tt*sini))/b0)
        arg1 = double( d1)
        arg2 =double( d1 + dd2 + dd3)
        arg3 =double( db2 + db3)

        f2d=double(exp(-arg2)*(1+tt*tt))
        f2b=double(exp(-arg1)*(bdrat*exp(-arg3)) * $
                   (1+tt*tt))
     endelse
  endelse
  y1=y
  tplint,rad,rotvel,nr,y1,res1
  noemer=double(y*y+(tt*sini+z/cosi)^2.D0)
  if (noemer eq 0.D0) then result=0.D0 else $
	   result=double(res1*sini*ABS(y)/SQRT(noemer))
  test=0.D0
  if (f2d lt 1.D-50) then f2d=0.D
  vadd=dblarr(nv)
  for i=0,nv-1 do begin
     vel=double(vlow+(i)*vscale)
     vadd[i]=double(EXP(-5.D-1*(((vel-result)/dveldisp)^2.D0)))
  endfor
  test=0.D0
  for i=0,nv-1 do begin
     test=test+vadd[i]
  endfor
  vadd[*]=vadd[*]*f2d/test
  hulp=dblarr(nv)
  for i=0,nv-1 do begin
     vel=double(vlow+(i)*vscale)
     hulp[i]=double(EXP(-5.D-1*(((vel-result)/buldisp)^2.)))
  endfor
  test=0.D0
  for i=0,nv-1 do begin
     test=test+hulp[i]
  endfor
  vadd[*]=double(vadd[*]+((hulp[*]*f2b)/test))
  f2=double(f2d)
end

Pro dustint,dustinte,t,NN,xint,yint,y2,boun,y,z,cosi,sini
  rad=double(SQRT(y*y+(z/cosi+TAN(t)*sini)*(z/cosi+TAN(t)*sini)))
  if (((rad ge boun[0]) and (rad le boun[1])) or $
      ((rad ge boun[2]) and (rad le boun[3]))) then begin
     tplint,xint,yint,NN,t,dustinte
     if (dustinte lt 0.) then dustinte=0.D0
  endif else dustinte=0.D0
end

Pro tplint,XA,YA,N,X,Y
  KLO=0.D0
  KHI=double(N-1.D0)
  WHILE (KHI-KLO GT 1.D0) DO BEGIN
     K=fix((KHI+KLO)/2.D0)
     IF(XA(K) GT X) THEN KHI=K ELSE KLO=K
  ENDWHILE
  Y=double(YA(KLO)+(X-XA(KLO))*(YA(KHI)-YA(KLO))/(XA(KHI)-XA(KLO)))
end

Pro dimake,NN,td,csttau,Rdi,y,z,sini,cosi,xint,yint,trRdi
xint=dblarr(NN)
yint=dblarr(NN)
for i=0.d0,NN-1 do begin
   Sum2=0.d0
   xst=double(-!dpi/2.D0+(i)*!dpi/NN)
   xend=double(xst+!dpi/NN)
   Gauleg,xst,xend,XX,W,16

   for k=0.d0,15.D0 do begin
      if (XX[k] le (-!dpi/2.D0+1.D-4)) or $
         (XX[k] ge (!dpi/2.D0-1.D-4)) or $
         (sqrt(y*y+(z/cosi+TAN(XX[k])*sini)*(z/cosi+TAN(XX[k])*sini)) gt trRdi) then begin
         Sum2=double(Sum2+0.D0)
      endif else begin
         tt=double(TAN(XX[k]))
         f2=double(EXP(-(SQRT(y*y+(z/cosi+tt*sini)*(z/cosi+tt*sini)))/Rdi-ABS(tt*cosi/td))*(1.D0+tt*tt))

         Sum2=double(Sum2+W[k]*f2)
     endelse
   endfor
   xint[i]=double(xst)
   yint[i]=double(Sum2*csttau)
endfor

for i=NN-2,0,-1 do begin
   yint[i]=double(yint[i]+yint[i+1])
endfor
for i=0,NN-2 do begin
   xint[i]=double((xint[i]+xint[i+1])/2.D0)
endfor

end

Pro Gauleg,X1,X2,X,W,N
  X=dblarr(N)
  W=dblarr(N)
  M=(N+1)/2
  XM=0.5D0*(X2+X1)
  XL=0.5D0*(X2-X1)
  EPS=3E-14
  For i=0.d0,m-1 do begin
     Z=double(COS(!dpi*((I+1.d0)-0.25D0)/(N+0.5D0)))
     fromhere:
     P1=1.D0
     P2=0.D0
     FOR J=0.d0,N-1.d0 Do begin
        P3=double(P2)
        P2=double(P1)
        P1=double(((2.D0*(J+1.d0)-1.d0)*Z*P2-(J)*P3)/(J+1.d0))
     ENDFOR
     PP=double(N*(Z*P1-P2)/(Z*Z-1.D0))
     Z1=double(Z)
     Z=double(Z1-P1/PP)
     IF double(ABS(Z-Z1)) GT EPS then GOTO,fromhere
     X[I]=double(XM-XL*Z)
     X[N-I-1]=double(XM+XL*Z)
     W[I]=double(2.D0*XL/((1.D0-Z*Z)*PP*PP))
     W[N-1-I]=W[I]
  endfor
end
