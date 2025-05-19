pro sfit_rectify_v4,check=check,red=red,blue=blue,chk=chk,clobber=clobber,debug=debug

forward_function mapfitw
  
;;v4, july 2015, dramatically redoes weights in surface fitting section

  ;;march 2015
;;trying to simplify a little; i found a way to smooth out the map fit
;;part, but there's a lot of kludge up front removing orders,
;;interpolating over them. I'm wondering if that's causing issues.


;;;new bit trying out: do a simple fit to some of the orders, then
;;;interpolate for the broad line winged orders. 


;this program does a couple of things- it fits the surface of the extracted orders 
;to try to remove the blaze function. it also then fits the order/blaze to try to mitigate
;edge/polynomial effects. Both are written out to file.
;;;the new additions in 2014 include: 
;adding fake orders to make the surface fitting a little easier (shallower gradient
;arbitrarity- works for red or blue data, any binning.

;;red keyword will start on the first red spectrum.

;window 0- shows surface of extracted orders
;window 2- shows surface of extracted orders (white) and the smoothed
;          order surface (orange) and the surface
;          fit to the smoothed order surface (green)
;window 1- res_sf (fit surface to smoothed surface) in blue, bfunc in
;          green (top panel). orders divided by each of those,
;          respectively (bottom panel)
;window 3- order joining code. white is the big array, orange the next
;          order being appended; top in green is the new big array.

   
if not windowavailable(0) then window,0,xsize=1000,ysize=800
if not windowavailable(1) then window,1,xsize=2000,ysize=800
if not windowavailable(2) then window,2,xsize=800,ysize=800
if not windowavailable(3) then window,3,xsize=800,ysize=800


loadct,39
;50 is blue, 150 is green, 200 is orangey yellow, 250 is red

;find the files
files=file_search('Extract/Obj_m*')

;;;;IMPORTANT: move the fake traces to an 'old' folder so you
;;;;don't fit them, too.
if keyword_set(red) then begin
spots=strpos(files,'_mr')
start=min(where(spots ge 0))
endif else begin
start=0
endelse

if keyword_set(blue) then begin
spots=where(strpos(files,'_mb') gt 0)
endpt=max(where(spots ge 0))
endif else begin
endpt=n_elements(files)-1
endelse

start=0 ;51
endpt=start ;53
for f=start,endpt do begin

filenamein=files[f]
;read it in
star=mrdfits(filenamein,1)
starname=star[0].field
filenum=strmid(star[0].img_fil,strpos(star[0].img_fil,'_')+2,5)
fileout='FMfiles/'+starname+'_'+filenum+'_surface.fits.gz'

flag=file_test(fileout)
if flag eq 1 and not keyword_set(clobber) then goto,nextfile
print,'Loop ',f,', object: ',files[f]

;;yikes, don't do super low SNR objects..
if mean(star.box_fx) lt 100 then print,'Low counts- skipping this file!'
if mean(star.box_fx) lt 100 then goto,nextfile

;initializations..
sz=size(star)
nrows=star[0].nrow
nplaces=n_elements(star[10].box_fx)
norders=sz[1]
maxes=make_array(norders)
for p=0,norders-1 do maxes[p]=max(star[p].box_fx)
norders=n_elements(where(maxes gt 0))

star2=star

meds=make_array(sz[1])
for l=0,sz[1]-1 do meds[l]=median(star[l].box_wv[0:nrows-1])

;;ok, quick ad-hoc fix for bad blue orders:
if max(star.box_wv) lt 7000 then begin
;;hg and hd orders really really bad.
xarr=findgen(nplaces)

indx=closest(meds,4955)
x1indx=closest(star[indx].box_wv,4930) & x2indx=closest(star[indx].box_wv,4915)
spots=where(xarr ge x1indx and xarr le x2indx,complement=nots)
star2[indx].box_fx[spots]=interpol(star[indx].box_fx[nots],xarr[nots],xarr[spots])

indx3=closest(meds,3950)
x1indx3=closest(star[indx3].box_wv,3980) & x2indx3=closest(star[indx3].box_wv,3955)
spots=where(xarr ge x1indx3 and xarr le x2indx3,complement=nots)
star2[indx3].box_fx[spots]=interpol(star[indx3].box_fx[nots],xarr[nots],xarr[spots])
x3indx3=closest(star[indx3].box_wv,3938) & x4indx3=closest(star[indx3].box_wv,3930)
spots=where(xarr ge x3indx3 and xarr le x4indx3,complement=nots)
star2[indx3].box_fx[spots]=interpol(star[indx3].box_fx[nots],xarr[nots],xarr[spots])

indx4=closest(meds,5018)
x1indx4=closest(star[indx4].box_wv,5022) & x2indx4=closest(star[indx4].box_wv,5013)
spots=where(xarr ge x1indx4 and xarr le x2indx4,complement=nots)
star2[indx4].box_fx[spots]=interpol(star[indx4].box_fx[nots],xarr[nots],xarr[spots])

indx5=closest(meds,4810)
x1indx5=closest(star[indx5].box_wv,4863)+10 & x2indx5=closest(star[indx5].box_wv,4850)
spots=where(xarr ge x1indx5 and xarr le x2indx5,complement=nots)
star2[indx5].box_fx[spots]=interpol(star[indx5].box_fx[nots],xarr[nots],xarr[spots])

indx6=closest(meds,3915)
x1indx6=closest(star[indx6].box_wv,3901.5) & x2indx6=closest(star[indx6].box_wv,3880)
spots=where(xarr ge x1indx6 and xarr le x2indx6,complement=nots)
star2[indx6].box_fx[spots]=interpol(star[indx6].box_fx[nots],xarr[nots],xarr[spots])
x3indx6=closest(star[indx6].box_wv,3940) & x4indx6=closest(star[indx6].box_wv,3925)
spots=where(xarr ge x3indx6 and xarr le x4indx6,complement=nots)
star2[indx6].box_fx[spots]=interpol(star[indx6].box_fx[nots],xarr[nots],xarr[spots])

indx7=closest(meds,4100)
x1indx7=closest(star[indx7].box_wv,3985) & x2indx7=closest(star[indx7].box_wv,3965)
spots=where(xarr ge x1indx7 and xarr le x2indx7,complement=nots)
star2[indx7].box_fx[spots]=interpol(star[indx7].box_fx[nots],xarr[nots],xarr[spots])

indx8=closest(meds,3870)
x1indx8=closest(star[indx8].box_wv,3900) & x2indx8=closest(star[indx8].box_wv,3875)
spots=where(xarr ge x1indx8 and xarr le x2indx8,complement=nots)
star2[indx8].box_fx[spots]=interpol(star[indx8].box_fx[nots],xarr[nots],xarr[spots])

indx9=closest(meds,4250)
x1indx9=closest(star[indx9].box_wv,4238) & x2indx9=closest(star[indx9].box_wv,4228)
spots=where(xarr ge x1indx9 and xarr le x2indx9,complement=nots)
star2[indx9].box_fx[spots]=interpol(star[indx9].box_fx[nots],xarr[nots],xarr[spots])

indx10=closest(meds,4570)
x1indx10=closest(star[indx10],4586) & x2indx10=closest(star[indx10],4579)
spots=where(xarr ge x1indx10 and xarr le x2indx10,complement=nots)
star2[indx10].box_fx[spots]=interpol(star[indx10].box_fx[nots],xarr[nots],xarr[spots])
x3indx10=closest(star[indx10].box_wv,4562) & x4indx10=closest(star[indx10].box_wv,4546)
spots=where(xarr ge x3indx10 and xarr le x4indx10,complement=nots)
star2[indx10].box_fx[spots]=interpol(star[indx10].box_fx[nots],xarr[nots],xarr[spots])

;interp over full orders
indx0=closest(meds,4860)
spot=closest(star[indx0].box_wv,4885)
for p=spot,nrows-1 do star2[indx0].box_fx[p]=interpol([star2[indx0-1].box_fx[p],star2[indx0+1].box_fx[p]],[indx0-1,indx0+1],indx0)
indx1=closest(meds,4340) 
for p=0,nrows-1 do star2[indx1].box_fx[p]=interpol([star[indx1-1].box_fx[p],star[indx1+1].box_fx[p]],[indx1-1,indx1+1],indx1)
indx2=closest(meds,4100)
for p=0,nrows-1 do star2[indx2].box_fx[p]=interpol([star[indx2-1].box_fx[p],star[indx2+1].box_fx[p]],[indx2-1,indx2+1],indx2)

;;orders 21-25 are all quite bad. would it even work to interpolate over that many orders? 
indices=[21,22,23,24,25]
others=[17,26,27,28]
for p=0,nrows-1 do star2[indices].box_fx[p]=interpol(star2[others].box_fx[p],others,indices)

;last ditch attempt to fix order 19...
temp=make_array(nrows)
ind1=[18,20]
for p=0,nrows-1 do temp[p]=mean(star2[ind1].box_fx[p])
spot=closest(star2[19].box_wv,3950)
star2[19].box_fx[0:spot]=temp[0:spot]
;same thing with edge of order next to Hbeta.. get rid of that damned
;artifact once and for all
temp=make_array(nrows)
ind1=[2,4]
for p=0,nrows-1 do temp[p]=mean(star2[ind1].box_fx[p])
spot=closest(star2[3].box_wv,4830)
star2[3].box_fx[0:spot]=temp[0:spot]
;what about Hgamma? it has two adjacent orders to fix up.. not sure i
;can get that.


endif


;;;mmk, doing similar stuff on the red side. 
;;things don't really go too far to hell until 7600A. the last
;;3 orders are really bad..
if max(star.box_wv) gt 7000 then begin
;;some O lines..
indx4=closest(meds,6350)
x1indx4=closest(star[indx4].box_wv,6355) & x2indx4=closest(star[indx4].box_wv,6340)
xarr=findgen(n_elements(star[10].box_fx))
spots=where(xarr ge x1indx4 and xarr le x2indx4,complement=nots)
for p=x1indx4,x2indx4 do star2[indx4].box_fx[x1indx4:x2indx4]=interpol(star[indx4].box_fx[nots],xarr[nots],xarr[spots])

x3indx4=closest(star[indx4].box_wv,6377) & x4indx4=closest(star[indx4].box_wv,6360)
spots=where(xarr ge x3indx4 and xarr le x4indx4,complement=nots)
for p=x3indx4,x4indx4 do star2[indx4].box_fx[x3indx4:x4indx4]=interpol(star[indx4].box_fx[nots],xarr[nots],xarr[spots])

indx00=closest(meds,7800)
x1indx00=closest(star[indx00].box_wv,7790) & x2indx00=closest(star[indx00].box_wv,7765)
spots=where(xarr ge x1indx00 and xarr le x2indx00,complement=nots)
for p=x1indx00,x2indx00 do star2[indx00].box_fx[x1indx00:x2indx00]=interpol(star[indx00].box_fx[nots],xarr[nots],xarr[spots])
x3indx00=closest(star[indx00].box_wv,7722) & x4indx00=closest(star[indx00].box_wv,7708)
spots=where(xarr ge x3indx00 and xarr le x4indx00,complement=nots)
for p=x3indx00,x4indx00 do star2[indx00].box_fx[x3indx00:x4indx00]=interpol(star[indx00].box_fx[nots],xarr[nots],xarr[spots])


indx0=closest(meds,7600)
indx1=closest(meds,6620) ;this is halpha
indx2=indx1-1 ;ha is on the edge of this order, too.
indx3=closest(meds,6850)
for p=0,nrows-1 do star2[indx0].box_fx[p]=mean([star[indx0-1].box_fx[p],star[indx0+1].box_fx[p]])
for p=nrows/2,nrows-1 do star2[indx1].box_fx[p]=mean([star[indx1-2].box_fx[p],star[indx1+1].box_fx[p]])
for p=0,nrows-1 do star2[indx2].box_fx[p]=mean([star[indx2-1].box_fx[p],star[indx2+2].box_fx[p]])
for p=0,nrows-1 do star2[indx3].box_fx[p]=mean([star[indx3-1].box_fx[p],star[indx3+1].box_fx[p]])

;ok, if it's the red side, use the trace instead of the blaze
;for the last few orders- get those wiggles out from the bg sub
;alright, now, for the red side- do the last few orders (28-33) a
;different way.
fakeflat=mrdfits('Extract/old/Obj_mr0301.fits.gz',1)
for p=28,33 do begin
if fakeflat[0].nrow ne star[0].nrow then begin
tmp=fakeflat[p].box_fx
fakeflat[p].box_fx=0.0
fakeflat[p].box_fx[0:star[0].nrow-1]=congrid(tmp[0:fakeflat[0].nrow-1],star[0].nrow)
endif

;divide obj by trace
tr=smooth(fakeflat[p].box_fx[0:nrows-1],40)
repl=where(tr lt 0)
if repl[0] ge 0 then tr[repl]=0.0
temp=star[p].box_fx[0:nrows-1]/tr
;nan trap
junk=where(finite(temp) eq 0)
if junk[0] ge 0 then temp[junk]=0.0
xarr=star[p].box_wv[0:nrows-1]
;fit with some kind of robust, low order poly. cut out order
;edges, or you'll get that 'this is too weird!' error.
;;res=robust_poly_fit(xarr[60:4000],temp[60:4000],3,out)
;;ybl=res[0]+res[1]*xarr+res[2]*xarr^2.0+res[3]*xarr^3.0 
;;;;upped the order of the fit, since these orders were consistently
;;;;high before
res=robust_poly_fit(xarr[100:nrows-100],temp[100:nrows-100],4,out)
ybl=res[0]+res[1]*xarr+res[2]*xarr^2.0+res[3]*xarr^3.0 +res[4]*xarr^4.0
;write trace/poly as the final blaze function
star2[p].box_fx[0:nrows-1]=smooth(fakeflat[p].box_fx[0:nrows-1],40)*ybl
;;also need to use this in lieu of whatever the fitting stuff puts in
;;here.
endfor
endif


;;;;how about removing the edge pixels?
edgeclip=50 ;pix
;;this cuts a few steps out..
surf00=make_array([norders,nrows-2.0*edgeclip])
surf0=make_array([norders,nrows-2.0*edgeclip])
for p=0,norders-1 do surf0[p,*]=star2[p].box_fx[edgeclip:nrows-1-edgeclip]
for p=0,norders-1 do surf00[p,*]=star[p].box_fx[edgeclip:nrows-1-edgeclip]
spots=where(surf0 lt 0)
if spots[0] ge 0 then surf0[spots]=0.0  ;;apparently this is a huge problem on the red side.
nrows_orig=nrows
nrows=nrows-2.0*edgeclip
surf=congrid(surf00,norders*2-1,nrows,/cubic);/interp)

smoot=25 ;smoothing length for initial median smoothing of orders
smoo=filter_image(surf,fwhm_g=[1,smoot],/all) ;(surf,smoo=3,/med)
smoo[0,*]=smooth(surf[0,*],50)
if max(star.box_wv) ge 7000 then begin
for p=0,5 do smoo[p*2,*]=smooth(surf0[p,*],50)
endif

wset,0
image_cont_uv,surf,/noc
wait,0.3
image_cont_uv,smoo,/noc

;look for old surfaces, use them if they exist.. i think in many
;cases, they're better :( but the jo_simple back then was worse.
oldfil=strmid(fileout,0,strpos(fileout,'/')+1)+'old/'+strmid(fileout,strpos(fileout,'/')+1,strpos(fileout,'.gz'))
flag=file_test(oldfil)
rejoldsurflag=0
if keyword_set(clobber) then rejoldsurflag=1


rej_old_surf:
if rejoldsurflag eq 1 then flag=0

if flag eq 1 then begin
oldsurf=mrdfits(oldfil,0)
finalblaze=make_array([sz[1],nrows_orig])
for s=0,sz[1]-1 do finalblaze[s,*]=oldsurf[s,0:nrows_orig-1] ;[*,s]
res_sf=congrid(finalblaze,norders*2-1,nrows_orig,/inter)
bfunc=finalblaze
endif else begin
print,'surface fitting 1'

;oh dear. weights on or off produces the same result :/
surf_patch,smoo,res_sf,surf,star,wgts_out,/weights,check=check

;;;;continue here! it looks like surf_patch is doing a good job.
s4=size(res_sf)
bfunc=make_array([norders,nplaces])
;don't just make one with the original orders, do the full
;stretched out thing
;bfuncfull=res_sf*0.0
;;spline to each order to remove hitches between fit boxes
for i=0,s4[1]-1,2 do begin
x=dindgen(nrows)+edgeclip
xfull=dindgen(nrows_orig)
;;;do NOT fit the whole thing- 0's at the end will mess it up badly.
;x1=100 & x2=nrows-100
x1=10 & x2=max(where(res_sf[i,*] gt 0))-10
if nrows gt 2000 then factr=nrows/16. else factr=nrows/12.
bkpts=[congrid(x[x1:x2],s4[2]/factr),x[x2]-10]
bkpts[0]=bkpts[0]+10
temp=0.0 ;;needs to be re-initialized each time?
if nrows_orig gt 3000 then bsp_ordr=5 else bsp_ordr=5
temp=bspline_iterfit(x[x1:x2],res_sf[i,x1:x2],nord=bsp_ordr,fullbkpt=bkpts,upper=2,lower=2)
out=bspline_valu(xfull,temp)
;;;should put in a failsafe, just in case the spline doesn't work.. 
;if nrows lt 2000 then begin
;test=robust_poly_fit(x[x1:x2],res_sf[i,x1:x2],5)
;out=xfull*test[1]+(xfull^2.)*test[2]+(xfull^3.)*test[3]+(xfull^4.)*test[4]+(xfull^5.)*test[5]
;endif
bfunc[i/2,0:nrows_orig-1]=out
;bfuncfull[i,0:nrows_orig-1]=out
if keyword_set(chk) then begin
wset,2
erase

plot,x,res_sf[i,*]
oplot,x[x1:x2],res_sf[i,x1:x2],col=220
vline,bkpts,li=3
oplot,xfull,out,col=120

endif

endfor

wset,2
surface,surf,zrange=[0,1.1*max(res_sf)],zstyle=1,noclip=0,xrange=[0,s4[1]-1],yrange=[0,s4[2]]
surface,smoo,/noerase,col=220,zrange=[0,1.1*max(res_sf)],zstyle=1,noclip=0,xrange=[0,s4[1]-1],yrange=[0,s4[2]]
factor,s4[2],rbnf,ns,/quiet
xarr=congrid(findgen(s4[2]),s4[2]/max(rbnf))
surface,rebin(res_sf,s4[1],s4[2]/max(rbnf)),findgen(s4[1]),xarr,/noerase,col=140,zrange=[0,1.1*max(res_sf)],zstyle=1,xrange=[0,s4[1]-1],yrange=[0,s4[2]]

finalblaze=bfunc
for p=1,norders-1 do begin
finalblaze[p,*]=bfunc[p,*]

;if nrows le 2000 then begin
;;for lower binnings, going to fit to figure out how the orders should
;;be shifted/scaled.
;xarr=findgen(s4[1])
;res=robust_poly_fit(xarr,surf[*,s4[2]/2],5)
;tmp=res[0]+xarr*res[1]+xarr^2.*res[2]+xarr^3.*res[3]+xarr^4.*res[4]+xarr^5.*res[5]
;finalblaze[p,*]=finalblaze[p,*]*(tmp[p*2]/finalblaze[p,s4[2]/2])
;endif


;;;;add in a little shifting code to triage issues with the surface
 ;;;;being systematically too low/high in some places. find the median,
;;;;and if it's >=1.2, shift. otherwise, probably leave it alone..
tmp=star[p].box_fx[edgeclip:nrows_orig-1-edgeclip]/finalblaze[p,edgeclip:nrows_orig-1-edgeclip]
print,p,'; ',median(tmp)
if median(tmp) ge 1.02 or median(tmp) le 0.977 then begin
print,'fixing normalization a little bit for order '
if max(star.box_wv) lt 6000 then begin
if p eq 18 then finalblaze[p,*]=finalblaze[p,*]*median(tmp)
endif else begin
if p le 16 then finalblaze[p,*]=finalblaze[p,*]*median(tmp)
endelse

endif

;red side, replace all the fit stuff with the extracted trace as
;scaled earlier.
if max(star.box_wv) gt 7000 then begin
if p ge 28 then begin
finalblaze[p,*]=star2[p].box_fx;[0:nrows-1]
endif
endif

endfor
endelse


;;;;;;;quick test. try smoothing perpendicular to dispersion direction
for w=0,(size(finalblaze))(2)-1 do finalblaze[*,w]=smooth(finalblaze[*,w],3)
;;;;;

;preview
wset,1
erase
!p.multi=[0,1,2]

plot,[0,0],xrange=[min2(star.box_wv),max(star.box_wv)],xstyle=1,yrange=[0,1.0*max(res_sf)],ystyle=1  
for i=0,norders-1 do oplot,star[i].box_wv[edgeclip:nrows_orig-1-edgeclip],star[i].box_fx[edgeclip:nrows_orig-1-edgeclip]
for i=0,norders-1 do oplot,star[i].box_wv[edgeclip:nrows_orig-1-edgeclip],res_sf[i*2,*],col=80
for i=0,norders-1 do oplot,star[i].box_wv[edgeclip:nrows_orig-1-edgeclip],bfunc[i,edgeclip:nrows_orig-1-edgeclip],col=180
for i=0,norders-1 do oplot,star[i].box_wv[edgeclip:nrows_orig-1-edgeclip],finalblaze[i,edgeclip:nrows_orig-1-edgeclip],col=200
for i=0,norders-1 do xyouts,median(star[i].box_wv),1e4,strtrim(string(i),2),chars=1.1
legend_leg,box=0,['Object','Surface fit orders','Spline fit orders','Finalblaze'],col=[255,80,180,200],psym=[4,4,4,4]

!p.multi=[1,1,2]
plot,[0,0],xrange=[min2(star.box_wv),max(star.box_wv)],xstyle=1,yrange=[0,3],ystyle=1  
for i=0,norders-1 do oplot,star[i].box_wv[edgeclip:nrows_orig-1-edgeclip],1.0+(star[i].box_fx[edgeclip:nrows_orig-1-edgeclip]/res_sf[i*2,*]),col=80
for i=0,norders-1 do oplot,star[i].box_wv[edgeclip:nrows_orig-1-edgeclip],(star[i].box_fx[edgeclip:nrows_orig-1-edgeclip]/finalblaze[i,edgeclip:nrows_orig-1-edgeclip]),col=200
!p.multi=[0,1,0]


if rejoldsurflag eq 0 then begin
read,yesno,prompt='Keep old surface fit? 0=yes, 1=no, redo: '
if yesno eq 1 then rejoldsurflag=1
if yesno eq 1 then goto,rej_old_surf
endif

if keyword_set(debug) then stop
;write the fit surface to file
mwrfits,finalblaze,fileout,/create

if finite(max(finalblaze)) eq 0 then print,'Blaze fit issue- stopping.'
if finite(max(finalblaze)) eq 0 then stop
;add em up
jo_simple,star,finalblaze,files[f],chk=chk

nextfile:
endfor

print,'All done! ',strtrim(string(f),2),' files processed.'
end



pro surf_patch,surf_in,fit_out,surf_orig,objstruct,weights_out,weights=weights,check=check,pass2=pass2

;;ok, going to write the 'moving patch' code to scoot a box along the
;;surface, fitting smaller bits of it at a time, and then medianing them.
s=size(surf_in)
nx=s[1]
ny=s[2]
;ny=28 ;;try not including the last few orders in the fitting, since i overwrite those later anyway

;;;consider setting the widths based on how broad the lines are/are
;not... 

;box widths in x and y- x is #orders direction, y is #pix
widthx=nx/15         ;4 orders
;;so, i want to divide it up into about 16 bins in pixel space-
;;that's 250 pixels wide for 4096, and 79 pixels wide for 1365 binning
widthy=fix(ny/16.)
print,'Moving surface size, Orders x rows: ',widthx,'x',widthy

;number of steps in each direction; want overlap, so div/2
dx=nx/(widthx/2)
dy=ny/(widthy/2)

xoff=indgen(dx)*(widthx/2)
yoff=indgen(dy)*(widthy/2)

;make structure to store the results
blank=make_array([nx,ny],/double)
temp={seg_fit:blank}
fits=replicate(temp,dx*(dy-1))

s2=size(blank)

big_fit=blank

;;;ok! really, really need a careful treatment of the weights. I could
;;;make the weights an enormous array.. with the weights good if the
;;;obj/smooth is within a stdev or something, and then 0 if
;;;it's an em line.
weights=blank
sigs=blank

if keyword_set(weights) then begin
print,'entering weight calculation loop'
    for p=0,s[1]-1 do begin
tmp=robust_poly_fit(findgen((size(surf_orig[p,*]))(2)),surf_orig[p,*],6,yfit)
sigs[p,*]=abs(1.0-surf_orig[p,*]/yfit)
weights[p,*]=1.0/sigs[p,*]

;if max(objstruct.box_wv) lt 7000 then begin
;;;want to give less weight to hbeta adjacent order- the wings extend to that order, but are being fit out
;weights[5:7,0:ny/2]=1e-3
;endif
;stop

           if keyword_set(check) then begin
              wset,2
              erase
              multiplot,[1,2]
              xval=objstruct[p/2].box_wv[0:s[2]-1]
              plot,xval,surf_orig[p,*] & oplot,objstruct[p/2].box_wv[0:s[2]-1],surf_in[p,*],col=220
              multiplot
;              plot,xval,tmp,yrange=[-0.5,2],ps=3,ystyle=1
;              oplot,xval,yfit,col=220
;              oplot,xval,tmp2,col=120
;              hline,nsig
;              oplot,xval[keep],tmp[keep],ps=3,col=80
;              if toss[0] ge 0 then oplot,xval[toss],tmp[toss],ps=3,col=250
;              print,der_snr(tmp),' ',p
plot,xval, surf_orig[p,*], ps=3,/ynoz
oploterror,surf_orig[p,*],sigs[p,*]

       ;       wait,0.1   
       hak
              multiplot,/reset
           endif
    endfor                  ;end loop on orders
endif ;end the weights if loop

tmp=finite(weights)
weights[where(tmp eq 0)]=0.0

weights_out=weights

counter=0
;;problem- last bits may run off the chip
for i=0,dx-1 do begin
for j=0,dy-2 do begin

   if xoff[i]+widthx gt s[1]-1 then begin
      x1=xoff[i]
      x2=s[1]-1
   endif else begin
      x1=xoff[i]
      x2=xoff[i]+widthx
   endelse

   if yoff[j]+widthy gt s2[2]-1 then begin
      y1=yoff[j]
      y2=s2[2]-1
   endif else begin
      y1=yoff[j]
      y2=yoff[j]+widthy
   endelse

;;ok! this works pretty well.
fit_seg=surf_in[x1:x2,y1:y2]
;print,x1,x2,y1,y2

if max(objstruct.box_wv) ge 7000 then weights[0,*]=1e-4 
;if max(objstruct.box_wv) ge 7000 then maporder=4 else
maporder=3 ;blue side, 4 for red. doesn't seem to change much.
;maporder=4

;;alright, so how does mapfitw use the weights..did i alter it to take
;;the weights but use them as errors? that would be ok. yes-
;;that's what I did. it's a subroutine called regress that
;;doesn't take weights anymore, so i put in 1/weights within mapfitw_aa

;if keyword_set(weights) then fit=mapfitw_aa(fit_seg,maporder,weights[x1:x2,y1:y2]) else fit=rob_mapfit(fit_seg,3)
if keyword_set(weights) then fit=mapfitw(fit_seg,maporder,sigs[x1:x2,y1:y2]) else fit=rob_mapfit(fit_seg,3)
fits[counter].seg_fit[x1:x2,y1:y2]=fit

;;;what if i intentionally loopped off the ends a little bit, to get
;;;rid of edge effects? Just an idea...
counter=counter+1
;print,'which segment # was being fit? ',counter
;if counter eq 40 or counter eq 10 or counter eq 20 or counter eq 30 then hak
endfor
endfor

;;;;;ANOTHER LOOP HERE--- BIGGER SURFACES, WEIGHT THEM IN COMBINE;;;;;;
;;fits=replicate(temp,dx*(dy-1))
fits2=replicate(temp,dx*(dy-1)+6)
fits2[0:counter-1].seg_fit=fits[0:counter-1].seg_fit

segxs=[0,(s[1]-1)/4,(s[1]-1)/2,(s[1]-1)*3/4,s[1]-1]
segys=[0,(s[2]-1)/4,(s[2]-1)/2,(s[2]-1)*3/4,s[2]-1]
for i=0,2 do begin
for j=0,1 do begin
;doing this a little bit klugely, but that's ok
seg1=surf_in[segxs[i]:segxs[i+1],segys[j]:segys[j+1]]
;;;edited to use the weights here too..
if keyword_set(weights) then begin
wts2=weights[segxs[i]:segxs[i+1],segys[j]:segys[j+1]]
fits2[counter+i].seg_fit[segxs[i]:segxs[i+1],segys[j]:segys[j+1]]=mapfitw_aa(seg1,3,wts2)
endif else begin
fits2[counter+i].seg_fit[segxs[i]:segxs[i+1],segys[j]:segys[j+1]]=rob_mapfit(seg1,4)
endelse
endfor
endfor


s3=size(blank)
npix=s3[4]

for p=0l,npix-1 do begin
;now, for each pixel position, take the median of all the overlapping fits.
   vals=fits2.seg_fit[p]
   spots=where(vals gt 0)
      if n_elements(spots) eq 1 then begin
         if spots eq -1 then medfitval=0 else begin
            medfitval=vals[spots]
            big_fit[p]=medfitval
         endelse
      endif else begin
         medfitval=median(vals[spots])
         resistant_mean,vals[spots],2.0,meanfitval
;big_fit[p]=medfitval
         big_fit[p]=meanfitval
      endelse
endfor


;did a division, so now I have some NaNs to take care of.
find_bad=finite(big_fit)
spots=where(find_bad eq 0)
s3=size(spots)
if s3[0] ne 0 then big_fit[spots]=0.0

;;;;;;alright, here's some code if you actually want to plot up
;;;;;;the surface fit segments across each order
if keyword_set(check) then begin
wset,3
erase
!p.multi=[0,1,0]

;;ok, this makes arrays of order number to match up with the segment fits
indx_arr=big_fit
szz=size(big_fit)
npixx=szz[2]
for p=0,szz[1]-1 do indx_arr[p,*]=replicate(p,npixx)
rindx=reform(indx_arr,szz[4])

for q=0,szz[1]-1 do begin
   print,'Order = ',q
   vals=0
   pixs=0
erase
   plot,surf_orig[q,*],yrange=[0,max(fits2.seg_fit[where(rindx eq q)])]
;;;what in the world is this?
   for i=0,6+counter-1 do begin
      values=fits2[i].seg_fit[where(rindx eq q)]
      keep=where(values gt 0)
      if keep[0] ne -1 then begin
         values=values[keep]
         oplot,keep,values,col=80
         vals=[vals,values]
         pixs=[pixs,keep]
      endif
   endfor
oplot,big_fit[q,*],col=220
wait,0.05
;hak
endfor

endif

;ok, now, return big_fit
fit_out=big_fit

end


