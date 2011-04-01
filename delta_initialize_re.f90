subroutine delta_initialize_re(x,y,ny,nsd,adist,snvolume,maxconn,rksh,ninf,inf)
!===========================================
real(8) x(nsd) ! point to be interpolated (RKPM center)
real(8) y(nsd,ny) ! possible surrounding nodes (back round mesh)
integer ny ! number of surrounding nodes
real(8) adist(nsd,ny)
integer nsd ! nsd
real(8) snvolume(ny)
integer maxconn ! max nodes in the support region
real(8) rksh(maxconn) ! rkpm coefficient for point x
integer ninf
integer inf(maxconn) ! node index in support region
!===========================================
real(8) y_tmp(nsd)
real(8) a(nsd)

rksh(:)=0.0d0

ninf=0
inf(:)=0

call getinf_re(inf,ninf,x,y,adist,ny,nsd,maxconn)

if (nsd == 2 ) then
		call correct2d(b,bd,x,y,adist,snvolume,1,inf,ninf,maxconn)
	elseif (nsd ==3 ) then
		call correct3d(b,bd,x,y,adist,snvolume,1,inf,ninf,maxconn)
end if


do n=1,ninf
	nnum=inf(n)
	do isd=1,nsd
		y_tmp(isd)=y(isd,nnum)
		a(isd)=adist(isd,nnum)
	end do
		if (nsd ==2) then
		call RKPMshape2d(shp,b,bd,x,y_tmp,a,snvolume(nnum))
		elseif (nsd ==3) then
		call RKPMshape2d(shp,b,bd,x,y_tmp,a,snvolume(nnum))
		end if
	rksh(n)=shp
end do
return

end subroutine
