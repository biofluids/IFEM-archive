subroutine getnormal(normal,ie,ieface,ien,x)
	! tetrahedral elements

	use fluid_variables
      	implicit none


      	real* 8 normal(nsdpad)
      	real* 8 x(3,nn)
	integer rngface(neface,ne),ien(nen,ne)
	integer ie,ieface,i

      	real* 8 xnode(nsdpad,4)
      	integer inface,node
      	real* 8 x12,x13,x14,y12,y13,y14,z12,z13,z14
      	real* 8 nnorm
	integer irng, irnn, inl
	
      	do inface = 1,nnface
        	node =ien(mapping(ieface,inface,etype),ie)
		
		xnode(1:3,inface) = x(1:3,node)
	end do


	x12 = xnode(1,2) - xnode(1,1)
        x13 = xnode(1,3) - xnode(1,1)
        y12 = xnode(2,2) - xnode(2,1)
        y13 = xnode(2,3) - xnode(2,1)
        z12 = xnode(3,2) - xnode(3,1)
        z13 = xnode(3,3) - xnode(3,1)




      	normal(1) = z12*y13 - y12*z13
      	normal(2) = x12*z13 - z12*x13
      	normal(3) = y12*x13 - x12*y13


	do i =1,3

	if ((xnode(1,i)>0).and.(xnode(2,i)>0)) then
		if (normal(1)<0) then
			normal(1)=-normal(1)
		end if
		if (normal(2)<0) then
                       normal(2)=-normal(2)
              end if

	end if
	
	if ((xnode(1,i)>0).and.(xnode(2,i)<0)) then
               if (normal(1)<0) then
                     normal(1)=-normal(1)
             end if
		if (normal(2)>0) then
                       normal(2)=-normal(2)
              end if

        end if

	if ((xnode(1,i)<0).and.(xnode(2,i)>0)) then
               if (normal(1)>0) then
                      normal(1)=-normal(1)
             end if
		if (normal(2)<0) then
                       normal(2)=-normal(2)
              end if

       end if

	if ((xnode(1,i)<0).and.(xnode(2,i)<0)) then
               if (normal(1)>0) then
                      normal(1)=-normal(1)
             end if
		if (normal(2)>0) then
                       normal(2)=-normal(2)
              end if
       end if
        

	enddo


      	nnorm = sqrt(normal(1)*normal(1)+normal(2)*normal(2)+normal(3)*normal(3))
      	normal(1:nsd) = normal(1:nsd)/nnorm
      return
      end
