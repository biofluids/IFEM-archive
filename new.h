
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	integer nn_on,nn_al,mmon
	integer node_on(1),node_al(1)
	pointer (node_onptr,node_on),(node_alptr,node_al)
	common  nn_on,nn_al,mmon,node_onptr,node_alptr

	integer on(1),onfrom(1),al(1),alfrom(1)
	integer destin(1),source(1)
	pointer (onptr,on),(onfromptr,onfrom)
	pointer (alptr,al),(alfromptr,alfrom)
	pointer (destinptr,destin),(sourceptr,source)
	common  onptr,onfromptr,alptr,alfromptr,destinptr,sourceptr

	real*8  bufon(1),bufal(1)
	pointer (bufonptr,bufon),(bufalptr,bufal)
	common  bufonptr,bufalptr
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
