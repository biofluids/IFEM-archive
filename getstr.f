	subroutine GETSTR (str)
	character str

	int istat,myid,numproc,ierr;
	
	ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid)
	ierr = MPI_Comm_size(MPI_COMM_WORLD, &numproc)

	if (myid .eq.0)

	ierr = MPI_Bcast(cstr,ilen,MPI_CHAR,0,MPI_COMM_WORLD)
	ierr = MPI_Barrier(MPI_COMM_WORLD)
	return	
