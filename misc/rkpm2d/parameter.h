c
c*********************************************************
c
c      parameter.h
c
c      Definition of limits of all working arrays
c
c      The parameter.h is for adaptive algorithm
c
c
c*********************************************************
c
       integer maxNumnp
       integer maxRF
       integer maxDof
                
       integer maxInt		! the max num of gauss integraion point in 1D
       integer maxIntElem	! the max num of gauss points in each 2D elem
       
c      integer max       
       integer maxGP
       
c      integer mnpeb
       integer maxEssbcX	! the max num of essbc along X or Y direction
       integer maxDispbcX
       integer maxVelbcX
       integer maxEssbcY	! the max num of essbc along X or Y direction
       integer maxDispbcY
       integer maxVelbcY
       integer maxEssbc
       integer maxTract
        
       integer mxndt
       integer mnsch
       
       integer maxNode		! the max num of nodes per element
       integer maxElem
       integer maxTF

       parameter (maxNode=4) 
c
       parameter (maxElem   = 5780)
       parameter (maxNumnp  = 6532)
       parameter (maxRF     = 2)
       parameter (maxDof    = maxNumnp + 3*maxRF)
       parameter (maxInt=2)
       parameter (maxIntElem=maxInt*maxInt)
c
       parameter (maxGP=maxElem*maxIntElem)
c       
       parameter (maxDispbcX = 85)
       parameter (maxDispbcY = 50)
       parameter (maxVelbcX  = 50)
       parameter (maxVelbcY  = 110)
       parameter (maxEssbcX  = maxDispbcX+maxVelbcX)
       parameter (maxEssbcY  = maxDispbcY+maxVelbcY)
c       parameter (maxEssbc   = max0(maxEssbcX,maxEssbcY))
       parameter (maxEssbc   = maxEssbcX+maxEssbcY)
       parameter (maxTract = 2)
       parameter (mxndt    = 2)
       parameter (mnsch    = 33)
       parameter (maxTF    = 4)
c
