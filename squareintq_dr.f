cc This software is being released under a modified FreeBSD license
cc (see COPYING in home directory). 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       
c       Quadratures and interpolation for smooth functions on a square
c
c       All quadratures are accurate to 34 digits
c       All weights are positive and inside the square
c
c       Input:
c
c       n - the degree of the quadrature (must not exceed 20)
c
c       Output:
c
c       fort.11 - quadrature nodes and weights
c       fort.12 - quadrature nodes in gnuplot compatible format
c       fort.14 - square as contructed in gnuplot compatible format
c       
c       in gnuplot:
c
c       plot 'fort.14' w l, 'fort.12'
c
c
c
        implicit real *8 (a-h,o-z)
        real *8 rnodes(2,1000),weights(1000)
        real *8 rints(10000),z(2),pols(10000)
        real *8 work(1000000)
c
c
        call prini(6,13)
c
c       SET ALL PARAMETERS
c
        PRINT *, 'ENTER mmax (interpolation, 1..20)'
        READ *, mmax
c
        call prinf('mmax=*',mmax,1)
c
        PRINT *, 'ENTER nmax (quadrature)'
        READ *, nmax
c
        call prinf('nmax=*',nmax,1)
c
c

c
c       ... retrieve the quadrature rule 
c
        call squareintq(mmax,rnodes,weights,numnodes)
c
        call prinf('nummodes=*',numnodes,1)
        call prin2('rnodes=*',rnodes,2*numnodes)       
        call prin2('weights=*',weights,numnodes)       
c
        d=0
        do i=1,numnodes
        d=d+weights(i)
        enddo
c
        call prin2('sum of weights=*',d,1)
c        
c       ... plot the square
c
        write(14,*) -1,-1
        write(14,*) +1,-1
        write(14,*) +1,+1
        write(14,*) -1,+1
        write(14,*) -1,-1
        write(14,*)
c
c
c       ... dump the nodes into a file
c       
        write(11,*) numnodes
        do i=1,numnodes
        write(11,1020) rnodes(1,i),rnodes(2,i),weights(i)
        write(12,*) rnodes(1,i),rnodes(2,i)
        enddo
 1010   format(i2)
 1020   format(4(e21.15,2x))
c

c
c       construct the matrix of values of the orthogonal polynomials
c       at the user-provided nodes        
c
        npols=(nmax+1)*(nmax+2)/2
c
c
        do 1200 j=1,npols
        rints(j)=0
 1200 continue
c
c
        do 2400 i=1,numnodes
c
        z(1)=rnodes(1,i)
        z(2)=rnodes(2,i)
c       
        call lege2eva(nmax,z,pols,work)
c
        do 2200 j=1,npols
c
        rints(j)=rints(j)+weights(i)*pols(j)

 2200 continue

 2400 continue
c
c
        call prin2('rints=*',rints,npols)
c
        area=4.0d0
c
        d=0
        d=(rints(1)-sqrt(area))**2
        do i=2,npols
        d=d+rints(i)**2
        enddo
c
        d=sqrt(d)/npols
c
        call prin2('error=*',d,1)
c
        ifinside = 1
        ifpositive = 1
        do i = 1,numnodes
        if( abs(rnodes(1,i)) .ge. 1.0d0 ) ifinside = 0
        if( abs(rnodes(2,i)) .ge. 1.0d0 ) ifinside = 0
        if( weights(i) .le. 0.0d0 ) ifpositive = 0
        enddo
        call prinf('ifinside=*',ifinside,1)
        call prinf('ifpositive=*',ifpositive,1)
        
        stop
        end
c
c
c
