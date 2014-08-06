cc Copyright (C) 2009: Zydrunas Gimbutas and Hong Xiao
cc Contact: Zydrunas Gimbutas <gimbutas@cims.nyu.edu>
cc          Hong Xiao <hxiao@ucdavis.edu>
cc 
cc This program is free software; you can redistribute it and/or modify 
cc it under the terms of the GNU General Public License as published by 
cc the Free Software Foundation; either version 2 of the License, or 
cc (at your option) any later version.  This program is distributed in 
cc the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
cc even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
cc PARTICULAR PURPOSE.  See the GNU General Public License for more 
cc details. You should have received a copy of the GNU General Public 
cc License along with this program; 
cc if not, see <http://www.gnu.org/licenses/>.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c       Orthogonal polynomials on the standard triangle
c       
c       (0,2/sqrt(3)), (-1,-1/sqrt(3)), (1,-1/sqrt(3))
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       This file contains a set of subroutines for the handling 
c       of the evaluation of orthogonal polynomials on the standard 
c       triangle.  It contains 2 user-callable subroutines.
c
c   ortho2eva - evaluates at the user-supplied point z
c       a collection of polynomials (of x,y) orthogonal on the
c       standard triangle
c
c   ortho2eva3 - evaluates at the user-supplied point z a collection
c       of polynomials (of x,y) orthogonal on the standard triangle,
c       together with their derivatives with respect to x and y.

c       
c
      subroutine ortho2eva(mmax,z,pols,w)
c
c     This subroutine evaluates at the user-supplied point z
c     a collection of polynomials (of x,y) orthogonal on the
c     standard triangle.
c
c     The "standard" triangle is the triangle with the
c     vertices
c     
c     (0,2/sqrt(3)), (-1,-1/sqrt(3)), (1,-1/sqrt(3)).       (1)
c     
c     The polynomials evaluated by this subroutine are all 
c     orthogonal polynomials up to order mmax, arranged in the 
c     increasing order. 
c
c     This subroutine is based on the Koornwinder's representation
c     of the orthogonal polynomials on the right triangle
c
c     (-1,-1), (-1,1), (1,-1)                                       (2)
c
c     given by
c
c     K_mn(x,y) = P_m ((2*x+1+y)/(1-y)) ((1-y)/2)^m P_n^{2m+1,0} (y)
c
c     where P_m are the Legendre polynomials or order m
c     and P_n^{2m+1,0} are the Jacobi polynomials of order n with
c     the parameters alpha=2*m+1 and beta=0.
c
c                   Input parameters:
c
c  mmax - the maximum order to which the polynomials are to be evaluated;
c  z - the location in R^2 where the polynomials are to be evaluated;
c       normally, expected to be inside (including boundary) the 
c       standard triangle (1) above.
c  w - the work array. must be at least mmax+1+(mmax+1)**2 long
c     
c                   Output parameters:
c
c  pols - the orthogonal polynomials evaluated at the point z 
c       ( (mmax+1)*(mmax+2)/2 of them things)
c
c       . . . allocate memory for the work arrays
c
      implicit real *8 (a-h,o-z)
      dimension z(2), pols(1), w(1)
      data init/1/
      save
c
        if( init .eq. 1 ) then
        c0=1/sqrt(3.0d0)*sqrt(sqrt(3.0d0))
        c1=sqrt(2.0d0)*sqrt(sqrt(3.0d0))
        c2=sqrt(2.0d0)*sqrt(sqrt(3.0d0))
        init=0
        endif
c
        if( mmax .eq. 0 ) then
        pols(1)=c0
        return
        endif
c
        if( mmax .eq. 1 ) then
        pols(1)=c0
        pols(2)=z(1)*c1
        pols(3)=z(2)*c2
        return
        endif
c
      iw1=1
      iw2=iw1+mmax+1
      iw3=iw2+(mmax+1)**2
      call ortho2eva0(mmax,z,pols,w(iw1),w(iw2))
      return
      end
c
c
c
c
c
      subroutine ortho2eva0(mmax,z,pols,f1,f2)
c
c     evaluate the orthonormal polynomials on the triangle 
c
      implicit real *8 (a-h,o-z)
      dimension z(2), pols(1), f1(1), f2(1)
      data init/1/
      save
c
      if( init.eq.1 ) then
         zero=0
         sqrt2=dsqrt(2.0d0)
         sqrt3=dsqrt(3.0d0)
         init=0
         r11=-1.0d0/3.0d0
         r12=-1.0d0/sqrt3
         r21=1.0d0/9.0d0*(-sqrt3)*sqrt3
         r22=1.0d0/9.0d0*(6.0d0)*sqrt3
      endif
c
      a=z(1)
      b=z(2)
c
c     ... map the standard triangle to the right
c     triangle with the vertices (-1,-1), (1,-1), (-1,1)
c
      x=r11+r12*b+a
      y=r21+r22*b
c
c     ... evaluate the Koornwinder's polynomials 
c     the via three terms recursion
c
      par1=(2*x+1+y)/2
      par2=(1-y)/2
      call klegeypols(par1,par2,mmax,f1)
      do 1200 m=0,mmax
         par1=2*m+1
         call kjacopols(y,par1,zero,(mmax-m),f2(1+m*(mmax+1)))
 1200 continue
	
      kk=0
      do 2200 m=0,mmax
         do 2000 n=0,m
         kk=kk+1
c
c     ... evaluate the polynomial (m-n, n)
c
         pols(kk)=f1(m-n+1)*f2(n+1+(m-n)*(mmax+1))
c
c     ... and normalize it
c
         scale=dsqrt((1.0d0+(m-n)+n)*(1.0d0+(m-n)+(m-n))/sqrt3) 
cccc         scale=1
         pols(kk)=pols(kk)*scale
 2000	 continue
 2200 continue
      return
      end
c
c
c
c
c
      subroutine ortho2eva3(mmax,z,pols,dersx,dersy,w)
c
c     This subroutine evaluates at the user-supplied point z
c     a collection of polynomials (of x,y) orthogonal on the
c     standard triangle, together with their derivatives with 
c     respect to x and y.
c
c     The "standard" triangle is the triangle with the
c     vertices
c     
c     (0,2/sqrt(3)), (-1,-1/sqrt(3), (1,-1/sqrt(3)).       (1)
c     
c     The polynomials evaluated by this subroutine are all 
c     orthogonal polynomials up to order mmax, arranged in the 
c     increasing order. 
c
c     This subroutine is based on the Koornwinder's representation
c     of the orthogonal polynomials on the right triangle
c
c     (-1,-1), (-1,1), (1,-1)                                       (2)
c
c     given by
c
c     K_mn(x,y) = P_m ((2*x+1+y)/(1-y)) ((1-y)/2)^m P_n^{2m+1,0} (y)
c
c     where P_m are the Legendre polynomials or order m
c     and P_n^{2m+1,0} are the Jacobi polynomials of order n with
c     the parameters alpha=2*m+1 and beta=0.
c
c                   Input parameters:
c
c  mmax - the maximum order to which the polynomials are to be evaluated;
c  z - the location in R^2 where the polynomials are to be evaluated;
c       normally, expected to be inside (including boundary) the 
c       standard triangle (1) above.
c  w - the work array. must be at least mmax+1+(mmax+1)**2 long
c     
c                   Output parameters:
c
c  pols - the orthogonal polynomials evaluated at the point z 
c       ( (mmax+1)*(mmax+2)/2 of them things)
c  dersx - the derivatives with respect to x of the polynomials 
c       returned in array pols
c  dersy - the derivatives with respect to y of the polynomials 
c       returned in array pols
c
c       . . . allocate memory for the work arrays
c
      implicit real *8 (a-h,o-z)
      dimension z(2),pols(1),dersx(1),dersy(1),w(1)
      data init/1/
      save
c
        if( init .eq. 1 ) then
        c0=1/sqrt(3.0d0)*sqrt(sqrt(3.0d0))
        c1=sqrt(2.0d0)*sqrt(sqrt(3.0d0))
        c2=sqrt(2.0d0)*sqrt(sqrt(3.0d0))
        init=0
        endif
c
        if( mmax .eq. 0 ) then
        pols(1)=c0
        dersx(1)=0
        dersy(1)=0
        return
        endif
c
        if( mmax .eq. 1 ) then
        pols(1)=c0
        pols(2)=z(1)*c1
        pols(3)=z(2)*c2
        dersx(1)=0
        dersx(2)=c1
        dersx(3)=0
        dersy(1)=0
        dersy(2)=0
        dersy(3)=c2
        return
        endif
c
      iw1=1
      iw2=iw1+mmax+1
      iw3=iw2+(mmax+1)**2
      iw4=iw3+mmax+1
      iw5=iw4+mmax+1
      iw6=iw5+(mmax+1)**2
      iw7=iw6+(mmax+1)**2
      call ortho2eva30(mmax,z,pols,dersx,dersy,
     $     w(iw1),w(iw2),w(iw3),w(iw4),w(iw5),w(iw6))
      return
      end
c
c
c
c
c
      subroutine ortho2eva30(mmax,z,pols,dersx,dersy,
     $     f1,f2,f3,f4,f5,f6)
c
c     evaluate the orthonormal polynomials on the triangle, 
c     together with their derivatives with respect to x and y.
c
      implicit real *8 (a-h,o-z)
      dimension z(2), pols(1), dersx(1), dersy(1)
      dimension f1(1), f2(1)
      dimension f3(1), f4(1), f5(1), f6(1)
      data init/1/
      save
c
      if( init.eq.1 ) then
         zero=0
         sqrt2=dsqrt(2.0d0)
         sqrt3=dsqrt(3.0d0)
         init=0
         r11=-1.0d0/3.0d0
         r12=-1.0d0/sqrt3
         r21=1.0d0/9.0d0*(-sqrt3)*sqrt3
         r22=1.0d0/9.0d0*(6.0d0)*sqrt3
      endif
c
      x=z(1)
      y=z(2)
c
c     ... map the standard triangle to the right
c     triangle with the vertices (-1,-1), (1,-1), (-1,1)
c
      a=r11+r12*y+x
      b=r21+r22*y

c
c     ... evaluate the Koornwinder's polynomials 
c     the via three terms recursion
c
      par1=(2*a+1+b)/2
      par2=(1-b)/2
      call klegeypols3(par1,par2,mmax,f1,f3,f4)

      do 1200 m=0,mmax
         par1=2*m+1
         call kjacopols2(
     $        b,par1,zero,(mmax-m),f2(1+m*(mmax+1)),f5(1+m*(mmax+1)))
 1200 continue
	
      kk=0
      do 2200 m=0,mmax
         do 2000 n=0,m
         kk=kk+1
c
c     ... evaluate the polynomial (m-n, n), and their derivatives
c     with respect to x,y
c
         pols(kk)=f1(m-n+1)*f2(n+1+(m-n)*(mmax+1))
c
         dersx(kk)=
     $        f3(m-n+1)*f2(n+1+(m-n)*(mmax+1))
c
         dersy(kk)=
     $        f1(m-n+1)*f5(n+1+(m-n)*(mmax+1))*r22 +
     $        f3(m-n+1)*f2(n+1+(m-n)*(mmax+1))*(r12+r22/2) +
     $        f4(m-n+1)*f2(n+1+(m-n)*(mmax+1))*(-r22/2) 

c
c     ... and normalize it
c
         scale=dsqrt((1.0d0+(m-n)+n)*(1.0d0+(m-n)+(m-n))/sqrt3) 
cccc         scale=1
         pols(kk)=pols(kk)*scale
         dersx(kk)=dersx(kk)*scale
         dersy(kk)=dersy(kk)*scale
 2000	 continue
 2200 continue


      return
      end
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c        this is the end of the debugging code and the beginning 
c        of the evaluation of scaled Legendre polynomials
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
        subroutine klegeypols(x,y,n,pols)
        implicit real *8 (a-h,o-z)
        dimension pols(1)
c
c       evaluate a sequence of scaled Legendre polynomials 
c       P_n(x/y) y^n, with the parameter y \in [0..1]
c
c       ...
c
        pkm1=1
        pk=x
c
        pk=1
        pkp1=x
c
c
c        if n=0 or n=1 - exit
c
        if(n .ge. 2) goto 1200
        pols(1)=1
        if(n .eq. 0) return
c
        pols(2)=x
        return
 1200 continue
c
        pols(1)=1
        pols(2)=x
c
c       n is greater than 2. conduct recursion
c
        do 2000 k=1,n-1
c
        pkm1=pk
        pk=pkp1
        pkp1=( (2*k+1)*x*pk-k*pkm1*y**2 )/(k+1)
        pols(k+2)=pkp1
 2000 continue
c
        return        
        end
c
c
c
c
c
        subroutine klegeypols3(x,y,n,pols,dersx,dersy)
        implicit real *8 (a-h,o-z)
        dimension pols(1),dersx(1),dersy(1)
c
c     evaluate a sequence of scaled Legendre polynomials P_n(x/y) y^n,
c     with the parameter y \in [0..1], together with their derivatives
c     with respect to the parameters x and y.
c
c       ... if n=0 or n=1 - exit
c
        if(n .ge. 2) goto 1200
        pols(1)=1
        dersx(1)=0
        dersy(1)=0
        if(n .eq. 0) return
c
        pols(2)=x
        dersx(2)=1
        dersy(2)=0
        return
 1200 continue
c
c       n is greater than 1. conduct recursion
c
        pols(1)=1
        dersx(1)=0
        dersy(1)=0
        pols(2)=x
        dersx(2)=1
        dersy(2)=0
c
        pkm1=1
        pk=x
        dkm1=0
        dk=1
        ykm1=0
        yk=0
c
        pk=1
        pkp1=x
        dk=0
        dkp1=1
        yk=0
        ykp1=0
c
        do 2000 k=1,n-1
c
        pkm1=pk
        pk=pkp1
        dkm1=dk
        dk=dkp1
        ykm1=yk
        yk=ykp1
c
        pkp1=( (2*k+1)*x*pk-k*pkm1*y*y )/(k+1)
        dkp1=( (2*k+1)*(x*dk+pk)-k*dkm1*y*y )/(k+1)
        ykp1=( (2*k+1)*(x*yk)-k*(pkm1*2*y+ykm1*y*y) )/(k+1)
        pols(k+2)=pkp1
        dersx(k+2)=dkp1
        dersy(k+2)=ykp1
 2000 continue
c
        return
        end
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c        this is the end of the debugging code and the beginning 
c        of the evaluation of Jacobi polynomials
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
        subroutine kjacopols(x,a,b,n,pols)
        implicit real *8 (a-h,o-z)
        dimension pols(1)
c
c       evaluates a bunch of Jacobi polynomials 
c       at the user-provided point x
c
c       ... if n=0 or n=1 - exit
c
        if(n .ge. 2) goto 1200
        pols(1)=1
        if(n .eq. 0) return
c
        pols(2)=(a/2-b/2)+(1+a/2+b/2)*x
        return
 1200 continue
c
        pols(1)=1
        pols(2)=(a/2-b/2)+(1+a/2+b/2)*x
c
c       n is greater than 2. conduct recursion
c
        pkm1=1
        pk=(a/2-b/2)+(1+a/2+b/2)*x
        pk=1
        pkp1=(a/2-b/2)+(1+a/2+b/2)*x
c
        do 2000 k=2,n
c
        pkm1=pk
        pk=pkp1
        alpha=(2*k+a+b-1)*(a**2-b**2+(2*k+a+b-2)*(2*k+a+b)*x)
        beta=2*(k+a-1)*(k+b-1)*(2*k+a+b)
        pkp1=(alpha*pk-beta*pkm1)/(2*k*(k+a+b)*(2*k+a+b-2))
        pols(k+1)=pkp1
 2000 continue
c
        return
        end
c
c
c
c
c
        subroutine kjacopols2(x,a,b,n,pols,ders)
        implicit real *8 (a-h,o-z)
        dimension pols(1), ders(1)
c
c       evaluates a bunch of Jacobi polynomials (together
c       with their derivatives) at the user-provided point x
c
c       ... if n=0 or n=1 - exit
c
        if(n .ge. 2) goto 1200
        pols(1)=1
        ders(1)=0
        if(n .eq. 0) return
c
        pols(2)=(a/2-b/2)+(1+a/2+b/2)*x
        ders(2)=(1+a/2+b/2)
        return
 1200 continue
c
        pols(1)=1
        ders(1)=0
        pols(2)=(a/2-b/2)+(1+a/2+b/2)*x
        ders(2)=(1+a/2+b/2)
c
c       n is greater than 2. conduct recursion
c
        pkm1=1
        dkm1=0
        pk=(a/2-b/2)+(1+a/2+b/2)*x
        dk=(1+a/2+b/2)
c
        pk=1
        dk=0
        pkp1=(a/2-b/2)+(1+a/2+b/2)*x
        dkp1=(1+a/2+b/2)
c
        do 2000 k=2,n
c
        pkm1=pk
        pk=pkp1
        dkm1=dk
        dk=dkp1
        alpha1=(2*k+a+b-1)*(a**2-b**2)
        alpha2=(2*k+a+b-1)*((2*k+a+b-2)*(2*k+a+b))
        beta=2*(k+a-1)*(k+b-1)*(2*k+a+b)
        gamma=(2*k*(k+a+b)*(2*k+a+b-2))
        pkp1=((alpha1+alpha2*x)*pk-beta*pkm1)/gamma
        dkp1=((alpha1+alpha2*x)*dk-beta*dkm1+alpha2*pk)/gamma
        pols(k+1)=pkp1
        ders(k+1)=dkp1
 2000 continue
c
        return
        end
