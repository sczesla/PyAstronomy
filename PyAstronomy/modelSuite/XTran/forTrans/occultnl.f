      subroutine occultnl(rl,c1,c2,c3,c4,b0,mulimb0,mulimbf,nb)
C Please cite Mandel & Agol (2002) if making use of this routine.
      implicit none
      integer i,j,nb,nr,i1,i2,nmax
      parameter (nmax=513)
      real*8 mulimbf(nb,5),pi,c1,c2,c3,c4,rl,bt0(nb),b0(nb),mulimb0(nb),
     &       mulimb(nb),mulimbp(nb),dt,t(nmax),th(nmax),r(nmax),sig,
     &       mulimb1(nb),mulimbhalf(nb),mulimb3half(nb),mulimb2(nb),
     &       sig1,sig2,omega,dmumax,fac,mu(nb),f1,f2
      pi=acos(-1.d0)
C  This routine uses the results for a uniform source to
C  compute the lightcurve for a limb-darkened source
C  (5-1-02 notes)
C Input:
C   rl        radius of the lens   in units of the source radius
C   c1-c4     limb-darkening coefficients
C   b0        impact parameter normalized to source radius
C Output:
C  mulimb0 limb-darkened magnification
C  mulimbf lightcurves for each component
C  
C  First, make grid in radius:
C  Call magnification of uniform source:
      call occultuniform(b0,rl,mulimb0,nb)
      i1=nb
      i2=1
      fac=0.d0
      do i=1,nb
        bt0(i)=b0(i)
        mulimbf(i,1)=1.d0
        mulimbf(i,2)=0.8d0
        mulimbf(i,3)=2.d0/3.d0
        mulimbf(i,4)=4.d0/7.d0
        mulimbf(i,5)=0.5d0
        mulimb(i)=mulimb0(i)
        if(mulimb0(i).ne.1.d0) then
          i1=min(i1,i)
          i2=max(i2,i)
        endif
        fac=max(fac,abs(mulimb0(i)-1.d0))
      enddo
C print,rl
      omega=4.*((1.d0-c1-c2-c3-c4)/4.+c1/5.+c2/6.+c3/7.+c4/8.)
      nr=2
      dmumax=1.d0
c      write(6,*) 'i1,i2 ',i1,i2
      do while (dmumax.gt.fac*1.d-3)
        do i=i1,i2
          mulimbp(i)=mulimb(i)
        enddo
        nr=nr*2
c        write(6,*) 'nr ',nr
        dt=0.5d0*pi/dble(nr)
        do j=1,nr+1
          t(j) =dt*dble(j-1)
          th(j)=t(j)+0.5d0*dt
          r(j)=sin(t(j))
        enddo
        sig=sqrt(cos(th(nr)))
        do i=i1,i2
          mulimbhalf(i) =sig**3*mulimb0(i)/(1.d0-r(nr))
          mulimb1(i)    =sig**4*mulimb0(i)/(1.d0-r(nr))
          mulimb3half(i)=sig**5*mulimb0(i)/(1.d0-r(nr))
          mulimb2(i)    =sig**6*mulimb0(i)/(1.d0-r(nr))
        enddo
        do j=2,nr
          do i=1,nb
            b0(i)=bt0(i)/r(j)
          enddo
C  Calculate uniform magnification at intermediate radii:
          call occultuniform(b0,rl/r(j),mu,nb)
C  Equation (29):
          sig1=sqrt(cos(th(j-1)))
          sig2=sqrt(cos(th(j)))
          dmumax=0.d0
          do i=i1,i2
            f1=r(j)*r(j)*mu(i)/(r(j)-r(j-1))
            f2=r(j)*r(j)*mu(i)/(r(j+1)-r(j))
            mulimbhalf(i) =mulimbhalf(i) +f1*sig1**3-f2*sig2**3
            mulimb1(i)    =mulimb1(i)    +f1*sig1**4-f2*sig2**4
            mulimb3half(i)=mulimb3half(i)+f1*sig1**5-f2*sig2**5
            mulimb2(i)    =mulimb2(i)    +f1*sig1**6-f2*sig2**6
            mulimb(i)=((1.d0-c1-c2-c3-c4)*mulimb0(i)+c1*mulimbhalf(i)*dt
     &        +c2*mulimb1(i)*dt+c3*mulimb3half(i)*dt+c4*mulimb2(i)*dt)
     &        /omega
            if(mulimb(i)+mulimbp(i).ne.0.d0) then 
              dmumax=max(dmumax,abs(mulimb(i)-mulimbp(i))/(mulimb(i)+
     &               mulimbp(i)))
            endif
          enddo
        enddo
      enddo
      do i=i1,i2
        mulimbf(i,1)=mulimb0(i)
        mulimbf(i,2)=mulimbhalf(i)*dt
        mulimbf(i,3)=mulimb1(i)*dt
        mulimbf(i,4)=mulimb3half(i)*dt
        mulimbf(i,5)=mulimb2(i)*dt
        mulimb0(i)=mulimb(i)
      enddo
      do i=1,nb
        b0(i)=bt0(i)
      enddo
      return
      end
      
      subroutine occultuniform(b0,w,muo1,nb)
      implicit none
      integer i,nb
      real*8 muo1(nb),w,b0(nb),z,pi,lambdae,kap0,kap1
      if(abs(w-0.5d0).lt.1.d-3) w=0.5d0
      pi=acos(-1.d0)
C  This routine computes the lightcurve for occultation
C  of a uniform source without microlensing  (Mandel & Agol 2002).
C Input:
C 
C  rs   radius of the source (set to unity)
C  b0   impact parameter in units of rs
C  w    occulting star size in units of rs
C 
C Output:
C  muo1 fraction of flux at each b0 for a uniform source
C 
C  Now, compute pure occultation curve:
      do i=1,nb
C  substitute z=b0(i) to shorten expressions
        z=b0(i)
C  the source is unocculted:
C  Table 3, I.
        if(z.ge.1.d0+w) then
          muo1(i)=1.d0
          goto 1
        endif
C  the  source is completely occulted:
C  Table 3, II.
        if(w.ge.1.d0.and.z.le.w-1.d0) then
          muo1(i)=0.d0
          goto 1
        endif
C  the source is partly occulted and the occulting object crosses the limb:
C  Equation (26):
        if(z.ge.abs(1.d0-w).and.z.le.1.d0+w) then
          kap1=acos(min((1.d0-w*w+z*z)/2.d0/z,1.d0))
          kap0=acos(min((w*w+z*z-1.d0)/2.d0/w/z,1.d0))
          lambdae=w*w*kap0+kap1
          lambdae=(lambdae-0.5d0*sqrt(max(4.d0*z*z-(1.d0+z*z-w*w)**2,
     &            0.d0)))/pi
          muo1(i)=1.d0-lambdae
        endif
C  the occulting object transits the source star (but doesn't
C  completely cover it):
        if(z.le.1.d0-w) muo1(i)=1.d0-w*w
 1      continue
      enddo
C muo1=1.d0-lambdae
      return
      end
