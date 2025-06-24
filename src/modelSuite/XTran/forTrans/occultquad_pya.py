"""
    Python implementation of the routines below

    subroutine occultquad(z0,u1,u2,p,muo1,mu0,nz)
C  This routine computes the lightcurve for occultation
C  of a quadratically limb-darkened source without microlensing.
C  Please cite Mandel & Agol (2002) if you make use of this routine
C  in your research.  Please report errors or bugs to agol@tapir.caltech.edu
      implicit none
      integer i,nz
      double precision z0(nz),u1,u2,p,muo1(nz),mu0(nz),
     &       mu(nz),lambdad(nz),etad(nz),lambdae(nz),lam,
     &       pi,x1,x2,x3,z,omega,kap0,kap1,q,Kk,Ek,Pk,n,ellec,ellk,rj
      if(abs(p-0.5d0).lt.1.d-3) p=0.5d0
C
C Input:
C
C rs   radius of the source (set to unity)
C z0   impact parameter in units of rs
C p    occulting star size in units of rs
C u1   linear    limb-darkening coefficient (gamma_1 in paper)
C u2   quadratic limb-darkening coefficient (gamma_2 in paper)
C
C Output:
C
C muo1 fraction of flux at each z0 for a limb-darkened source
C mu0  fraction of flux at each z0 for a uniform source
C
C Limb darkening has the form:
C  I(r)=[1-u1*(1-sqrt(1-(r/rs)^2))-u2*(1-sqrt(1-(r/rs)^2))^2]/(1-u1/3-u2/6)/pi
C 
C To use this routine
C
C Now, compute pure occultation curve:
      omega=1.d0-u1/3.d0-u2/6.d0
      pi=acos(-1.d0)
C Loop over each impact parameter:
      do i=1,nz
C substitute z=z0(i) to shorten expressions
        z=z0(i)
        x1=(p-z)**2
        x2=(p+z)**2
        x3=p**2-z**2
C the source is unocculted:
C Table 3, I.
        if(z.ge.1.d0+p) then
          lambdad(i)=0.d0
          etad(i)=0.d0
          lambdae(i)=0.d0
          goto 10
        endif
C the  source is completely occulted:
C Table 3, II.
        if(p.ge.1.d0.and.z.le.p-1.d0) then
          lambdad(i)=1.d0
          etad(i)=1.d0
          lambdae(i)=1.d0
          goto 10
        endif
C the source is partly occulted and the occulting object crosses the limb:
C Equation (26):
        if(z.ge.abs(1.d0-p).and.z.le.1.d0+p) then
          kap1=acos(min((1.d0-p*p+z*z)/2.d0/z,1.d0))
          kap0=acos(min((p*p+z*z-1.d0)/2.d0/p/z,1.d0))
          lambdae(i)=p*p*kap0+kap1
          lambdae(i)=(lambdae(i)-0.5d0*sqrt(max(4.d0*z*z-
     &               (1.d0+z*z-p*p)**2,0.d0)))/pi
        endif
C the occulting object transits the source star (but doesn't
C completely cover it):
        if(z.le.1.d0-p) lambdae(i)=p*p
C the edge of the occulting star lies at the origin- special 
C expressions in this case:
        if(abs(z-p).lt.1.d-4*(z+p)) then
C Table 3, Case V.:
          if(z.ge.0.5d0) then
            lam=0.5d0*pi
            q=0.5d0/p
            Kk=ellk(q)
            Ek=ellec(q)
C Equation 34: lambda_3
            lambdad(i)=1.d0/3.d0+16.d0*p/9.d0/pi*(2.d0*p*p-1.d0)*Ek-
     &                 (32.d0*p**4-20.d0*p*p+3.d0)/9.d0/pi/p*Kk
C Equation 34: eta_1
            etad(i)=1.d0/2.d0/pi*(kap1+p*p*(p*p+2.d0*z*z)*kap0-
     &              (1.d0+5.d0*p*p+z*z)/4.d0*sqrt((1.d0-x1)*(x2-1.d0)))
            if(p.eq.0.5d0) then
C Case VIII: p=1/2, z=1/2
              lambdad(i)=1.d0/3.d0-4.d0/pi/9.d0
              etad(i)=3.d0/32.d0
            endif
            goto 10
          else
C Table 3, Case VI.:
            lam=0.5d0*pi
            q=2.d0*p
            Kk=ellk(q)
            Ek=ellec(q)
C Equation 34: lambda_4
            lambdad(i)=1.d0/3.d0+2.d0/9.d0/pi*(4.d0*(2.d0*p*p-1.d0)*Ek+
     &                 (1.d0-4.d0*p*p)*Kk)
C Equation 34: eta_2
            etad(i)=p*p/2.d0*(p*p+2.d0*z*z)
            goto 10
          endif
        endif
C the occulting star partly occults the source and crosses the limb:
C Table 3, Case III:
        if((z.gt.0.5d0+abs(p-0.5d0).and.z.lt.1.d0+p).or.(p.gt.0.5d0.
     &      and.z.gt.abs(1.d0-p)*1.0001d0.and.z.lt.p)) then
          lam=0.5d0*pi
          q=sqrt((1.d0-(p-z)**2)/4.d0/z/p)
          Kk=ellk(q)
          Ek=ellec(q)
          n=1.d0/x1-1.d0
          Pk=Kk-n/3.d0*rj(0.d0,1.d0-q*q,1.d0,1.d0+n)
C Equation 34, lambda_1:
          lambdad(i)=1.d0/9.d0/pi/sqrt(p*z)*(((1.d0-x2)*(2.d0*x2+
     &        x1-3.d0)-3.d0*x3*(x2-2.d0))*Kk+4.d0*p*z*(z*z+
     &        7.d0*p*p-4.d0)*Ek-3.d0*x3/x1*Pk)
          if(z.lt.p) lambdad(i)=lambdad(i)+2.d0/3.d0
C Equation 34, eta_1:
          etad(i)=1.d0/2.d0/pi*(kap1+p*p*(p*p+2.d0*z*z)*kap0-
     &          (1.d0+5.d0*p*p+z*z)/4.d0*sqrt((1.d0-x1)*(x2-1.d0)))
          goto 10
        endif
C the occulting star transits the source:
C Table 3, Case IV.:
        if(p.le.1.d0.and.z.le.(1.d0-p)*1.0001d0) then
          lam=0.5d0*pi
          q=sqrt((x2-x1)/(1.d0-x1))
          Kk=ellk(q)
          Ek=ellec(q)
          n=x2/x1-1.d0
          Pk=Kk-n/3.d0*rj(0.d0,1.d0-q*q,1.d0,1.d0+n)
C Equation 34, lambda_2:
          lambdad(i)=2.d0/9.d0/pi/sqrt(1.d0-x1)*((1.d0-5.d0*z*z+p*p+
     &         x3*x3)*Kk+(1.d0-x1)*(z*z+7.d0*p*p-4.d0)*Ek-3.d0*x3/x1*Pk)
          if(z.lt.p) lambdad(i)=lambdad(i)+2.d0/3.d0
          if(abs(p+z-1.d0).le.1.d-4) then
            lambdad(i)=2/3.d0/pi*acos(1.d0-2.d0*p)-4.d0/9.d0/pi*
     &            sqrt(p*(1.d0-p))*(3.d0+2.d0*p-8.d0*p*p)
          endif
C Equation 34, eta_2:
          etad(i)=p*p/2.d0*(p*p+2.d0*z*z)
        endif
 10     continue
C Now, using equation (33):
        muo1(i)=1.d0-((1.d0-u1-2.d0*u2)*lambdae(i)+(u1+2.d0*u2)*
     &      lambdad(i)+u2*etad(i))/omega
C Equation 25:
        mu0(i)=1.d0-lambdae(i)
      enddo
      return
      end

      FUNCTION rc(x,y)
      REAL*8 rc,x,y,ERRTOL,TINY,SQRTNY,BIG,TNBG,COMP1,COMP2,THIRD,C1,C2,
     *C3,C4
      PARAMETER (ERRTOL=.04d0,TINY=1.69d-38,SQRTNY=1.3d-19,BIG=3.d37,
     *TNBG=TINY*BIG,COMP1=2.236d0/SQRTNY,COMP2=TNBG*TNBG/25.d0,
     *THIRD=1.d0/3.d0,C1=.3d0,C2=1.d0/7.d0,C3=.375d0,C4=9.d0/22.d0)
      REAL*8 alamb,ave,s,w,xt,yt
      if(x.lt.0..or.y.eq.0..or.(x+abs(y)).lt.TINY.or.(x+
     *abs(y)).gt.BIG.or.(y.lt.-COMP1.and.x.gt.0..and.x.lt.COMP2))pause 
     *'invalid arguments in rc'
      if(y.gt.0.d0)then
        xt=x
        yt=y
        w=1.
      else
        xt=x-y
        yt=-y
        w=sqrt(x)/sqrt(xt)
      endif
1     continue
        alamb=2.d0*sqrt(xt)*sqrt(yt)+yt
        xt=.25d0*(xt+alamb)
        yt=.25d0*(yt+alamb)
        ave=THIRD*(xt+yt+yt)
        s=(yt-ave)/ave
      if(abs(s).gt.ERRTOL)goto 1
      rc=w*(1.d0+s*s*(C1+s*(C2+s*(C3+s*C4))))/sqrt(ave)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software

      FUNCTION rj(x,y,z,p)
      REAL*8 rj,p,x,y,z,ERRTOL,TINY,BIG,C1,C2,C3,C4,C5,C6,C7,C8
      PARAMETER (ERRTOL=.05d0,TINY=2.5d-13,BIG=9.d11,C1=3.d0/14.d0,
     *C2=1.d0/3.d0,C3=3.d0/22.d0,C4=3.d0/26.d0,C5=.75d0*C3,
     *C6=1.5d0*C4,C7=.5d0*C2,C8=C3+C3)
CU    USES rc,rf
      REAL*8 a,alamb,alpha,ave,b,beta,delp,delx,dely,delz,ea,eb,ec,ed,ee,
     *fac,pt,rcx,rho,sqrtx,sqrty,sqrtz,sum,tau,xt,yt,zt,rc,rf
      if(min(x,y,z).lt.0..or.min(x+y,x+z,y+z,abs(p)).lt.TINY.or.max(x,y,
     *z,abs(p)).gt.BIG)pause 'invalid arguments in rj'
      sum=0.d0
      fac=1.d0
      if(p.gt.0.d0)then
        xt=x
        yt=y
        zt=z
        pt=p
      else
        xt=min(x,y,z)
        zt=max(x,y,z)
        yt=x+y+z-xt-zt
        a=1.d0/(yt-p)
        b=a*(zt-yt)*(yt-xt)
        pt=yt+b
        rho=xt*zt/yt
        tau=p*pt/yt
        rcx=rc(rho,tau)
      endif
1     continue
        sqrtx=sqrt(xt)
        sqrty=sqrt(yt)
        sqrtz=sqrt(zt)
        alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
        alpha=(pt*(sqrtx+sqrty+sqrtz)+sqrtx*sqrty*sqrtz)**2
        beta=pt*(pt+alamb)**2
        sum=sum+fac*rc(alpha,beta)
        fac=.25d0*fac
        xt=.25d0*(xt+alamb)
        yt=.25d0*(yt+alamb)
        zt=.25d0*(zt+alamb)
        pt=.25d0*(pt+alamb)
        ave=.2d0*(xt+yt+zt+pt+pt)
        delx=(ave-xt)/ave
        dely=(ave-yt)/ave
        delz=(ave-zt)/ave
        delp=(ave-pt)/ave
      if(max(abs(delx),abs(dely),abs(delz),abs(delp)).gt.ERRTOL)goto 1
      ea=delx*(dely+delz)+dely*delz
      eb=delx*dely*delz
      ec=delp**2
      ed=ea-3.d0*ec
      ee=eb+2.d0*delp*(ea-ec)
      rj=3.d0*sum+fac*(1.d0+ed*(-C1+C5*ed-C6*ee)+eb*(C7+delp*
     *(-C8+delp*C4))+delp*ea*(C2-delp*C3)-C2*delp*ec)/(ave*sqrt(ave))
      if (p.le.0.d0) rj=a*(b*rj+3.d0*(rcx-rf(xt,yt,zt)))
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software

      function ellec(k)
      implicit none
      double precision k,m1,a1,a2,a3,a4,b1,b2,b3,b4,ee1,ee2,ellec
C Computes polynomial approximation for the complete elliptic
C integral of the second kind (Hasting's approximation):
      m1=1.d0-k*k
      a1=0.44325141463d0
      a2=0.06260601220d0
      a3=0.04757383546d0
      a4=0.01736506451d0
      b1=0.24998368310d0
      b2=0.09200180037d0
      b3=0.04069697526d0
      b4=0.00526449639d0
      ee1=1.d0+m1*(a1+m1*(a2+m1*(a3+m1*a4)))
      ee2=m1*(b1+m1*(b2+m1*(b3+m1*b4)))*log(1.d0/m1)
      ellec=ee1+ee2
      return
      end

      function ellk(k)
      implicit none
      double precision a0,a1,a2,a3,a4,b0,b1,b2,b3,b4,ellk,
     &       ek1,ek2,k,m1
C Computes polynomial approximation for the complete elliptic
C integral of the first kind (Hasting's approximation):
      m1=1.d0-k*k
      a0=1.38629436112d0
      a1=0.09666344259d0
      a2=0.03590092383d0
      a3=0.03742563713d0
      a4=0.01451196212d0
      b0=0.5d0
      b1=0.12498593597d0
      b2=0.06880248576d0
      b3=0.03328355346d0
      b4=0.00441787012d0
      ek1=a0+m1*(a1+m1*(a2+m1*(a3+m1*a4)))
      ek2=(b0+m1*(b1+m1*(b2+m1*(b3+m1*b4))))*log(m1)
      ellk=ek1-ek2
      return
      end

      FUNCTION rf(x,y,z)
      REAL*8 rf,x,y,z,ERRTOL,TINY,BIG,THIRD,C1,C2,C3,C4
      PARAMETER (ERRTOL=.08d0,TINY=1.5d-38,BIG=3.d37,THIRD=1.d0/3.d0,
     *C1=1.d0/24.d0,C2=.1d0,C3=3.d0/44.d0,C4=1.d0/14.d0)
      REAL*8 alamb,ave,delx,dely,delz,e2,e3,sqrtx,sqrty,sqrtz,xt,yt,zt
      if(min(x,y,z).lt.0.d0.or.min(x+y,x+z,y+z).lt.TINY.or.max(x,y,
     *z).gt.BIG)pause 'invalid arguments in rf'
      xt=x
      yt=y
      zt=z
1     continue
        sqrtx=sqrt(xt)
        sqrty=sqrt(yt)
        sqrtz=sqrt(zt)
        alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
        xt=.25d0*(xt+alamb)
        yt=.25d0*(yt+alamb)
        zt=.25d0*(zt+alamb)
        ave=THIRD*(xt+yt+zt)
        delx=(ave-xt)/ave
        dely=(ave-yt)/ave
        delz=(ave-zt)/ave
      if(max(abs(delx),abs(dely),abs(delz)).gt.ERRTOL)goto 1
      e2=delx*dely-delz**2
      e3=delx*dely*delz
      rf=(1.d0+(C1*e2-C2-C3*e3)*e2+C4*e3)/sqrt(ave)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software 0NL&WR2.
"""

import numpy as np
from numpy import sqrt, log
from PyAstronomy.pyaC import pyaErrors as PE

class OccultQuadPy:

    def occultquad(self,z0,u1,u2,p,nz=None):
        """
        Occultquad routine
        
        Python implementation of code provided by Mandel and Agol (2002).
        Please cite Mandel & Agol (2002) if you make use of this routine.
        The original authors are not responsible for the python
        re-implmentation.
        
        Parameters
        ----------
        z0 : array
            Impact parameter in units of stellar radius
        u1 : int
            Linear limb-darkening coefficient (gamma_1 in paper)
        u2 : int
            Quadratic limb-darkening coefficient (gamma_2 in paper)   
        p : int
            Occulting star (or planet) size in units of stellar radius
        nz : int, optional
            Number of values in z0. Determined here from z0, parameter
            introduced to maintain API.
        
        Returns
        -------
        muo1 : array
            Fraction of flux at each z0 (impact parameter) for a limb-darkened source
        mu0 : array
            Fraction of flux at each z0 (impact parameter) for a uniform source
        """
        # Determine number of impact parameter values (ignore keyword altogether)
        nz = len(z0)

        if abs(p-0.5) < 1e-3:
            p = 0.5

        omega = 1 - u1/3 - u2/6
        pi = np.pi
        lambdad = np.zeros(nz)
        etad = np.zeros(nz)
        lambdae = np.zeros(nz)
# C Loop over each impact parameter:
        for i in range(nz):
# C substitute z=z0(i) to shorten expressions
            z = z0[i]
            x1 = (p-z)**2
            x2 = (p+z)**2
            x3 = p**2-z**2
# C the source is unocculted:
# C Table 3, I.
            if z >= 1.0+p:
                lambdad[i] = 0
                etad[i] = 0
                lambdae[i] = 0
                continue
            elif (p >= 1.0) and (z <= p-1.0):
# C the  source is completely occulted:
# C Table 3, II.
                lambdad[i] = 1
                etad[i] = 1
                lambdae[i] = 1
                continue          
            elif (z >= abs(1-p)) and (z <= 1+p):
# C the source is partly occulted and the occulting object crosses the limb:
# C Equation (26):
                kap1 = np.arccos(min((1-p*p+z*z)/2/z,1))
                kap0 = np.arccos(min((p*p+z*z-1)/2/p/z,1))
                lambdae[i] = p*p*kap0+kap1
                lambdae[i] = (lambdae[i]-0.5*sqrt(max(4*z*z-(1+z*z-p*p)**2,0)))/pi

# C the occulting object transits the source star (but doesn't
# C completely cover it):
            if z <= 1-p:
                lambdae[i] = p*p
# C the edge of the occulting star lies at the origin- special 
# C expressions in this case:
            if abs(z-p) < 1e-4*(z+p):
# C Table 3, Case V.:
                if z >= 0.5:
                    lam = 0.5*pi
                    q = 0.5/p
                    Kk = self.ellk(q)
                    Ek = self.ellec(q)
# C Equation 34: lambda_3
                    lambdad[i] = 1/3+16*p/9/pi*(2*p*p-1)*Ek- \
                        (32*p**4-20*p*p+3)/9/pi/p*Kk
# C Equation 34: eta_1
                    etad[i] = 1/2/pi*(kap1+p*p*(p*p+2*z*z)*kap0- \
                                      (1+5*p*p+z*z)/4*sqrt((1-x1)*(x2-1)))
                    if p == 0.5:
# C Case VIII: p=1/2, z=1/2
                        lambdad[i] = 1/3-4/pi/9
                        etad[i] = 3/32
                    continue
                else:
# C Table 3, Case VI.:
                    lam = 0.5*pi
                    q = 2*p
                    Kk = self.ellk(q)
                    Ek = self.ellec(q)
# C Equation 34: lambda_4
                    lambdad[i] = 1/3+2/9/pi*(4*(2*p*p-1)*Ek + (1-4*p*p)*Kk)
# C Equation 34: eta_2
                    etad[i] = p*p/2*(p*p+2*z*z)
                    continue

# C the occulting star partly occults the source and crosses the limb:
# C Table 3, Case III:
            if ((z > 0.5+abs(p-0.5)) and (z < 1+p)) or ((p > 0.5) and (z > abs(1-p)*1.000001) and (z < p)):
                lam = 0.5*pi
                q = sqrt((1-(p-z)**2)/4/z/p)
                Kk = self.ellk(q)
                Ek = self.ellec(q)
                n = 1/x1-1
                Pk = Kk-n/3*self.rj(0,1-q*q,1,1+n)
# C Equation 34, lambda_1:
                lambdad[i] = 1/9/pi/sqrt(p*z)*(((1-x2)*(2*x2+ \
                                x1-3)-3*x3*(x2-2))*Kk+4*p*z*(z*z+ \
                                7*p*p-4)*Ek-3*x3/x1*Pk)
                if z < p:
                    lambdad[i] = lambdad[i]+2/3
# C Equation 34, eta_1:
                etad[i] = 1/2/pi*(kap1+p*p*(p*p+2*z*z)*kap0- \
                            (1+5*p*p+z*z)/4*sqrt((1-x1)*(x2-1)))
                continue
# C the occulting star transits the source:
# C Table 3, Case IV.:
            if (p <= 1) and (z <= (1-p)*1.000001):
                lam = 0.5*pi
                q = sqrt((x2-x1)/(1-x1))
                Kk = self.ellk(q)
                Ek = self.ellec(q)
                n = x2/x1-1
                Pk = Kk-n/3*self.rj(0,1-q*q,1,1+n)
# C Equation 34, lambda_2:
                lambdad[i] = 2/9/pi/sqrt(1-x1)*((1-5*z*z+p*p+ \
                                x3*x3)*Kk+(1-x1)*(z*z+7*p*p-4)*Ek-3*x3/x1*Pk)
                if z < p:
                    lambdad[i] = lambdad[i]+2/3
                if abs(p+z-1) <= 1e-4:
                    lambdad[i] = 2/3/pi*np.arccos(1-2*p)-4/9/pi*sqrt(p*(1-p))*(3+2*p-8*p*p)
# C Equation 34, eta_2:
                etad[i] = p*p/2*(p*p+2*z*z)
        
 # 10     continue
# C Now, using equation (33):
        
        muo1 = 1-((1-u1-2*u2)*lambdae+(u1+2*u2)*lambdad+u2*etad)/omega
# C Equation 25:
        mu0 = 1-lambdae
        return muo1, mu0


    def rf(self,x,y,z, ERRTOL=0.08, TINY=1.5e-38, BIG=3e37):
        if (min(x,y,z) < 0) or (min(x+y,x+z,y+z) < TINY) or (max(x,y,z) > BIG):
            raise(ValueError('invalid arguments in rf'))
        
        xt, yt, zt = x, y, z

        while True:
            sqrtx = sqrt(xt)
            sqrty = sqrt(yt)
            sqrtz = sqrt(zt)
            alamb = sqrtx*(sqrty+sqrtz) + sqrty*sqrtz
            xt = 0.25*(xt+alamb)
            yt = 0.25*(yt+alamb)
            zt = 0.25*(zt+alamb)
            ave = (xt+yt+zt)/3
            delx = (ave-xt)/ave
            dely = (ave-yt)/ave
            delz = (ave-zt)/ave
            if(max(abs(delx),abs(dely),abs(delz)) <= ERRTOL):
                break
     
        C1 = 1/24
        C2 = 0.1
        C3 = 3/44
        C4 = 1/14
     
        e2 = delx*dely-delz**2
        e3 = delx*dely*delz
        rf = (1.0+(C1*e2-C2-C3*e3)*e2+C4*e3)/sqrt(ave)
        return rf
      
    def rj(self,x,y,z,p,ERRTOL=0.05,TINY=2.5e-13,BIG=9.e11):

        if min(x,y,z) < 0:
            raise(PE.PyAValError(f"min(x,y,z) < 0 violated (x,y,z = {x}, {y}, {z})", \
                              where="rj in OccultQuadPy"))
        if min(x+y,x+z,y+z,abs(p)) < TINY:
            raise(PE.PyAValError(f"min(x+y,x+z,y+z,abs(p)) < TINY violated (x+y,x+z,y+z,abs(p) = {x+y}, {x+z}, {y+z}, {abs(p)})", \
                              where="rj in OccultQuadPy"))     
        if max(x,y,z,abs(p)) > BIG:
            raise(PE.PyAValError(f"max(x,y,z,abs(p) > BIG violated (x,y,z,abs(p) = {x}, {y}, {z}, {abs(p)})", \
                              where="rj in OccultQuadPy"))             
            
        
        summe = 0.0
        fac = 1.0
        if(p > 0.0):
            xt, yt, zt, pt = x, y, z, p
        else:
            xt = min(x,y,z)
            zt = max(x,y,z)
            yt = x+y+z-xt-zt
            a = 1.0/(yt-p)
            b = a*(zt-yt)*(yt-xt)
            pt = yt+b
            rho = xt*zt/yt
            tau = p*pt/yt
            rcx = self.rc(rho,tau)

        while True:
            sqrtx, sqrty, sqrtz = sqrt(xt), sqrt(yt), sqrt(zt)
            alamb = sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
            alpha = (pt*(sqrtx+sqrty+sqrtz)+sqrtx*sqrty*sqrtz)**2
            beta = pt*(pt+alamb)**2
            summe += fac*self.rc(alpha,beta)
            fac = 0.25*fac
            xt = 0.25*(xt+alamb)
            yt = 0.25*(yt+alamb)
            zt = 0.25*(zt+alamb)
            pt = 0.25*(pt+alamb)
            ave = 0.2*(xt+yt+zt+pt+pt)
            delx=(ave-xt)/ave
            dely=(ave-yt)/ave
            delz=(ave-zt)/ave
            delp=(ave-pt)/ave
            if max(abs(delx),abs(dely),abs(delz),abs(delp)) <= ERRTOL:
                break
         
        ea = delx*(dely+delz)+dely*delz
        eb = delx*dely*delz 
        ec = delp**2
        ed = ea-3.0*ec
        ee = eb+2.0*delp*(ea-ec)
        
        C1=3/14
        C2=1/3
        C3=3/22
        C4=3/26
        C5=0.75*C3
        C6=1.5*C4
        C7=0.5*C2
        C8=C3+C3
        
        rj = 3.0*summe+fac*(1.0+ed*(-C1+C5*ed-C6*ee)+eb*(C7+delp*(-C8+delp*C4))+delp*ea*(C2-delp*C3)-C2*delp*ec)/(ave*sqrt(ave))
        if (p < 0):
            rj=a*(b*rj+3.0*(rcx-self.rf(xt,yt,zt)))
        return rj
      
    def rc(self,x,y,SQRTNY=1.3e-19,ERRTOL=0.04,TINY=1.69e-38,BIG=3e37):
        TNBG = TINY * BIG
        COMP1 = 2.236/SQRTNY
        COMP2 = TNBG*TNBG/25.0
        if (x < 0) or (y == 0) or (x+abs(y) < TINY) or (x+abs(y) > BIG) \
            or ( (y < -COMP1) and (x > 0) and (x < COMP2) ):
                raise(ValueError('invalid arguments in rc'))
        if y > 0:
            xt, yt, w = x, y, 1.
        else:
            xt = x-y
            yt = -y
            w = sqrt(x)/sqrt(xt)

        while True:
            alamb = 2.0*sqrt(xt)*sqrt(yt)+yt
            xt = 0.25*(xt+alamb)
            yt = 0.25*(yt+alamb)
            ave = (xt+yt+yt)/3
            s = (yt-ave)/ave
            if abs(s) <= ERRTOL:
                break
            
        C1 = 0.3
        C2 = 1.0/7.0
        C3 = 0.375
        C4 = 9.0/22
            
        rc = w*(1.0+s*s*(C1+s*(C2+s*(C3+s*C4))))/sqrt(ave)
        return rc

      
    def ellk(self, k):
        m1=1.0-k*k
        a0=1.38629436112
        a1=0.09666344259
        a2=0.03590092383
        a3=0.03742563713
        a4=0.01451196212
        b0=0.5
        b1=0.12498593597
        b2=0.06880248576
        b3=0.03328355346
        b4=0.00441787012
        ek1 = a0+m1*(a1+m1*(a2+m1*(a3+m1*a4)))
        ek2 = (b0+m1*(b1+m1*(b2+m1*(b3+m1*b4))))*log(m1)
        ellk = ek1-ek2
        return ellk
    
    def ellec(self, k):
        m1=1.0-k*k
        a1=0.44325141463
        a2=0.06260601220
        a3=0.04757383546
        a4=0.01736506451
        b1=0.24998368310
        b2=0.09200180037
        b3=0.04069697526
        b4=0.00526449639
        ee1 = 1.0+m1*(a1+m1*(a2+m1*(a3+m1*a4)))
        ee2 = m1*(b1+m1*(b2+m1*(b3+m1*b4)))*log(1.0/m1)
        ellec = ee1+ee2
        return ellec


      