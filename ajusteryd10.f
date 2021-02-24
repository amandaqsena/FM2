      program minimo

      implicit double precision (a-h,o-z)
      common /const/de,ef,x0,we

      dimension xi(12,12)

c     x0 is the equilibrium distance
      write(*,*) 'entre com o valor de Re'
      read(*,*) x0
C
      write(*,*) 'entre o numero de parametros a serem minimizados'
      read(*,*) npar 
      call driver(npar,xi)
      end

      subroutine driver(nknots1,xi)
      implicit double precision (a-h,o-z)

      dimension a(12),xi(nknots1,nknots1)

      common /idata/ndat
      common /const/de,ef,x0,we
c     common /pesos/wght

      include 'sysm.cm'
      character qid*20
c     call pregs
c     call secnd(io)
      nknots = nknots1
      icc = 0
      write(*,*) 'entre com o numero de distancias e energias'
      read(*,*) ndat
      write(*,*) 'entre com o nome do arquivo das energias'
      read(*,1) qid
c     write(*,*)'valor de qid eh igual a', qid
      write(*,*) 'entre com o numero: 0(aleatorio),1(file)'
      read(*,*) ifnots
c      write(*,*)'valor de ifnots eh iqaul a', ifnots
      write(*,*) 'entre com a tolerancia'
      read(*,*) tol
c      write(*,*) 'tol=',tol
1     format(a20)
      ln = index(qid,' ') - 1
      open(1,file=qid(1:ln),status='unknown')
      read(1,*) (xo(j),y(j),j=1,ndat)
      close(1)

      if(ifnots.eq.1) then
         open(10,file='par10.in',status='unknown')
         write(*,*)'nknots',nknots
         read(10,*) (a(i),i=1,nknots)
         close(10)
      endif
      do i = 1, nknots
         xi(i,i) = 1.d0
      enddo
      call powell(a,xi,nknots,nknots,tol,iter,fret)
      u = 0.5*float(ndat-nknots)
      open(4,file='finalryd10',status='unknown')
           write(4,40) (a(i),i=1,nknots)
           write(4,40) x0
40    format(E18.12)
      close(4)
      end


      function func(p)
      implicit double precision (a-h,o-z)

c     the include sysm.cm contains common /engine/xo(120),y(120)

      include 'sysm.cm'

      dimension p(12),c(12)

      common /idata/ndat
      common /pesos/wght
      common /const/de,ef,x0,we
C
      sum2 = 0.0d0
C
      c1=abs(p(1))
      c2=p(2)
      c3=p(3)
      c4=p(4)
      c5=p(5)
      c6=p(6)
      c7=p(7)
      c8=p(8)
      c9=p(9)
      c10=p(10)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC    Dissociation energy
      de=0.1020554739d0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      ef=0.0d0

      do l=1,ndat
      x = xo(l)-x0
      eta=exp(-c1*x)
      XU=-de*(1.0d0+c1*x+c2*x**2+c3*x**3+c4*x**4+
     >c5*x**5+c6*x**6+c7*x**7+c8*x**8+c9*x**9+c10*x**10)*eta
      fx=XU+ef
      sum2 = sum2 + ((fx-y(l)))**2.0d0
      enddo
      dat = ndat
      func =dsqrt(sum2/dat)
      write(*,*)'dqm=',func
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FUNCTION brent(ax,bx,cx,f,tol,xmin)
      implicit double precision (a-h,o-z)
      INTEGER ITMAX
      double precision brent,ax,bx,cx,tol,xmin,f,CGOLD,ZEPS
      EXTERNAL f
      PARAMETER (ITMAX=100,CGOLD=.3819660,ZEPS=1.0e-10)
      INTEGER iter
      double precision a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2
     +,u,v,w,x,xm
      a=min(ax,cx)
      b=max(ax,cx)
      v=bx
      w=v
      x=v
      e=0.
      fx=f(x)
      fv=fx
      fw=fx
      do 11 iter=1,ITMAX
        xm=0.5*(a+b)
        tol1=tol*abs(x)+ZEPS
        tol2=2.*tol1
        if(abs(x-xm).le.(tol2-.5*(b-a))) goto 3
        if(abs(e).gt.tol1) then
          r=(x-w)*(fx-fv)
          q=(x-v)*(fx-fw)
          p=(x-v)*q-(x-w)*r
          q=2.*(q-r)
          if(q.gt.0.) p=-p
          q=abs(q)
          etemp=e
          e=d
          if(abs(p).ge.abs(.5*q*etemp).or.p.le.q*(a-x).or.p.ge.q*(b-x)) 
     *goto 1
          d=p/q
          u=x+d
          if(u-a.lt.tol2 .or. b-u.lt.tol2) d=sign(tol1,xm-x)
          goto 2
        endif
1       if(x.ge.xm) then
          e=a-x
        else
          e=b-x
        endif
        d=CGOLD*e
2       if(abs(d).ge.tol1) then
          u=x+d
        else
          u=x+sign(tol1,d)
        endif
        fu=f(u)
        if(fu.le.fx) then
          if(u.ge.x) then
            a=x
          else
            b=x
          endif
          v=w
          fv=fw
          w=x
          fw=fx
          x=u
          fx=fu
        else
          if(u.lt.x) then
            a=u
          else
            b=u
          endif
          if(fu.le.fw .or. w.eq.x) then
            v=w
            fv=fw
            w=u
            fw=fu
          else if(fu.le.fv .or. v.eq.x .or. v.eq.w) then
            v=u
            fv=fu
          endif
        endif
11    continue
      pause 'brent exceed maximum iterations'
3     xmin=x
      brent=fx
      return
      END
      FUNCTION f1dim(x)
      implicit double precision (a-h,o-z)
      INTEGER NMAX
      double precision f1dim,func,x
      PARAMETER (NMAX=100)
CU    USES func
      INTEGER j,ncom
      double precision pcom(NMAX),xicom(NMAX),xt(NMAX)
      COMMON /f1com/ pcom,xicom,ncom
      do 11 j=1,ncom
        xt(j)=pcom(j)+x*xicom(j)
11    continue
      f1dim=func(xt)
      return
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE linmin(p,xi,n,fret)
      implicit double precision (a-h,o-z)
      INTEGER n,NMAX
      double precision fret,p(n),xi(n),TOL
      PARAMETER (NMAX=100,TOL=1.e-4)
CU    USES brent,f1dim,mnbrak
      INTEGER j,ncom
      double precision ax,bx,fa,fb,fx,xmin,xx,pcom(NMAX),xicom(NMAX)
     +,brent
      COMMON /f1com/ pcom,xicom,ncom
      EXTERNAL f1dim
      ncom=n
      do 11 j=1,n
        pcom(j)=p(j)
        xicom(j)=xi(j)
11    continue
      ax=0.
      xx=1.
      call mnbrak(ax,xx,bx,fa,fx,fb,f1dim)
      fret=brent(ax,xx,bx,f1dim,TOL,xmin)
      do 12 j=1,n
        xi(j)=xmin*xi(j)
        p(j)=p(j)+xi(j)
12    continue
      return
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE mnbrak(ax,bx,cx,fa,fb,fc,func)
      implicit double precision (a-h,o-z)
      double precision ax,bx,cx,fa,fb,fc,func,GOLD,GLIMIT,TINY
      EXTERNAL func
      PARAMETER (GOLD=1.618034, GLIMIT=100., TINY=1.e-20)
      double precision dum,fu,q,r,u,ulim
      fa=func(ax)
      fb=func(bx)
      if(fb.gt.fa)then
        dum=ax
        ax=bx
        bx=dum
        dum=fb
        fb=fa
        fa=dum
      endif
      cx=bx+GOLD*(bx-ax)
      fc=func(cx)
1     if(fb.ge.fc)then
        r=(bx-ax)*(fb-fc)
        q=(bx-cx)*(fb-fa)
        u=bx-((bx-cx)*q-(bx-ax)*r)/(2.*sign(max(abs(q-r),TINY),q-r))
        ulim=bx+GLIMIT*(cx-bx)
        if((bx-u)*(u-cx).gt.0.)then
          fu=func(u)
          if(fu.lt.fc)then
            ax=bx
            fa=fb
            bx=u
            fb=fu
            return
          else if(fu.gt.fb)then
            cx=u
            fc=fu
            return
          endif
          u=cx+GOLD*(cx-bx)
          fu=func(u)
        else if((cx-u)*(u-ulim).gt.0.)then
          fu=func(u)
          if(fu.lt.fc)then
            bx=cx
            cx=u
            u=cx+GOLD*(cx-bx)
            fb=fc
            fc=fu
            fu=func(u)
          endif
        else if((u-ulim)*(ulim-cx).ge.0.)then
          u=ulim
          fu=func(u)
        else
          u=cx+GOLD*(cx-bx)
          fu=func(u)
        endif
        ax=bx
        bx=cx
        cx=u
        fa=fb
        fb=fc
        fc=fu
        goto 1
      endif
      return
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE powell(p,xi,n,np,ftol,iter,fret)
      implicit double precision (a-h,o-z)
      INTEGER iter,n,np,NMAX,ITMAX
c     REAL fret,ftol,p(np),xi(np,np),func
      double precision fret,ftol,p(np),xi(np,np),func
      EXTERNAL func
c+      PARAMETER (NMAX=100,ITMAX=200)
      PARAMETER (NMAX=500000,ITMAX=500000)
CU    USES func,linmin
      INTEGER i,ibig,j
      double precision del,fp,fptt,t,pt(NMAX),ptt(NMAX),xit(NMAX)
      fret=func(p)
      do 11 j=1,n
        pt(j)=p(j)
11    continue
      iter=0
1     iter=iter+1
      fp=fret
      ibig=0
      del=0.
      do 13 i=1,n
        do 12 j=1,n
          xit(j)=xi(j,i)
12      continue
        fptt=fret
        call linmin(p,xit,n,fret)
        if(abs(fptt-fret).gt.del)then
          del=abs(fptt-fret)
          ibig=i
        endif
13    continue
      write(*,*) fret
      if(2.*abs(fp-fret).le.ftol*(abs(fp)+abs(fret)))return
      if(iter.eq.ITMAX) pause 'powell exceeding maximum iterations'
      do 14 j=1,n
        ptt(j)=2.*p(j)-pt(j)
        xit(j)=p(j)-pt(j)
        pt(j)=p(j)
14    continue
      fptt=func(ptt)
      if(fptt.ge.fp)goto 1
      t=2.*(fp-2.*fret+fptt)*(fp-fret-del)**2-del*(fp-fptt)**2
      if(t.ge.0.)goto 1
      call linmin(p,xit,n,fret)
      do 15 j=1,n
        xi(j,ibig)=xi(j,n)
        xi(j,n)=xit(j)
15    continue
      goto 1
      END
