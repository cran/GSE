C***********************************************************************
        subroutine wctrsc(x,w1,n,p,xbar,sdv,mvcode)
C Centers and scales the data matrix so that the observed data in every
C column have mean zero and variance one. If a column has zero variance
C or less than 2 observations then the data are centered (set equal 
C to zero)
        integer n,p
        double precision count
C real x(n,p),mvcode was replaced by line below
        double precision x(n,p),mvcode,w1(n)
        double precision xbar(p),sdv(p),sum1,sum2
        do 10 j=1,p
           sum1=0
           sum2=0
           count=0
           do 5 i=1,n
              if(x(i,j).ne.mvcode) then
                 count=count+w1(i)
                 sum1=sum1+w1(i)*x(i,j)
                 sum2=sum2+w1(i)*x(i,j)**2.
              endif
5           continue
           if(count.gt.0) then
              xbar(j)=sum1/count
              sdv(j)=sqrt((sum2-(sum1**2.)/count)/count)
              do 7 i=1,n
                 if(x(i,j).ne.mvcode) x(i,j)=x(i,j)-xbar(j)
7             continue
              if(sdv(j).gt.0.) then
                 do 9 i=1,n
                    if(x(i,j).ne.mvcode) x(i,j)=x(i,j)/sdv(j)
9                continue
               else 
                 sdv(j)=1.
              endif
           else
              sdv(j)=1.
           endif
10        continue
        return
        end
C***********************************************************************
        subroutine mkpsi(p,psi)
C Generates a symmetric matrix of integers indicating the linear
C position in packed storage of the matrix elements
        integer p,psi(0:p,0:p),posn
        posn=0
        do 10 j=0,p
           posn=posn+1
           psi(j,j)=posn
           do 5 k=j+1,p
              posn=posn+1
              psi(j,k)=posn
              psi(k,j)=posn
5          continue
10      continue
        return
        end
C***********************************************************************
        subroutine swp(d,theta,pivot,p,psi,submat,dir)
C Performs sweep on a symmetric matrix in packed storage.
C Sweeps on pivot position. Sweeps only the (0:submat,0:submat)
C submatrix.
C If dir=1, performs ordinary sweep. If dir=-1, performs reverse sweep.
        integer d,p,pivot,psi(0:p,0:p),submat,dir
        double precision theta(d),a,b,c
        a=theta(psi(pivot,pivot))
        theta(psi(pivot,pivot))=-1./a
        do 10 j=0,submat
           if(j.ne.pivot) theta(psi(j,pivot))=theta(psi(j,pivot))/a*dir
10      continue
        do 30 i=0,submat
           do 20 j=i,submat
              if((i.ne.pivot).and.(j.ne.pivot))then
                 b=theta(psi(i,pivot))
                 c=theta(psi(j,pivot))
                 theta(psi(i,j))=theta(psi(i,j))-a*b*c
              endif
20         continue
30      continue
        return
        end
C***********************************************************************
        subroutine initn(d,theta)
C Initializes theta
        integer d
        double precision theta(d)
        theta(1)=1.
        do 1 i=2,d
           theta(i)=0.
1       continue
        return
        end
C***********************************************************************
        subroutine stvaln(d,theta,p,psi)
C Gets starting value of theta: mu=0 and sigma=I
        integer d,p,psi(0:p,0:p)
        double precision theta(d)
        call initn(d,theta)
        theta(1)=-1.
        do 5 j=1,p
           theta(psi(j,j))=1.
5       continue
        return
        end
C***********************************************************************
        subroutine wtobsn(d,tobs,p,psi,n,w1n,w2n,x,w1,
     /     npatt,r,mdpst,nmdp,oc)
C Tabulates the known part of the sscp matrix for all missingness
C patterns.
        integer d,p,psi(0:p,0:p),n,npatt,r(npatt,p),nmdp(npatt)
        integer mdpst(npatt),oc(p),noc,patt
C real x(n,p) was changed to the line below
        double precision x(n,p),w1(n)
        double precision tobs(d)
        double precision w1n,w2n
        call initn(d,tobs)
        do 40 patt=1,npatt
           call gtoc(p,npatt,r,patt,oc,noc,p)
           do 30 i=mdpst(patt),mdpst(patt)+nmdp(patt)-1
              do 20 j=1,noc
                 tobs(psi(0,oc(j)))=tobs(psi(0,oc(j)))+
     /                w1(i)*x(i,oc(j))/w1n
                 do 10 k=j,noc
                    tobs(psi(oc(j),oc(k)))=tobs(psi(oc(j),oc(k)))+
     /                   w1(i)*x(i,oc(j))*x(i,oc(k))/w2n
10               continue
20            continue
30         continue
40      continue
        return
        end
C************************************************************************
        subroutine wemn(d,theta,t,tobs,p,psi,n,w1n,w2n,x,w1,w2,npatt,r,
     /     mdpst,nmdp,oc,mc,c)
C Performs one step of em. Theta must be in sweep(0) condition. 
C After execution, the new parameter value is contained in t, and
C theta is left swept on columns corresponding to observed variables
C in the first missingness pattern. 
        integer d,p,psi(0:p,0:p),n,npatt,r(npatt,p),nmdp(npatt)
        integer mdpst(npatt),oc(p),noc,mc(p),nmc,patt
C real x(n,p) was changed to the line below
        double precision x(n,p),w1(n),w2(n)
        double precision theta(d),t(d),tobs(d),c(p)
C w1n is the average weight
        double precision w1n,w2n
C copy the part corresponding to observed components
        do 1 i=1,d
           t(i)=tobs(i)
1       continue
        do 200 patt=1,npatt
           call swpobs(d,theta,p,psi,npatt,r,patt)
           call gtmc(p,npatt,r,patt,mc,nmc,p)
           call gtoc(p,npatt,r,patt,oc,noc,p)
           do 150 i=mdpst(patt),mdpst(patt)+nmdp(patt)-1
C compute c()=E(xmis_i|xobs_i,theta)
              do 50 j=1,nmc
                 c(mc(j))=theta(psi(0,mc(j)))
                 do 40 k=1,noc
                    c(mc(j))=c(mc(j))+
     /                      theta(psi(oc(k),mc(j)))*x(i,oc(k))
40               continue
50            continue
              do 100 j=1,nmc
C update location vector (corresponding to missing part)
                 t(psi(0,mc(j)))=t(psi(0,mc(j)))+w1(i)*c(mc(j))/w1n
C covariance of missing vs observed: cross-product
                 do 70 k=1,noc
                    t(psi(oc(k),mc(j)))=t(psi(oc(k),mc(j)))+
     /                   w1(i)*x(i,oc(k))*c(mc(j))/w2n
70               continue
C covariance of missing: cross-product + Ci
                 do 80 k=j,nmc
                    t(psi(mc(k),mc(j)))=t(psi(mc(k),mc(j)))+
     /                w1(i)*c(mc(k))*c(mc(j))/w2n+
     /                      w2(i)*theta(psi(mc(k),mc(j)))/w2n
80               continue
100           continue
150           continue
200     continue
        do 210 i=2,d
           t(i)=t(i)/dfloat(n)
210     continue
C 2008-10-30 by Mike. Add (mean(w)-1)*mu*mu' to XX' matrix
C 2010-05-06 Because what SWP(0) does is basically 1/n*X*X'-mu*mu'
C and we need this correction to accomodate discrepacies in the 
C denominators of mu (w1n) and Sigma (w2n).
        do 250 j=1,p
           do 230 k=j,p
              t(psi(j,k)) = t(psi(j,k)) +
     /             (1-w1n/w2n)*t(psi(0,j))*t(psi(0,k))
230           continue
250        continue
        call swp(d,t,0,p,psi,p,1)
        return
        end
C************************************************************************
        subroutine lobsn(d,theta,t,p,psi,n,x,npatt,r,mdpst,nmdp,oc,
     /    c,loglik)
C Evaluates observed-data loglikelihood at theta
        integer d,p,psi(0:p,0:p),n,npatt,r(npatt,p),nmdp(npatt)
        integer mdpst(npatt),oc(p),noc,patt
C real x(n,p) was changed to the line below
        double precision x(n,p)
        double precision theta(d),t(d),c(p),loglik,logdet,trace
        loglik=dfloat(0)
        logdet=dfloat(0)  
        do 5 j=1,p
           c(j)=theta(psi(0,j))
5       continue
        do 200 patt=1,npatt
           call initn(d,t)
           do 10 j=1,p
              if((r(patt,j).eq.1).and.(theta(psi(j,j)).gt.0.))then
                 logdet=logdet+log(theta(psi(j,j)))
                 call swp(d,theta,j,p,psi,p,1)
              elseif((r(patt,j).eq.0).and.(theta(psi(j,j)).lt.0.))then
                 call swp(d,theta,j,p,psi,p,-1)
                 logdet=logdet-log(theta(psi(j,j)))
              endif
10         continue
           call gtoc(p,npatt,r,patt,oc,noc,p)
           do 100 i=mdpst(patt),mdpst(patt)+nmdp(patt)-1
              do 30 j=1,noc
                 t(psi(0,oc(j)))=x(i,oc(j))-c(oc(j))
30            continue
              do 50 j=1,noc
                 do 40 k=j,noc
                    t(psi(oc(j),oc(k)))=t(psi(oc(j),oc(k)))+
     /                 t(psi(0,oc(j)))*t(psi(0,oc(k)))
40               continue
50            continue
100        continue
           trace=dfloat(0)
           do 120 j=1,noc
              do 110 k=1,noc
                 trace=trace-(theta(psi(oc(j),oc(k)))*
     /              t(psi(oc(j),oc(k))))
110           continue
120        continue
           loglik=loglik-(dfloat(nmdp(patt))*logdet/2.)-(trace/2.)
200     continue
        return
        end
C************************************************************************
        subroutine swpobs(d,theta,p,psi,npatt,r,patt)
C Sweeps theta to condition on the observed variables
        integer d,p,psi(0:p,0:p),npatt,r(npatt,p),patt
        double precision theta(d)
        do 10 j=1,p
           if((r(patt,j).eq.1).and.(theta(psi(j,j)).gt.0.))then
              call swp(d,theta,j,p,psi,p,1)
           elseif((r(patt,j).eq.0).and.(theta(psi(j,j)).lt.0.))then
              call swp(d,theta,j,p,psi,p,-1)
           endif
10      continue
        return
        end
C************************************************************************
        subroutine gtmc(p,npatt,r,patt,mc,nmc,last)
C Finds the column numbers of the missing variables, and stores them
C in the first nmc elements of mc. Does not go beyond column=last.
        integer p,npatt,r(npatt,p),patt,mc(p),nmc,last
        nmc=0
        do 10 j=1,last
           if(r(patt,j).eq.0)then
              nmc=nmc+1
              mc(nmc)=j
           endif
10      continue
        return
        end
C************************************************************************
        subroutine gtoc(p,npatt,r,patt,oc,noc,last)
C Finds the column numbers of the observed variables, and stores them
C in the first noc elements of oc. Does not go beyond column=last.
        integer p,npatt,r(npatt,p),patt,oc(p),noc,last
        noc=0
        do 10 j=1,last
           if(r(patt,j).eq.1)then
              noc=noc+1
              oc(noc)=j
           endif
10      continue
        return
        end
C***********************************************************************
        subroutine sjn(p,npatt,r,sj)
C computes s_j, the number of the last missingness pattern for which
C the jth variable needs to be imputed to complete the monotone pattern
        integer p,npatt,r(npatt,p),sj(p),patt,tmp
        do 30 j=1,p
           patt=npatt+1
 10        continue
           patt=patt-1
           if(patt.ge.1)then
              if(r(patt,j).eq.0)goto 10
           endif
           sj(j)=patt
 30     continue
        tmp=sj(p)
        do 40 j=p-1,1,-1
           sj(j)=max0(sj(j),tmp)
           tmp=sj(j)
 40     continue
        return
        end
C***********************************************************************
        subroutine nmons(p,npatt,nmdp,sj,nmon)
C Computes the number of observations in (Xobs,Xmis*) for each column
C (called N_j in Figure 6.11 of book)
        integer p,npatt,nmdp(npatt),sj(p),nmon(p),patt
        do 40 j=1,p
           nmon(j)=0
           do 20 patt=1,sj(j)
              nmon(j)=nmon(j)+nmdp(patt)
20         continue
40      continue
        return
        end
C***********************************************************************
        subroutine lasts(p,npatt,sj,last)
C Finds last variable in each missingness pattern
C to complete a monotone pattern
        integer p,npatt,sj(p),last(npatt),patt,start
        do 50 j=p,1,-1
           if(j.eq.p)then
             start=1
           else
              start=sj(j+1)+1
           endif
           do 40 patt=start,sj(j)
              last(patt)=j
40         continue
50      continue
        return
        end
C***********************************************************************
        subroutine layers(p,sj,layer,nlayer)
C Finds layers for observed parts of the sufficient statistics
        integer p,sj(p),layer(p),nlayer
        nlayer=0
        do 10 j=p,1,-1
           if(j.eq.p)then
              if(sj(j).gt.0) nlayer=nlayer+1
           else
              if(sj(j).gt.sj(j+1)) nlayer=nlayer+1
           endif
           layer(j)=nlayer
10      continue
        return
        end
C***********************************************************************
