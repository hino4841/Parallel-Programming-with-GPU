* Short description of the code *************************************************************************
*  
*  A set of "Number" vortons is generated randomly within the box [0,1]*[0,1]*[0,1]
*  Each vorton has a strength "Omega_v" and a radius "Radius"
*  Strength is a random vector with components r_x,r_y,r_z.  Each component is a random number in the range [-0.5,0.5]
*  
*  The initially generated vortons move according to the Computational Vortex Method procedure "Ntime" time steps
*  The time step "Delta_t" is constant 
*  All vortons have the same radius "Radius"
*
*  Output:
*  After simulation we print the distribution of vortons and their strength
*  The velocity and vorticity fields in the plane x=0.5, 0<y<1 , 0<z<1  are also printed
*  The total time of simulations is also printed
*
*  Written by Nikolai Kornev on 14 January 2013
*  Last change 24 May 2017
*********************************************************************************************************

******************************************* Description of variables ************************************

      real,allocatable:: 
     &Vortex(:,:),                                                 ! coordinates of vortons
     &Omega_v(:,:),                                                ! strength vector of vortons
     &VortexN(:,:),                                                ! coordinates of vortons (Working field)
     &Omega_vN(:,:),                                               ! strength vector of vortons  (Working field)
     &Sigma(:)                                                     ! radius of vortons  (Working field)
      
	Real StatisticalMoments(4)	               ! New Line 1

      open(774,file='Velocities.dat')

      open(773,file='MaxValue.dat')
*******************************************  Input data   Begin  ***********************************************

	Ntime=10000	                               ! Number of time steps
	Delta_t=0.01					   ! Time step
	Radius=0.1								   ! Radius of vorton (the same for all)
	Number=1000								   ! number of vortons
      V_mean=1.0
 
      Do i_order=1,4						! New Line	2
      StatisticalMoments(i_order)=0.000	! New Line	3
	End Do								! New Line	4

*******************************************  Input data  End  ***********************************************



* description of variables **********************************************************************************
      allocate(Vortex(Number,3),VortexN(Number,3),Sigma(Number),
	&Omega_v(Number,3),Omega_vN(Number,3),stat=ierr)

* description of variables ***********************************************************************************
      

* generation of initial distribution of vortons Begin

      Do ivorton=1,Number

	call Random(xxx)		   ! Any generator of random numbers Uniform within the range [0,1]
	call Random(yyy)
	call Random(zzz)
	Vortex(ivorton,1)=xxx
	Vortex(ivorton,2)=yyy
	Vortex(ivorton,3)=zzz

 	call Random(Omega_x)
	call Random(Omega_y)
	call Random(Omega_z)
	Omega_v(ivorton,1)=(Omega_x-0.5)
	Omega_v(ivorton,2)=(Omega_y-0.5)
	Omega_v(ivorton,3)=(Omega_z-0.5)
	Sigma  (ivorton)  = Radius
	End Do
* generation of initial distribution of vortons End

        time0=0
        call system_clock(jcount1,jcount_rate1,jcount_max1)
        time0=(jcount1+0.0)/(jcount_rate1+0.0)
 
      do itime=1,Ntime                               ! loop over time steps
     	write(6,*)itime,Amagni,Energy,Ncout

      Do ivorton=1,Number		 ! this loop can be done parallel
*
*       calculation of the velocity and tensor S_ij induced at vorton with number "ivorton"
	  Vxc=V_mean
	  Vyc=0.0
	  Vzc=0.0
	      dvxdxmov=0.0
	      dvxdymov=0.0
	      dvxdzmov=0.0 	 
	      dvydxmov=0.0 
	      dvydymov=0.0
	      dvydzmov=0.0	
	      dvzdxmov=0.0
	      dvzdymov=0.0
	      dvzdzmov=0.0
	    Do induced=1,Number	                          ! this loop can be done parallel
	      vxx=Vortex(ivorton,1)-Vortex(induced,1)
	      vyy=Vortex(ivorton,2)-Vortex(induced,2)
	      vzz=Vortex(ivorton,3)-Vortex(induced,3)
            radiika=vxx*vxx+vyy*vyy+vzz*vzz
	      t1=vyy*	Omega_v(induced,3)- vzz*Omega_v(induced,2)
	      t2=vzz*	Omega_v(induced,1)- vxx*Omega_v(induced,3)
	      t3=vxx*	Omega_v(induced,2)- vyy*Omega_v(induced,1)
	Om22P=3.1416/Sigma(induced)/Sigma(induced)/2.0

	      ssss=Exp(-radiika*Om22P)
*        Here are the velocities
            Vxc           = Vxc+ssss*t1
            Vyc           = Vyc+ssss*t2
            Vzc           = Vzc+ssss*t3

*       Here are the strain rate tensor S_ij

		  dssss_dr=(-Om22P)*ssss

	      dvxdxmov=dssss_dr*vxx*t1                         +dvxdxmov
	      dvxdymov=dssss_dr*vyy*t1+Omega_v(induced,3)*ssss +dvxdymov 
	      dvxdzmov=dssss_dr*vzz*t1-Omega_v(induced,2)*ssss +dvxdzmov 
	 
	      dvydxmov=dssss_dr*vxx*t2-Omega_v(induced,3)*ssss +dvydxmov 
	      dvydymov=dssss_dr*vyy*t2                         +dvydymov 
	      dvydzmov=dssss_dr*vzz*t2+Omega_v(induced,1)*ssss +dvydzmov
	
	      dvzdxmov=dssss_dr*vxx*t3+Omega_v(induced,2)*ssss +dvzdxmov 
	      dvzdymov=dssss_dr*vyy*t3-Omega_v(induced,1)*ssss +dvzdymov 
	      dvzdzmov=dssss_dr*vzz*t3                         +dvzdzmov 
	    End do

*        new coordinates of the vorton "ivorton"	calculated using explicit Euler method
        VortexN(ivorton,1)=Vortex(ivorton,1)+Delta_t*Vxc
        VortexN(ivorton,2)=Vortex(ivorton,2)+Delta_t*Vyc
        VortexN(ivorton,3)=Vortex(ivorton,3)+Delta_t*Vzc

*        new strengths of the vorton "ivorton"  calculated using explicit Euler method
	  domxdt=dvxdxmov*Omega_v(ivorton,1)+dvxdymov*Omega_v(ivorton,2)+
	&         dvxdzmov*Omega_v(ivorton,3)
	  domydt=dvydxmov*Omega_v(ivorton,1)+dvydymov*Omega_v(ivorton,2)+
	&         dvydzmov*Omega_v(ivorton,3)	
        domzdt=dvzdxmov*Omega_v(ivorton,1)+dvzdymov*Omega_v(ivorton,2)+
	&         dvzdzmov*Omega_v(ivorton,3)
	  Omega_vN(ivorton,1)=Omega_v(ivorton,1)+domxdt*Delta_t
	  Omega_vN(ivorton,2)=Omega_v(ivorton,2)+domydt*Delta_t
	  Omega_vN(ivorton,3)=Omega_v(ivorton,3)+domzdt*Delta_t
      End do

* mapping to the cube back------------------------------------------------------------------
      Ncout=0
	    Do ivorton=1,Number
	Replace=0.0
	do kkk=1,3
      If(VortexN(ivorton,kkk).lt.0.0)Replace=1.0
      If(VortexN(ivorton,kkk).gt.1.0)Replace=1.0
	end do
	If(Replace.eq.1.0)then
      Ncout=Ncout+1
	call Random(xxx)		   
	call Random(yyy)
	call Random(zzz)
	VortexN(ivorton,1)=xxx
	VortexN(ivorton,2)=yyy
	VortexN(ivorton,3)=zzz

 	call Random(Omega_x)
	call Random(Omega_y)
	call Random(Omega_z)
	Omega_vN(ivorton,1)=(Omega_x-0.5)
	Omega_vN(ivorton,2)=(Omega_y-0.5)
	Omega_vN(ivorton,3)=(Omega_z-0.5)
	Sigma  (ivorton)  = Radius
	End if
	    End do
* mapping to the cube back------------------------------------------------------------------

* the old parameters became new ones
      Amagni=0.0
	    Do ivorton=1,Number
          Vortex(ivorton,1)=VortexN(ivorton,1)
          Vortex(ivorton,2)=VortexN(ivorton,2)
          Vortex(ivorton,3)=VortexN(ivorton,3)
          Amagnit_old=Sqrt(
	&Omega_v(ivorton,1)**2+Omega_v(ivorton,2)**2+
     &Omega_v(ivorton,3)**2)
	    Omega_v(ivorton,1)=Omega_vN(ivorton,1)			 
	    Omega_v(ivorton,2)=Omega_vN(ivorton,2)			 
	    Omega_v(ivorton,3)=Omega_vN(ivorton,3)			 
          Amagnit_new=Sqrt(
	&Omega_v(ivorton,1)**2+Omega_v(ivorton,2)**2+
     &Omega_v(ivorton,3)**2)
	Sigma(ivorton)=Sigma(ivorton)*Sqrt(Amagnit_old/Amagnit_new)  

      If(Amagnit_new.ge.Amagni)then
	Amagni=Amagnit_new
	Energy=(Amagnit_new**2)*(Sigma(ivorton)**5)        
	Speed_max=Amagnit_new*Sigma(ivorton)
	Sigmas=Sigma(ivorton)
	end if
	    End do


* we print the maximum vorticity energy velocity radius
	write(773,773)itime*Delta_t,Amagni,Energy,Speed_max,Sigmas
773   format(2f12.4,f12.7,2f12.4)


*eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee
*       calculation of the velocity induced at the centre of the computational domain
          Vx=0.0																							 
	    Do induced=1,Number	                          
	      vxx=0.5-Vortex(induced,1)																		 
	      vyy=0.5-Vortex(induced,2)																		
	      vzz=0.5-Vortex(induced,3)																		 
            radiika=vxx*vxx+vyy*vyy+vzz*vzz																 
	      t1=vyy*	Omega_v(induced,3)- vzz*Omega_v(induced,2)
	Om22P=3.1416/Sigma(induced)/Sigma(induced)/2.0											 
	      ssss=Exp(-radiika*Om22P)																	 
		 Vx=Vx+ssss*t1																	             
	    End do																							 

	write(774,773)itime*Delta_t,Vx

      Do ier=1,4						                                                                     
      StatisticalMoments(ier)=StatisticalMoments(ier)+Vx**ier	                                             
	End Do							             	     												 
*eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee


	End do                                            ! end of time loop

*      From this it is just printing results------------------------------------------------------
  	close(773)
  	close(774)

      call system_clock(jcount1,jcount_rate1,jcount_max1)                                  
      time1=(jcount1+0.0)/(jcount_rate1+0.0)                                                                                                    
                                                         
      write(6,*)' It took',(time1-time0)/60.,' minutes'                                    


* This part is not important This is just printing
* Print Print	 Begin
      open(1,file='coordinates.dat')
      open(2,file='strengths.dat')
	Do ivorton=1,Number
      write(1,*)Vortex(ivorton,1),Vortex(ivorton,2),Vortex(ivorton,3)
      write(2,*)Omega_v(ivorton,1),Omega_v(ivorton,2),Omega_v(ivorton,3)
	End do
	close(1)
	close(2)
 
      open(1,file='sliceVelocity.dat')
      open(2,file='sliceVorticity.dat')
      NNN=100
	Delta=1./(NNN-1.)
	xxx=0.5
	do iy=1,NNN
	do iz=1,NNN
	yyy=(iy-1.)*Delta
	zzz=(iz-1.)*Delta

 	  Vxc=0.0
	  Vyc=0.0
	  Vzc=0.0
	      dvxdymov=0.0
	      dvxdzmov=0.0 	 
	      dvydxmov=0.0 
	      dvydzmov=0.0	
	      dvzdxmov=0.0
	      dvzdymov=0.0

	    Do induced=1,Number
	      vxx=xxx-Vortex(induced,1)
	      vyy=yyy-Vortex(induced,2)
	      vzz=zzz-Vortex(induced,3)
            radiika=vxx*vxx+vyy*vyy+vzz*vzz
	      t1=vyy*	Omega_v(induced,3)- vzz*Omega_v(induced,2)
	      t2=vzz*	Omega_v(induced,1)- vxx*Omega_v(induced,3)
	      t3=vxx*	Omega_v(induced,2)- vyy*Omega_v(induced,1)
	Om22P=3.1416/Sigma(induced)/Sigma(induced)/2.0
	      ssss=Exp(-radiika*Om22P/2.)
*        Here are the velocities
            Vxc           = Vxc+ssss*t1
            Vyc           = Vyc+ssss*t2
            Vzc           = Vzc+ssss*t3

*       Here are the strain rate tensor S_ij
		  dssss_dr=(-Om22P)*Exp(-radiika*Om22P/2.)

	      dvxdymov=dssss_dr*vyy*t1+Omega_v(induced,3)*ssss +dvxdymov 
	      dvxdzmov=dssss_dr*vzz*t1-Omega_v(induced,2)*ssss +dvxdzmov 
	 
	      dvydxmov=dssss_dr*vxx*t2-Omega_v(induced,3)*ssss +dvydxmov 
	      dvydzmov=dssss_dr*vzz*t2+Omega_v(induced,1)*ssss +dvydzmov
	
	      dvzdxmov=dssss_dr*vxx*t3+Omega_v(induced,2)*ssss +dvzdxmov 
	      dvzdymov=dssss_dr*vyy*t3-Omega_v(induced,1)*ssss +dvzdymov 
		  
	    end do 
	write(1,1)yyy,zzz,Vxc,Vyc,Vzc
	write(2,1)yyy,zzz,dvzdymov-dvydzmov,dvxdzmov-dvzdxmov,
	&dvydxmov-dvxdymov
	end do
	end do
	close(1)
	close(2)
1     format(5f12.4)

	Umean=StatisticalMoments(1)/Ntime						                    
	Variance=Sqrt(StatisticalMoments(2)/Ntime)				                    
	Skew=(StatisticalMoments(3)/Ntime)/(Variance**3)		                    
	Curt=(StatisticalMoments(4)/Ntime)/(Variance**4)		                    
      write(6,*)Umean,Variance,Skew,Curt	                                        
           
* Print Print	 End
6677  continue
      deallocate(Vortex,VortexN,Omega_v,Omega_vN,Sigma,stat=ierr)
               
			 call frank(Ntime)  


	end 

              subroutine frank(Ntime)  
              real f_mix(120000)

      open(1,file='Velocities.dat')

	sred=0.0
	do ipa=1,Ntime
	read(1,*)pind,f_mix(ipa)
	sred=sred+f_mix(ipa)
  	end do
	sred=sred/(Ntime+0.000)
	write(6,*)sred
	close(1)

* Berechnung von Schwankungen
 	do ipa=1,Ntime
 	f_mix(ipa)=f_mix(ipa)-sred
  	end do

              pi2N =2.*3.1416/Ntime
              kshag=3
              Nbasa=kshag
              K_min=0
              K_max=Ntime/2-1
              Call Spectrum(Ntime,f_mix,K_min,K_max,Pi2N,kshag,Nbasa)
              return
	end

       SUBROUTINE Spectrum(N,x,K_min,K_max,Pi2N,kshag,Nbasa)
       REAL x(N),sever(0:K_max)
        open(23,file='spectrum.dat')
        open(13,file='spectrum_sm.dat')        
        do iii=0,K_max
        sever(iii)=0.0
        end do
        k=K_min
 1      continue
        summa=0.
        summa1=0.
        do j=1,N
        Summa=Summa+cos(pi2N*(j-1)*k)*x(j)   
        Summa1=Summa1+sin(pi2N*(j-1)*k)*x(j)
        end do
        sever(k)=(summa*summa+summa1*summa1)*2.0/N
        write(23,*)k,sever(k)
        k=k+1
        if(k.le.K_max)go to 1
        close(1)
* averaging     	
        do k=K_min+Nbasa+1,K_max-Nbasa,kshag
        sevv=0.0
        N_basa=0
        do kel=-Nbasa,Nbasa,kshag
        sevv=sevv+sever(k+kel)
        N_basa=N_basa+1
        end do
        write(13,*)k,sevv/(N_basa)
        end do  
        close(13)
        close(23)
        RETURN
        END
