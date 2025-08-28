        program translacional
        implicit none       
        integer nhis,i,npart, N_total
        parameter(nhis=250)
        real ri,gri,t,r,gr,dx,Lado
        real xi,x,rho1,sex,dx_ent
        real cutoff,t_cm,sex_cm
        read(*,*) Lado
        open(1,file="gr_cm.dat")
        open(2,file="tr_cm.dat")
        open(3,file="ent_cm.dat")

        t_cm = 0
        sex_cm = 0
        
        N_total = 2000
        rho1 = N_total/((2.d0*Lado)**3.d0)
        cutoff=(Lado/2.0)*rho1**(1./3.)

        read(1,*)ri,gri
        xi = ri*rho1**(1.0/3.0)
        do i=2,nhis
           read(1,*)r,gr
           x = r*rho1**(1.0/3.0)

           dx_ent = r - ri
           dx = x - xi
        
           if(gr.gt.0.D0) sex_cm = sex_cm+(gr*log(gr)-gr+1)*(r**2.D0)*dx_ent
             ri = r
           if(r.le.cutoff) t_cm = t_cm + abs(gr-1.0)*dx
        
           xi = x
           ri = r

        enddo
 
       sex_cm = -2.D0*3.14159265*rho1*sex_cm

        write(2,*) rho1, t_cm
        write(3,*) rho1, sex_cm
        
        
        endprogram
        
