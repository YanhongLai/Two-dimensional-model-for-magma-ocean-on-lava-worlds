clc;clear;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %Set Parameters
  t0=0; nt0 = 1;
  dt = 2000;     
  dt_mom = dt;  num_dt = 1; num2_dt = 10;
  dt_temp = dt*num_dt; 
  dt_eta = dt*num2_dt;
  
  nt1 = 250E5;
  dtsave = 2E3;          %save once every 10 timesteps
  dt_output=5E4;
  nt = nt1/dtsave + 1;   %+1 plus the initial value
  
  kappa = 1E3;  
  kh = kappa;    kz = 1E-4;   
  kh_GM = kappa;
  kz_convec = 1E1;
  
  vis_liq = 1E4;   vis_sol = 1E8;
  Ah_liq = vis_liq;    Ah_sol = vis_sol; 
  Az_liq = 1E2;    Az_sol = vis_sol; 
  
  
  T_liq = 2000;     T_sol = 1700;
 
  tau = 10*86400;
  bottomDragLinear = 1E3;
  
  %%%%Initial values
     path1='Snapshot20000000.nc';
     t0=200E5; 
     nt=(nt1-t0)/dtsave+1;
     nt0=t0/dtsave+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Set grids and resolutions

  nx = 101;
  x1 = -18000000.0; x2 = 18000000.0;
  dx = ( x2 - x1 ) / ( nx - 1 );
  x = ( linspace ( x1, x2, nx ) )';
  xa=x./x2*110;
  x22=x2/110*90;
  
  %   nz=101;
  %   z1 = 0.0; z2 = -2000.0;   
  %   dz = ( z2 - z1 ) / ( nz - 1 ); %dz is negative
  %   z = ( linspace ( z1, z2, nz ) )';

  %%%%%%%%%------------5 km deep-------vertical-----
%  delRc= [10, 10, 10, 20, 20, 20, 30, 30, 30, 40, 40, 40, 50, 50, 50, 60, 60, 70, 70, 80, 80, 90, 90, ...
%    100, 100, 100, 100, 100, 100, 150, 150, 150, 200, 200, 200, 200, 200, 300, 300, 300, 300, 300, 400];
%  delRc= [20, 20, 20, 20, 20, 20, 30, 30, 40, 40, 40, 50, 50, 50, 60, 60, 70, 70, 80, 80, 90, 90, ...
%    100, 100, 100, 100, 100, 100, 150, 150, 150, 200, 200, 200, 200, 200, 300, 300, 300, 300, 300, 400];
  delRc= [40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,...
    40, 40, 40, 40, 40, 60, 60, 60, 60, 60, 100, 100, 100, 200, 200, 300, 300, 300, 300, 400, 400, 400, 400];

  delRc=-delRc;
  nz1=max(size(delRc));
  nz=nz1+1;
  z=zeros(nz,1);
for i=2:nz
    z(i)=z(i-1)+delRc(i-1);
end
% H1=sum(delRc);

nxz=nx*nz;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Initialization
    % Setup numbering
  num = 1; Number=zeros(nz,nx);
  for i=1:nz
      for j=1:nx
          Number(i,j)=num;
          num=num+1;
      end
  end

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%Determine the initial state
 
 Tsurf = 3000;   Tbot = 50; 
 Twest = 50;   Teast = 50;

 Tbot2=T_liq;
  Tvert = zeros ( nz, 1);Tx = zeros ( nx, 1);
 for i=1:nz
     Tvert(i)=Tsurf-((Tsurf-Tbot2)/(nz-1))*(i-1); 
 end
%  for j=1:nx
%      Tx(j)=Tbot+((Tsurf-Tbot)*cos((abs(x(j))/x2)*pi/2)); 
%  end

for j=1:nx
     if xa(j)>=-85
         Tc = Tbot+((Tsurf-Tbot)*((cos((abs(x(j))/x22)*pi/2)^(1/4))));
         xc = x(j);
         xac = xa(j);
         break;
     end
 end    
   for j=1:nx
     if -85<=xa(j) && xa(j)<=85
         Tx(j) = Tbot+((Tsurf-Tbot)*((cos((abs(x(j))/x22)*pi/2)^(1/4))));
     elseif xa(j)>=110 || xa(j)<=-110
        Tx(j) = Tbot;
     elseif xa(j)>85 && xa(j)<110
        Tx(j) = Tbot+(Tc-Tbot)/(-xac-110)*(xa(j)-110);
     elseif xa(j)>-110 && xa(j)<-85
        Tx(j) = Tbot+(Tc-Tbot)/(xac+110)*(xa(j)+110);
     end     
   end
 
   %Steady wind stress
 taux=zeros(nx,1); tau0=0.0;
% for j=1:nx
%    taux(j)=tau0*sin(x(j)/x2*pi);
% end
% fid=fopen('taux_nx101_100Earth.bin','rb');
% fid=fopen('taux_nx51_SiO.bin','rb');
  fid=fopen('taux_nx101_Kepler-10b_CD1E-2.bin','rb');
  taux = fread(fid,[nx 1],'double');
  fclose(fid);
% taux=taux';

 usurf = 1;    ubot = 0;
 uwest = 0;    ueast = 0; 

T0 = zeros ( nx*nz, 1 );  
for i=1:nz
    for j=1:nx
        ii = Number(i,j);
%         T0(ii,1) = Tx(j);
        T0(ii,1) = Tbot;
     end
end
  
 T2=T0;  u2=zeros(nxz,1);    w2=zeros(nxz,1);
 eta2=zeros(nx,1);
 kh2=zeros(nxz,1);     kh2(:,1)=kh;
 kz2=zeros(nxz,1);     kz2(:,1)=kz;
 
   % 1D to 2D
 u3=zeros(nz,nx,nt);w3=zeros(nz,nx,nt);T3=zeros(nz,nx,nt);
 kh3=zeros(nz,nx,nt);kz3=zeros(nz,nx,nt);
 Ah3=zeros(nz,nx,nt);Az3=zeros(nz,nx,nt);
 rho3=zeros(nz,nx,nt);    pres3=zeros(nz,nx,nt); 
 eta1=zeros(nx,nt); 
 pot_rho3=zeros(nz,nx,nt);
 ue3=zeros(nz,nx,nt); we3=zeros(nz,nx,nt);
 situ_temp3=zeros(nz,nx,nt);
 
  %store the initial value
  for i=1:nz
      for j=1:nx
          ii=Number(i,j);
          T3(i,j,1)=T2(ii,1);         
      end
  end
  
 if nt0~=1
    u21=ncread(path1,'u');  w21=ncread(path1,'w');
    T21=ncread(path1,'T');  eta21=ncread(path1,'eta');
    pres21=ncread(path1,'pres');    rho21=ncread(path1,'rho');
    Ah21=ncread(path1,'Ah');    Az21=ncread(path1,'Az');
    
    for i=1:nz
        for j=1:nx
            ii=Number(i,j);
            u2(ii,1)=u21(i,j);  w2(ii,1)=w21(i,j);
            T2(ii,1)=T21(i,j);  
            
            u3(i,j,1)=u21(i,j);  w3(i,j,1)=w21(i,j);  T3(i,j,1)=T21(i,j);
            Ah3(i,j,1)=Ah21(i,j);   Az3(i,j,1)=Az21(i,j);
            rho3(i,j,1)=rho21(i,j);    pres3(i,j,1)=pres21(i,j); 
            eta1(j,1)=eta21(j,1); 
        end
    end
    
    eta2(:)=eta21(:); 
 end

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %Start calculating
ut1=u2; wt1=w2; etat=eta2;
tt=1;
 for k=(t0+1):nt1
        Tt1=T2;  %%potential temperature
        kht=kh2; kzt=kz2;
        eta2=zeros(nx,1);   %values of the current timestep
        
        %%%% conversion from potential temp. to in situ temp.
        alpha = 2E-4; cp = 1800; g = 22;
%         [situ_temp]=conversion_theta2temp(nx, nz, dz, alpha, cp, g, Tt1); 
        [situ_temp]=conversion_theta2temp_delR(nx, nz, delRc, alpha, cp, g, Tt1);
        
        pref=0; deltarho=0.1;
    	[rho2] = SilicatesDensity_modified2( nx, nz, z, pref, T_sol, T_liq, g, situ_temp, deltarho);
        rho1t=rho2;
        rhos=zeros(nx,1);   rhos(1:end)=rho2(1:nx,1);  
        
        pres2=zeros(nxz,1);     %values of the last timestep using the density of last timestep
        p1=zeros(nx,1);   p2=zeros(nx,1);
        for i=1:nz
            for j=1:nx
                ii=Number(i,j);
                if i==1
                    p1(j)=rhos(j)*g*etat(j,1); p2(j)=rhos(j)*g*etat(j,1)+abs(rho2(ii)*g*0.5*delRc(i));
                else
                    p1(j)=p2(j);
                    p2(j)=p2(j)+abs(rho2(ii)*g*delRc(i-1));
                end
                pres2(ii)=0.5*(p1(j)+p2(j));                
            end
        end
        pres1t=pres2; 
        
        %%%%%%%%%%%%%%%%%%%%%%%
        % Visocisty, using the temp. of the last timestep
        Aht=zeros(nxz,1);   Azt=zeros(nxz,1); 
        for i=1:nz
            for j=1:nx
                ii = Number(i,j);
                a=5.3;    %tanh(a)=1;
                b=atanh(Ah_liq/Ah_sol);
                if situ_temp(ii)>=T_liq
                    Aht(ii,1)=Ah_liq*exp(-10/situ_temp(ii)*(situ_temp(ii)-T_liq));
                elseif situ_temp(ii)< T_sol
                    Aht(ii,1)=Ah_sol*exp(-0.5/situ_temp(ii)*(situ_temp(ii)-T_sol));
                elseif (T_sol<=situ_temp(ii)) && (situ_temp(ii)<T_liq)
                    Aht(ii,1)=Ah_sol*tanh(a-(a-b)*(situ_temp(ii)-T_sol)/(T_liq-T_sol));
                end
            end
        end
        %Vertical viscosity
        for i=1:nz
            for j=1:nx
                ii = Number(i,j);
                a=5.3;    %tanh(a)=1;
                b=atanh(Az_liq/Az_sol);
                if situ_temp(ii)>=T_liq
                    Azt(ii,1)=Az_liq*exp(-10/situ_temp(ii)*(situ_temp(ii)-T_liq));
                elseif situ_temp(ii)< T_sol
                    Azt(ii,1)=Az_sol*exp(-0.5/situ_temp(ii)*(situ_temp(ii)-T_sol)); 
                elseif (T_sol<=situ_temp(ii)) && (situ_temp(ii)<T_liq)
                    Azt(ii,1)=Az_sol*tanh(a-(a-b)*(situ_temp(ii)-T_sol)/(T_liq-T_sol));
                end
            end
        end
       
        %%%%%%%%%%%%%%%%%%%%%%%% momentum equation loop
       for ti=1:num_dt

%         [u2] = AdvectionDiffusion2Dforu_implicit_bottomDrag( nx, nz, dt, dx, dz, x, z, Aht, Azt, ut1, wt1, pres1t, rho1t, taux, bottomDragLinear);
        [u2] = AdvectionDiffusion2Dforu_delR( nx, nz, dt, dx, delRc, Aht, Azt, ut1, wt1, pres1t, rho1t, taux, bottomDragLinear);
        u4=zeros(nz,nx);            % 1D to 2D
        for i=1:nz
            for j=1:nx
                 ii=Number(i,j);
                 u4(i,j)=u2(ii,1);
            end
        end
        
        %%%%%%%%%%%%%%%%%%Solve SSH
        %periodic boundaries
        hdiv_inter=zeros(nx,1);
        for j=1:nx
            for i=2:nz
                ii=Number(i,j);
                if j==1
                    ii2=Number(i,nx-1);
                    hdiv_inter(j)=hdiv_inter(j)-dt/2/dx*(u2(ii+1,1)-u2(ii2,1))*delRc(i-1);
                elseif j==nx
                    ii2=Number(i,2);
                    hdiv_inter(j)=hdiv_inter(j)-dt/2/dx*(u2(ii2,1)-u2(ii-1,1))*delRc(i-1);
                else
                    hdiv_inter(j)=hdiv_inter(j)-dt/2/dx*(u2(ii+1,1)-u2(ii-1,1))*delRc(i-1);
                end
            end
            
            eta2(j,1)=etat(j,1) - hdiv_inter(j);
        end

        %%%%%%%%%%%%%%%Solution of w
        %w=\partial \eta /\partial t at the surface
        w4=zeros(nz,nx);
        for i=1:nz
            for j=1:nx 
                if i==nz
                    w4(i,j)=0;
                elseif i==1
                    w4(1,j)=(eta2(j,1)-etat(j,1))/dt;  
                elseif j==1 
                    w4(i,j)=w4(i-1,j)-(z(i)-z(i-1))/2/dx*(u4(i,j+1)-u4(i,nx-1));
                elseif j==nx
                    w4(i,j)=w4(i-1,j)-(z(i)-z(i-1))/2/dx*(u4(i,2)-u4(i,j-1));   
                else
                    w4(i,j)=w4(i-1,j)-(z(i)-z(i-1))/2/dx*(u4(i,j+1)-u4(i,j-1));   
                end
            end
        end
                
        %2D to 1D
         for i=1:nz
             for j=1:nx
                ii = Number(i,j);
                w2(ii,1) = w4(i,j);
             end
         end
         
         ut1 = u2; wt1 = w2; etat = eta2;
       end
       
       %%%%%%%%%%Gent-McWilliams parameterization
       pot_rho=PotentiaDensity(nx, nz, T_sol, T_liq, Tt1, deltarho);       
       for i=2:nz
           for j=1:nx
               ii=Number(i,j);
               ii1=Number(i-1,j); 
               if pot_rho(ii)<pot_rho(ii1)
                   kzt(ii1)=kz_convec;
               end
           end
       end
       
        %%%%%%%%%%%%%%%%%%%%%%
%         [T2]=AdvectionDiffusion2D_FT_ConstK_periodic( nx, nz, dt_temp, dx, dz, x, z, kht, kzt, Tx, Tbot, Twest, Teast, ut1, wt1, Tt1,tau);
        [T2] = AdvectionDiffusion2DforT_delR( nx, nz, dt, dx, delRc, kht, kzt, Tx, ut1, wt1, Tt1,tau);
        T4=zeros(nz,nx);            % 1D to 2D
        for i=1:nz
            for j=1:nx
                 ii=Number(i,j);
                 T4(i,j)=T2(ii,1);
            end
        end

        kh4=zeros(nz,nx);   kz4=zeros(nz,nx);
        Ah4=zeros(nz,nx);   Az4=zeros(nz,nx);
        rho4=zeros(nz,nx);  pres4=zeros(nz,nx);
        pot_rho4=zeros(nz,nx);
        situ_temp4=zeros(nz,nx);
        
        for i=1:nz
            for j=1:nx
                 ii=Number(i,j);
                 kh4(i,j)=kht(ii,1);    %values of the last time step
                 kz4(i,j)=kzt(ii,1);
                 Ah4(i,j)=Aht(ii,1);
                 Az4(i,j)=Azt(ii,1);
                 
                 rho4(i,j)=rho2(ii,1);
                 pres4(i,j)=pres2(ii,1);
                 pot_rho4(i,j)=pot_rho(ii,1);
                 
                 situ_temp4(i,j)=situ_temp(ii);
            end
        end
        
        %%%%%%%%%%%%%%%%%save data but not for all
        if rem(k,dtsave)==0
            tt=tt+1;
            u3(:,:,tt)=u4(:,:);     w3(:,:,tt)=w4(:,:);     T3(:,:,tt)=T4(:,:);
            kh3(:,:,tt)=kh4(:,:);    kz3(:,:,tt)=kz4(:,:);
            Ah3(:,:,tt)=Ah4(:,:);   Az3(:,:,tt)=Az4(:,:);
            rho3(:,:,tt)=rho4(:,:);     pres3(:,:,tt)=pres4(:,:);
            eta1(:,tt)=eta2(:);
            pot_rho3(:,:,tt)=pot_rho4(:,:);
            situ_temp3(:,:,tt)=situ_temp4(:,:);
        end
        
       %%%%%%%%%%%%%%%%%Output many times
       tt2=k*dt;
         if rem(k,dt_output)==0
             path=['Snapshot', num2str(k), '.nc'];
             nccreate(path,'u','Dimensions', {'z',nz,'x',nx},'FillValue','disable');
             nccreate(path,'T','Dimensions', {'z',nz,'x',nx},'FillValue','disable');
             nccreate(path,'w','Dimensions', {'z',nz,'x',nx},'FillValue','disable');
             nccreate(path,'rho','Dimensions', {'z',nz,'x',nx},'FillValue','disable');
             nccreate(path,'pot_rho','Dimensions', {'z',nz,'x',nx},'FillValue','disable');
             nccreate(path,'pres','Dimensions', {'z',nz,'x',nx},'FillValue','disable');
             nccreate(path,'Ah','Dimensions', {'z',nz,'x',nx},'FillValue','disable');
             nccreate(path,'Az','Dimensions', {'z',nz,'x',nx},'FillValue','disable');
             nccreate(path,'eta','Dimensions', {'x',nx},'FillValue','disable');
             nccreate(path,'x','Dimensions',{'x',nx});
             nccreate(path,'z','Dimensions',{'z',nz});
             nccreate(path,'ST','Dimensions', {'z',nz,'x',nx},'FillValue','disable');
             
             ncwrite(path,'x',x);    ncwrite(path,'z',z);
             ncwrite(path,'u',u4);   ncwrite(path,'T',T4);  ncwrite(path,'w',w4);
             ncwrite(path,'rho',rho4);  ncwrite(path,'pres',pres4);
             ncwrite(path,'pot_rho',pot_rho4);
             ncwrite(path,'Ah',Ah4);    ncwrite(path,'Az',Az4);
             ncwrite(path,'eta',eta2);
             ncwrite(path,'ST',situ_temp4);
         end
 end
  
deltaT=1;
ti=0:deltaT:(nt-1);  time=ti.*dt_temp;  time=time.*dtsave; 
time(:)=time(:)+t0;

  %%
path='twod_H5km_delRc2_bd1E3_Keplertaux2_cd1E-2_melt2000K_kh1000kz1E-4_Az1E2_restart5.nc';
nccreate(path,'u','Dimensions', {'z',nz,'x',nx,'time',nt},'FillValue','disable');
nccreate(path,'T','Dimensions', {'z',nz,'x',nx,'time',nt},'FillValue','disable');
nccreate(path,'w','Dimensions', {'z',nz,'x',nx},'FillValue','disable');
nccreate(path,'rho','Dimensions', {'z',nz,'x',nx},'FillValue','disable');
nccreate(path,'pot_rho','Dimensions', {'z',nz,'x',nx},'FillValue','disable');
nccreate(path,'pres','Dimensions', {'z',nz,'x',nx},'FillValue','disable');
nccreate(path,'Ah','Dimensions', {'z',nz,'x',nx},'FillValue','disable');
nccreate(path,'Az','Dimensions', {'z',nz,'x',nx},'FillValue','disable');
nccreate(path,'eta','Dimensions', {'x',nx,'time',nt},'FillValue','disable');
nccreate(path,'ST','Dimensions', {'z',nz,'x',nx,'time',nt},'FillValue','disable');

% 
nccreate(path,'x','Dimensions',{'x',nx});
nccreate(path,'z','Dimensions',{'z',nz});
nccreate(path,'time','Dimensions',{'time',nt});
% 
ncwrite(path,'x',x);
ncwrite(path,'z',z);
ncwrite(path,'time',time);
% 
ncwrite(path,'u',u3);
ncwrite(path,'T',T3);
ncwrite(path,'ST',situ_temp3);
ncwrite(path,'w',squeeze(w3(:,:,end)));
ncwrite(path,'rho',squeeze(rho3(:,:,end)));
ncwrite(path,'pot_rho',squeeze(pot_rho3(:,:,end)));
ncwrite(path,'pres',squeeze(pres3(:,:,end)));
ncwrite(path,'Ah',squeeze(Ah3(:,:,end)));
ncwrite(path,'Az',squeeze(Az3(:,:,end)));
ncwrite(path,'eta',eta1);


