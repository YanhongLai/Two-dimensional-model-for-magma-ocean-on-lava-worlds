    function [rho] = SilicatesDensity_modified2( nx, nz, z, pref, T_sol, T_liq, g, temp0, deltarho)
    
    %默认pref为0，所以这里的计算没有涉及
    rho = zeros ( nx*nz, 1 );
    pres = zeros ( nx*nz, 1 );
    
    %令表面的密度为参考密度，pres=0, temp=2000 K，这个温度一般高于融化温度
    %密度和温度的系数
    drhodT=-1.86E-4;
    rho0_2000=2.673;    temp00=2000;    pres0=0; 

    %二次拟合cftool得到压强的系数
    pa=-0.006667;   pb=0.1022;  pc=2.823;%压强的单位是GPa，密度是g/cm3。后面再转换
    
    % density at the surface, 用来算压强
    rho_liq=rho0_2000+drhodT*(T_liq-temp00);
    rho_sol=rho_liq*(1+deltarho);
    T_transit=T_liq-T_sol; 

    a=5.3;
%     a=50;
    T_surf=zeros(nx,1); T_surf(:)=temp0(1:nx,1);
    rho_surf=zeros(nx,1);
    for j=1:nx
        if T_surf(j)>T_liq
            rho_surf(j)=rho0_2000+drhodT*(T_surf(j)-temp00);  %
        elseif T_surf(j)<T_sol
            rho_surf(j)=rho_sol+drhodT*(T_surf(j)-T_sol);   
        else
            c=rho_liq/rho_sol;
            b=atanh(c);
            rho_surf(j) = rho_sol*tanh(a-(a-b)*(T_surf(j)-T_sol)/T_transit);
%             rho_surf(j) = rhom+rhodif*tanh(a*((T_surf(j)-Tm))/(T_liq-Tm));
        end
    end
    
    rho_solid=zeros(nz*nx,1 );
    
  % Setup numbering
  num = 1; Number=zeros(nz,nx);
  for i=1:nz
      for j=1:nx
          Number(i,j)=num;
          num=num+1;
      end
  end
  
  %%%%%%%%%%%%%%%% calculate density using eos
    for i=1:nz
        for j=1:nx
            ii = Number(i,j);
            pres(ii,1)=abs(1000*rho_surf(j)*g*z(i));  %计算密度的压强只是估算，给的不是准确的
            pres(ii,1)=1E-9*pres(ii,1); %单位换算为GPa
            
            rho_solid(ii,1) = rho_sol-pa*pres(ii,1)*pres(ii,1)+pb*pres(ii,1);
            c=(rho_liq-pa*pres(ii,1)*pres(ii,1)+pb*pres(ii,1))/(rho_solid(ii,1));
            b=atanh(c);
            
            if temp0(ii,1)<T_sol
                rho(ii,1) = rho_solid(ii)+drhodT*(temp0(ii,1)-T_sol);
            elseif temp0(ii,1)>T_liq
                rho(ii,1) = rho0_2000+drhodT*(temp0(ii,1)-temp00)-pa*pres(ii,1)*pres(ii,1)+pb*pres(ii,1);
            else
                rho(ii,1) = rho_solid(ii,1)*tanh(a-(a-b)*(temp0(ii,1)-T_sol)/T_transit);
            end
        end
    end
    
    rho=1000.*rho;
    
    end