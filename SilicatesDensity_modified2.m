    function [rho] = SilicatesDensity_modified2( nx, nz, z, pref, T_sol, T_liq, g, temp0, deltarho)
    
    %Ĭ��prefΪ0����������ļ���û���漰
    rho = zeros ( nx*nz, 1 );
    pres = zeros ( nx*nz, 1 );
    
    %�������ܶ�Ϊ�ο��ܶȣ�pres=0, temp=2000 K������¶�һ������ڻ��¶�
    %�ܶȺ��¶ȵ�ϵ��
    drhodT=-1.86E-4;
    rho0_2000=2.673;    temp00=2000;    pres0=0; 

    %�������cftool�õ�ѹǿ��ϵ��
    pa=-0.006667;   pb=0.1022;  pc=2.823;%ѹǿ�ĵ�λ��GPa���ܶ���g/cm3��������ת��
    
    % density at the surface, ������ѹǿ
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
            pres(ii,1)=abs(1000*rho_surf(j)*g*z(i));  %�����ܶȵ�ѹǿֻ�ǹ��㣬���Ĳ���׼ȷ��
            pres(ii,1)=1E-9*pres(ii,1); %��λ����ΪGPa
            
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