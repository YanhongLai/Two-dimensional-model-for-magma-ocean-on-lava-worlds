function pot_rho=PotentiaDensity(nx, nz, T_sol, T_liq, Tt1, deltarho)  

  pot_rho=zeros(nx*nz,1);
    
    drhodT=-1.86E-4;
    rho0_2000=2.673;    temp00=2000;    
    rho_liq=rho0_2000+drhodT*(T_liq-temp00);
    rho_sol=rho_liq*(1+deltarho);
    T_transit=T_liq-T_sol; 
    
    a=5.3;
      % Setup numbering
  num = 1; Number=zeros(nz,nx);
  for i=1:nz
      for j=1:nx
          Number(i,j)=num;
          num=num+1;
      end
  end
  
   for i=1:nz
        for j=1:nx
            ii = Number(i,j);
            rho_solid = rho_sol;
            c=(rho_liq)/(rho_solid);
            b=atanh(c);
            
            if Tt1(ii,1)<T_sol
                pot_rho(ii,1) = rho_solid+drhodT*(Tt1(ii,1)-T_sol);
            elseif Tt1(ii,1)>T_liq
                pot_rho(ii,1) = rho0_2000+drhodT*(Tt1(ii,1)-temp00);
            else
                pot_rho(ii,1) = rho_solid*tanh(a-(a-b)*(Tt1(ii,1)-T_sol)/T_transit);
            end
        end
    end
    
    pot_rho=1000.*pot_rho;
    
end