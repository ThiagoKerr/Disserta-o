clc; clear all; close all; format long g; tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rotina para gerar valores de coeficientes polinomiais para a Colocação por mínimos quadrados %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Entrada de dados %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Importando os pontos de amostragem da Anomalia de gravidade residual (mGal)
Dados=load('Pontos_GO100.txt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Preparação para o cálculo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Lat=Dados(:,1); Long=Dados(:,2); h=Dados(:,3);
AG_res=Dados(:,4);

% Transformação para o Sistema Geodésico Local com origem no centro da área - Elipsoide GRS80
Lat0=mean(Lat); Long0=mean(Long); h0=min(h);
grs80 = referenceEllipsoid('Geodetic Reference System 1980');
[e,n,u] = geodetic2enu(Lat,Long,h,Lat0,Long0,h0,grs80);

clear Dados_terrestre Dados_aereo Dados_grade Dados_total_3D Long Lat h Lat0 Long0 grs80

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Cáclculo dos coeficientes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[covariance,covdist]=spatialCov(e,n,AG_res);
s=covdist'; % m

% Ajustamento pelo MMQ - covariance = a0 + a1*s + a2*s^2 + a3*s^3 + a4*s^4 + a5*s^5 + a6*s^6 - melhor estimativa com grau 11
Lb=covariance;  % mGal^2
P=eye(length(covariance),length(covariance));
% A=[ones(length(Lb),1) s s.^2  s.^3  s.^4  s.^5  s.^6   s.^7   s.^8   s.^9   s.^10   s.^11];
% A=[ones(length(Lb),1) s s.^2  s.^3  s.^4  s.^5  s.^6   s.^7   s.^8];
% A=[ones(length(Lb),1) s s.^2  s.^3  s.^4  s.^5  s.^6   s.^7];
A=[ones(length(Lb),1) s s.^2  s.^3  s.^4  s.^5  s.^6];

Xa=inv(A'*P*A)*A'*P*Lb;
% a0=Xa(1); a1=Xa(2); a2=Xa(3); a3=Xa(4); a4=Xa(5); a5=Xa(6); a6=Xa(7);  a7=Xa(8); a8=Xa(9);  a9=Xa(10); a10=Xa(11);  a11=Xa(12);
% a0=Xa(1); a1=Xa(2); a2=Xa(3); a3=Xa(4); a4=Xa(5); a5=Xa(6); a6=Xa(7);  a7=Xa(8); a8=Xa(9);
% a0=Xa(1); a1=Xa(2); a2=Xa(3); a3=Xa(4); a4=Xa(5); a5=Xa(6); a6=Xa(7); a7=Xa(8);
a0=Xa(1); a1=Xa(2); a2=Xa(3); a3=Xa(4); a4=Xa(5); a5=Xa(6); a6=Xa(7);

% covariance_ajust = a0 + a1*s + a2*s.^2 + a3*s.^3 + a4*s.^4 + a5*s.^5 + a6*s.^6 + a7*s.^7 + a8*s.^8 + a9*s.^9 + a10*s.^10 + a11*s.^11;
% covariance_ajust = a0 + a1*s + a2*s.^2 + a3*s.^3 + a4*s.^4 + a5*s.^5 + a6*s.^6 + a7*s.^7 + a8*s.^8;
% covariance_ajust = a0 + a1*s + a2*s.^2 + a3*s.^3 + a4*s.^4 + a5*s.^5 + a6*s.^6 + a7*s.^7;
covariance_ajust = a0 + a1*s + a2*s.^2 + a3*s.^3 + a4*s.^4 + a5*s.^5 + a6*s.^6;
R2=1 - sum((covariance - covariance_ajust).^2)/sum((covariance - mean(covariance)).^2)

figure(1);
plot(s,covariance_ajust,'DisplayName','Função covariância ajustada - polinomial','LineWidth',1.5);   hold on
plot(s,covariance,'DisplayName','Função covariância empírica','LineWidth',1.5); hold off
str=['R^2 = ' num2str(R2)];    dim = [.6 .45 .3 .3];
annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',14,'BackgroundColor','white');
grid on; grid minor;  lgd = legend; lgd.NumColumns = 1;
xlabel('Distância (m)','FontName','Times New Roman','FontSize',14)
ylabel('Covariância (mGal^2)','FontName','Times New Roman','FontSize',14)
set(gca,'FontName','Times New Roman','FontSize',14)
print(gcf, 'F:\DISSERTACAO\PPTE_SegundoTeste\GRADE\Funções covariância\Resultados\GO100\R2_Poly_GO100.jpeg', '-djpeg', '-r300');
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions which fit the covariance function parameters to the data
% by Jack (2024). griddataLSC (https://www.mathworks.com/matlabcentral/fileexchange/57342-griddatalsc),
% MATLAB Central File Exchange. Retrieved January 13, 2024.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A function to find the 2-D Empirical Covariance of the input data
function [covariance,covdist]=spatialCov(X,Y,G)
smax=sqrt((max(X)-min(X))^2+(max(Y)-min(Y))^2); % maximum range of the data
% smax=50000; % maximum range of the data
ds=sqrt((2*pi*(smax/2)^2)/(length(X)));         % average spatial density
covariance=zeros(round(smax/ds)+2,1);
ncov=zeros(round(smax/ds)+2,1);
for i = 1:length(G)
        for j = 1:length(G)
            if j~=i
          xx=(X(i)-X(j));
          yy=(Y(i)-Y(j)); 
          r = sqrt(xx.^2+yy.^2);
          ir = round((r/ds))+1;
          if r<smax
            covariance(ir) = covariance(ir) + G(i)*G(j) ;
            ncov(ir) = ncov(ir)+1;
          end
            end
        end
end
for i = 1:length(covariance)
  if ncov(i)~=0
      covariance(i) = covariance(i)/ncov(i);
  end
end  
covdist=(0:length(covariance)-1)*ds;
end