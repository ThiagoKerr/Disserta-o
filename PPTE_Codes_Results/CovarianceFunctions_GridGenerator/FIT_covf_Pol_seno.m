clc; clear all; close all; format long g; tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rotina para gerar valores de coeficientes de função polinomial senoides para a Colocação por mínimos quadrados %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Entrada de dados %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Importando os pontos de amostragem da Anomalia de gravidade residual (mGal)
Dados=load('Pontos_GO200.txt');

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
s=covdist'/108000*pi/180; % rad

% Ajustamento pelo MMQ - covariance = a0 + a1*sen(s) + a2*sen(2*s) + a3*sen(3*s) + a4*sen(4*s) ...
%                                        melhor ajuste = grau 4
Lb=covariance;
P=eye(length(covariance),length(covariance));
% A=[ones(length(Lb),1) sin(s) sin(2*s) sin(3*s)];
A=[ones(length(Lb),1) sin(s) sin(2*s) sin(3*s) sin(4*s)];

Xa=inv(A'*P*A)*A'*P*Lb;
% a0=Xa(1); a1=Xa(2); a2=Xa(3); a3=Xa(4);
a0=Xa(1); a1=Xa(2); a2=Xa(3); a3=Xa(4); a4=Xa(5);

% covariance_ajust = a0 + a1*sin(s) + a2*sin(2*s) + a3*sin(3*s);
covariance_ajust = a0 + a1*sin(s) + a2*sin(2*s) + a3*sin(3*s) + a4*sin(4*s);
R2=1 - sum((covariance - covariance_ajust).^2)/sum((covariance - mean(covariance)).^2)

figure(1);
plot(covdist',covariance_ajust,'DisplayName','Função covariância ajustada - polinomial de seno','LineWidth',1.5);   hold on
plot(covdist',covariance,'DisplayName','Função covariância empírica','LineWidth',1.5); hold off
str=['R^2 = ' num2str(R2)];    dim = [.6 .45 .3 .3];
annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',14,'BackgroundColor','white');
grid on; grid minor;  lgd = legend; lgd.NumColumns = 1;
xlabel('Distância (m)','FontName','Times New Roman','FontSize',14)
ylabel('Covariância (mGal^2)','FontName','Times New Roman','FontSize',14)
set(gca,'FontName','Times New Roman','FontSize',14)
print(gcf, 'F:\DISSERTACAO\PPTE_Metodologia\GRADE\Funções covariância\Resultados\R2_Polinomial_GO200.jpeg', '-djpeg', '-r300');
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