clc; clear all; format long g; tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Rotina para gerar valores de C0, D e T para a Colocação por mínimos quadrados %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Entrada de dados %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Importando os pontos de amostragem da Anomalia de gravidade residual (mGal)
Dados_terrestre=load('TERRESTRE_MGG300_ERTM.txt');
Dados_aereo=load('AEREO_MGG300_ERTM.txt');
desvio_g_terrestre=1;     % mGal
desvio_g_aereo=2;         % mGal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Preparação para o cálculo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Dados_total_3D=[Dados_terrestre; Dados_aereo];
% Dados_total_3D=Dados_terrestre;

Lat=Dados_total_3D(:,1); Long=Dados_total_3D(:,2); h=Dados_total_3D(:,3);
AG_res=Dados_total_3D(:,4);
N=[desvio_g_terrestre*ones(length(Dados_terrestre),1); desvio_g_aereo*ones(length(Dados_aereo),1)];

% Transformação para o Sistema Geodésico Local com origem no centro da área - Elipsoide GRS80
Lat0=mean(Lat); Long0=mean(Long); h0=min(h);
grs80 = referenceEllipsoid('Geodetic Reference System 1980');
[e,n,u] = geodetic2enu(Lat,Long,h,Lat0,Long0,h0,grs80);

clear Dados_terrestre Dados_aereo Dados_grade Dados_total_3D Long Lat h Lat0 Long0 grs80

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Cáclculo de C0, D e T %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figplease='covfigure'; % ou 'no'
[covariance,covdist]=spatialCov(e,n,AG_res);
[C0,D,T,func_cov_ajust]=logcovarfit(e,n,AG_res,covariance,covdist,figplease)
covdist=covdist';
R2=1 - sum((covariance - func_cov_ajust).^2)/sum((covariance - mean(covariance)).^2)
toc

figure(1);
plot(covdist,func_cov_ajust,'DisplayName','Função covariância ajustada - Log3D','LineWidth',1.5);   hold on
plot(covdist,covariance,'DisplayName','Função covariância empírica','LineWidth',1.5); hold off
str=['R^2 = ' num2str(R2)];    dim = [.7 .45 .3 .3];
annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',14,'BackgroundColor','white');
grid on; grid minor;  lgd = legend; lgd.NumColumns = 1;
xlabel('Distância (m)','FontName','Times New Roman','FontSize',14)
ylabel('Covariância (mGal^2)','FontName','Times New Roman','FontSize',14)
set(gca,'FontName','Times New Roman','FontSize',14)

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A function to fit the logarithmic covariance function parameters.
% Disclaimer, the function is being fitted in two dimensional space.
function [C0,Dbest,Tbest,func_cov_ajust]=logcovarfit(X,Y,G,covariance,covdist,figplease)
Dmax=sqrt((max(X)-min(X))^2+(max(Y)-min(Y))^2);
Tmax=sqrt((max(X)-min(X))^2+(max(Y)-min(Y))^2);
% Dmax=50000;
% Tmax=50000;
Step=sqrt((2*pi*(Tmax/2)^2)/(length(X)));
C0=std(G)^2;
s=(covdist).^2;
Fit=[];
fiterr=999999999999999999999999999;
for D=Step:Step:Dmax
 for T=Step:Step:Tmax
     d11=D;
     d12=D+T;
     d13=D+2*T;
     d14=D+3*T;
     z1 = d11;
     z2 = d12;
     z3 = d13;
     z4 = d14;
     fc = C0/log(((d12.^3)/(d13.^3))*d14/d11);
     r1 = sqrt(s + z1.^2);
     r2 = sqrt(s + z2.^2);
     r3 = sqrt(s + z3.^2);
     r4 = sqrt(s + z4.^2);
     covp = fc*log((((z2+r2).^3)./((z3+r3).^3)).*(z4+r4)./(z1+r1));
     Fit=[Fit;sum((covp-covariance').^2),D,T];
  if sum((covp-covariance').^2)<fiterr
     fiterr=sum((covp-covariance').^2);
     Dbest=D;
     Tbest=T;
     if strcmp(figplease,'covfigure')==1
     figure(10)
     clf
     hold on
     plot(covdist,covariance(1:end),'o')
     plot(covdist,covp,'g')
     xlabel('Distance')
     ylabel('Covariance')
     title('Empirical covariance in blue, fitted model in green.')
     drawnow
     end
  func_cov_ajust=covp;
  end
 end
end
func_cov_ajust=func_cov_ajust';
end