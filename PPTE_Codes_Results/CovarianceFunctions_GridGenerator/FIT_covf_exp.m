clc; clear all; close all; format long g; tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Rotina para gerar valores de C0, D e T para a Colocação por mínimos quadrados %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Entrada de dados %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Importando os pontos de amostragem da Anomalia de gravidade residual (mGal)
Dados=load('Pontos_GO100.txt');
desvio_Dg=1;     % mGal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Preparação para o cálculo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Lat=Dados(:,1); Long=Dados(:,2); h=Dados(:,3);
DG_res=Dados(:,4);
N=desvio_Dg*ones(length(Dados),1);

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
[covariance,covdist]=spatialCov(e,n,DG_res);
[C0,D,func_cov_ajust]=expcovarfit(e,n,DG_res,covariance,covdist,'no');
covdist=covdist';
R2=1 - sum((covariance - func_cov_ajust').^2)/sum((covariance - mean(covariance)).^2)
toc

figure(1);
plot(covdist,func_cov_ajust,'DisplayName','Função covariância ajustada - exponencial','LineWidth',1.5);   hold on
plot(covdist,covariance,'DisplayName','Função covariância empírica','LineWidth',1.5); hold off
str=['R^2 = ' num2str(R2)];    dim = [.6 .45 .3 .3];
annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',14,'BackgroundColor','white');
grid on; grid minor;  lgd = legend; lgd.NumColumns = 1;
xlabel('Distância (m)','FontName','Times New Roman','FontSize',14)
ylabel('Covariância (mGal^2)','FontName','Times New Roman','FontSize',14)
set(gca,'FontName','Times New Roman','FontSize',14)
print(gcf, 'F:\DISSERTACAO\PPTE_SegundoTeste\GRADE\Funções covariância\Resultados\GO100\R2_Exponencial_GO100.jpeg', '-djpeg', '-r300');
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
%% A function to fit the exponential covariance function parameters.
function [C0,Dbest,func_cov_ajust]=expcovarfit(X,Y,G,covariance,covdist,figplease)
Dmax=sqrt((max(X)-min(X))^2+(max(Y)-min(Y))^2);
Step=sqrt((2*pi*(Dmax/2)^2)/(length(X)));
C0=std(G)^2;
s=(covdist);
fiterr=999999999999999999999999999;
for D=Step:Step:Dmax
     covp = C0*exp(-s/D);
  if sum((covp-covariance').^2)<fiterr
     fiterr=sum((covp-covariance').^2);
     Dbest=D;
     if strcmp(figplease,'covfigure')==1
     figure(10)
     clf
     hold on
     plot(covdist,covariance(1:end),'o')
     plot(covdist,covp,'g')
     xlabel('Distance')
     ylabel('Covariance')
     title('Empirical covariance in blue, fitted model in green.')
     end
  func_cov_ajust=covp;
  end
end
end