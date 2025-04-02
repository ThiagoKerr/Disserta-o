clc; clear all; format long g; tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Rotina para gerar grade usando Colocação por mínimos quadrados usando função logarítmica 3D %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Importando os pontos de amostragem da Anomalia de gravidade residual (mGal)
Dados_terrestre=load('TERRESTRE_MGG300_ERTM.txt');
Dados_aereo=load('AEREO_MGG300_ERTM.txt');
desvio_g_terrestre=1;     % mGal
desvio_g_aereo=2;         % mGal

% Importando os pontos da grade interpolada de Anomalia de gravidade residual (mGal)
Dados_grade = load('GRADE_INTERPOLAR.txt');
resolucao_graus=0.0166666666667;    % resolução em graus decimais

% Dados de entrada
C0=82.3287097396916e-5;  % (m/s2)^2
D=7376.09721992145;         % m
T=18440.2430498036;        % m

% Tamanho dos blocos de interpolação em graus
d_lat=0.25;
d_long=0.25;

% Offsets em graus para o efeito de borda nas interpoalções
offset=0.6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Preparação para o cálculo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dados das estações
Dados_total_3D=[Dados_terrestre; Dados_aereo];
% Dados_total_3D=Dados_terrestre;

Latitudes=Dados_total_3D(:,1); Longitudes=Dados_total_3D(:,2); Altitudes=Dados_total_3D(:,3);
AG_residuais=Dados_total_3D(:,4)*1D-5; % m/s2
Desvios_g=[desvio_g_terrestre*ones(length(Dados_terrestre),1); desvio_g_aereo*ones(length(Dados_aereo),1)]*1D-5; % m/s2
% Desvios_g=desvio_g_terrestre*ones(length(AG_res),1)*1D-5; % m/s2
Latitudes_interp=Dados_grade(:,1); Longitudes_interp=Dados_grade(:,2); Altitudes_interp=Dados_grade(:,3);

clear Dados_terrestre Dados_aereo Dados_grade Dados_total_3D

% Resolvendo por blocos de interpolação
blocos_lat=ceil((max(Latitudes_interp)-min(Latitudes_interp))/d_lat);      % quantidade de blocos em lat
blocos_long=ceil((max(Longitudes_interp)-min(Longitudes_interp))/d_long);  % quantidade de blocos em long
k=0;

for i=1:blocos_lat
for j=1:blocos_long

Lat=Latitudes; Long=Longitudes; h=Altitudes; AG_res=AG_residuais; N=Desvios_g;
Lati=Latitudes_interp; Longi=Longitudes_interp; hi=Altitudes_interp;

min_lat=min(Lati)+(i-1)*d_lat;
max_lat=min(Lati)+(i)*d_lat;
min_long=min(Longi)+(j-1)*d_long;
max_long=min(Longi)+(j)*d_long;

Lat(Lat<min_lat-offset)=NaN;       % latitude mínima
if i==blocos_lat
Lat(Lat>max_lat+offset)=NaN;       % latitude máxima chegando no limite da área
else
Lat(Lat>=max_lat+offset)=NaN;      % latitude máxima
end
index_lat=~isnan(Lat); index_lat=double(index_lat); index_lat(index_lat==0)=NaN; % criando a máscara para exclusão
Long=Long.*index_lat;  % excluindo das longitudes os ponto fora dos limites de latitude
Long(Long<min_long-offset)=NaN;       % longitude mínima
if j==blocos_long
Long(Long>max_long+offset)=NaN;       % longitude máxima chegando no limite da área
else
Long(Long>=max_long+offset)=NaN;      % longitude máxima
end
index_long=~isnan(Long); index_long=double(index_long); index_long(index_long==0)=NaN; % criando a máscara para exclusão
Lat=Lat.*index_long;   % excluindo das latitudes os pontos fora dos limites de longitude

AG_res=AG_res.*index_long;  % Filtrando os dados para a região de interpolação
h=h.*index_long;            % Filtrando os dados para a região de interpolação
N=N.*index_long;            % Filtrando os dados para a região de interpolação

AG_res(isnan(AG_res))=[];   % eliminando os NaN
h(isnan(h))=[];             % eliminando os NaN
N(isnan(N))=[];             % eliminando os NaN
Lat(isnan(Lat))=[];         % eliminando os NaN
Long(isnan(Long))=[];       % eliminando os NaN

% Dados da grade de interpolação
Lati(Lati<min_lat)=NaN;       % latitude mínima
if i==blocos_lat
Lati(Lati>max_lat)=NaN;       % latitude máxima chegando no limite da área
else
Lati(Lati>=max_lat)=NaN;      % latitude máxima
end
index_lati=~isnan(Lati); index_lati=double(index_lati); index_lati(index_lati==0)=NaN;
Longi=Longi.*index_lati;
Longi(Longi<min_long)=NaN;       % longitude mínima
if j==blocos_long
Longi(Longi>max_long)=NaN;        % longitude máxima chegando no limite da área
else
Longi(Longi>=max_long)=NaN;      % longitude máxima
end
index_longi=~isnan(Longi); index_longi=double(index_longi); index_longi(index_longi==0)=NaN;

Lati=Lati.*index_longi;     % Filtrando os dados para a região de interpolação
hi=hi.*index_longi;         % Filtrando os dados para a região de interpolação

hi(isnan(hi))=[];           % eliminando os NaN
Lati(isnan(Lati))=[];       % eliminando os NaN
Longi(isnan(Longi))=[];     % eliminando os NaN

% Transformação para o Sistema Geodésico Local com origem no centro da área - Elipsoide GRS80

Lat0=mean(Lat); Long0=mean(Long); h0=min(h);
grs80 = referenceEllipsoid('Geodetic Reference System 1980');
[e,n,u] = geodetic2enu(Lat,Long,h,Lat0,Long0,h0,grs80);
[ei,ni,ui] = geodetic2enu(Lati,Longi,hi,Lat0,Long0,h0,grs80);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Colocação 3D e gridagem %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Gouti,Varsi]=LSClogarithmic(ei,ni,hi,e,n,h,C0,D,T,N,AG_res);

aloca=(k+1):(max(k))+length(Gouti); aloca=aloca';
k=max(aloca);
Resultado(aloca,:)=[Longi Lati Gouti*1D5 aloca];

end
end

toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Exportar em formato tif %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Resultado=sortrows(Resultado);
nlin=length(unique(Latitudes_interp));
ncol=length(unique(Longitudes_interp));

Imagem=reshape(Resultado(:,3),nlin,ncol);

lat_S_grade=min(Latitudes_interp)-(resolucao_graus/2);
lat_N_grade=max(Latitudes_interp)+(resolucao_graus/2);
lon_W_grade=min(Longitudes_interp)-(resolucao_graus/2);
lon_E_grade=max(Longitudes_interp)+(resolucao_graus/2);
Arquivo=georasterref('RasterSize',size(Imagem),'LatitudeLimits',[lat_S_grade,lat_N_grade],'LongitudeLimits',[lon_W_grade,lon_E_grade]);
geotiffwrite('Grade_DG_res_3DLSC_ERTM_V1.tif',(Imagem),Arquivo)

load handel
sound(y,Fs)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions which fit the covariance function parameters to the data
% by Jack (2024). griddataLSC (https://www.mathworks.com/matlabcentral/fileexchange/57342-griddatalsc),
% MATLAB Central File Exchange. Retrieved January 13, 2024.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to perform the least squares collocation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Least squares collocation using logarithmic covariance function
function [SolG,Vars]=LSClogarithmic(Xi,Yi,Zi,X,Y,Z,C0,D,T,N,G)
a=[1,-3,3,-1];
Czz=zeros(length(X),length(X));
Csz=zeros(length(Xi),length(X));
Css=zeros(length(Xi),length(Xi));
for i=0:length(a)-1
Di=D+i*T;
Zm=Z*ones(size(Z))'+ones(size(Z))*Z'+Di;
s2=(X*ones(size(X))'-ones(size(X))*X').^2+(Y*ones(size(Y))'-ones(size(Y))*Y').^2;
r=sqrt(s2+Zm.^2);
Czz=Czz+a(i+1)*log(Zm+r);
Zmi=Zi*ones(size(Z))'+ones(size(Zi))*Z'+Di;
s2i=(Xi*ones(size(X))'-ones(size(Xi))*X').^2+(Yi*ones(size(Y))'-ones(size(Yi))*Y').^2;
ri=sqrt(s2i+Zmi.^2);
Csz=Csz+a(i+1)*log(Zmi+ri);
Zmii=Zi*ones(size(Zi))'+ones(size(Zi))*Zi'+Di;
s2ii=(Xi*ones(size(Xi))'-ones(size(Xi))*Xi').^2+(Yi*ones(size(Yi))'-ones(size(Yi))*Yi').^2;
rii=sqrt(s2ii+Zmii.^2);
Css=Css+a(i+1)*log(Zmii+rii);
end
f=C0/log((D+T)^3*(D+3*T)/((D)*(D+2*T)^3));
Czz=-f*Czz;
Csz=-f*Csz;
LF=((Czz+diag(N))\eye(size(Czz)));
SolG=Csz*(LF*G);
Css=-f*Css;
Vars=diag(Css-Csz*LF*Csz');
end