clc; clear all; close all; format long g; tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Rotina para gerar grade usando Colocação por mínimos quadrados usando função logarítmica 3D %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Importando os pontos de amostragem da Anomalia de gravidade residual (mGal)
Dados=load('Pontos_GO100.txt');
desvio_g=1;     % mGal

% Importando os pontos da grade interpolada de Anomalia de gravidade residual (mGal)
Dados_grade = load('GRADE_INTERPOLAR_elip.txt');
resolucao_graus=0.05;    % resolução em graus decimais

% Dados de entrada - coeficientes polinomiais
% Xa=[1456.93194903678;
% 0.0137061702790406;
% -2.90376799894151e-07;
% 1.99819976255458e-12;
% -7.86800482621596e-18;
% 2.02127135928127e-23;
% -3.54601159811198e-29;
% 4.25471341465806e-35;
% -3.40265643739201e-41;
% 1.72085909909805e-47;
% -4.95290219913017e-54;
% 6.16212216453955e-61];

Xa=[323.777196295906
-0.00854950602630977
7.753653280279e-08
-3.32736810051934e-13
7.56888811585078e-19
-9.34520628564655e-25
5.88543026604383e-31
-1.4699848318701e-37];

% Xa=[1888.72796638329;
% -0.0102105499521405;
% 3.85146063110817e-08;
% -1.07382194281858e-13;
% 1.51237834754826e-19;
% -9.67981673312113e-26;
% 2.26966397870036e-32];

% Tamanho dos blocos de interpolação em graus
d_lat=0.5;
d_long=0.5;

% Offsets em graus para o efeito de borda nas interpoalções
offset=0.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Preparação para o cálculo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dados das estações
Latitudes=Dados(:,1); Longitudes=Dados(:,2); Altitudes=Dados(:,3);
AG_residuais=Dados(:,4)*1D-5; % m/s2
Desvios_g=desvio_g*ones(length(Dados),1)*1D-5; % m/s2
Latitudes_interp=Dados_grade(:,1); Longitudes_interp=Dados_grade(:,2); Altitudes_interp=Dados_grade(:,3);

clear Dados_grade Dados

% Resolvendo por blocos de interpolação
blocos_lat=ceil((max(Latitudes_interp)-min(Latitudes_interp))/d_lat);      % quantidade de blocos em lat
blocos_long=ceil((max(Longitudes_interp)-min(Longitudes_interp))/d_long);  % quantidade de blocos em long
k=0;
% for i=1:length(Latitudes_interp)
%     geoid = geoidheight(Latitudes_interp(i),Longitudes_interp(i));
%     Altitudes_interp(i) = Altitudes_interp(i)+geoid;
% end



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
Gouti=LSCPoly(ei,ni,e,n,Xa,N,AG_res);

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
geotiffwrite('Grade_Polinomial_GO100_Elip.tif',(Imagem),Arquivo)

load handel
sound(y,Fs)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions which fit the covariance function parameters to the data
% by Tiago Lima Rodrigues with base in the Jack (2024)'s griddataLSC function (https://www.mathworks.com/matlabcentral/fileexchange/57342-griddatalsc),
% Developed July 16, 2024.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to perform the least squares collocation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Least squares collocation using polynomial covariance function
function [SolG]=LSCPoly(Xi,Yi,X,Y,Xa,N,G)
% a0=Xa(1); a1=Xa(2); a2=Xa(3); a3=Xa(4); a4=Xa(5); a5=Xa(6); a6=Xa(7);  a7=Xa(8); a8=Xa(9);  a9=Xa(10); a10=Xa(11);  a11=Xa(12);
a0=Xa(1); a1=Xa(2); a2=Xa(3); a3=Xa(4); a4=Xa(5); a5=Xa(6); a6=Xa(7); a7=Xa(8);
% a0=Xa(1); a1=Xa(2); a2=Xa(3); a3=Xa(4); a4=Xa(5); a5=Xa(6); a6=Xa(7);
s2=(X*ones(size(X))'-ones(size(X))*X').^2+(Y*ones(size(Y))'-ones(size(Y))*Y').^2;
s=sqrt(s2);
% Czz=(a0+a1*s+a2*s.^2+a3*s.^3+a4*s.^4+a5*s.^5+a6*s.^6+a7*s.^7+a8*s.^8+a9*s.^9+a10*s.^10+a11*s.^11)*1D-5;
Czz=(a0+a1*s+a2*s.^2+a3*s.^3+a4*s.^4+a5*s.^5+a6*s.^6+a7*s.^7)*1D-5;
% Czz=(a0+a1*s+a2*s.^2+a3*s.^3+a4*s.^4+a5*s.^5+a6*s.^6)*1D-5;
s2i=(Xi*ones(size(X))'-ones(size(Xi))*X').^2+(Yi*ones(size(Y))'-ones(size(Yi))*Y').^2;
si=sqrt(s2i);
Csz=(a0+a1*si+a2*si.^2+a3*si.^3+a4*si.^4+a5*si.^5+a6*si.^6+a7*si.^7)*1D-5;
% Csz=(a0+a1*si+a2*si.^2+a3*si.^3+a4*si.^4+a5*si.^5+a6*si.^6)*1D-5;
LF=(Czz+diag(N))\eye(size(Czz));
SolG=Csz*(LF*G);
end