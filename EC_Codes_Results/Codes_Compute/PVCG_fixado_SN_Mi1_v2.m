clc; clear all; format long g; close all; tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Entrada de dados %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Importando as grades 
Dg_res_MDE = double(imread('Grade_DG_res_3DLSC_ERTM_V1.tif'))*1D-5;  % distúrbios de gravidade residuais em m/s2 com resolução espacial igual do MDE
h = double(imread('./imagens/MDE_COLORADO_RECORTADO_1m_Mi0_ELIPSOIDAL.tif'));  % MDE detalhado com altitude elipsoidal
H = double(imread('./imagens/MDE_COLORADO_RECORTADO_1m_Mi0.tif'));  % MDE detalhado com altitude normal
Eta_res_Mi0_MDE = double(imread('Eta_residual_WG300_Mi0_ERTM_V1.tif'));  % anomalias de altura residuais do termo Mi0 em m com resolução espacial igual do MDE

resolucao_graus=0.016666667;          % Resolução do MDE em °
raio_integracao_graus=0.8;            % Raio de integração em ° para cálculo do termo Mi1

% Coordenadas geodésicas do canto superior esquerdo do primeiro pixel da imagem no canto NW - superior esquerdo - depois passando para o centro geométrico do pixel
% Distúrbios de gravidade residuais e MDE
latN_Dg=40.0083333333333471-(resolucao_graus/2);     % latN
lonW_Dg=-110.0083333333333542+(resolucao_graus/2);    % longW

% Anomalias de altura residuais com o termo Mi0
latN_Eta_Mi0=39.2083333173333557-(resolucao_graus/2);     % latN  
lonW_Eta_Mi0=-109.2083333173333557+(resolucao_graus/2);    % longW

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Preparação para o cálculo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Grades de coordenadas geodésicas da área de integração para a grade de distúrbios de gravidade residuais e MDE
[lin_Dg,col_Dg]=size(Dg_res_MDE);

latS_Dg=(latN_Dg-(lin_Dg-1)*resolucao_graus);
lonE_Dg=(lonW_Dg+(col_Dg-1)*resolucao_graus);

grade_lat_Dg=(latN_Dg:-resolucao_graus:latS_Dg)'.*ones(1,col_Dg)*pi/180;  % rad
grade_long_Dg=(lonW_Dg:resolucao_graus:lonE_Dg).*ones(lin_Dg,1)*pi/180;   % rad

% Localizando a posição da grade de anomalias de altura com o termo Mi0 na grade de lat/long dos distúrbios de gravidade residuais
Closest_Value_lat=interp1(grade_lat_Dg(:,1), grade_lat_Dg(:,1), latN_Eta_Mi0*pi/180, 'nearest');
find_lat=find(grade_lat_Dg(:,1)==Closest_Value_lat);
Closest_Value_long=interp1(grade_long_Dg(1,:), grade_long_Dg(1,:), lonW_Eta_Mi0*pi/180, 'nearest');
find_long=find(grade_long_Dg(1,:)==Closest_Value_long);

% Grades de coordenadas geodésicas da área de cálculo
[lin_Eta_Mi0,col_Eta_Mi0]=size(Eta_res_Mi0_MDE);

grade_lat_calc=grade_lat_Dg(find_lat:find_lat+lin_Eta_Mi0-1, find_lat:find_lat+col_Eta_Mi0-1);         % descontando os raios de integração da grade de área de integração
grade_long_calc=grade_long_Dg(find_long:find_long+lin_Eta_Mi0-1 , find_long:find_long+col_Eta_Mi0-1);  % descontando os raios de integração da grade de área de integração

% Distúrbios de gravidade residuais e H dos pontos de cálculo P
gradeDgP_calc=Dg_res_MDE(find_lat:find_lat+lin_Eta_Mi0-1, find_lat:find_lat+col_Eta_Mi0-1);            % descontando os raios de integração da grade de área de integração;
HP_calc=H(find_lat:find_lat+lin_Eta_Mi0-1, find_lat:find_lat+col_Eta_Mi0-1);                           % descontando os raios de integração da grade de área de integração;
hP_calc=h(find_lat:find_lat+lin_Eta_Mi0-1, find_lat:find_lat+col_Eta_Mi0-1);                           % descontando os raios de integração da grade de área de integração;

% Coordenadas em forma de vetor, cálculo da gravidade normal no teluroide, MDE e anomalia de altura residual com o termo Mi0 da grade de cálculo
% Pontos de cálculo
latP=grade_lat_calc(:);      longP=grade_long_calc(:);        DgP=gradeDgP_calc(:);     Eta_res_Mi0P=Eta_res_Mi0_MDE(:);      HP=HP_calc(:);      hP=hP_calc(:);
gama0=(9.7803267715*(1+0.0052790414*(sin(latP).^2)+0.0000232718*(sin(latP).^4)+0.0000001262*(sin(latP).^6)+0.0000000007*(sin(latP).^8)));  % m/s2 - GRS80
gama=gama0.*(1-2*HP.*(1+0.00335281068118-2*0.00335281068118*(sin(latP).^2)+0.00344978600308)/6378137+(3*HP.^2)/6378137^2); % m/s2 - GRS80
sin_latP=sin(latP);

% Pontos de integração
latQ=grade_lat_Dg(:);       longQ=grade_long_Dg(:);     DgQ=Dg_res_MDE(:);      hQ=h(:);
sin_latQ=sin(latQ);     cos_latQ=cos(latQ);

% Constantes
res_rad=resolucao_graus*pi/180;                % Resolução em rad
raio_integ_rad=raio_integracao_graus*pi/180;   % Raio de integração em rad
R=6371000; R_2=2*R;                            % Raio da esfera de mesmo volume do Elpsoide em metros e 2*R

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Solução do PVCG gravimétrico fixado %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Mi1=zeros(length(latP),1); % alocando espaço para ganhar velocidade

% % Áreas dividido pelo cosseno da latitude do pontos de integração
% AQ_div_cos_phi=abs((R^2)*res_rad*(sin(latQ+(res_rad/2))-sin(latQ-(res_rad/2))))./latQ;
% Áreas
AQ=abs((R^2)*res_rad*(sin(latQ+(res_rad/2))-sin(latQ-(res_rad/2))));


for p=1:length(latP)

% Filtrando os pontos de integração dentro da área de integração
cos_psi=sin_latQ*sin_latP(p)+cos_latQ*cos(latP(p)).*cos(longQ-longP(p));
psi=real(acos(cos_psi));

psi(psi>raio_integ_rad)=NaN; index=~isnan(psi); index=double(index); index(index==0)=NaN;
psi(isnan(psi))=[];
cos_psi_integracao=cos_psi.*index; cos_psi_integracao(isnan(cos_psi_integracao))=[];
hQ_integracao=hQ.*index; hQ_integracao(isnan(hQ_integracao))=[];
DgQ_integracao=DgQ.*index; DgQ_integracao(isnan(DgQ_integracao))=[];
AQ_integracao=AQ.*index; AQ_integracao(isnan(AQ_integracao))=[];

l=R_2*sin(psi/2);

% Termo Mi1
Mi1(p)=sum(AQ_integracao.*(((hQ_integracao-hP(p))./(l.^3)).*DgQ_integracao),'omitnan'); % m/s2

end

Mi1=1/(2*pi)*Mi1;  % m/s2

toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Exportação dos resultados %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Salvando o arquivo .tif
Mi1_grade=reshape(Mi1,size(Eta_res_Mi0_MDE));
lat_S_grade=min(min(grade_lat_calc))*180/pi-(resolucao_graus/2);
lat_N_grade=max(max(grade_lat_calc))*180/pi+(resolucao_graus/2);
lon_W_grade=min(min(grade_long_calc))*180/pi-(resolucao_graus/2);
lon_E_grade=max(max(grade_long_calc))*180/pi+(resolucao_graus/2);
Mi1_arquivo=georasterref('RasterSize',size(Mi1_grade),'LatitudeLimits',[lat_S_grade,lat_N_grade],'LongitudeLimits',[lon_W_grade,lon_E_grade]);
geotiffwrite('Mi1_WG300_ERTM_V1_1m.tif',flipud(Mi1_grade),Mi1_arquivo)

load handel
sound(y,Fs)
