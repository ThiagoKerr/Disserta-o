clc; clear all; format long g; tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Entrada de dados %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Importando as grades
Dg_res=double(imread('Grade_Polinomial_GO200.tif'))*1D-5;  % distúrbios de gravidade residuais em m2/s2 com resolução espacial normal
H=double(imread('MDE.tif'));        % MDE na resoução original dele

resolucao_graus=0.04999999999999991257;;         % Resolução espacial normal em °
resolucao_graus_MDE=0.0002777777777778145862;   % Resolução espacial normal em ° do MDE
raio_integracao_graus=0.7;          % Raio de integração em °

% Coordenadas geodésicas do canto superior esquerdo do primeiro pixel da imagem no canto NW - superior esquerdo - depois passando para o centro geométrico do pixel
latN=-20.0949000000000026-(resolucao_graus/2);     % latN
lonW=-53.4334999999999951+(resolucao_graus/2);    % longW

latN_MDE=-19.3998611111216306-(resolucao_graus_MDE/2);      % latN
lonW_MDE=-54.2001388888722317+(resolucao_graus_MDE/2);    % longW

grau_modificacao=200; % Grau de modificação do núcleo da função de Hotine

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Preparação para o cálculo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Grades de coordenadas geodésicas da área de integração para a grade de anomalias residuais e MDE
[lin,col]=size(Dg_res);
[lin_MDE,col_MDE]=size(H);

latS=(latN-(lin-1)*resolucao_graus);
lonE=(lonW+(col-1)*resolucao_graus);
latS_MDE=(latN_MDE-(lin_MDE-1)*resolucao_graus_MDE);
lonE_MDE=(lonW_MDE+(col_MDE-1)*resolucao_graus_MDE);

grade_lat=(latN:-resolucao_graus:latS)'.*ones(1,col)*pi/180; % rad
grade_long=(lonW:resolucao_graus:lonE).*ones(lin,1)*pi/180;  % rad
grade_lat_MDE=(latN_MDE:-resolucao_graus_MDE:latS_MDE)'.*ones(1,col_MDE)*pi/180; % rad
grade_long_MDE=(lonW_MDE:resolucao_graus_MDE:lonE_MDE).*ones(lin_MDE,1)*pi/180;  % rad

% Grade de coordenadas geodésicas da área de cálculo
raio_pixels=ceil(raio_integracao_graus/resolucao_graus); % raio em pixels arredondado para cima para dar inteiro

grade_lat_calc=grade_lat(raio_pixels+1:lin-raio_pixels , raio_pixels+1:col-raio_pixels);    % descontando os raios de integração da grade de área de integração
grade_long_calc=grade_long(raio_pixels+1:lin-raio_pixels , raio_pixels+1:col-raio_pixels);  % descontando os raios de integração da grade de área de integração

% Distúrbio de gravidade dos pontos de cálculo P

gradeDgP_calc=Dg_res(raio_pixels+1:lin-raio_pixels , raio_pixels+1:col-raio_pixels);    % descontando os raios de integração da grade de área de integração;

% Coordenadas em forma de vetor e cálculo da gravidade normal da grade de cálculo
% Pontos de cálculo
latP=grade_lat_calc(:);      longP=grade_long_calc(:);        DgP=gradeDgP_calc(:);        
HP=interp2(grade_long_MDE,grade_lat_MDE,H,longP,latP,'cubic');    % interpolação para obtenção do valor de HN do MDE nas posições de cálculo
sin_latP=sin(latP);
gama0=(9.7803267715*(1+0.0052790414*(sin_latP.^2)+0.0000232718*(sin_latP.^4)+0.0000001262*(sin_latP.^6)+0.0000000007*(sin_latP.^8)));  % m/s2
gama=gama0.*(1-2*HP.*(1+0.00335281068118-2*0.00335281068118*(sin(latP).^2)+0.00344978600308)/6378137+(3*HP.^2)/6378137^2); % m/s2 - GRS80

% Pontos de integração
latQ=grade_lat(:);       longQ=grade_long(:);     DgQ=Dg_res(:);
sin_latQ=sin(latQ);      cos_latQ=cos(latQ);

% Constantes
res_rad=resolucao_graus*pi/180;              % Resolução em rad
raio_integ_rad=raio_integracao_graus*pi/180; % raio de integração em rad
R=6371000;                                   % Raio da esfera de mesmo volume do Elpsoide em metros
res_rad_div_2=(res_rad/2);
cte=4*pi*R*gama;
n_mod=2:grau_modificacao;
cte_mod=(2*n_mod+1)./(n_mod+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Solução do PVCG gravimétrico fixado %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Zeta_res=zeros(length(latP),1); % alocando espaço para ganhar velocidade

% Áreas
Ak=abs((R^2)*res_rad*(sin(latQ+res_rad_div_2)-sin(latQ-res_rad_div_2)));

% Contribuição da zona mais interna
Zetai=(R./gama).*sqrt((cos(latP)*res_rad^2)/pi).*DgP;

for p=1:length(latP)

% Filtrando os pontos de integração dentro da área de integração
cos_psi=sin_latQ*sin_latP(p)+cos_latQ*cos(latP(p)).*cos(longQ-longP(p));
psi=real(acos(cos_psi));

psi(psi>raio_integ_rad)=NaN; index=~isnan(psi); index=double(index); index(index==0)=NaN;
psi(isnan(psi))=[];
cos_psi_integracao=cos_psi.*index; cos_psi_integracao(isnan(cos_psi_integracao))=[];
latQ_integracao=latQ.*index; latQ_integracao(isnan(latQ_integracao))=[];
longQ_integracao=longQ.*index; longQ_integracao(isnan(longQ_integracao))=[];
DgQ_integracao=DgQ.*index; DgQ_integracao(isnan(DgQ_integracao))=[];
Ak_integracao=Ak.*index; Ak_integracao(isnan(Ak_integracao))=[];

% Detectando a localização do ponto de cálculo nos vetores de coordenadas
L1=find((latP(p)-latQ_integracao)==0);
L2=find((longP(p)-longQ_integracao)==0);
pos_p=intersect(L1,L2);

% Função de Hotine com modificação de Wang e Gore (1969)
s2=sin((latP(p)-latQ_integracao)./2).*sin((latP(p)-latQ_integracao)./2)+sin((longP(p)-longQ_integracao)./2).*sin((longP(p)-longQ_integracao)./2).*cos(latP(p)).*cos(latQ_integracao);
s=s2.^0.5;
H_psi=1./s-log(1+(1./s));
H_psi(pos_p)=NaN;   % desconsiderando H-psi do ponto de cálculo (infinito)

Pn = legendreN(grau_modificacao,cos_psi_integracao); % Polinômios de Legendre de grau máximo n e x = cos(psi)
H_psi_SH=cte_mod*Pn(:,3:end)';                       % começando de 3 para dispensar os termos de grau 0 e 1

% Anomalia de altura residual
Zeta_res(p)=sum(((H_psi-H_psi_SH').*Ak_integracao.*DgQ_integracao)/cte(p),'omitnan')+Zetai(p);

end

toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Exportação dos resultados %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transformando em formato de matriz
Zeta_res_grade=reshape(Zeta_res,size(grade_lat_calc));

% Salvando o arquivo .tif
lat_S_grade=min(min(grade_lat_calc))*180/pi-(resolucao_graus/2);
lat_N_grade=max(max(grade_lat_calc))*180/pi+(resolucao_graus/2);
lon_W_grade=min(min(grade_long_calc))*180/pi-(resolucao_graus/2);
lon_E_grade=max(max(grade_long_calc))*180/pi+(resolucao_graus/2);
Mi0_arquivo=georasterref('RasterSize',size(Zeta_res_grade),'LatitudeLimits',[lat_S_grade,lat_N_grade],'LongitudeLimits',[lon_W_grade,lon_E_grade]);
geotiffwrite('Zeta_residual_Mi0_GO200_WG200.tif',flipud(Zeta_res_grade),Mi0_arquivo)

load handel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Função rápida dos Polinômios de Legendre de grau n e ordem zero %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Pn = legendreN(nmax,t)
Pn = zeros(length(t), nmax+1);
Pn(:, 1) = 1;    
Pn(:, 2) = t; 
for n = 2:nmax 
    Pn(:, n+1) = (1 - n) / n * Pn(:, n - 1) + (2 * n - 1) / n * t .* Pn(:, n);
end
end