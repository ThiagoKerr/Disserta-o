clc; clear all; format long g; tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Entrada de dados %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Importando as grades
Dg_res=double(imread('Grade_DG_res_3DLSC_ERTM_V1.tif'))*1D-5;         % distúrbios de gravidade de Molodenskii residuais em m2/s2 com resolução espacial normal
H=double(imread('./imagens/MDE_COLORADO_RECORTADO.tif'));        % MDE na resoução original dele

resolucao_graus=0.016666667;      % Resolução espacial normal em °
resolucao_graus_MDE=0.00083333333;   % Resolução espacial normal em ° do MDE

% Coordenadas geodésicas do canto superior esquerdo do primeiro pixel da imagem no canto NW - superior esquerdo - depois passando para o centro geométrico do pixel
latN=40.0083333333333471-(resolucao_graus/2);     % latN
lonW=-110.0083333333333542+(resolucao_graus/2);    % longW

latN_MDE=41-(resolucao_graus_MDE/2);      % latN
lonW_MDE=-111+(resolucao_graus_MDE/2);    % longW

% Importando estações de análise
estacao_GSVS=load('Pontos_GSVS.txt');

% Importando Zeta XGM 2016 - 300 (m),	Termo de grau zero IHRS (m),	Zeta RTM (m) - EARTH2014+ERTM2160,	Zeta NAVD88 (m) e Zeta médio 12 soluções (m)

Zeta_outros=load('Dados_ZETA_total.txt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Preparação para o cálculo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Varredura de variáveis
i=1;  % índice para cada raio de integração
j=1;  % índice para cada grau de modificação
graus_de_modificacao=30:30:300;
raios=0:0.1:4;
Grafico_GSVS17=zeros(length(raios),length(graus_de_modificacao));
Grafico_media12=zeros(length(raios),length(graus_de_modificacao));

for ppp=30:30:300   % grau de modificação
for pp=0:0.1:4      % raio de integração

grau_modificacao=ppp; % Grau de modificação do núcleo da função de Stokes
raio_integracao_graus=pp;       % Raio de integração em °

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Grades de coordenadas geodésicas da área de integração para a grade de distúrbios residuais e MDE
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

% Pontos de cálculo
latP=deg2rad(estacao_GSVS(:,1));    longP=deg2rad(estacao_GSVS(:,2));

% Distúrbio de gravidade e altitude normal dos pontos de cálculo P
DgP=interp2(grade_long,grade_lat,Dg_res,longP,latP,'cubic');
HP=interp2(grade_long_MDE,grade_lat_MDE,H,longP,latP,'cubic');    % interpolação para obtenção do valor de HN do MDE nas posições de cálculo

% Gravidade normal no Teluroide dos pontos de cálculo P
sin_latP=sin(latP);
gama0=(9.7803267715*(1+0.0052790414*(sin_latP.^2)+0.0000232718*(sin_latP.^4)+0.0000001262*(sin_latP.^6)+0.0000000007*(sin_latP.^8)));  % m/s2
gama=gama0.*(1-2*HP.*(1+0.00335281068118-2*0.00335281068118*(sin_latP.^2)+0.00344978600308)/6378137+(3*HP.^2)/6378137^2); % m/s2 - GRS80

% Pontos de integração
latQ=grade_lat(:);       longQ=grade_long(:);     DgQ=Dg_res(:);
sin_latQ=sin(latQ);     cos_latQ=cos(latQ);

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

Eta_res=zeros(length(latP),1); % alocando espaço para ganhar velocidade

% �?reas
Ak=abs((R^2)*res_rad*(sin(latQ+res_rad_div_2)-sin(latQ-res_rad_div_2)));

% Contribuição da zona mais interna
Etai=(R./gama).*sqrt((cos(latP)*res_rad^2)/pi).*DgP;

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
L1=find(abs((latP(p)-latQ_integracao))==min(abs((latP(p)-latQ_integracao))));
L2=find(abs((longP(p)-longQ_integracao))==min(abs((longP(p)-longQ_integracao))));
pos_p=intersect(L1,L2);

% Função de Hotine com modificação de Wang e Gore (1969)
s2=sin((latP(p)-latQ_integracao)./2).*sin((latP(p)-latQ_integracao)./2)+sin((longP(p)-longQ_integracao)./2).*sin((longP(p)-longQ_integracao)./2).*cos(latP(p)).*cos(latQ_integracao);
s=s2.^0.5;
H_psi=1./s-log(1+(1./s));
H_psi(pos_p)=NaN;   % desconsiderando S-psi do ponto de cálculo (infinito)

Pn = legendreN(grau_modificacao,cos_psi_integracao); % Polinômios de Legendre de grau máximo n e x = cos(psi)
H_psi_SH=cte_mod*Pn(:,3:end)';                       % começando de 3 para dispensar os termos de grau 0 e 1

% Anomalia de altura residual
% Eta_res(p)=sum((fator_F_psi.*(S_psi-S_psi_SH').*Ak_integracao.*AgQ_integracao)/cte(p),'omitnan')+Etai(p);
Eta_res(p)=sum(((H_psi-H_psi_SH').*Ak_integracao.*DgQ_integracao)/cte(p),'omitnan')+Etai(p);

end

% Restore, Discrepâncias e Desvios-padrão para cada raio de integração e grau de modificação
Zeta=Eta_res+Zeta_outros(:,1)+Zeta_outros(:,2)+Zeta_outros(:,3);  % Zeta_residual + Zeta_MGG + Zeta_termo de grau zero + Zeta_RTM
D_Zeta_GSVS17=Zeta-Zeta_outros(:,4);   % discrepância em relação ao GSVS17
D_Zeta_media12=Zeta-Zeta_outros(:,5);   % discrepância em relação à média das 12 soluções

Desv_pad_D_Zeta_GSVS17=std(D_Zeta_GSVS17)*100; % cm
Desv_pad_D_Zeta_media12=std(D_Zeta_media12)*100; % cm
Grafico_GSVS17(i,j)=Desv_pad_D_Zeta_GSVS17;
Grafico_media12(i,j)=Desv_pad_D_Zeta_media12;
i=i+1;     % índice i atualizado para cada raio de integração

end
i=1;
j=j+1;     % índice j atualizado para cada grau de modificação
end
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Geração dos gráficos %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1);
plot(raios,Grafico_GSVS17(:,1),'DisplayName','WG 30','LineWidth',1.5);   hold on
plot(raios,Grafico_GSVS17(:,2),'DisplayName','WG 60','LineWidth',1.5);
plot(raios,Grafico_GSVS17(:,3),'DisplayName','WG 90','LineWidth',1.5);
plot(raios,Grafico_GSVS17(:,4),'DisplayName','WG 120','LineWidth',1.5);
plot(raios,Grafico_GSVS17(:,5),'DisplayName','WG 150','LineWidth',1.5);
plot(raios,Grafico_GSVS17(:,6),'DisplayName','WG 180','LineWidth',1.5);
plot(raios,Grafico_GSVS17(:,7),'DisplayName','WG 210','LineWidth',1.5);
plot(raios,Grafico_GSVS17(:,8),'DisplayName','WG 240','LineWidth',1.5);
plot(raios,Grafico_GSVS17(:,9),'DisplayName','WG 270','LineWidth',1.5);
plot(raios,Grafico_GSVS17(:,10),'DisplayName','WG 300','LineWidth',1.5);  hold off
title('Standard deviations from GNSS-levelling comparisons');
grid on; grid minor; lgd = legend; lgd.NumColumns = 5;
xlabel('Cap size (�)','FontName','Times New Roman','FontSize',14)
ylabel('Stardard Deviation (cm)','FontName','Times New Roman','FontSize',14)
set(gca,'FontName','Times New Roman','FontSize',14)

figure(2);
plot(raios,Grafico_media12(:,1),'DisplayName','WG 30','LineWidth',1.5);   hold on
plot(raios,Grafico_media12(:,2),'DisplayName','WG 60','LineWidth',1.5);
plot(raios,Grafico_media12(:,3),'DisplayName','WG 90','LineWidth',1.5);
plot(raios,Grafico_media12(:,4),'DisplayName','WG 120','LineWidth',1.5);
plot(raios,Grafico_media12(:,5),'DisplayName','WG 150','LineWidth',1.5);
plot(raios,Grafico_media12(:,6),'DisplayName','WG 180','LineWidth',1.5);
plot(raios,Grafico_media12(:,7),'DisplayName','WG 210','LineWidth',1.5);
plot(raios,Grafico_media12(:,8),'DisplayName','WG 240','LineWidth',1.5);
plot(raios,Grafico_media12(:,9),'DisplayName','WG 270','LineWidth',1.5);
plot(raios,Grafico_media12(:,10),'DisplayName','WG 300','LineWidth',1.5);  hold off
title('Standard deviations from mean solution comparisons');
grid on; grid minor; lgd = legend; lgd.NumColumns = 5;
xlabel('Cap size (�)','FontName','Times New Roman','FontSize',14)
ylabel('Stardard Deviation (cm)','FontName','Times New Roman','FontSize',14)
set(gca,'FontName','Times New Roman','FontSize',14)

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