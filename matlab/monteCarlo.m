%% 2022, BP - vyvoj digitalniho odometru na platforme arduino - metoda Monte Carlo pro vypocet teoreticke presnosti
% Vojtech Remes

clc; format long g; clear all;

% nazev vystupniho souboru
nazev_vystupu = "MC3D_fin.txt";

% delka kterou ma vozidlo v rozboru presnosti ujet
delka_k_ujeti = 100;

% velikost souboru pro ktery je pocitana vyberova smerodatna odchylka
MC = 10000;

% polomery kruznic po kterych se vozidlo pohybuje
r_1_vec = [5, 10, 50, 250, 1000, 2000, 3500, Inf];
% r_1_vec = Inf;
% delka mezi koleckem a zavesem - T(ouchPoint)-M(ountPoint)
    d = 0.1;         
sig_d = 0.001;
% delka od zavesu k referencnimu bodu vozidla - M(ountPoint)-R(eferencePoint)
    m = 0.0;        
sig_m = 0.001;
% delka od zavesu k bodu kde normala k trajektorii prochazi stredem otaceni
% M(ountPoint)-K(olmyPoint)
    t = 0.4;        
sig_t = 0.001;
% polomer odvalovaciho kolecka
    r_k = 0.05;     
sig_r_k = 0.0001;
% enkoder pro urceni ujete delky
    N_d = 360;      % rozliseni delkoveho enkoderu
    n_d = 11;       % urcuje delku integracniho kroku 'k'
sig_n_d = 1;
    k   = r_k * 2 * pi() * n_d / N_d;
% rozliseni smeroveho enkoderu
    N_s = 4000;
sig_n_s = 1;
% naklon v podelnem smeru -- ovlivni ujetou delku
    zeta_x = 0;        
sig_zx     = pi()/180*(1/60);
% naklon v pricnem smeru -- ovlivni polomer zataceni
    zeta_y = 0;
sig_zy     = pi()/180*(1/60);
% uhel od predozadni osy vozidla ke spojnici M(ountPoint)-R(eferencePoint)
gama = 0;               

% format vystupni zpravy
vypis_hlavicka = "=== === === r = %5d m === === ===\n velikost souboru = %7d\n ujeta delka = %7.2f m\n";
vypis_rozmery = " d      = %5.3f  +- %5.3f m\n t      = %5.3f  +- %5.3f m\n m      = %5.3f  +- %5.3f m\n r_k    = %6.4f +- %6.4f m\n";
vypis_naklony = " zeta_x = %5.3f  +- %5.3f °\n zeta_y = %5.3f  +- %5.3f °\n\n";
vypis_cteni =   " N_d = %5d  n_d = %2d    +- %3.1f\n N_s = %5d  n_s = %5.2f +- %3.1f\n\n";
vypis_sour =    " --- TEORETICKE\n X = %8.3f +- %5.3f m   Y = %8.3f +- %5.3f m\n";
vypis_real =    " --- PRUMERNE VYPOCTENE\n X = %8.3f m            Y = %8.3f m\n=== === === === === === === === ===\n\n\n";
vypis = vypis_hlavicka + vypis_rozmery + vypis_naklony + vypis_cteni + vypis_sour + vypis_real;

for index1 = 1:length(r_1_vec)
    r = r_1_vec(index1);
[n_st, n_sm] = vypocitej_cteniSmerove3D(N_s, d, t, r, zeta_y);
    theta_teoreticka = n_st * 2 * pi() / N_s;
    theta_merena =     n_sm * 2 * pi() / N_s;

%% SOURADNICOVE VYPOCTY
% souradnice pocatecniho bodu a smernik 
sigT0_T1 = 0;
X_T1 = 0;
Y_T1 = 0;
n = round( delka_k_ujeti/k );

% vypocet teoreticke polohy R(eferencePoint) na konci
[~, ~, blb, X_Rnt, Y_Rnt] = odometr3D(X_T1, Y_T1, sigT0_T1, zeta_x, zeta_y, d, t, m, gama, k, theta_teoreticka, n);

% Monte Carlo vypocet
X_TnMC(1:MC) = 0; Y_TnMC(1:MC) = 0;
sig_TnMC(1:MC) = sigT0_T1;
X_RnMC(1:MC) = 0; Y_RnMC(1:MC) = 0;
ujeta_delka = 0;
for index2 = 1:MC
    d_now = d + sig_d * (2 * rand - 1);
    t_now = t + sig_t * (2 * rand - 1);
    m_now = m + sig_m * (2 * rand - 1);
    rk_now = r_k + sig_r_k * (2 * rand - 1);
    zetax_now = zeta_x + sig_zx * (2 * rand - 1);
    zetay_now = zeta_y + sig_zy * (2 * rand - 1);
    ujeta_delka = 0;
    while (ujeta_delka < delka_k_ujeti)
        nd_now = round( n_d  + sig_n_d * (2 * rand - 1) );
        ns_now = round( n_st + sig_n_s * (2 * rand - 1) );
        k_now     = 2 * pi() * nd_now / N_d * rk_now;
        theta_now = 2 * pi() * ns_now / N_s;
        ujeta_delka = ujeta_delka + 2 * pi() * n_d / N_d * r_k;
        [X_TnMC(index2), Y_TnMC(index2), sig_TnMC(index2), X_RnMC(index2), Y_RnMC(index2)] = ...
            odometr3D(X_TnMC(index2), Y_TnMC(index2), sig_TnMC(index2), zetax_now, zetay_now, d_now, t_now, m_now, 0, k_now, theta_now, 1);
    end
end
[sdXR, mXR] = std(X_RnMC, 1, 'all');
[sdYR, mYR] = std(Y_RnMC, 1, 'all');

fid = fopen(nazev_vystupu, "a+");
fprintf(fid, vypis, r, MC, ujeta_delka, d, sig_d, t, sig_t, m, sig_m, r_k, sig_r_k,...
    180/pi()*(zeta_x), 180/pi()*(sig_zx), 180/pi()*(zeta_y), 180/pi()*(sig_zy), N_d, n_d, sig_n_d, N_s, n_st, sig_n_s,...
    X_Rnt, sdXR, Y_Rnt, sdYR, mXR, mYR);
fclose(fid);

end

%%
function [X_Tn, Y_Tn, sigTn_1_Tn, X_Rn, Y_Rn] = odometr3D(X_T0, Y_T0, sig_0, zeta_x, zeta_y, d, t, m, gama, k, theta, n)
% vypocita souradnice XY TouchPointu a ReferencePointu pri pohybu po
% elementu krivky v delce 'n'*'k'
% VSTUP X_T0, Y_T0  souradnice bodu T pred zacatkem pohybu
%       sig_0       smernik odometru v case 0
%       zeta_x      okamzity podelny naklon
%       zeta_y      okamzity pricny naklon
%       d           vodorovna delka mezi TouchPointem a MountPointem
%       t           vodorovna delka mezi MountPointem a bodem kde
%                   normala k trajektorii prochazi stredem otaceni
%       m           vodorovna delka mezi MountPointem a ReferencePointem
%       gama        uhel od predozadni osy vozidla a spojnici
%                   MountPointu a ReferencePointu
%       k           delka integracniho kroku (zavisi na polomeru kolecka a
%                   poctu pulzu delkoveho enkoderu)
%       theta       okamzity (mereny) uhel natoceni
%       n           pocet kroku 'k', kdy se odometr pohybuje po elementu krivky

X_Ti = X_T0; Y_Ti = Y_T0; sigTi_1_Ti = sig_0;
for i = 1:n
    fi = 2 * asin( k * sin(theta) / ( 2*( d * cos(theta) + t) * cos(zeta_y) ) );
    sigTi_1_Ti = sigTi_1_Ti + fi;
    dXTi_1_T = k * cos(zeta_x) * cos( sigTi_1_Ti );
    dYTi_1_T = k * cos(zeta_x) * sin( sigTi_1_Ti );
    X_Ti = X_Ti + dXTi_1_T;
    Y_Ti = Y_Ti + dYTi_1_T;
end
X_Tn = X_Ti; Y_Tn = Y_Ti; sigTn_1_Tn = sigTi_1_Ti;

X_Mn = X_Tn + d * cos(zeta_x) * cos( sigTn_1_Tn + fi/2 );
Y_Mn = Y_Tn + d * cos(zeta_x) * sin( sigTn_1_Tn + fi/2 );

X_Rn = X_Mn + m * cos( sigTn_1_Tn + fi/2 + theta + gama );
Y_Rn = Y_Mn + m * sin( sigTn_1_Tn + fi/2 + theta + gama );
end


function [n_s_teoreticke, n_s_merene] = vypocitej_cteniSmerove3D(N_s, d, t, r, zeta_y)
% vypocita teoreticky okamzity uhel natoceni a 
%              mereny okamzity uhel natoceni zaokrouhleny na cely pocet pulzu
% VSTUP N_s     rozliseni smeroveho enkoderu
%       d       vodorovna delka mezi TouchPointem a MountPointem
%       t       vodorovna delka mezi MountPointem a bodem kde
%               normala k trajektorii prochazi stredem otaceni
%       r       polomer zataceni (trajektorie)
%       zeta_y  okamzity pricny naklon


n_s = N_s * t / ( 2* pi() * r );
theta = 2 * pi() * n_s / N_s;

n_s0 = n_s + 2;     % pro aspon jeden pruchod while cyklem
while abs( n_s - n_s0 ) > 0.5
        n_s0 = n_s;
        theta0 = theta;
        theta = asin( ( d * cos(theta0) + t ) / r * cos(zeta_y) );
        n_s  = N_s * theta / ( 2 * pi() );
end
n_s_teoreticke = n_s;
n_s_merene = round(n_s);


end
