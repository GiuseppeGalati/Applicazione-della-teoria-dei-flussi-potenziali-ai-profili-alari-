%% Descrizione
%  Il codice applica la teoria dei flussi potenziali per definire 
%  il flusso e le azioni dinamiche su un profilo alare definito dall'utente
%  attraverso la teoria di Witoszynski. I parametri del campo sono definiti
%  dall'utente tramite input (ad eccezione di V_infinito, da definire
%  preliminarmente). Analogamente, il profilo in esame è definito
%  dall'utente mediante trasformazione conforme definita da Muller, che
%  permette di ottenere profili di Joukowski. Per ottenere risultati
%  migliori si consiglia di annullare l'angolo di incidenza.



%% Definizione preliminare costanti 
clearvars
Vinf = 7;                       
alfa_deg = input('Angolo di incidenza: ');
alfa = alfa_deg*pi/180;
beta_deg = input('Angolo primo asse: ');                  
beta = beta_deg*pi/180;
gamma_deg = input('Angolo secondo asse: ');                  
gamma = gamma_deg*pi/180;
r0 = input('Raggio circonferenza primitiva: ');                       
mi = r0^2*(2*pi*Vinf);
p = input('Parametro p di Muller: ');
Z_Muller = r0*exp(i*beta)-p*exp(i*gamma);   % Parametro z0 di Muller                             
angle = linspace (0,2*pi,200);  
rv = linspace (1, 8, 200);
tetav = linspace(0,2*pi,200);
[r,teta] = meshgrid(rv,tetav);              % Costruzione mesh
xv = linspace(-2,2,200); 
yv = xv;
[x,y] = meshgrid(xv,yv);                
z = x + 1i*y;                               % Piano complesso
gap = 0.09;                     % se la soluzione mostra curve indesiderate
                                % provare ad agire su questo parametro
xM = 0;                         
yM = 0;                         % Coordinate centro della circonferenza
M = xM + 1i*yM;                 % Localizzazione centro della circonferenza
rM = sqrt(xM^2 + yM^2);
tetaM = atan2(yM,xM);

for rn = 1:length(r)
    for tetan = 1:length(teta)
        if r(rn,tetan) <= r0+ gap
            r(rn,tetan) = NaN;      % Eliminazione punti interni al cerchio
            teta(rn,tetan) = NaN; 
        end             
    end
end

for xn = 1:length(x)
    for yn = 1:length(y)
        if abs(z(xn,yn)-M) <=  r0 - gap 
            z(xn,yn) = NaN;  
        end       
    end
end

%% Trasformazione di Joukowski
xc = r0*cos(angle)+xM;
yc = r0*sin(angle)+yM;
zc = xc + 1i.*yc; 
Z_prof = zc + (((p.*exp(i*gamma)).^2)./(zc+Z_Muller));
Xprofilo = real(Z_prof); 
Yprofilo = imag(Z_prof); 

%% Trasformazione del campo con Joukowski

figure
axis equal 
axis ([-5 5 -5 5])
hold on
xv = linspace(-2,2,200); 
yv = xv;
[x,y] = meshgrid(xv,yv);
z = x + 1i*y;
Z = z + (((p.*exp(i*gamma)).^2)./(z + Z_Muller));   %Trasformazione secondo
                                                    %formula di Muller
for xn = 1:length(x)
    for yn = 1:length(y)
        if abs(z(xn,yn)-M) <=  r0 - gap 
            z(xn,yn) = NaN;  
        end       
    end
end

%% Definizione del campo e plot curve di livello per funzione di corrente
W = - Vinf*exp(-1i*alfa).*(z-M) - (mi./(2*pi.*(z-M)))*exp(1i*alfa)- 4*1i*Vinf*r0*sin(beta).*((z-M).^(1/2) - ((r0^(1/2))*exp(1i*beta/2)))./((z-M).^(1/2) + ((r0^(1/2))*exp(1i*beta/2))); %L'inserimento di M permette di spostare il cilindo in un punto qualsiasi del piano di Gauss
psi = imag (W);
contour(real(Z),imag(Z),psi,[2:2:10]);              
xv = linspace(-8,8,200); 
yv = xv;
[x,y] = meshgrid(xv,yv);
z = x + 1i*y;
Z = z + (((p.*exp(i*gamma)).^2)./(z + Z_Muller));
for xn = 1:length(x)
    for yn = 1:length(y)
        if abs(z(xn,yn)-M) <=  r0 - gap 
            z(xn,yn) = NaN;  
        end       
    end
end
W = - Vinf*exp(-1i*alfa).*(z-M) - (mi./(2*pi.*(z-M)))*exp(1i*alfa)- 4*1i*Vinf*r0*sin(beta).*((z-M).^(1/2) - ((r0^(1/2))*exp(1i*beta/2)))./((z-M).^(1/2) + ((r0^(1/2))*exp(1i*beta/2))); %L'inserimento di M permette di spostare il cilindo in un punto qualsiasi del piano di Gauss
psi = imag (W);
contour(real(Z),imag(Z),psi,[0 0]);
P = (r.^2 + rM^2 - (2.*r.*rM.*cos(teta-tetaM))).^0.5;
phi = atan2((r.*sin(teta-beta)-rM.*sin(tetaM-beta)),(r.*cos(teta-beta)-rM.*cos(tetaM-beta)))./2;
psiwito = - 4*Vinf*r0*sin(beta).*((P-r0)./(P+r0+(2.*r0^0.5.*P.^(1/2).*cos(phi)))).*cos(alfa);         
psi0 =  - Vinf.*((r.*sin(teta-alfa) -rM*sin(tetaM-alfa))-((r0^2.*(r.*sin(teta-alfa) -rM*sin(tetaM-alfa)))./(r.^2 +rM^2-(2.*r.*rM.*cos(teta-tetaM))))) + psiwito;
x = r.*cos(teta);
y = r.*sin(teta);
z = x + 1i*y;
Z = z + (((p.*exp(i*gamma)).^2)./(z + Z_Muller));
contour (real(Z),imag(Z),psi0,[-100:2:100]);
contour (real(Z),imag(Z),psi0,[0 0]);
fill (Xprofilo,Yprofilo,'k');

%% Calcolo Azioni Dinamiche
clear z;
z0 = Z_Muller;                          
rho = 1.225;                            %Densità del fluido    
dr = 0.001;                             %modificare il parametro per 
                                        %variare l'accuratezza del
                                        %risultato 
disp('Immettere scelta del caso')
disp('---')
disp('scegliere tra "W" (confronto con Witoszynski) o "M" (confronto con Muller)')
disp('---')
choise = input('Immettere scelta: ');

%% Calcolo azioni dinamiche con formule analitiche su PROFILI RETTILINEI
P_Wito = 8.*rho.*Vinf.^2.*r0.*sin(beta).*cos(beta).^2;
R_Wito = 8.*rho.*Vinf.^2.*r0.*sin(beta).^2.*cos(beta);
cd_Muller = 0.48.*sin(beta).^2;
cl_Wito = 4.*sin(beta).*cos(beta).^2;
cl_Muller = 4.*sin(beta);
cd_Wito = 4.*sin(beta).^2.*cos(beta);
R_Muller = cd_Muller.*0.5.*rho.*Vinf^2.*4.*r0;

        %% Calcolo di Lift e Drag tramite Blasius
        % si definisce una curva g lungo la quale effettuare l'integrazione
        % si osserva che dr deve essere diverso da 0 affinchè il calcolo
        % sia effettuabile
        
switch choise
    
    case 'W'

        %%Definizione integrale
        g = @(teta) (r0+dr)*cos(teta) + (r0+dr)*1i*sin(teta);
        gprime = @(teta) -(r0+dr)*sin(teta) + (r0+dr)*1i*cos(teta);
        Int_Blasius = @(z) (((exp(-1i*alfa)-exp(1i*alfa).*(r0.^2)./(z.^2))+ exp(1i*alfa).*(4.*i.*r0.^(3/2).*sin(beta).*exp(i.*beta.*0.5)./((z.^0.5).*((z.^0.5)+(r0.^0.5).*exp(i.*beta.*0.5)).^2))).^2) .*(((z+z0).^2)./((z+z0).^2-(p.^2).*exp(i.*2.*gamma)));
        R_attrito = rho.*Vinf.^2.*r0.*(9/4).*pi.*i.*exp(-i.*beta).*sin(beta).^2;
        Forces = -0.5.*rho.*(Vinf.^2).*integral(@(t) Int_Blasius(g(t)).*gprime(t),0,2*pi)-R_attrito;
        Lift = real(Forces); 
        Drag = -imag (Forces); 
        %% Calcolo dei coefficienti

        cl = Lift./(4.*r0.*0.5.*rho.*Vinf.^2)
        cd = Drag./(4.*r0.*0.5.*rho.*Vinf.^2)

    case 'M'
        
        g = @(teta) (r0+dr)*cos(teta) + (r0+dr)*1i*sin(teta);
        gprime = @(teta) -(r0+dr)*sin(teta) + (r0+dr)*1i*cos(teta);
        Int_Blasius = @(z) (((exp(-1i*alfa)-exp(1i*alfa).*(r0.^2)./(z.^2))+ exp(1i*alfa).*(4.*i.*r0.^(3/2).*sin(beta).*exp(i.*beta.*0.5)./((z.^0.5).*((z.^0.5)+(r0.^0.5).*exp(i.*beta.*0.5)).^2))).^2) .*(((z+z0).^2)./((z+z0).^2-(p.^2).*exp(i.*2.*gamma)));
        R_attrito = rho.*Vinf.^2.*r0.*(9/4).*pi.*i.*exp(-i.*beta).*sin(beta).^2;
        Forces = -0.5.*rho.*(Vinf.^2).*integral(@(t) Int_Blasius(g(t)).*gprime(t),0,2*pi);
        Lift = real(Forces); 
        Drag = -imag (Forces); 
            
        %% Calcolo dei coefficienti

        cl = Lift./(4.*r0.*0.5.*rho.*Vinf.^2)
        cd = Drag./(4.*r0.*0.5.*rho.*Vinf.^2)
    
end

%% Diritti
%  Questo codice è stato prodotto da Giuseppe Galati, matricola 216291
%  nell'ambito dell' elaborato finale per il corso di Laurea Triennale
%  in Ingegneria Meccanica presso il Politecnico di Torino,
%  da titolo "Applicazione della teoria dei flussi
%  potenziali ai profili alari". Come tale, l'elaborato e tutti i codici
%  prodotti sono da considerarsi opere tutelate da diritto d'autore,
%  incluso quest'ultimo, che pur non essendo riportato nell'elaborato
%  finale in appendice, rappresenta un codice riassuntivo del lavoro svolto
%
% Copyright 2018 Giuseppe Galati

