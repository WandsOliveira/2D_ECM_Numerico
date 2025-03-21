% Parâmetros geométricos (em metros)
r1 = 53e-3/2;          % Raio do tubo interno
R = 150e-3/2;            % Raio da câmara externa
r3 = R;
L = 62.5e-3;
d1 = 2*r1;           % Diâmetro do tubo interno

% Propriedades do ar (a 20°C)
rho = 1.21;          % Densidade do ar (kg/m³)
c = 343;             % Velocidade do som no ar (m/s)
nu = 1.5e-5;         % Viscosidade cinemática do ar (m²/s)

% Vetor de frequência
df = 10;
freq = 1:df:10000;          % Frequência de análise (Hz)
omega = 2*pi*freq;      % Frequência angular (rad/s)
k0 = omega/c;        % Número de onda
STL = zeros(1,length(freq));
% Cálculo da porosidade (Eq. C.1)
% Número de modos a considerar
N = 2;              % Número de modos para a análise

