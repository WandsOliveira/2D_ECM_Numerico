% ECM_PURO
clc; close all; clear
% Implementação do método de matriz de transferência 2D para MCPR
digits(50) 
% Importar parâmetros geométricos e físicos
Data;

%% ========== REGIÃO A ==========
Kr_A = get_radial_wavenumber_A(r1,N+1,2);
Kx_A = get_axial_wavenumber(Kr_A,k0,freq,1,'A',N+1);
%% ======== Regiao B =======
Kr_B = get_radial_wavenumber_A(r3,N+1,2);
Kx_B = get_axial_wavenumber(Kr_B,k0,freq,2,'B',N+1);
%% ======== Regiao C =======
Kr_C = get_radial_wavenumber_A(r1,N+1,2);
Kx_C = get_axial_wavenumber(Kr_C,k0,freq,2,'C',N+1);
%% Inicializar matriz A e vetor b
% X  = [A-,B+,B-,C+]
A = zeros(4*(N+1), 4*(N+1));
b = zeros(4*(N+1), 1);
%% etapa: Criar a matriz A, que sao os valores das integrais
%% e o vetor b, que é o vetor dos resultados
%% por fim resolver Ax = b, para obter os coeficientes da onda
STL = zeros(1,length(freq));
for j = 1:length(freq)
    %% Reiniciar A e b para cada frequencia
A = zeros(4*(N+1), 4*(N+1));
b = zeros(4*(N+1), 1);
% X  = [A-,B+,B-,C+]
%% ======== Regiao I ============
aux = 1;
for i = 1:N+1
    aux_col = (aux-1)*4; % Poe o valor no elemento correto da matriz A
    limite_superior = (aux / (N+1)) * r1; % Exemplo para r_m,p1
    limite_inferior = 0;
    % Coeficiente A- posicao: 1
    A(i,1 + aux_col) = calcIntegral(Kr_A(aux),limite_superior,'J0') - calcIntegral(Kr_A(aux),limite_inferior,'J0');
    
    % Coeficiente B+ posicao: 3
    A(i,2 + aux_col) = -(calcIntegral(Kr_B(aux),limite_superior,'J0') - calcIntegral(Kr_B(aux),limite_inferior,'J0'));
    
    % Coeficiente B- posicao: 4
    A(i,3 + aux_col) = -(calcIntegral(Kr_B(aux),limite_superior,'J0') - calcIntegral(Kr_B(aux),limite_inferior,'J0'));
    
    % Termo conhecido (A_+) vai para b
%     if i == 1
        b(i) = -(calcIntegral(Kr_A(aux),limite_superior,'J0') - calcIntegral(Kr_A(aux),limite_inferior,'J0')); % A0+ = 1
%     end
    aux = aux + 1;
end
%% ======== Regiao II ============
aux = 1;
for i = (N+1)+1 : 2*(N+1)
    r_mu = (aux / (N+1)) * r3;
    if r_mu <= r1
        aux_col = (aux-1)*4; % Poe o valor no elemento correto da matriz A
        limite_superior = (aux / (N+1)) * r3; % Exemplo para r_m,p1
        limite_inferior = 0;
        % Coeficiente A- posicao: 1
        A(i,1 + aux_col) = Kx_A(j,aux)*(calcIntegral(Kr_A(aux),limite_superior,'J0') - calcIntegral(Kr_A(aux),limite_inferior,'J0'));
        
        % Coeficiente B+ posicao: 
        A(i,2 + aux_col) = Kx_B(j,aux)*(calcIntegral(Kr_B(aux),limite_superior,'J0') - calcIntegral(Kr_B(aux),limite_inferior,'J0'));
        
        % Coeficiente B- posicao: 
        A(i,3 + aux_col) = -Kx_B(j,aux)*(calcIntegral(Kr_B(aux),limite_superior,'J0') - calcIntegral(Kr_B(aux),limite_inferior,'J0'));
        
        % Termo conhecido (A_+) vai para b
        if i == (N+1)+1
            b(i) = Kx_A(j,aux)*(calcIntegral(Kr_A(aux),limite_superior,'J0') - calcIntegral(Kr_A(aux),limite_inferior,'J0')); % A0+ = 1
        end
    else
        aux_col = (aux-1)*4; % Poe o valor no elemento correto da matriz A
        limite_superior = (aux / (N+1)) * r3; % Exemplo para r_m,p1
        limite_inferior = 0   ;    
        % Coeficiente A- posicao: 1
        A(i,1 + aux_col) = Kx_A(j,aux)*(calcIntegral(Kr_A(aux),r1,'J0') - calcIntegral(Kr_A(aux),0,'J0'));
        
        % Coeficiente B+ posicao: 
        A(i,2 + aux_col) = Kx_B(j,aux)*(calcIntegral(Kr_B(aux),limite_superior,'J0') - calcIntegral(Kr_B(aux),limite_inferior,'J0'));
        
        % Coeficiente B- posicao: 
        A(i,3 + aux_col) = -Kx_B(j,aux)*(calcIntegral(Kr_B(aux),limite_superior,'J0') - calcIntegral(Kr_B(aux),limite_inferior,'J0'));
        
        % Termo conhecido (A_+) vai para b
        if i == (N+1)+1
            b(i) = Kx_A(j,aux)*(calcIntegral(Kr_A(aux),r1,'J0') - calcIntegral(Kr_A(aux),0,'J0')); % A0+ = 1
        end
    end
    aux = aux + 1;
end
%% ======== Regiao III ============
aux = 1;
for i = 2*(N+1)+1 : 3*(N+1)
    aux_col = (aux-1)*4; % Poe o valor no elemento correto da matriz A
    limite_superior = (aux / (N+1)) * r1; % Exemplo para r_m,p1
    limite_inferior = 0;
    % Coeficiente C+ posicao: 1
    A(i,4 + aux_col) = calcIntegral(Kr_C(aux),limite_superior,'J0') - calcIntegral(Kr_C(aux),limite_inferior,'J0');
    
    % Coeficiente B+ posicao: 3
    A(i,2 + aux_col) = -exp(-1i*Kx_B(j,aux)*L)*(calcIntegral(Kr_B(aux),limite_superior,'J0') - calcIntegral(Kr_B(aux),limite_inferior,'J0'));
    
    % Coeficiente B- posicao: 4
    A(i,3 + aux_col) = -exp(1i*Kx_B(j,aux)*L)*(calcIntegral(Kr_B(aux),limite_superior,'J0') - calcIntegral(Kr_B(aux),limite_inferior,'J0'));
         
    aux = aux + 1;
end
%% ======== Regiao IV ============
aux = 1;
for i = 3*(N+1)+1 : 4*(N+1)
    r_mu = (aux / (N+1)) * r3;
    if r_mu <= r1
        aux_col = (aux-1)*4; % Poe o valor no elemento correto da matriz A
        limite_superior = (aux / (N+1)) * r3; % Exemplo para r_m,p1
        limite_inferior = 0;
        % Coeficiente C+ posicao: 4
        A(i,4 + aux_col) = -Kx_C(j,aux)*(calcIntegral(Kr_C(aux),limite_superior,'J0') - calcIntegral(Kr_C(aux),limite_inferior,'J0'));
        
        % Coeficiente B+ posicao: 
        A(i,2 + aux_col) = exp(-1i*Kx_B(j,aux)*L)*Kx_B(j,aux)*(calcIntegral(Kr_B(aux),limite_superior,'J0') - calcIntegral(Kr_B(aux),limite_inferior,'J0'));
        
        % Coeficiente B- posicao: 
        A(i,3 + aux_col) = -exp(1i*Kx_B(j,aux)*L)*Kx_B(j,aux)*(calcIntegral(Kr_B(aux),limite_superior,'J0') - calcIntegral(Kr_B(aux),limite_inferior,'J0'));
        
    else
        aux_col = (aux-1)*4; % Poe o valor no elemento correto da matriz A
        limite_superior = (aux / (N+1)) * r3; % Exemplo para r_m,p1
        limite_inferior = 0;
        % Coeficiente C+ posicao: 4
        A(i,4 + aux_col) = -Kx_C(j,aux)*(calcIntegral(Kr_C(aux),r1,'J0') - calcIntegral(Kr_C(aux),0,'J0'));
        
        % Coeficiente B+ posicao: 
        A(i,2 + aux_col) = exp(-1i*Kx_B(j,aux)*L)*Kx_B(j,aux)*(calcIntegral(Kr_B(aux),limite_superior,'J0') - calcIntegral(Kr_B(aux),limite_inferior,'J0'));
        
        % Coeficiente B- posicao: 
        A(i,3 + aux_col) = -exp(1i*Kx_B(j,aux)*L)*Kx_B(j,aux)*(calcIntegral(Kr_B(aux),limite_superior,'J0') - calcIntegral(Kr_B(aux),limite_inferior,'J0'));

    end
    aux = aux + 1;
end
%% ==== Resolvendo o sistema linear ============
X = A \ b;
%% Calcular a STL para cada frequência
STL(j) = -20 * log10(abs(X(4)));
% fprintf('A-: %d C+: %d C-: %d B+: %d \n', ...
%         abs(X(1)), abs(X(2)), abs(X(3)), abs(X(4)));


end

% Plotar a Perda de Transmissão
figure('Name','Perda de Transmissão')
hold on
plot(freq, STL, 'LineWidth', 1.5, 'DisplayName', sprintf('STL - Modo : %d', N));
xlabel('Frequência (Hz)')
ylabel('STL (dB)')
grid on
grid minor
legend show

% Função para calcular a integral definida
function integral_result = calcIntegral(k_r, limite_superior, tipo)
    % Verifica o tipo de função de Bessel (J0 ou Y0)
    if strcmp(tipo, 'J0')
        if k_r == 0
            integral_result = limite_superior^2 / 2;
        else
            integral_result = besselj(1, k_r * limite_superior) * limite_superior / k_r;
        end
    elseif strcmp(tipo, 'Y0')
        if k_r == 0
            error('k_r não pode ser zero para Y0.');
        else
            integral_result = (limite_superior * bessely(1, k_r * limite_superior) / k_r) + 2 / (pi * k_r^2);
        end
    else
        error('Tipo de função inválido. Use "J0" ou "Y0".');
    end
end



