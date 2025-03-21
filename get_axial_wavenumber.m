function Kx = get_axial_wavenumber(Kr,k0,freq,type,region,N)

            
% Ajusta Kr para ser uma matriz coluna
%     Kr = Kr(:);
    
    % Inicializa a matriz Kx
    Kx = zeros(length(freq),size(Kr,2));
    
    % Calcula Kx para cada frequência
    for i = 1:length(freq)
        for m = 1:N
                   % Calcula Kx dependendo da relação entre k0 e Kr
            if k0(i) >= Kr(m)
                Kx(i,m) = sqrt(k0(i)^2 - Kr(m)^2); % Número de onda real
            else
                Kx(i,m) = -1i * sqrt(Kr(m)^2 - k0(i)^2); % Número de onda imaginário
            end

        end
    end


    if type == 1
        % Define o gráfico para a parte real e imaginária do número de onda axial
    figure('Position', [100, 100, 1200, 600]);
    
    % Loop para plotar as partes real e imaginária para cada modo

            % Plota a parte real
            plot(freq, real(Kx),'LineWidth', 1.5);
            hold on;
            
            % Plota a parte imaginária com linha tracejada
            plot(freq, imag(Kx), '--','LineWidth', 1.5);

        
        % Adiciona rótulos e título
        xlabel('Frequência (Hz)', 'FontSize', 14);
        ylabel('Número de Onda Axial', 'FontSize', 14);
        title(sprintf('Número de Onda Axial : %s', num2str(region)), 'FontSize', 16);
        
        % Exibe a legenda e ativa a grade
        legend('show');
        grid on;
        
        % Mantém o gráfico aberto
        hold off;

    end

    


end