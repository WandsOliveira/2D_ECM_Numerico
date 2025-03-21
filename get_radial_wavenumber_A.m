function Kr_A = get_radial_wavenumber_A(r1, N,type)

    % Define a função anônima que representa J1(r1*x)
    f = @(x) besselj(1, r1*x);

    % Passo pequeno para derivada numérica
    h = 1e-6;

    % Inicializa vetor para armazenar as raízes
    Kr_A = zeros(1, N);
    x0 = 0; % Chute inicial para a primeira raiz
%        x0 = pi/(r1);
    for i = 1:N
        x = x0; % Define o ponto inicial

        % Método de Newton-Raphson
        while true
            % Derivada numérica central: f'(x) ≈ (f(x + h) - f(x - h)) / (2h)
            df = (f(x + h) - f(x - h)) / (2 * h);

            % Atualização do método de Newton-Raphson
            x_novo = x - f(x) / df;

            % Critério de convergência
            if abs(x_novo - x) < 1e-6
                break;
            end

            x = x_novo;
        end

        Kr_A(i) = x; % Armazena a raiz encontrada

        % Atualiza o chute inicial para a próxima raiz
        x0 = x + pi / r1;
    end
    
    
    % Plot da função e das raízes
    x = linspace(0, Kr_A(end) + 20, 1000); % Intervalo para o plot
    y = f(x); % Valores da função
    if type == 1
        figure;
        plot(x, y, 'b', 'LineWidth', 1.5); % Plot da função
        hold on;
        scatter(Kr_A, f(Kr_A), 'ro', 'filled'); % Raízes como círculos vermelhos
        title('Function J1(r1*x) and its Roots');
        xlabel('x');
        ylabel('J1(r1*x)');
        grid on;
        legend('J1(r1*x)', 'Roots');
        yline(0, '--k'); % Linha horizontal em y=0
        hold off;
    end
end
