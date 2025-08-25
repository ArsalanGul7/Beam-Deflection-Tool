%% Beam Deflection Calculator
% Author:    Syed Arsalan Gul
% Date:      19/08/2025
% About:     This script calculates deflection, shear force, and bending moment
%            diagrams for beams with different boundary conditions and load cases.

%% ------------------ INPUT SECTION ------------------
Beam.L     = 5;            % Beam Length [m]
Beam.E     = 200e9;        % Young's Modulus [Pa]
Beam.I     = 8e-6;         % Second Moment of Area [m^4]

Beam.support  = 'simply';  % Options: 'cantilever', 'simply', 'fixedfixed'
Beam.loadType = 'point';   % Options: 'point', 'udl'

% Point Load Parameters
Beam.P = 10e3;             % Load Magnitude [N]
Beam.a = 2.5;              % Distance from left end [m]

% UDL Parameters
Beam.w = 5e3;              % UDL intensity [N/m]

%% ------------ COMPUTATIONAL ANALYSIS --------------
n = 200;                                % Discretization pointsMaximum Deflection	1.6276
x = linspace(0, Beam.L, n);

% Initialize
y = zeros(size(x));
M = zeros(size(x));
V = zeros(size(x));

switch lower(Beam.loadType)
    case 'point'
        switch lower(Beam.support)
            case 'cantilever'
                % Deflection due to point load at a distance 'a'
                y = -(Beam.P * x.^2 / (6*Beam.E*Beam.I)) .* (3*Beam.a - x) .* (x<=Beam.a) ...
                    - (Beam.P * Beam.a^2 / (6*Beam.E*Beam.I)) .* (3*x - Beam.a) .* (x>Beam.a);
                M = -(Beam.P * (Beam.a - x)) .* (x<=Beam.a);
                V = -Beam.P * (x<=Beam.a);

            case 'simply'
                % Reaction forces
                R1 = Beam.P * (Beam.L - Beam.a) / Beam.L;
                R2 = Beam.P * Beam.a / Beam.L;

                % Moment and shear
                M = R1 .* x - Beam.P .* max(0, x - Beam.a);
                V = R1 - Beam.P .* (x >= Beam.a);

                % Deflection (approx)
                y = (R1 .* x .* (Beam.L.^2 - x.^2) / (6*Beam.E*Beam.I)) ...
                    - (Beam.P .* max(0, x-Beam.a).^3 / (6*Beam.E*Beam.I));

            case 'fixedfixed'
                if abs(Beam.a - Beam.L/2) > 1e-6
                    error('Fixed–fixed + off-center point load requires FEM (not included here).');
                end
                y = -(Beam.P * x .* (Beam.L^3 - 2*Beam.L*x.^2 + x.^3)) / (48*Beam.E*Beam.I);
                M = Beam.P*Beam.L/8 * (1 - 2*x/Beam.L);
                V = [repmat(Beam.P/2,1,ceil(n/2)) repmat(-Beam.P/2,1,floor(n/2))];
        end

    case 'udl'
        switch lower(Beam.support)
            case 'cantilever'
                y = -(Beam.w * x.^2 / (24*Beam.E*Beam.I)) .* (6*Beam.L^2 - 4*Beam.L*x + x.^2);
                M = -Beam.w/2 .* (Beam.L - x).^2;
                V = -Beam.w * (Beam.L - x);

            case 'simply'
                R = Beam.w*Beam.L/2;
                M = R.*x - (Beam.w.*x.^2)/2;
                V = R - Beam.w.*x;
                y = (Beam.w .* x .* (Beam.L^3 - 2*Beam.L*x.^2 + x.^3)) / (24*Beam.E*Beam.I);

            case 'fixedfixed'
                y = -(Beam.w * x .* (Beam.L^3 - 2*Beam.L*x.^2 + x.^3)) / (24*Beam.E*Beam.I);
                M = Beam.w*Beam.L^2/12 * (1 - 6*(x/Beam.L).*(1-x/Beam.L));
                V = gradient(M, x(2)-x(1));
        end
end

%% ---------------- OUTPUT SECTION ------------------
fprintf('>> Beam Deflection Calculator Results:\n');
fprintf('--------------------------------------\n');
fprintf('Beam Length             : %.2f m\n', Beam.L);
fprintf('Support Type            : %s\n', Beam.support);
fprintf('Load Type               : %s\n\n', Beam.loadType);

fprintf('--- MAX RESULTS ---\n');
fprintf('Maximum Deflection      : %.3f mm\n', min(y)*1e3);
fprintf('Maximum Moment          : %.2f kN·m\n', max(abs(M))/1e3);
fprintf('Maximum Shear           : %.2f kN\n', max(abs(V))/1e3);
fprintf('--------------------------------------\n\n');

%% ---------------- PLOTTING SECTION ------------------
figure;
plot(x, y*1e3, 'LineWidth', 2);
xlabel('x [m]'); ylabel('Deflection [mm]');
title('Deflection Curve'); grid on;

figure;
plot(x, M/1e3, 'LineWidth', 2);
xlabel('x [m]'); ylabel('Moment [kN·m]');
title('Bending Moment Diagram'); grid on;

figure;
plot(x, V/1e3, 'LineWidth', 2);
xlabel('x [m]'); ylabel('Shear Force [kN]');
title('Shear Force Diagram'); grid on;

%% ---------------- SAVE RESULTS ------------------
save('Beam_Deflection_Calculator', 'Beam', 'x', 'y', 'M', 'V');