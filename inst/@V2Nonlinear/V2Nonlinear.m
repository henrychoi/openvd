classdef V2Nonlinear < VehicleSimple
    % V2Nonlinear Nonlinear simple vehicle model.
    %
    % It inherits properties from VehicleSimple.

    methods
        % Constructor
        function self = V2Nonlinear()
            self.mF0 = 700;
            self.mR0 = 600;
            self.IT = 10000;
            self.lT = 3.5;
            self.nF = 2;
            self.nR = 2;
            self.wT = 2;
            self.muy = .8;
            self.deltaf = 0;
            self.Fxf = 0;
            self.Fxr = 0;
        end

        %% Model
        % Função com as equações de estado do modelo
        function dx = Model(self, t, states,tspan)
            % Data
            m = self.mT;
            I = self.IT;
            a = self.a;
            b = self.b;
            nF = self.nF;
            nR = self.nR;
            muy = self.muy;


            g = 9.81;                 % Gravity [m/s^2]

            FzF = self.mF0 * g;       % Vertical load @ F [N]
            FzR = self.mR0 * g;       % Vertical load @ R [N]

            X = states(1);
            Y = states(2);
            PSI = states(3);
            Vd = states(4);
            Vq = states(5);
            dPSI = states(6);

            if isa(self.deltaf,'function_handle') % steering algorithm
                DELTA = self.deltaf(states, t);
            elseif length(self.deltaf)>1
                DELTA = interp1(tspan,self.deltaf,t);
            else
                DELTA = self.deltaf;
            end

            % wheel slip angles
            ALPHAF = atan2(Vq + a * dPSI, Vd) - DELTA;
            ALPHAR = atan2(Vq - b * dPSI, Vd);

            if isa(self.Fxf,'function_handle') % front wheel drive
                FxF = self.Fxf(states, t);
            elseif length(self.Fxf)>1
                FxF = interp1(tspan,self.Fxf,t);
            else
                FxF = self.Fxf;
            end

            if isa(self.Fxr,'function_handle') % rear wheel drive
                FxR = self.Fxr(states, t);
            elseif length(self.Fxr)>1
                FxR = interp1(tspan,self.Fxr,t);
            else
                FxR = self.Fxr;
            end

            % Lateral forces on EACH wheel
            FyF = nF * self.tire.Characteristic(ALPHAF, FzF/nF, muy);
            FyR = nR * self.tire.Characteristic(ALPHAR, FzR/nR, muy);

            % Equations of motion
            dx(1,1) = Vd * cos(PSI) - Vq * sin(PSI); % dx
            dx(2,1) = Vd * sin(PSI) + Vq * cos(PSI); % dy
            dx(3,1) = dPSI; % w
            % generalized forces Qd, Qq, Q\psi
            dx(4,1) = (FxF * cos(DELTA) - FyF * sin(DELTA) + FxR)/m;
            dx(5,1) = (FxF * sin(DELTA) + FyF * cos(DELTA) + FyR)/m;
            dx(6,1) = (FxF * a * sin(DELTA) + FyF * a * cos(DELTA) - FyR * b) / I;
        end
    end
end

%% See Also
%
% <../../index.html Home>
%
