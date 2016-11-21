classdef VehicleArticulatedLinear < VehicleDynamicsLateral.VehicleArticulated
    % VehicleArticulatedLinear Linear articulated vehicle model.
    %
    % It inherits properties from VehicleArticulated.

    methods
        % Constructor
        function self = VehicleArticulatedLinear()
            self.mF0 = 5200;
            self.mR0 = 2400;
            self.mF = 6000;
            self.mR = 10000;
            self.mM = 17000;
            self.IT = 46000;
            self.IS = 450000;
            self.lT = 3.5;
            self.lS = 7.7;
            self.c = -0.3;
            self.nF = 2;
            self.nR = 4;
            self.nM = 8;
            self.wT = 2.6;
            self.wS = 2.4;
            self.muy = 0.3;
            self.deltaf = 0;
            self.Fxf = 0;
            self.Fxr = 0;
            self.Fxm = 0;
        end

        function dx = Model(self,t,estados,tspan)
            % Vehicle parameters
            mT = self.mT;
            mS = self.mS;
            IT = self.IT;
            IS = self.IS;
            a = self.a;
            b = self.b;
            c = self.c;
            d = self.d;
            e = self.e;
            nF = self.nF;
            nR = self.nR;
            nM = self.nM;
            g = 9.81;

            % Steering angle
            deltaf = interp1(tspan,self.deltaf,t);

            % Vertical forces
            FzF = self.mF * g;
            FzR = self.mR * g;
            FzM = self.mM * g;
            muy = self.muy;

            v0 = 20;

            % State
            X = estados(1,1);
            Y = estados(2,1);
            PSI     = estados(3,1);
            PHI     = estados(4,1);
            V       = estados(5,1);
            ALPHAT  = estados(6,1);
            dPSI    = estados(7,1);
            dPHI    = estados(8,1);

            % Slip angles - linear
            ALPHAF = ALPHAT + a/v0*dPSI - deltaf;
            ALPHAR = ALPHAT - b/v0*dPSI;
            ALPHAM = ALPHAT + PHI - (dPSI*(b + c + d + e))/v0 + (dPHI*(d + e))/v0;

            % Longitudinal forces
            if isa(self.Fxf,'function_handle')
                FxF = self.Fxf([X;Y;PSI;PHI;V;ALPHAT;dPSI;dPHI],t);
            elseif length(self.Fxf)>1
                FxF = interp1(tspan,self.Fxf,t);
            else
                FxF = self.Fxf;
            end

            if isa(self.Fxr,'function_handle')
                FxR = self.Fxr([X;Y;PSI;PHI;V;ALPHAT;dPSI;dPHI],t);
            elseif length(self.Fxr)>1
                FxR = interp1(tspan,self.Fxr,t);
            else
                FxR = self.Fxr;
            end

            if isa(self.Fxm,'function_handle')
                FxM = self.Fxm([X;Y;PSI;PHI;V;ALPHAT;dPSI;dPHI],t);
            elseif length(self.Fxm)>1
                FxM = interp1(tspan,self.Fxm,t);
            else
                FxM = self.Fxm;
            end


            % Lateral forces - Characteristic curve
            FyF = nF*self.tire.Characteristic(ALPHAF,FzF/nF,muy);
            FyR = nR*self.tire.Characteristic(ALPHAR,FzR/nR,muy);
            FyM = nM*self.tire.Characteristic(ALPHAM,FzM/nM,muy);

            A = [ 0  0   0  0  1   0                       0  0;...
                  0  0  v0  0  0  v0                       0  0;...
                  0  0   0  0  0   0                       1  0;...
                  0  0   0  0  0   0                       0  1;...
                  0  0   0  0  0   0                       0  0;...
                  0  0   0  0  0   0           (-v0*(mS + mT))  0;...
                  0  0   0  0  0   0  (mS*(d*v0 + v0*(b + c)))  0;...
                  0  0   0  0  0   0                (-d*mS*v0)  0];

            B = [ 0  0  0  0  0   0                0;...
                  0  0  0  0  0   0                0;...
                  0  0  0  0  0   0                0;...
                  0  0  0  0  0   0                0;...
                  0  1  1  1  0   0                0;...
                  0  0  0  0  1   1                1;...
                  0  0  0  0  a  -b  (- b - c - d - e);...
                  0  0  0  0  0   0            (d + e)];

            vetEst = [X ; Y ; PSI ; PHI ; V ; ALPHAT ; dPSI ; dPHI];
            vetEnt = [deltaf ; FxF ; FxR ; FxM ; FyF ; FyR ; FyM];

            % Integrator output
            dx = A*vetEst + B*vetEnt;
        end

        function [value,isterminal,direction] = velocity(~,~,estados)
            % If the velocity is less than 0.1m/s the integrator stops.
            % The MassMatrix is singular when the velocity is 0 m/s.
            value = estados(5,1) - 0.1;
            isterminal = 1;
            direction = -1;
        end

        function E = MassMatrix(self,~,~)
            % Vehicle parameters
            mT = self.mT;
            mS = self.mS;
            IT = self.IT;
            IS = self.IS;
            a = self.a;
            b = self.b;
            c = self.c;
            d = self.d;
            e = self.e;
            deltaf = self.deltaf;
            nF = self.nF;
            nR = self.nR;
            nM = self.nM;

            g = 9.81;

            v0 = 20;

            % Mass matrix
            E = [ 1  0  0  0        0                        0                                               0                            0;...
                  0  1  0  0        0                        0                                               0                            0;...
                  0  0  1  0        0                        0                                               0                            0;...
                  0  0  0  1        0                        0                                               0                            0;...
                  0  0  0  0  (mS + mT)                      0                                               0                            0;...
                  0  0  0  0        0             (v0*(mS + mT))                                 (-mS*(b + c + d))                         (d*mS);...
                  0  0  0  0        0  ( -mS*v0*(b + c + d))            (IS + IT + mS*(b + c + d)^2)  (- IS - mS*(d^2 + (b + c)*d));...
                  0  0  0  0        0                  (d*mS*v0)                     (- IS - mS*(d^2 + (b + c)*d))                  (mS*d^2 + IS)];

        end
    end
end

%% See Also
%
% <../../../index.html index>
%
