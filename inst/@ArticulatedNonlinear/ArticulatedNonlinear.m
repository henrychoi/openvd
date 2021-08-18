classdef ArticulatedNonlinear < VehicleArticulated
    % Nonlinear articulated vehicle model using the d-q frame
    %
    methods
        % Constructor
        function self = ArticulatedNonlinear()
            self.mF0 = 5200; % Mass over the front axle [kg]
            self.mR0 = 2400;  % Mass over the rear axle [kg]
            self.mF = 6000; % Mass over the front axle when the semitrailer is coupled 
            self.mR = 10000; % mass over rear wheel
            self.mM = 17000; % mass over trailer axle
            self.IT = 46000; % tracotr Izz
            self.IS = 450000; % trailer Izz
            self.lT = 3.5; % tractor wheelbase
            self.lS = 7.7; % A-M
            self.c = -0.3; % 5th wheel distance from rear axle
            self.nF = 2; % N front
            self.nR = 4; % N rear
            self.nM = 8; % N trailer wheels
            self.wT = 2.6; % tractor track width
            self.wS = 2.4; % trailer track width
            self.muy = 0.3;
            self.deltaf = 0; % steering angle
            self.Fxf = 0; % axial forces
            self.Fxr = 0;
            self.Fxm = 0;
        end

        function dx = Model(self,t, states,tspan)
            % Vehicle parameters
            mT = self.mT;
            mS = self.mS;
            % IT = self.IT;
            % IS = self.IS;
            a = self.a;
            b = self.b;
            c = self.c;
            d = self.d;
            e = self.e;
            nF = self.nF;
            nR = self.nR;
            nM = self.nM;

            g = 9.81;

            FzF = self.mF * g;
            FzR = self.mR * g;
            FzM = self.mM * g;
            muy = self.muy;


            % States
            X = states(1,1);
            Y = states(2,1);
            PSI = states(3,1); % tractor heading
            PHI = states(4,1); % tractor heading - trailer heading 
            Vd = states(5,1);
            Vq = states(6,1);
            dPSI = states(7,1); % tractor yaw rate
            dPHI = states(8,1); % trailer yaw rate - trailer heading

            % Steering input
            if isa(self.deltaf,'function_handle')
                deltaf = self.deltaf(states(:,1), t);
            elseif length(self.deltaf)>1
                deltaf = interp1(tspan,self.deltaf,t);
            else
                deltaf = self.deltaf;
            end

            % Slip angles are determined from the Vq_i/Vd_i (and the steering angle)
            ALPHAF = atan2(Vq + a * dPSI, Vd) - deltaf;
            ALPHAR = atan2(Vq - b * dPSI, Vd);
            % \vec{Vm} = [∂Pc/∂d, ∂Pc/∂q, ∂Pc/∂\psi, ∂Pc/∂\phi] • [Vd; Vq; \dot\psi; \dot\phi]
            % ∂Pc/∂d • Vd = \hat{d} Vd 
            % ∂Pc/∂q • Vq = \hat{q} Vq
            % L = d + e
            % ∂Pc/∂𝜓 • w = {-L sin(ø) \hat{d} - (L cos(ø) + b + c) \hat{q} } w
            % ∂Pc/∂ø • \dot ø = { L sin(ø) \hat{d} +  L cos(ø)     \hat{q} } \dot ø
            % \vec{Vm} is then = {Vd - L sin(ø) w + L sin(ø) w} \hat{d}
            %                  + {Vq - (L cos(ø) + b + c) w + L cos(ø) \dot{ø}}} \hat{q}
            % from which we can calculate the angle of the \vec{Vm} = atan(Vm_q, Vm_d)
            ALPHAM = atan2(Vq - ((d+e)*cos(PHI) + b + c) * dPSI + (d+e)*cos(PHI) * dPHI, ...
                           Va - ((d+e)*sin(PHI)        ) * dPSI + (d+e)*sin(PHI) * dPHI);

            % front wheel drive/drive
            if isa(self.Fxf,'function_handle')
                FxF = self.Fxf(states(:,1), t);
            elseif length(self.Fxf)>1
                FxF = interp1(tspan,self.Fxf,t);
            else
                FxF = self.Fxf;
            end

            % rear wheel drive
            if isa(self.Fxr,'function_handle')
                FxR = self.Fxr(states(:,1), t);
            elseif length(self.Fxr)>1
                FxR = interp1(tspan,self.Fxr,t);
            else
                FxR = self.Fxr;
            end

            % tractor drive/drag
            if isa(self.Fxm,'function_handle')
                FxM = self.Fxm(states(:,1), t);
            elseif length(self.Fxm)>1
                FxM = interp1(tspan,self.Fxm,t);
            else
                FxM = self.Fxm;
            end

            % Lateral forces on EACH wheel
            FyF = nF * self.tire.Characteristic(ALPHAF, FzF/nF, muy);
            FyR = nR * self.tire.Characteristic(ALPHAR, FzR/nR, muy);
            FyM = nM * self.tire.Characteristic(ALPHAM, FzM/nM, muy);

            % d(p)/dt = vel
            dx(1,1) = Vd * cos(PSI) - Vq * sin(PSI); % dx
            dx(2,1) = Vd * sin(PSI) + Vq * cos(PSI); % dy
            dx(3,1) = dPSI; % w
            dx(4,1) = dPHI; % dPSI

            % M • d(\vec V)/dt = \vec Q_{ext} - \vec Q_{cen}, where
            % M is a generalized inertia matrix
            % Q_{ext} is generalized external force, Q_{cen} is the generalized centripetal force
            Qext = [ FxF*cos(deltaf) - FyF*sin(deltaf) + FxR + FxM*cos(PHI) + FyM*sin(PHI); ...
                    FxF*sin(deltaf) + FyF*cos(deltaf) + FyR - FxM*sin(PHI) + FyM*cos(PHI); ...
                    a*FxF*sin(deltaf) + a*FyF*cos(deltaf) - b*FyR + (b+c)*FxM*sin(PHI) ...
                        - (d+e + (b+c)*cos(PHI)) * FyM; ...
                    (d+e) * FyM ];
            
            dx(5,1) =  + FxM*cos(PSI - PHI) - FyF*sin(PSI + deltaf) - FyR*sin(PSI) - FyM*sin(PSI - PHI) - mS*(b+c)*dPSI^2*cos(PSI) - mS*d*(dPSI - dPHI)^2*cos(PSI - PHI) + (mT + mS)*VT*sin(PSI+ALPHAT)*dPSI;

            FxF*sin(PSI + deltaf) + FxR*sin(PSI) + FxM*sin(PSI - PHI) + FyF*cos(PSI + deltaf) + FyR*cos(PSI) + FyM*cos(PSI - PHI) - mS*(b+c)*dPSI^2*sin(PSI) - mS*d*(dPSI - dPHI)^2*sin(PSI - PHI) - (mT + mS)*VT*cos(PSI+ALPHAT)*dPSI;...
            FxF*a*sin(deltaf) + FxM*(b + c)*sin(PHI) + FyF*a*cos(deltaf) - FyR*b - FyM*((b+c)*cos(PHI) + (d+e)) - mS*(b+c)*d*(dPSI - dPHI)^2*sin(PHI) + mS*(b+c)*d*dPSI^2*sin(PHI) + mS*((b+c)*VT*cos(ALPHAT) + d*VT*cos(ALPHAT + PHI))*dPSI;...
            FyM*(d + e) - mS*(b+c)*d*dPSI^2*sin(PHI) - mS*d*VT*cos(ALPHAT + PHI)*dPSI;
        end

        function [value,isterminal,direction] = velocity(~,~,states)
            % If the velocity is less than 0.1m/s the integrator stops.
            % The MassMatrix is singular when the velocity is 0 m/s.
            value = states(5,1) - 0.1;
            isterminal = 1;
            direction = -1;
        end

        function M = MassMatrix(self,~,states)
            % Vehicle Parameters
            mT = self.mT;
            mS = self.mS;
            IT = self.IT;
            IS = self.IS;
            % a = self.a;
            b = self.b;
            c = self.c;
            d = self.d;
            % e = self.e;
            % deltaf = self.deltaf;
            % nF = self.nF;
            % nR = self.nR;
            % nM = self.nM;

            % g = 9.81;

            % FzF = self.mF * g;
            % FzR = self.mR * g;
            % FzM = self.mM * g;
            % muy = self.muy;


            % States
            PSI = states(3,1);
            PHI = states(4,1);
            VT = states(5,1);
            ALPHAT = states(6,1);
            % dPSI = states(7,1);
            % dPHI = states(8,1);

            % Matriz de massa
            M55 = (mT + mS)*cos(PSI + ALPHAT);
            M56 = -(mT + mS)*VT*sin(PSI + ALPHAT);
            M57 = mS*( (b+c)*sin(PSI) + d*sin(PSI - PHI) );
            M58 = -mS*d*sin(PSI - PHI);
            M65 = (mT + mS)*sin(PSI + ALPHAT);
            M66 = (mT + mS)*VT*cos(PSI + ALPHAT);
            M67 = -mS*( (b+c)*cos(PSI) + d*cos(PSI - PHI) );
            M68 = mS*d*cos(PSI - PHI);
            M75 = -mS*( (b+c)*sin(ALPHAT) + d*sin(ALPHAT + PHI) );
            M76 = -mS*( (b+c)*VT*cos(ALPHAT) + d*VT*cos(ALPHAT + PHI) );
            M77 = mS*( (b+c)^2 + 2*(b+c)*d*cos(PHI) + d^2 ) + IT + IS;
            M78 = -( mS*( (b+c)*d*cos(PHI) + d^2 ) + IS);
            M85 = mS*d*sin(ALPHAT + PHI);
            M86 = mS*d*VT*cos(ALPHAT + PHI);
            M87 = - (mS*(d^2 + (b+c)*d*cos(PHI)) + IS);
            M88 = mS*d^2 + IS;

            M = [   1 0 0 0  0   0   0   0 ;...
                    0 1 0 0  0   0   0   0 ;...
                    0 0 1 0  0   0   0   0 ;...
                    0 0 0 1  0   0   0   0 ;...
                    0 0 0 0 M55 M56 M57 M58 ;...
                    0 0 0 0 M65 M66 M67 M68 ;...
                    0 0 0 0 M75 M76 M77 M78 ;...
                    0 0 0 0 M85 M86 M87 M88 ];
        end
    end
end

%% See Also
%
% <../../index.html Home>
%
