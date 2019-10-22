  function [A,delta] = harmonicAddition(Avals,deltavals)
            % Calculates Harmonic Addition Theorem:
            % Acos(x-delta) = sum_i Avals(i)cos(x-deltavals(i))
            % Params
            %   Avals     : vector of A_i
            %   deltavals : vector of delta_i (same length as Avals)
            % Returns
            %   A         : scalar A value
            %   delta     : scalar delta value

            bx = sum(Avals.*cos(deltavals));
            by = sum(Avals.*sin(deltavals));
            A = sqrt(bx^2+by^2);
            delta = atan2(by,bx);
          
        end
