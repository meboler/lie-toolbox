classdef SE2 < lie.LieGroupBase
    %SE2 The space of 2D rigid transformation matrices
    %   SE(2) is the space of 2D rigid motions T:
    %       T = [ R , t ]
    %           [ 0 , 1 ]
    %       where
    %           R \in SO(2)
    %           t \in R^2
    %   The Lie algebra of SE(2), se(2), is the space of 3x3 matrices of the
    %   form:
    %       \tau^ = [ [\theta]_x , \rho ]
    %               [     0      ,   0  ]
    %   The tangent vectors of SE(2) are of the form:
    %       \tau = [  \rho  ]
    %              [ \theta ]
    %       where
    %           \rho \in R^2
    %           \theta \in R^1
    
    properties (Constant)
        dof = 3
        size = [3, 3]
    end

    methods (Static)
        
        function xi = hat(x)
            arguments (Input)
                x (3, 1) double {mustBeReal, mustBeFinite}
            end
            arguments (Output)
                xi (3, 3) double {mustBeReal, mustBeFinite}
            end
            import lie.SO2;
            rho = x(1:2);
            theta = x(3);
            xi = [SO2.hat(theta), rho; ...
                      zeros(1,2),   0];
        end

        function x = vee(xi)
            arguments (Input)
                xi (3, 3) double {mustBeReal, mustBeFinite}
            end
            arguments (Output)
                x (3, 1) double {mustBeReal, mustBeFinite}
            end
            import lie.SO2;
            rho = xi(1:2, 3);
            theta = SO2.vee(xi(1:2, 1:2));
            x = [rho; theta];
        end

        function X = exp(xi)
            arguments (Input)
                xi (3, 3) double {mustBeReal, mustBeFinite}
            end
            arguments (Output)
                X (3, 3) double {mustBeReal, mustBeFinite}
            end
            import lie.SO2;
            import lie.SE2;
            x = SE2.vee(xi);
            rho = x(1:2);
            theta = x(3);
            R = SO2.bigexp(theta);
            V = vmatrix(theta);
            t = V * rho;
            X = [         R, t; ...
                 zeros(1,2), 1];
        end

        function xi = log(X)
            arguments (Input)
                X (3, 3) double {mustBeReal, mustBeFinite}
            end
            arguments (Output)
                xi (3, 3) double {mustBeReal, mustBeFinite}
            end
            import lie.SO2;
            R = X(1:2, 1:2);
            t = X(1:2, 3);
            theta_x = SO2.log(R);
            V = vmatrix(SO2.vee(theta_x));
            rho = V \ t;
            xi = [   theta_x, rho; ...
                  zeros(1,2),   0];
        end

        function X = bigexp(x)
            arguments (Input)
                x (3, 1) double {mustBeReal, mustBeFinite}
            end
            arguments (Output)
                X (3, 3) double {mustBeReal, mustBeFinite}
            end
            import lie.SO2;
            rho = x(1:2);
            theta = x(3);
            R = SO2.bigexp(theta);
            V = vmatrix(theta);
            t = V * rho;
            X = [         R, t; ...
                 zeros(1,2), 1];
        end

        function x = biglog(X)
            arguments (Input)
                X (3, 3) double {mustBeReal, mustBeFinite}
            end
            arguments (Output)
                x (3, 1) double {mustBeReal, mustBeFinite}
            end
            import lie.SO2;
            R = X(1:2, 1:2);
            t = X(1:2, 3);
            theta = SO2.biglog(R);
            V = vmatrix(theta);
            rho = V \ t;
            x = [rho; theta];
        end

    end

    methods (Static)
        
        function I = identity()
            arguments (Input)
            end
            arguments (Output)
                I (3, 3) double {mustBeReal, mustBeFinite}
            end
            I = eye(3);
        end

        function Xinv = inverse(X)
            arguments (Input)
                X (3, 3) double {mustBeReal, mustBeFinite}
            end
            arguments (Output)
                Xinv (3, 3) double {mustBeReal, mustBeFinite}
            end
            R = X(1:2, 1:2);
            t = X(1:2, 3);
            Xinv = [       R.', -R.' * t; ...
                    zeros(1,2),       1];
        end

    end

end

function V = vmatrix(theta)
    arguments (Input)
        theta (1, 1) double {mustBeReal, mustBeFinite}
    end
    arguments (Output)
        V (2, 2) double {mustBeReal, mustBeFinite}
    end
    import lie.SO2;
    if norm(theta) < 1e-10
        A = 1 - (theta^2 / 6) + (theta^4 / 120);
        B = (theta/2) - (theta^3 / 24) + (theta^5 / 720);
    else
        A = sin(theta) / theta;
        B = (1 - cos(theta)) / theta;
    end
    V = A * eye(2) + B * SO2.hat(1);
end

