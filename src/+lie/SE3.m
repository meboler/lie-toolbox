classdef SE3 < lie.LieGroupBase
    %SE2 The space of 3D rigid transformation matrices
    %   SE(3) is the space of 3D rigid motions T:
    %       T = [ R , t ]
    %           [ 0 , 1 ]
    %       where
    %           R \in SO(3)
    %           t \in R^3
    %   The Lie algebra of SE(3), se(3), is the space of 3x3 matrices of the
    %   form:
    %       \tau^ = [ [\theta]_x , \rho ]
    %               [     0      ,   0  ]
    %   The tangent vectors of SE(3) are of the form:
    %       \tau = [  \rho  ]
    %              [ \theta ]
    %       where
    %           \rho \in R^3
    %           \theta \in R^3
    
    properties (Constant)
        dof = 6
        size = [4, 4]
    end

    methods (Static)
        
        function xi = hat(x)
            arguments (Input)
                x (6, 1) double {mustBeReal, mustBeFinite}
            end
            arguments (Output)
                xi (4, 4) double {mustBeReal, mustBeFinite}
            end
            import lie.SO3;
            rho = x(1:3);
            theta = x(4:6);
            xi = [SO3.hat(theta), rho; ...
                      zeros(1,3),   0];
        end

        function x = vee(xi)
            arguments (Input)
                xi (4, 4) double {mustBeReal, mustBeFinite}
            end
            arguments (Output)
                x (6, 1) double {mustBeReal, mustBeFinite}
            end
            import lie.SO3;
            rho = xi(1:3, 4);
            theta = SO3.vee(xi(1:3, 1:3));
            x = [rho; theta];
        end

        function X = exp(xi)
            arguments (Input)
                xi (4, 4) double {mustBeReal, mustBeFinite}
            end
            arguments (Output)
                X (4, 4) double {mustBeReal, mustBeFinite}
            end
            import lie.SO3;
            import lie.SE3;
            x = SE3.vee(xi);
            rho = x(1:3);
            theta = x(4:6);
            R = SO3.bigexp(theta);
            V = vmatrix(theta);
            t = V * rho;
            X = [         R, t; ...
                 zeros(1,3), 1];
        end

        function xi = log(X)
            arguments (Input)
                X (4, 4) double {mustBeReal, mustBeFinite}
            end
            arguments (Output)
                xi (4, 4) double {mustBeReal, mustBeFinite}
            end
            import lie.SO3;
            R = X(1:3, 1:3);
            t = X(1:3, 4);
            theta_x = SO3.log(R);
            V = vmatrix(SO3.vee(theta_x));
            rho = V \ t;
            xi = [   theta_x, rho; ...
                  zeros(1,3),   0];
        end

        function X = bigexp(x)
            arguments (Input)
                x (6, 1) double {mustBeReal, mustBeFinite}
            end
            arguments (Output)
                X (4, 4) double {mustBeReal, mustBeFinite}
            end
            import lie.SO3;
            rho = x(1:3);
            theta = x(4:6);
            R = SO3.bigexp(theta);
            V = vmatrix(theta);
            t = V * rho;
            X = [         R, t; ...
                 zeros(1,3), 1];
        end

        function x = biglog(X)
            arguments (Input)
                X (4, 4) double {mustBeReal, mustBeFinite}
            end
            arguments (Output)
                x (6, 1) double {mustBeReal, mustBeFinite}
            end
            import lie.SO3;
            R = X(1:3, 1:3);
            t = X(1:3, 4);
            theta = SO3.biglog(R);
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
                I (4, 4) double {mustBeReal, mustBeFinite}
            end
            I = eye(4);
        end

        function Xinv = inverse(X)
            arguments (Input)
                X (4, 4) double {mustBeReal, mustBeFinite}
            end
            arguments (Output)
                Xinv (4, 4) double {mustBeReal, mustBeFinite}
            end
            R = X(1:3, 1:3);
            t = X(1:3, 4);
            Xinv = [       R.', -R.' * t; ...
                    zeros(1,3),       1];
        end

    end

end

function V = vmatrix(theta)
    arguments (Input)
        theta (3, 1) double {mustBeReal, mustBeFinite}
    end
    arguments (Output)
        V (3, 3) double {mustBeReal, mustBeFinite}
    end
    import lie.SO3;
    angle = norm(theta);
    if angle < 1e-10
        A = (1/2) - (angle^2 / 24) + (angle^4 / 720);
        B = (1/6) - (angle^2 / 120) + (angle^4 / 5040);
    else
        A = (1 - cos(angle)) / angle^2;
        B = (angle - sin(angle)) / angle^3;
    end
    theta_x = SO3.hat(theta);
    theta_x_2 = theta_x * theta_x;
    V = eye(3) + A * theta_x + B * theta_x_2;
end

