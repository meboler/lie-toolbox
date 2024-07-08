classdef SO3 < lie.LieGroupBase
    %SO3 The group of 3D rotation matrices
    %   SO(3) is the space of 3D rotation matrices R, which obey:
    %       * R^T * R = I
    %       * det(R) = 1
    %   The Lie algebra of SO(3), so(3), is the space of 3x3 skew-symmetric
    %   matrices [w]_x:
    %       [w]_x = [   0 , -wz ,  wy ]
    %               [  wz ,   0 , -wx ]
    %               [ -wy ,  wx ,   0 ]
    %       where
    %           w = [wx, wy, wz]^T \in R^3
    %   The tangent vectors of SO(3) are of the form:
    %       \tau = \theta
    %       where
    %           \theta \in R^3

    properties (Constant)
        dof = 3
        size = [3, 3]
    end

    methods (Static)
        
        function xi = hat(x)
            arguments (Input)
                x (3, 1) double
            end
            arguments (Output)
                xi (3, 3) double
            end
            xi = [    0, -x(3),  x(2); ...
                   x(3),     0, -x(1); ...
                  -x(2),  x(1),     0];
        end

        function x = vee(xi)
            arguments (Input)
                xi (3, 3) double
            end
            arguments (Output)
                x (3, 1) double
            end
            x = [xi(3,2); xi(1,3); xi(2,1)];
        end

        function X = exp(xi)
            arguments (Input)
                xi (3, 3) double
            end
            arguments (Output)
                X (3, 3) double
            end
            import lie.SO3;
            angle = norm(SO3.vee(xi));
            if angle < 1e-10
                % First-order Taylor approximation
                A = 1 - (angle^2 / factorial(3));
                B = 0.5 - (angle^2 / factorial(4));
            else
                % Closed-form Rodrigues solution
                A = sin(angle) / angle;
                B = (1 - cos(angle)) / (angle^2);
            end
            X = eye(3) + (A * xi) + (B * xi * xi);
        end

        function xi = log(X)
            arguments (Input)
                X (3, 3) double
            end
            arguments (Output)
                xi (3, 3) double
            end
            angle = acos((trace(X) - 1) / 2);
            if angle < 1e-8
                % First-order Taylor approximation
                A = 0.5 + (angle^2 / 12);
            else
                % Closed-form solution
                A = angle / (2 * sin(angle));
            end
            xi = A * (X - X.');
        end

        function X = bigexp(x)
            arguments (Input)
                x (3, 1) double
            end
            arguments (Output)
                X (3, 3) double
            end
            import lie.SO3;
            X = SO3.exp(SO3.hat(x));
        end

        function x = biglog(X)
            arguments (Input)
                X (3, 3) double
            end
            arguments (Output)
                x (3, 1) double
            end
            import lie.SO3;
            x = SO3.vee(SO3.log(X));
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
            Xinv = X.';
        end

    end

end

