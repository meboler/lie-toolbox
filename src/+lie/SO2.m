classdef SO2 < lie.LieGroupBase
    %SO2 The space of 2D rotation matrices
    %   SO(2) is the space of 2D rotation matrices R, which obey:
    %       R^T * R = I
    %       det(R) = 1
    %   The Lie algebra of SO(2), so(2), is the space of 2x2 skew-symmetrix
    %   matrices [s]_x:
    %       [s]_x = [  0 , -s]
    %               [  s ,  0]
    %       where
    %           s \in R^1
    %   The tangent vectors of SO(2) are of the form:
    %       \tau = \theta
    %       where
    %           \theta \in R^1

    properties (Constant)
        dof = 1
        size = [2, 2]
    end
    
    methods (Static)
        
        function xi = hat(x)
            %HAT Map a vector in R^n to an element of the Lie algebra
            %   $ hat(x) : R^n -> m $
            arguments (Input)
                x (1, 1) double {mustBeReal, mustBeFinite}
            end
            arguments (Output)
                xi (2, 2) double {mustBeReal, mustBeFinite}
            end
            xi = [ 0, -x;
                   x,  0];
        end

        function x  = vee(xi)
            %VEE Map an element of the Lie algebra to a vector in R^n
            %   $ vee(x) : m -> R^n $
            arguments (Input)
                xi (2, 2) double {mustBeReal, mustBeFinite}
            end
            arguments (Output)
                x (1, 1) double {mustBeReal, mustBeFinite}
            end
            x = xi(2,1);
        end

        function X = exp(xi)
            %EXP Map an element of the Lie algebra to an element of the Lie group
            %   $ exp(x) : m -> M $
            arguments (Input)
                xi (2, 2) double {mustBeReal, mustBeFinite}
            end
            arguments (Output)
                X (2, 2) double {mustBeReal, mustBeFinite}
            end
            import lie.SO2;
            angle = SO2.vee(xi);
            X = [ cos(angle), -sin(angle); ...
                  sin(angle),  cos(angle)];
        end
        
        function xi = log(X)
            %LOG Map an element of the Lie group to an element of the Lie algebra
            %   $ log(x) : M -> m $
            arguments (Input)
                X (2, 2) double {mustBeReal, mustBeFinite}
            end
            arguments (Output)
                xi (2, 2) double {mustBeReal, mustBeFinite}
            end
            import lie.SO2;
            angle = atan2(X(2,1), X(1,1));
            xi = SO2.hat(angle);
        end

        function X = bigexp(x)
            %BIGEXP Map a vector in R^n to an element of the Lie group
            %   $ Exp(x) : R^n -> M $
            arguments (Input)
                x (1, 1) double {mustBeReal, mustBeFinite}
            end
            arguments (Output)
                X (2, 2) double {mustBeReal, mustBeFinite}
            end
            X = [ cos(x), -sin(x); ...
                  sin(x),  cos(x)];
        end
        
        function x = biglog(X)
            %BIGLOG Map an element of the Lie group to a vector in R^n
            %   $ Log(x) : M -> R^n $
            arguments (Input)
                X (2, 2) double {mustBeReal, mustBeFinite}
            end
            arguments (Output)
                x (1, 1) double {mustBeReal, mustBeFinite}
            end
            x = atan2(X(2,1), X(1,1));
        end

    end 

    methods (Static)
        
        function I = identity()
            arguments (Input)
            end
            arguments (Output)
                I (2, 2) double {mustBeReal, mustBeFinite}
            end
            I = eye(2);
        end

        function Xinv = inverse(X)
            arguments (Input)
                X (2, 2) double {mustBeReal, mustBeFinite}
            end
            arguments (Output)
                Xinv (2, 2) double {mustBeReal, mustBeFinite}
            end
            Xinv = X.';
        end

    end
    
end