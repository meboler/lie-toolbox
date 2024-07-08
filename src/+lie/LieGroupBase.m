classdef (Abstract) LieGroupBase
    %LIEGROUPBASE Abstract base class for Lie groups

    properties (Abstract, Constant)
        %DOF The number of degrees of freedom of the group
        dof (1, 1) {mustBePositive, mustBeInteger}

        %SIZE The size of a group element
        size (1, 2) {mustBePositive, mustBeInteger}
    end

    methods (Abstract, Static)
        
        %HAT Map a vector in R^n to an element of the Lie algebra
        %   $ hat(x) : R^n -> m $
        xi = hat(x);

        %VEE Map an element of the Lie algebra to a vector in R^n
        %   $ vee(x) : m -> R^n $
        x = vee(xi);
        
        %EXP Map an element of the Lie algebra to an element of the Lie group
        %   $ exp(x) : m -> M $
        X = exp(xi);

        %LOG Map an element of the Lie group to an element of the Lie algebra
        %   $ log(x) : M -> m $
        xi = log(X);
        
        %BIGEXP Map a vector in R^n to an element of the Lie group
        %   $ Exp(x) : R^n -> M $
        X = bigexp(x);

        %BIGLOG Map an element of the Lie group to a vector in R^n
        %   $ Log(x) : M -> R^n $
        x = biglog(x);

    end

    methods (Abstract, Static)
        %IDENTITY Generate the identity element of the Lie group
        %   $ identity() : -> M $
        I = identity();

        %INVERSE Generate the inverse of an element of the Lie group
        %   $ inverse(X) : M -> M $
        Xinv = inverse(X);
    end

end

