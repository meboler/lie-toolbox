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
        %   $ hat(x) : \mathbb{R}^n \to \mathfrak{m} $
        xi = hat(x);

        %VEE Map an element of the Lie algebra to a vector in R^n
        %   $ vee(x) : \mathfrak{m} \to \mathbb{R}^n $
        x = vee(xi);
        
        %EXP Map an element of the Lie algebra to an element of the Lie group
        %   $ exp(x) : \mathfrak{m} \to \mathcal{M} $
        X = exp(xi);

        %LOG Map an element of the Lie group to an element of the Lie algebra
        %   $ log(x) : \mathcal{M} \to \mathfrak{m} $
        xi = log(X);
        
        %BIGEXP Map a vector in R^n to an element of the Lie group
        %   $ Exp(x) : \mathbb{R}^n \to \mathcal{M} $
        X = bigexp(x);

        %BIGLOG Map an element of the Lie group to a vector in R^n
        %   $ Log(x) : \mathcal{M} \to \mathbb{R}^n $
        x = biglog(x);

    end

    methods (Abstract, Static)
        %IDENTITY Generate the identity element of the Lie group
        %   $ identity() : \to \mathcal{M} $
        I = identity();

        %INVERSE Generate the inverse of an element of the Lie group
        %   $ inverse(X) : \mathcal{M} \to \mathcal{M} $
        Xinv = inverse(X);
    end

end

