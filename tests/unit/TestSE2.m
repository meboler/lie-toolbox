classdef TestSE2 < matlab.unittest.TestCase
    
    methods(TestClassSetup)
        % Shared setup for the entire test class
    end
    
    methods(TestMethodSetup)
        % Setup for each test
    end
    
    methods(Test)
        % Test methods
        function test_inverse(testCase)
            import lie.SE2;
            x = randn(3, 1);
            X = SE2.bigexp(x);
            Xinv = SE2.inverse(X);
            dX1 = Xinv * X;
            dX2 = X * Xinv;

            err1 = eye(3) - dX1;
            err2 = eye(3) - dX2;

            testCase.verifyLessThan(norm(err1), 1e-10);
            testCase.verifyLessThan(norm(err2), 1e-10);
        end

        function test_identity(testCase)
            import lie.SE2;
            x = randn(3, 1);
            X = SE2.bigexp(x);
            I = SE2.identity();

            err1 = X - X * I;
            err2 = X - I * X;

            testCase.verifyLessThan(norm(err1), 1e-10);
            testCase.verifyLessThan(norm(err2), 1e-10);
        end

        function test_associativity(testCase)
            import lie.SE2;
            X = SE2.bigexp(randn(3, 1));
            Y = SE2.bigexp(randn(3, 1));
            Z = SE2.bigexp(randn(3, 1));

            err = ( (X * Y) * Z) - ( X * (Y * Z) );
            
            testCase.verifyLessThan(norm(err), 1e-10);
        end
    end
    
end