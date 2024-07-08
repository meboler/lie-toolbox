classdef TestSO2 < matlab.unittest.TestCase
    
    methods(TestClassSetup)
        % Shared setup for the entire test class
    end
    
    methods(TestMethodSetup)
        % Setup for each test
    end
    
    methods(Test)
        % Test methods
        function test_inverse(testCase)
            import lie.SO2
            x = randn();
            X = SO2.bigexp(x);
            Xinv = SO2.inverse(X);
            dX1 = Xinv * X;
            dX2 = X * Xinv;

            err1 = eye(2) - dX1;
            err2 = eye(2) - dX2;

            testCase.verifyLessThan(norm(err1), 1e-10);
            testCase.verifyLessThan(norm(err2), 1e-10);
        end

        function test_identity(testCase)
            import lie.SO2;
            x = randn();
            X = SO2.bigexp(x);
            I = SO2.identity();

            err1 = X - X * I;
            err2 = X - I * X;

            testCase.verifyLessThan(norm(err1), 1e-10);
            testCase.verifyLessThan(norm(err2), 1e-10);
        end

        function test_associativity(testCase)
            import lie.SO2;
            X = SO2.bigexp(randn());
            Y = SO2.bigexp(randn());
            Z = SO2.bigexp(randn());

            err = ( (X * Y) * Z) - ( X * (Y * Z) );
            
            testCase.verifyLessThan(norm(err), 1e-10);
        end
    end
    
end