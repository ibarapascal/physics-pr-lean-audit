# Problem with StepSmoother

In MSL 3.2.1 the Modelica.Fluid.Dissipation.Utilities.Functions.General.Stepsmoother function is written as follows:
```mo
        function Stepsmoother "Continuous interpolation for x "

          extends Modelica.Icons.Function;
          input Real func "input value for that result = 100%";
          input Real nofunc "input value for that result = 0%";
          input Real x "input variable for continuous interpolation";
          output Real result "output value";

        protected
          Real m=Modelica.Constants.pi/(func - nofunc);
          Real b=-Modelica.Constants.pi/2 - m*nofunc;
          Real r_1=tan(m*x + b);

        algorithm
          result := if x >= 0.999999*(func - nofunc) + nofunc and func > nofunc or x
             <= 0.999999*(func - nofunc) + nofunc and nofunc > func then 1 else if x
             <= 0.000001*(func - nofunc) + nofunc and func > nofunc or x >= 0.000001*(
            func - nofunc) + nofunc and nofunc > func then 0 else ((0.5*(exp(r_1) - exp(
            -r_1))/(0.5*(exp(r_1) + exp(-r_1))) + 1)/2);
          annotation (
            derivative=Stepsmoother_der,
            Documentation(info="<html> ... </html>"),  smoothOrder=5);
        end Stepsmoother;
```
The problem is that the value of r_1 is computed unnecessarily, and in cases where it might introduce a NaN. In addition, the final clause of the if can be written much more efficiently.

These two problems can be easily corrected by re-writing the function as follows:
```mo
        function Stepsmoother "Continuous interpolation for x "

          extends Modelica.Icons.Function;
          input Real func "input value for that result = 100%";
          input Real nofunc "input value for that result = 0%";
          input Real x "input variable for continuous interpolation";
          output Real result "output value";

        protected
          Real m=Modelica.Constants.pi/(func - nofunc);
          Real b=-Modelica.Constants.pi/2 - m*nofunc;

        algorithm
          result := if x >= 0.999999*(func - nofunc) + nofunc and func > nofunc or x
             <= 0.999999*(func - nofunc) + nofunc and nofunc > func then 1 else if x
             <= 0.000001*(func - nofunc) + nofunc and func > nofunc or x >= 0.000001*(
            func - nofunc) + nofunc and nofunc > func then 0 else
            (1-Modelica.Math.tanh(Modelica.Math.tan(m*x + b)))/2);
          annotation (
            derivative=Stepsmoother_der,
            Documentation(info="<html> ... </html>"),  smoothOrder=5);
        end Stepsmoother;
```
Note that this problem is not present in Stepsmoother_der, the derivative of this function. In that case the efficient form is already in use.

---

You are a maintainer of the Modelica Standard Library (modelica/ModelicaStandardLibrary). Implement this request.
