// Judge stage (ii) Stepsmoother test harness (msl_1335): function value and derivative for dx=1 at each parameter point
model StepHarness
  parameter Real func = 0.75;
  parameter Real nofunc = 0.25;
  parameter Real x = 0.5;
  Real r;
  Real d;
equation
  r = Modelica.Fluid.Dissipation.Utilities.Functions.General.Stepsmoother(func, nofunc, x);
  d = Modelica.Fluid.Dissipation.Utilities.Functions.General.Stepsmoother_der(func, nofunc, x, 0, 0, 1);
end StepHarness;
