// Judge stage (ii) ideal-gas test harness (msl_1530): h-u, p/d and R_gas*T of SimpleAir BaseProperties and specificInternalEnergy at a given (p,T)
model GasHarness
  package Medium = Modelica.Media.Air.SimpleAir;
  parameter Real pp = 1e5;
  parameter Real TT = 300;
  Medium.BaseProperties med;
  Real hmu;
  Real RT;
  Real pd;
  Real du;
  Real u;
equation
  med.p = pp;
  med.T = TT;
  hmu = med.h - med.u;
  RT = Medium.R_gas*med.T;
  pd = med.p/med.d;
  du = med.u - Medium.specificInternalEnergy(med.state);
  u = med.u;
end GasHarness;
