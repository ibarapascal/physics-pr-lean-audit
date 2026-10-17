// Judge stage (ii) boiler test harness (shared by buildings_2725 / buildings_4775): steady-state BoilerPolynomial, constant efficiency a={eta}, heatPort connected to a fixed ambient temperature; flue-gas heat E computed from port quantities
model BoilerHarness
  package Medium = Buildings.Media.Water;
  parameter Real y_in = 0.5;
  parameter Real eta = 0.9;
  parameter Real T_in = 323.15;
  parameter Real T_amb = 288.15;
  parameter Real Q_nom = 3000;
  parameter Real m_rel = 1;
  parameter Real m_nom = Q_nom/20/4200;
  Buildings.Fluid.Boilers.BoilerPolynomial boi(
    redeclare package Medium = Medium,
    m_flow_nominal = m_nom,
    dp_nominal = 3000,
    Q_flow_nominal = Q_nom,
    fue = Buildings.Fluid.Data.Fuels.NaturalGasLowerHeatingValue(),
    effCur = Buildings.Fluid.Types.EfficiencyCurves.Constant,
    a = {eta},
    energyDynamics = Modelica.Fluid.Types.Dynamics.SteadyState,
    T_start = T_in);
  Buildings.Fluid.Sources.MassFlowSource_T sou(redeclare package Medium = Medium, m_flow = m_rel*m_nom, T = T_in, nPorts = 1);
  Buildings.Fluid.Sources.Boundary_pT sin(redeclare package Medium = Medium, nPorts = 1);
  Modelica.Thermal.HeatTransfer.Sources.FixedTemperature TAmb(T = T_amb);
  Modelica.Blocks.Sources.Constant yS(k = y_in);
  Real QFue = QFUE_EXPR;
  Real Qport = boi.heatPort.Q_flow;
  Real Hw = boi.port_a.m_flow*(boi.port_b.h_outflow - inStream(boi.port_a.h_outflow));
  Real E = QFue + Qport - Hw;
equation
  connect(sou.ports[1], boi.port_a);
  connect(boi.port_b, sin.ports[1]);
  connect(TAmb.port, boi.heatPort);
  connect(yS.y, boi.y);
end BoilerHarness;
