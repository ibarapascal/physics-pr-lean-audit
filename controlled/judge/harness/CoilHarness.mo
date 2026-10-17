// Judge stage (ii) wet-coil test harness (buildings_3027): steady-state WetCoilCounterFlow, mass-flow sources on both sides, both air-side boundaries at the same temperature and humidity (m2 may reverse); outputs net enthalpy inflow and net water-vapour inflow over the four ports (public quantities only)
model CoilHarness
  package Medium1 = Buildings.Media.Water;
  package Medium2 = Buildings.Media.Air;
  parameter Real T_w = 280.15;
  parameter Real T_a = 300.15;
  parameter Real X_w = 0.01;
  parameter Real m1_rel = 1;
  parameter Real m2_rel = 1;
  parameter Real m1_nom = 3.78;
  parameter Real m2_nom = 2.646;
  parameter Real UA = 4748;
  Buildings.Fluid.HeatExchangers.WetCoilCounterFlow hex(
    redeclare package Medium1 = Medium1,
    redeclare package Medium2 = Medium2,
    m1_flow_nominal = m1_nom,
    m2_flow_nominal = m2_nom,
    dp1_nominal = 3000,
    dp2_nominal = 200,
    allowFlowReversal1 = true,
    allowFlowReversal2 = true,
    UA_nominal = UA,
    energyDynamics = Modelica.Fluid.Types.Dynamics.SteadyState);
  Buildings.Fluid.Sources.MassFlowSource_T sou1(redeclare package Medium = Medium1, m_flow = m1_rel*m1_nom, T = T_w, nPorts = 1);
  Buildings.Fluid.Sources.Boundary_pT sin1(redeclare package Medium = Medium1, p = 300000, T = T_w, nPorts = 1);
  Buildings.Fluid.Sources.MassFlowSource_T sou2(redeclare package Medium = Medium2, m_flow = m2_rel*m2_nom, T = T_a, X = {X_w, 1 - X_w}, nPorts = 1);
  Buildings.Fluid.Sources.Boundary_pT sin2(redeclare package Medium = Medium2, p = 101325, T = T_a, X = {X_w, 1 - X_w}, nPorts = 1);
  Real H1 = hex.port_a1.m_flow*actualStream(hex.port_a1.h_outflow) + hex.port_b1.m_flow*actualStream(hex.port_b1.h_outflow);
  Real H2 = hex.port_a2.m_flow*actualStream(hex.port_a2.h_outflow) + hex.port_b2.m_flow*actualStream(hex.port_b2.h_outflow);
  Real mW2 = hex.port_a2.m_flow*actualStream(hex.port_a2.Xi_outflow[1]) + hex.port_b2.m_flow*actualStream(hex.port_b2.Xi_outflow[1]);
  Real mWat = hex.mWat_flow;
  Real Q1 = hex.Q1_flow;
  Real Q2 = hex.Q2_flow;
equation
  connect(sou1.ports[1], hex.port_a1);
  connect(hex.port_b1, sin1.ports[1]);
  connect(sou2.ports[1], hex.port_a2);
  connect(hex.port_b2, sin2.ports[1]);
end CoilHarness;
