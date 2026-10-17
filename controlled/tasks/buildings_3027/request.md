# The wet coil generates an unrealistic sensible heat exchange rate

In <code>Buildings.Fluid.HeatExchangers.WetCoilCounterFlow</code>, the sensible heat exchange rate is calculated by:
```
  Modelica.Units.SI.HeatFlowRate QSen2_flow=Q2_flow - QLat2_flow
    "Sensible heat input into air stream (negative if air is cooled)";
```
If <code>Q2_flow</code> is lower than <code>QLat2_flow</code>, then <code>QSen2_flow</code> is positive even <code>Q2_flow</code> and<code>QLat2_flow</code> are negative.
In this case, the coil is heating the air.

This error can be triggered when the water flow rate is lower while the air flow rate is still high (e.g., the air-side economizer is on). In this case, the <code>Q2_flow</code> is lower but <code>QLat2_flow</code> may be higher since it is not directly related to the water side condition, as shown below:

[image omitted]

---

You are a maintainer of the Modelica Buildings Library (lbl-srg/modelica-buildings). Implement this request.
