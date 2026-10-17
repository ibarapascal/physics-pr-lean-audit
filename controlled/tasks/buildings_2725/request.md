# Correct boiler heat flow with loss to the ambient

First discussed in comments in #2695.
In `Buildings.Fluid.Boilers`, the thermal conduction between the boiler and the ambient is considered (via `UAOve`), but this is factored into the calculation of the heat from the fuel to the working fluid `QWat_flow`.
This part of the calculation is also explained in the documentation.

---

You are a maintainer of the Modelica Buildings Library (lbl-srg/modelica-buildings). Implement this request.
