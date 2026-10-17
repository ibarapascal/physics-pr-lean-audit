# Boilers: efficiency as useful heat, and heat loss to the ambient

In `Buildings.Fluid.Boilers.BaseClasses.PartialBoiler` (used by `Buildings.Fluid.Boilers.BoilerPolynomial` and `Buildings.Fluid.Boilers.BoilerTable`), the heat loss of the boiler to the ambient through `UAOve` and `heatPort` is not correctly accounted for in the energy balance. Please revise the boiler models as follows:

- The efficiency `eta` and the parameter `eta_nominal` are to be interpreted as useful heat divided by fuel use, where the useful heat is the difference of the enthalpy streams of the working fluid between outlet and inlet. `Q_flow_nominal` is the nominal heating power (useful heat).
- The heat released by the fuel must also cover the conduction losses from the water through the boiler insulation into the ambient air.
- Add variables that report the heat released by the fuel (`QFue_flow`), the useful heat (`QUse_flow`), the heat loss to the ambient (`QAmb_flow`) and the heat exhausted by the flue gas (`QExh_flow`).
- With this change, `heatPort` must be connected, as otherwise the heat loss is not properly calculated. Translation should stop with an error message if `heatPort` is unconnected. Update the examples that use the boiler so that `heatPort` is connected to the ambient temperature.
- Check whether the steam boiler needs the same change.
- Clean up the code and update the documentation and the revision notes.

---

You are a maintainer of the Modelica Buildings Library (lbl-srg/modelica-buildings). Implement this request.
