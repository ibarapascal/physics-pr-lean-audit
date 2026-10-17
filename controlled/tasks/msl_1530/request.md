# PartialSimpleIdealGasMedium uses different offsets in BaseProperties and in specificInternalEnergy

When simulating the test model ModelicaTest.Media.TestAllProperties.SimpleAir an error occurs. The reason are different offsets in the medium model.

In PartialSimpleIdealGasMedium.BaseProperties:
```
h = specificEnthalpy_pTX(p,T, X);  // = cp_const*(T - T0);
u = h - R*T; // = (cp_const - R)*T - cp_const*T0
```
In PartialSimpleIdealGasMedium.specificInternalEnergy:
```
u = (cp_const - R)*(T - T0);
    // = (cp_const - R)*T - cp_const*T0 - (cp_const - R)*T0
```
Therefore in the first case the offset is
 cp_const*T0
and in the second case it is
 (cp_const - R)*T0

---

You are a maintainer of the Modelica Standard Library (modelica/ModelicaStandardLibrary). Implement this request.
