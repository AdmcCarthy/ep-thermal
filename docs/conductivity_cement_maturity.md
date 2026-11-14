# Conductivity, quartz cement and maturity

Three small models that turn a temperature or a rock description into a property: how well a rock conducts heat, how much quartz cement a burial history makes, and how mature the organic matter is. All three read packaged defaults (`conductivity.yaml`, `quartz_cement.yaml`, `maturity.yaml`) at the level **published, assumed valid for use**. They are to be checked against data when it exists; the checks are listed at the end.

## Conductivity (`ep_thermal.conductivity`)

The bulk conductivity of a rock is the geometric mean of its components weighted by volume fraction, `k = prod(k_i ^ f_i)` (Allen and Allen 2005, Table 9.2 for the values, Fig. 9.13a for the check). `bulk_conductivity` takes any set of components; `sandstone_conductivity` builds quartz, clay and water from a porosity and a clay share of the solids, which reproduces the three curves of the book's plot. `in_rock_range` says whether a value lies in the range the book gives for a rock type.

Limits.
- The values are for surface temperature. Conductivity falls with temperature, and no correction is made.
- The book's two tables disagree on dolomite rock (3.2 to 3.5 here, 3.2 to 5 in Table 2.1).
- Fractions that do not sum to one, or an unknown component, give NaN and an event. Nothing is renormalised.

## Quartz cement (`ep_thermal.quartz_cement`)

The rate of quartz precipitation is `r = a 10^(bT)` in mol per cm2 per second with `T` in degrees C, and zero below 80 degrees C (Walderhaug 1996). The cement made in a step is `M r A dt / rho`, where the quartz surface area `A` is that of spheres (`6 q / d`), reduced by any clay coat and falling in proportion to the porosity left. Within a step of constant temperature the solution is exact, so the answer does not depend on the step size.

`b` is not a single number. Across 19 Norwegian shelf wells the best fit clusters at 0.020 to 0.021 per degree C with a spread of 0.018 to 0.023, and the author puts much of the spread down to uncertain temperature histories. `quartz_cement_band` therefore runs the low, central and high value and returns all three. Report the band.

Stated assumptions: the area falls in proportion to the intergranular porosity left (the paper says "in proportion to the porosity filled"); porosity is lost to cement only, with no compaction inside the step; the temperature of a step is the mean of its ends. The value of `a` is quoted from Walderhaug (1994a), which has not been read.

## Maturity (`ep_thermal.maturity`)

EASY%Ro (Sweeney and Burnham 1990) treats vitrinite reflectance as 20 parallel first-order reactions with activation energies from 34 to 72 kcal/mol and a common frequency factor. The temperature integral of each reaction uses a rational approximation, the fraction reacted is a weighted sum, and `Ro = exp(-1.6 + 3.7 F)` per cent. At no conversion `Ro` is 0.20 per cent; at full conversion (F = 0.85) it is 4.69 per cent.

Departure from the open-source code the constants were taken from. That code replaces a zero heating rate by 1e-30, which gives a step with no temperature change no reaction at all. Here such a step is isothermal and uses the direct rate `A exp(-E / RT) dt`, so a long hold at constant temperature matures the rock. The two paths agree: a hold and a ramp of 0.01 degrees over the same time give the same reflectance to 0.5 per cent (tested).

Tests reproduce five check values recorded in the knowledge store (a ramp from 20 to 150 degrees C over 50 Myr gives 0.46 per cent at 85 and 1.00 at 150; 10 degrees C per Myr from 20 gives 0.36 at 70, 0.65 at 120 and 1.18 at 170). They come from an implementation of the same printed equations, so they check the arithmetic, not the model. The paper's own table has not been reproduced.

## To verify when more data are available

| What | Against |
| :--- | :--- |
| Conductivity values and ranges | Measured conductivity on core, or the gradient implied by a heat flow for a well of known lithology |
| The cement band | Cement and porosity measured in core of a well with a known temperature history |
| EASY%Ro constants | The printed table of Sweeney and Burnham (1990), then measured vitrinite reflectance |
