# Boeing 787 forward fuselage pressurization: problem sheet

This sheet defines one structural analysis. Two writers use it, and only
it, to write the same model in two different input formats. A fresh LLM
also gets it for the LLM runs. The sheet uses no words from either format.

The model is an idealized barrel of the 787 forward fuselage under cabin
pressure, with one passenger door cutout and one door-corner submodel.
Each value has a source in [section 9](#9-references), or the mark
**[assumption]** and a reason. Boeing does not publish the frame pitch,
the stringer pitch, the section shapes or the laminates.

Units: SI (N, m, Pa, kg). All lengths are in metres unless the table
gives other units.

---

## 1. Researched values

| Item | Value | Source |
|---|---|---|
| Fuselage exterior width, 787-8 | 5.77 m (18 ft 11 in) | [2] fig. 2.2.1 |
| Fuselage exterior height, 787-8 | 5.94 m (19 ft 6 in) | [2] fig. 2.2.1 |
| Door 2 (main entry, left and right) station from the nose, 787-8 | 15.32 m (50 ft 3 in) | [2] table 2.7.1 |
| Passenger exit type at door 2 | Type A | [1] section on exit arrangements |
| Type A exit minimum clear opening | 42 in × 72 in (1.067 m × 1.829 m), corner radii not more than 7 in (0.178 m) | [3] 25.807(a) |
| Passenger window | about 10.7 in × 18.4 in (0.272 m × 0.467 m) | [8] (secondary source, not in the model) |
| Maximum operating altitude | 43,100 ft (13,137 m) pressure altitude | [1] section 11 |
| Cabin altitude | 6,000 ft (1,829 m) | [10]; this sheet applies it at the maximum operating altitude **[assumption]** |
| Fuselage construction | one-piece CFRP barrels, made by fiber placement, with co-cured stringers | [6], [7] |
| Stringer shape | hat section, co-cured with the skin | [9] |
| Prepreg | Toray T800S fiber with 3900 series toughened epoxy | [6], [7] |
| Ultimate factor | 1.5 on limit load | [3] 25.303 |
| Pressure-only condition | 1.33 × maximum relief valve setting, for operation up to 45,000 ft | [3] 25.365(d) |

### 1.1 Cabin pressure differential

The International Standard Atmosphere [5] gives:

| Altitude | Static pressure |
|---|---|
| 6,000 ft (cabin) | 81,200 Pa (11.78 psi) |
| 43,100 ft (outside) | 16,158 Pa (2.34 psi) |
| Difference | 65,042 Pa (9.43 psi) |

### 1.2 Ply properties, T800S/3900 unidirectional tape, room temperature

| Property | Value | Source |
|---|---|---|
| E1 | 148 GPa | [4], 0° tensile modulus |
| E2 | 10.1 GPa | [4], 90° tensile modulus |
| G12 | 3.94 GPa | [4] |
| ν12 | 0.32 | **[assumption]** The data sheet does not give it. 0.32 is a usual value for carbon/epoxy tape. |
| Cured ply thickness | 0.191 mm | [4] |

[4] gives values at three conditions. This sheet uses the room-temperature,
ambient column (RTA). **[assumption]** Cabin pressure acts in all flight
phases, and the study needs one set of values only.

---

## 2. Geometry

### 2.1 Coordinates

- x is the fuselage axis, positive aft.
- z is up. y is positive to the left side of the airplane (the door side).
- θ is the angle in the y-z plane from the top (+z) to the left side (+y).
  A point on the skin is at y = R sin θ, z = R cos θ.
- s = R θ is the arc length from the top.

### 2.2 Barrel

| Item | Value | Reason |
|---|---|---|
| Cross-section | circle | **[assumption]** The real section is 5.77 m wide and 5.94 m high. A circle makes the hand checks exact. |
| Skin mid-surface radius R | 2.883 m | **[assumption]** Half the 5.77 m width, less half the skin thickness, rounded. |
| Barrel length | 12.192 m (20 frame pitches), from x = −6.096 to x = +6.096 | **[assumption]** The door is 10 frame pitches from each end, far from end effects. |
| Frame pitch | 0.6096 m (24 in) | **[assumption]** A round value for a wide-body fuselage. Not public. |
| Frame stations | x = ±(0.3048 + 0.6096 k), k = 0 to 9 | The door centre is mid-way between two frames, so no frame lies on the plane x = 0. |
| Stringer count | 90 around the full circle | **[assumption]** Gives a stringer pitch of 0.2013 m (about 8 in). Not public. |
| Stringer positions | θ = 2° + 4° k, k = 0 to 89 | No stringer lies on the plane y = 0. |

The barrel has no windows. **[assumption]** The study measures the model
text, and the door corner is the critical detail. Windows add many
cutouts of one shape, which give no new information.

### 2.3 Door cutout

| Item | Value | Reason |
|---|---|---|
| Side | left (+y) | [2] |
| Width, along x | 1.067 m (42 in) | the Type A minimum, [3] 25.807(a) |
| Height, along the arc s | 1.829 m (72 in) | the Type A minimum, [3] 25.807(a) |
| Corner radius | 0.178 m (7 in) | the Type A maximum, [3] 25.807(a) |
| Centre | x = 0, θ = 75° | **[assumption]** 15° above the horizontal puts the sill near the cabin floor. |
| Edges | x = ±0.5334, s = 2.8593 (upper) to 4.6883 (lower), which is θ = 56.83° to 93.17° | from the size and the centre: s = R × 75° ∓ 1.829 / 2 |

The skin cutout is equal to the clear opening. **[assumption]** The real
cutout is larger, but its size is not public.

The cutout removes the skin, the stringers and the frames inside it. It
cuts the frames at x = ±0.3048 and the stringers at θ = 58° to 90°
(9 stringers on the left side).

### 2.4 Skin pad-up around the door

A band of thicker skin goes around the cutout. **[assumption]** A real
door surround is a thick, stiff region, and its design is not public.

| Item | Value |
|---|---|
| Band | all skin within 0.150 m of the cutout edge, measured in the x-s plane |
| Laminate | the pad-up laminate in section 3.2 |
| Mid-surface | the same as the skin mid-surface, at radius R |

---

## 3. Materials and laminates

### 3.1 Plies

All plies are T800S/3900 tape with the properties in section 1.2. The
skin uses the plies directly.

Ply angles: 0° is along +x. +45° turns from +x toward the direction of
increasing θ, seen from outside the fuselage. The stacking sequence
starts at the outer surface.

### 3.2 Laminates

| Name | Stacking sequence, outer to inner | Plies | Thickness |
|---|---|---|---|
| Skin | [45/0/−45/90]3s | 24 | 4.584 mm |
| Pad-up | [45/0/−45/90]6s | 48 | 9.168 mm |
| Stringer wall | [45/0/−45/90]s | 8 | 1.528 mm |
| Frame wall | [45/0/−45/90]2s | 16 | 3.056 mm |

**[assumption]** All four are quasi-isotropic (25/50/25). [4] uses this
family for its notched laminate tests. The skin thickness is near the
4.0 mm to 4.2 mm average that public mass and area estimates give.

### 3.3 Stringer and frame material

The stringers and frames are beams of one homogeneous, isotropic material.
Its values are the in-plane values of the quasi-isotropic laminate, which
classical laminate theory gives from section 1.2:

| Property | Value |
|---|---|
| E | 55.85 GPa |
| ν | 0.326 |

**[assumption]** The stringer and frame walls are quasi-isotropic, so
this is exact for the membrane stiffness. Because the skin has the same
in-plane modulus, the axial check in section 7 has a closed form.

---

## 4. Stringers and frames

Each stringer and each frame is a beam on its line on the skin
mid-surface. The beam shares the skin points on that line. The section
centroid is offset inward, along the radius, by the value in the table.
Transverse shear deformation of the beams may be included or not. It has
a negligible effect here.

### 4.1 Section shapes **[assumption]**

- **Stringer:** a hat. Two flanges, each 20 mm wide, lie on the inner
  skin face. Two vertical webs, 25 mm apart (outer faces), are 35 mm
  high from the inner skin face. A crown joins the webs. All walls are
  1.528 mm thick.
- **Frame:** a J-section. An outer flange, 50 mm wide, lies on the inner
  skin face, centred on the web. A web is 150 mm high from the inner skin
  face. An inner flange, 40 mm wide, goes to one side of the web only.
  All walls are 3.056 mm thick. The frame is
  continuous over the stringers (the model has no stringer
  clearance holes).

### 4.2 Section properties

The beams use these values. The shapes above only show where the values
come from.

| Property | Stringer | Frame |
|---|---|---|
| Area A | 2.016e-4 m² | 7.148e-4 m² |
| Centroid offset inward from the skin mid-surface | 0.01750 m | 0.07415 m |
| I about the centroidal axis parallel to the skin (radial bending) | 3.639e-8 m⁴ | 2.237e-6 m⁴ |
| I about the centroidal axis normal to the skin (lateral bending) | 4.908e-8 m⁴ | 8.305e-8 m⁴ |
| Torsion constant J (open section) | 1.569e-10 m⁴ | 2.225e-9 m⁴ |

The frame centroid is 3.16 mm from the web line, on the side of the
inner flange. The model ignores this lateral offset. **[assumption]**

---

## 5. Model extent and symmetry

The structure and the loads are symmetric about the plane y = 0 and the
plane x = 0. The model is one quarter of the barrel:
0 ≤ x ≤ 6.096 and 0° ≤ θ ≤ 180°.

| Plane | Condition, for the skin and for the beams |
|---|---|
| x = 0 | displacement along x is zero; rotations about y and z are zero |
| y = 0 (θ = 0° and θ = 180°) | displacement along y is zero; rotations about x and z are zero |
| one point: x = 0, θ = 180° (bottom) | displacement along z is zero |

The last condition only removes the rigid-body motion along z. Its
reaction must be near zero (section 7, check C4).

The end x = 6.096 is free except for the end closure load.

---

## 6. Analysis

### 6.1 Analysis type

One linear static analysis, with small displacements and linear
materials. **[assumption]** Large-displacement effects increase near the
cutout, but a linear model is enough for a study of model text, and it
makes the hand checks exact.

### 6.2 Loads

All loads are at the ultimate pressure p_U.

| Item | Value |
|---|---|
| Limit pressure p_L | 65.0 kPa (the differential in section 1.1, rounded) |
| Ultimate pressure p_U | 1.5 × p_L = 97.5 kPa, from [3] 25.303 |

This sheet takes the maximum relief valve setting as equal to p_L.
**[assumption]** The setting is not public. Then [3] 25.365(d) gives
1.33 × 65.0 = 86.5 kPa, which is less than p_U, so 25.303 governs.

1. **Internal pressure:** p_U on the skin and pad-up, normal to the
   mid-surface and outward. The loaded area is the mid-surface area at
   radius R. No pressure acts in the cutout.
2. **End closure load:** at x = 6.096, a force along +x on the skin edge,
   uniform per unit arc length, equal to p_U R / 2 = 140.5 kN/m. The
   total for the quarter model is p_U π R² / 2 = 1273.0 kN.
3. **Door load:** the door is not in the model. A plug door transfers its
   pressure load to the cutout edges. Apply the door load on the aft half
   of the cutout edge in the quarter model: the straight aft edge and the
   two aft corner arcs. At each point, the load is outward along the local
   radius. Its value per unit arc length in s (not per unit edge length) is
   p_U × x_e(s). Here, x_e(s) is the x coordinate of the cutout edge at
   that s. On the straight aft edge, this is p_U × 0.5334 = 52.01 kN/m.
   The strip of the door from x = 0 to x_e(s) puts its pressure on the
   edge at the same s, so the resultant is equal to the pressure on the
   half door. The half door area is 0.9620 m², and the resultant is
   89.12 kN along +y and 23.88 kN along +z.

### 6.3 Mesh

| Region | Element type | Size |
|---|---|---|
| Skin and pad-up | quadrilateral shell, with triangles only where the cutout corners need them | about 0.050 m |
| Stringers and frames | beam, one element per skin element edge on the beam line | the same as the skin |

The element order is the choice of the writer. The writer states it in
the model.

### 6.4 Door-corner submodel

A second analysis gives a refined result at the upper aft door corner.

| Item | Value |
|---|---|
| Region | a square in the x-s plane, 0.600 m × 0.600 m, centred on the centre of the corner radius: x_c = 0.5334 − 0.178 = 0.3554, s_c = s_upper + 0.178 = 3.0373 m, where s_upper = 2.8593 m (section 2.3). The region contains the skin, the pad-up and the parts of the stringers and frames in it. |
| Idealization | the same as the global model: shells for the skin and pad-up, offset beams for the stringers and frames |
| Mesh | about 0.005 m on the cutout edge, and not more than 0.020 m elsewhere |
| Boundary | on the region boundary, the translations and rotations of the skin and the beam ends come from the global result |
| Loads | the same pressure and door edge load as the global model, inside the region |
| Material and laminates | the same as the global model |

---

## 7. Required outputs and built-in checks

### 7.1 Outputs

Global model:

1. The radial displacement over the skin, and its maximum value and point.
2. The skin hoop and axial membrane stress over the skin. Membrane stress
   is the laminate force per unit width divided by the laminate
   thickness.
3. The axial force in each stringer and the hoop force in each frame, at
   the maximum value for each beam.
4. The sum of the reaction forces in x, y and z.

Submodel:

5. The maximum principal membrane strain on the cutout edge, and its
   point.
6. The strain concentration factor K = (item 5) / (the hoop membrane
   strain at the check point P of the global model).
7. The margin against an ultimate strain allowable of 4,000 µε:
   MS = 4000 µε / (item 5) − 1. **[assumption]** The allowable is a usual
   value for damaged carbon/epoxy laminates. Boeing does not publish it.

### 7.2 Checks

The check point P is in plain skin, far from the door and the end:
x = 4.8768 (mid-way between two frames), θ = 140° (mid-way between two
stringers). The values below are at p_U, with t = 4.584 mm.

| Check | Quantity | Expected value | Pass when |
|---|---|---|---|
| C1 | Hoop membrane stress at P | between p_U R / (t + A_frame / frame pitch) = 48.8 MPa and p_U R / t = 61.3 MPa | inside the range; report the ratio to p_U R / t |
| C2 | Axial membrane stress at P | (p_U R / 2) / (t + A_stringer / stringer pitch) = 25.2 MPa | within 5 %; report the ratio to p_U R / (2t) = 30.7 MPa |
| C3 | Sum of x reactions | −1273.0 kN (equal and opposite to the end closure load) | within 0.1 % |
| C4 | Sum of z reactions | 0 | less than 0.1 % of 1273.0 kN |
| C5 | Sum of y reactions | −p_U × 2R × 6.096 = −3427.1 kN (the pressure resultant on the left half, including the door load) | within 0.1 % |
| C6 | Submodel boundary | the hoop membrane strain at the mid-point of each region edge that has skin (the aft and top edges; the forward and lower edge mid-points are in the cutout) agrees with the global model | within 5 % |

C1 is a range because the frames carry part of the hoop load near each
frame. In the middle of the bay, the skin carries most of it.

---

## 8. Request for the LLM runs

Each run gets this sheet and this one line, with the target format name in
place of `<format>`:

> Using only this problem sheet, write one complete `<format>` input model
> for the global analysis and the door-corner submodel in sections 2 to 7.

---

## 9. References

1. EASA, Type-Certificate Data Sheet EASA.IM.A.115, Boeing 787, issue 30,
   17 December 2025. <https://www.easa.europa.eu/en/downloads/7302/en>
2. Boeing, *787 Airplane Characteristics for Airport Planning*, D6-58333,
   revision Q, October 2025.
   <https://www.boeing.com/content/dam/boeing/v2/airports/acaps/787_ACAP_Rev_Q.pdf>
3. 14 CFR Part 25, sections 25.303 (factor of safety), 25.365 (pressurized
   compartment loads) and 25.807 (emergency exits).
   <https://www.law.cornell.edu/cfr/text/14/part-25>
4. Toray Composite Materials America, *3900 Prepreg System* data sheet,
   P2362W-19L (T800S-24K unidirectional) lamina properties.
   <https://www.toraycma.com/wp-content/uploads/3900-Prepreg-System.pdf>
5. NOAA, NASA and USAF, *U.S. Standard Atmosphere, 1976*, NASA-TM-X-74335.
   The equations below 20 km are the same as the ICAO standard atmosphere.
6. CompositesWorld, *Boeing 787 update*.
   <https://www.compositesworld.com/articles/boeing-787-update>
7. CompositesWorld, *Plant tour: Spirit AeroSystems, Wichita, KS* (787
   section 41). <https://www.compositesworld.com/articles/the-first-composite-fuselage-section-for-the-first-composite-commercial-jet>
8. Aviation A2Z, *This aircraft has the largest airplane windows in the
   world*, 19 August 2025.
   <https://aviationa2z.com/index.php/2025/08/19/this-aircraft-has-largest-airplane-windows-in-the-world/>
9. IntechOpen, *The Evolution of the Composite Fuselage: A Manufacturing
   Perspective*. <https://www.intechopen.com/chapters/64957>
10. Boeing, *787 Dreamliner By Design*.
    <https://www.boeing.com/commercial/787/by-design>
