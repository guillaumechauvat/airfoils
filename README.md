# airfoils
Airfoil meshing tools

## Example usage:

```py
import airfoils as af

# coordinates of the NACA 65(1)-412 from https://m-selig.ae.illinois.edu/ads/coord/naca651412.dat
x = [1, 0.95, 0.901, 0.851, 0.801, 0.751, 0.701, 0.651, 0.601, 0.55, 0.5, 0.45, 0.399, 0.349, 0.298, 0.249, 0.198, 0.148, 0.097, 0.072, 0.048, 0.023, 0.011, 0.006, 0.003, 0, 0.007, 0.009, 0.014, 0.027, 0.052, 0.078, 0.103, 0.152, 0.202, 0.252, 0.302, 0.351, 0.401, 0.45, 0.5, 0.55, 0.599, 0.649, 0.699, 0.749, 0.799, 0.849, 0.899, 0.95, 1]
y = [0, 0.01, 0.02, 0.03, 0.039, 0.048, 0.057, 0.064, 0.071, 0.076, 0.08, 0.081, 0.081, 0.08, 0.077, 0.072, 0.066, 0.057, 0.047, 0.04, 0.032, 0.022, 0.016, 0.012, 0.01, 0, -0.008, -0.01, -0.012, -0.015, -0.02, -0.023, -0.026, -0.03, -0.034, -0.036, -0.038, -0.039, -0.039, -0.038, -0.036, -0.032, -0.028, -0.023, -0.018, -0.013, -0.008, -0.003, 0.001, 0.003, 0]

# generate the interpolated airfoil
airfoil = af.Airfoil(x, y)
```

The airfoil is parametrised by `t ∈ [0, 1]`: `t = 0` and `t = 1` should be the trailing edge, while `t = 0.5` should be somewhere near the leading edge.
The value of `t` corresponding to a specific chord position can be calculated with `curvilinear_abscissa()`, giving a first approximation for `t` since there are usually multiple solutions:
```py
# suction side
t3s = airfoil.curvilinear_abscissa(0.3, 0.25)  # 0.34761588340824057

# pressure side
t3p = airfoil.curvilinear_abscissa(0.3, 0.75)  # 0.65480812484917
```

The coordinates for a given parameter `t` can be queried using `pos()`:
```py
# suction side
x3s, y3s = airfoil.pos(t3s)  # 0.3, 0.0771686

# pressure side
x3p, y3p = airfoil.pos(t3p)  # 0.3, -0.03793124
```

The normal vector can also be calculated:
```py
# suction side
nx3s, ny3s = airfoil.normal(t3s)  # -0.08282943442805163, 0.9965637384494426

# pressure side
nx3p, ny3p = airfoil.normal(t3p)  # -0.0349763290619618, -0.9993881410169672
```

Finally, given a point sufficiently close to the airfoil, the closest point on the airfoil can be found with `closest()`:
```py
t_closest, x_closest, y_closest = airfoil.closest(0.1, 0.04)  # 0.4479957791950058, 0.09823843, 0.04730137
```
