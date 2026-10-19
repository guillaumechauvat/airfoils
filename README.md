# airfoils
Parametrises an airfoil as a smooth spline given input coordinate points.
The coordinates must run from trailing edge to trailing edge.

## Example usage:

```py
import airfoil_spline as af

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

The leading edge position can be found directly via `le()`, which looks for a horizontal normal:
```py
tle = airfoil.le()
```

The coordinates for a given parameter `t` can be queried using `pos()`:
```py
# suction side
p3s = airfoil.pos(t3s)  # array([0.3      , 0.0771686])

# pressure side
p3p = airfoil.pos(t3p)  # array([ 0.3       , -0.03793124])

# or both at the same time:
p3 = airfoil.normal([t3s, t3p])  # array([[ 0.3       ,  0.0771686 ],
                                 #        [ 0.3       , -0.03793124]])
```

Normal vectors can also be calculated:
```py
n3 = airfoil.normal([t3s, t3p])  # array([[-0.08282943,  0.99656374],
                                 #        [-0.03497633, -0.99938814]])
```

Finally, given a point sufficiently close to the airfoil, the closest point on the airfoil can be found with `closest()`:
```py
t_closest = airfoil.closest(0.1, 0.04, 0.25)  # array([0.4479957791950058])
airfoil.pos(t_closest)  # array([[0.09823843, 0.04730137]])
# be carteful about starting point:
t_closest = airfoil.closest(0.2, 0.07, 0.5)  # array([0.5040716]) - wrong!
t_closest = airfoil.closest(0.2, 0.07, 0.25)  # array([0.39682072]) - ok
airfoil.pos(t_closest)  # array([[0.2005281 , 0.06637495]])

```
