## Runge-Kutta

Here we implement a Runge-Kutta method for solving a 
two-variable second-order differential equation
and use it to demonstrate chaotic motion in the plane:

![Double-pendulum](/images/double_pendulum.gif)

If we choose the initial angles to be vertical
to working precision, we get a very nice demonstration
of the inexactness of floating point representation.

![Double-pendulum](/images/double_pendulum_vertical.gif)

For the modeling, we simply write down the kinematic equations
and balance forces. To avoid some of the hairy algebra,
we use SymPy to solve for the second derivatives in terms
of the zeroth and first dervivatives. We then have a system
of first order differential equations in four variables which
is ripe for presentation to our Runge-Kutta solver.


```text
 ⎛                                                                                                        2                                        ⎞ 
 ⎜                                                                                             ⎛d        ⎞                                         ⎟ 
 ⎜                                                                l₁⋅m₂⋅sin(2⋅θ₁(t) - 2⋅θ₂(t))⋅⎜──(θ₁(t))⎟                                        2⎟ 
 ⎜                  g⋅m₂⋅sin(θ₁(t) - 2⋅θ₂(t))   g⋅m₂⋅sin(θ₁(t))                                ⎝dt       ⎠                             ⎛d        ⎞ ⎟ 
-⎜g⋅m₁⋅sin(θ₁(t)) + ───────────────────────── + ─────────────── + ───────────────────────────────────────── + l₂⋅m₂⋅sin(θ₁(t) - θ₂(t))⋅⎜──(θ₂(t))⎟ ⎟ 
 ⎝                              2                      2                              2                                                ⎝dt       ⎠ ⎠ 
─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
                                                           ⎛           2                    ⎞                                                        
                                                        l₁⋅⎝m₁ - m₂⋅cos (θ₁(t) - θ₂(t)) + m₂⎠                                                        


            ⎛                                                2⎞   ⎛                                                                        2⎞                   
            ⎜                                     ⎛d        ⎞ ⎟   ⎜                                                             ⎛d        ⎞ ⎟                   
- (m₁ + m₂)⋅⎜g⋅sin(θ₂(t)) - l₁⋅sin(θ₁(t) - θ₂(t))⋅⎜──(θ₁(t))⎟ ⎟ + ⎜g⋅m₁⋅sin(θ₁(t)) + g⋅m₂⋅sin(θ₁(t)) + l₂⋅m₂⋅sin(θ₁(t) - θ₂(t))⋅⎜──(θ₂(t))⎟ ⎟⋅cos(θ₁(t) - θ₂(t))
            ⎝                                     ⎝dt       ⎠ ⎠   ⎝                                                             ⎝dt       ⎠ ⎠                   
────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
                                                                ⎛           2                    ⎞                                                              
                                                             l₂⋅⎝m₁ - m₂⋅cos (θ₁(t) - θ₂(t)) + m₂⎠                                                              
```
