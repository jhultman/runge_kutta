import fire
from sympy import Function, solve, cos, sin, symbols, diff, init_printing, pprint

init_printing(use_latex=False)


def solve_symbolically():
    t = symbols('t')
    m1, m2, l1, l2, g = symbols('m1 m2 l1 l2 g')

    theta1 = Function('theta1')(t)
    theta2 = Function('theta2')(t)

    x1 = l1 * sin(theta1)
    y1 = -l1 * cos(theta1)

    x2 = x1 + l2 * sin(theta2)
    y2 = y1 - l2 * cos(theta2)

    d2x1 = x1.diff(t, 2)
    d2y1 = y1.diff(t, 2)

    d2x2 = x2.diff(t, 2)
    d2y2 = y2.diff(t, 2)

    d2theta1 = theta1.diff(t, 2)
    d2theta2 = theta2.diff(t, 2)

    LHS = -cos(theta1) * (m1 * d2x1 + m2 * d2x2)
    RHS = sin(theta1) * (m1 * d2y1 + m2 * d2y2 + m2 * g + m1 * g)
    eqn_1 = LHS - RHS

    LHS = cos(theta2) * m2 * d2x2
    RHS = -sin(theta2) * (m2 * g + m2 * d2y2)
    eqn_2 = LHS - RHS

    system = [eqn_1, eqn_2]
    soln = solve(system, d2theta1, d2theta2)

    pprint(soln[d2theta1])
    print('\n' * 1)
    pprint(soln[d2theta2])


if __name__ == '__main__':
    fire.Fire()
