import tqdm
import numpy as np
from numpy import sin, cos
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


class DoublePendulum:
    
    def __init__(self, m1=1.0, m2=1.0, l1=1.0, l2=1.0, g=9.81):
        self.m1 = m1
        self.m2 = m2
        self.l1 = l1
        self.l2 = l2
        self.g = g
    
    def q_prime(self, q, r, u, v, t):
        val = u
        return val

    def r_prime(self, q, r, u, v, t):
        val = v
        return val

    def u_prime(self, q, r, u, v, t):
        m1, m2, l1, l2, g = self.m1, self.m2, self.l1, self.l2, self.g
        num = 2 * g * m1 * sin(q) + g * m2 * sin(q) + g * m2 * sin(q - 2 * r)
        num += l1 * m2 * u ** 2 * sin(2 * (q - r)) + 2 * l2 * m2 * v ** 2 * sin(q - r)
        den = -2 * l1 * (m1 - m2 * cos(q - r) ** 2 + m2)
        val = num / den
        return val

    def v_prime(self, q, r, u, v, t):
        m1, m2, l1, l2, g = self.m1, self.m2, self.l1, self.l2, self.g
        num = -(m1 + m2) * (g * sin(r) - l1 * u ** 2 * sin(q - r))
        num += (cos(q - r)) * (g * m1 * sin(q) + g * m2 * sin(q) + l2 * m2 * v ** 2 * sin(q - r))
        den = l2 * (m1 - m2 * cos(q - r) ** 2 + m2)
        val = num / den
        return val
    
    def derivatives(self, q, r, u, v, t):
        args = (q, r, u, v, t)
        qruv_prime = np.array([
            self.q_prime(*args),
            self.r_prime(*args),
            self.u_prime(*args),
            self.v_prime(*args),
        ])
        return qruv_prime

    def x1(self, q, r):
        val = +self.l1 * sin(q)
        return val

    def y1(self, q, r):
        val = -self.l1 * cos(q)
        return val

    def x2(self, q, r):
        val = self.x1(q, r) + self.l2 * sin(r)
        return val

    def y2(self, q, r):
        val = self.y1(q, r) - self.l2 * cos(r)
        return val
    
    def positions(self, q, r):
        args = (q, r)
        vals = np.array([
            self.x1(*args),
            self.y1(*args),
            self.x2(*args),
            self.y2(*args),
        ])
        return vals


class RungeKuttaSimulator:
    
    def __init__(self, pendulum, q0=0, r0=0, u0=0, v0=0, t0=0, h=0.05, T=15):
        self.pendulum = pendulum
        self._init_state = np.array([q0, r0, u0, v0, t0])
        self.t = np.mgrid[t0 : T : h]
        self.h = h
    
    def _get_state(self, qruvti):
        q, r, u, v, t = qruvti
        state = [*self.pendulum.positions(q, r), q, r, t]
        return state
    
    def simulate(self):
        history = []
        qruvti = self._init_state
        
        h = self.h
        f = self.pendulum.derivatives
        
        print('Simulating...')
        for ti in tqdm.tqdm(self.t):
            history += [self._get_state(qruvti)]
            k1 = h * f(*qruvti)

            offset = 0.5 * np.array([*k1, h])
            k2 = h * f(*(qruvti + offset))

            offset = 0.5 * np.array([*k2, h])
            k3 = h * f(*(qruvti + offset))

            offset = 1.0 * np.array([*k3, h])
            k4 = h * f(*(qruvti + offset))

            offset = np.concatenate((1 / 6 * (k1 + 2 * k2 + 2 * k3 + k4), (h,)))
            qruvti = qruvti + offset
        
        history = np.array(history)
        return history


class Animator:
    '''Maintains state for FuncAnimation'''
    
    def __init__(self, state, title, incr):
        self.state = state
        self.incr = incr
        self.title = title
        self.fig, self.ax = plt.subplots()

    def _get_savepath(self):
        fname = self.title.lower().replace(' ', '_')
        savepath = f'../images/{fname}.gif'
        return savepath

    def _get_limits(self):
        min_vals = [-2, -2]
        max_vals = [2, 2]
        limits = np.vstack((min_vals, max_vals))
        return limits.T

    def _configure_axes(self):
        limits = self._get_limits()
        self.ax.set_title(self.title)
        self.ax.set_xlim(limits[0, :])
        self.ax.set_ylim(limits[1, :])
        self.ax.set_aspect('equal')
        self.ax.set_axis_off()

    def _init_ani(self):
        self._configure_axes()
        self.masses = self.ax.scatter([], [], s=25, c='black')
        self.arms, = self.ax.plot([], [], c='grey', linestyle='dotted')
        return self.masses, self.arms
    
    def _step(self, frame):
        index = self.incr * frame
        xy = self.state[index, :4].reshape(2, 2)
        xy = np.concatenate((np.zeros((1, 2)), xy), axis=0)
        self.masses.set_offsets(xy)
        self.arms.set_xdata(xy[:, 0])
        self.arms.set_ydata(xy[:, 1])
        return self.masses, self.arms
    
    def animate(self):
        num_points = len(self.state)
        frames = tqdm.trange(num_points // self.incr)
        savepath = self._get_savepath()

        print('Animating...')
        ani = FuncAnimation(
            fig=self.fig, 
            func=self._step, 
            init_func=self._init_ani,
            frames=frames,
            interval=50, 
            repeat=False, 
            blit=True,
        )

        ani.save(savepath, dpi=80, writer='imagemagick')
        

def main():
    q0 = np.deg2rad(180.0)
    r0 = np.deg2rad(180.0)

    pendulum = DoublePendulum()
    simulator = RungeKuttaSimulator(pendulum, q0=q0, r0=r0, T=25)
    history = simulator.simulate()

    title = 'Double Pendulum'
    animator = Animator(history, title, incr=1)
    animator.animate()


if __name__ == '__main__':
    main()
