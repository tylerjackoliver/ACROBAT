import numpy as np
import scipy.integrate as scii
import spiceypy as spice
import matplotlib.pyplot as plt
import copy
import IPython.core.debugger as ipdb

MERCURY_MU = 2.20329e4
SUN_MU = 1.327124400189e11
MERCURY_RADIUS = 2439.7
P_MUS = {'Venus': 3.248599e5, 'Jupiter barycenter': 1.266865349e8, 'Saturn barycenter': 3.79311879e7}
time_scale = np.sqrt(MERCURY_RADIUS ** 3/MERCURY_MU)
h_list = []
radius_list = []

def get_sun_vector(time: float):

    output, _ = spice.spkezr("Sun", time, "J2000", "NONE", "Mercury")
    return output


def force_function(t, y):

    # Get the derivative vector
    dy = np.zeros(y.shape)
    # Populate the first elements of dy with the last elements of y
    dy[:3] = y[3:]
    # Get the force on the s/c from Mercury
    dy[3:] -= MERCURY_MU/MERCURY_MU * y[:3] / ( np.linalg.norm(y[:3]) ** 3 )
    # Get the force on the s/c d/t the Sun
    sun_vector = get_sun_vector(634798080 + t * time_scale)
    sun_vector /= MERCURY_RADIUS
    dy[3:] -= SUN_MU/MERCURY_MU * sun_vector[:3] / ( np.linalg.norm(sun_vector[:3]) ** 3 )
    # Difference vector
    diff_vec = y - sun_vector
    dy[3:] -= SUN_MU/MERCURY_MU * diff_vec[:3] / ( np.linalg.norm(diff_vec[:3]) ** 3)

    for pname in ("Venus", "Jupiter barycenter", "Saturn barycenter"):
        p_vector, _ = spice.spkezr(pname, 634798080 + t * time_scale, "J2000", "NONE", "Mercury")
        p_vector /= MERCURY_RADIUS
        diff_vec = y - p_vector
        dy[3:] -= P_MUS[pname] / MERCURY_MU * (p_vector[:3] / (np.linalg.norm(p_vector[:3]) ** 3) + diff_vec[:3] / (np.linalg.norm(diff_vec[:3]) ** 3) )

    return dy


def event(t, y):
    #kepler_energy = np.linalg.norm(y[3:] / 2439.7 * np.sqrt(2439.7**3 / MERCURY_MU)) ** 2 / 2.0 - 1./np.linalg.norm(y[:3]/2439.7)
    #radius = np.linalg.norm(y[:3]/2439.7)
    radius = np.linalg.norm(y[:3])
    kepler_energy = np.linalg.norm(y[3:]) ** 2. / 2.0 - 1./np.linalg.norm(y[:3])
    h_list.append(kepler_energy)
    radius_list.append(radius)
    return int(not (kepler_energy > 0 and radius > 0.112e6 / 2439.7) or not radius < 1)

# Declare terminal
event.terminal = False


def main():

    # Original position of the s/c
    y0 = np.array([-20547.8,-61754.8,5061.37,0.727552,-0.218763,0.284494])
    t = 0.
    y0[:3] /= MERCURY_RADIUS
    y0[3:] *= time_scale / MERCURY_RADIUS
    spice.furnsh('spice/naif0012.tls')
    spice.furnsh('spice/de430.bsp')

    output = scii.solve_ivp(force_function, (-0., -100 * 86400. / time_scale), y0, atol=1e-14, rtol=1e-014, events=event)
    output1 = scii.solve_ivp(force_function, (0., 67.83 * 86400. / time_scale), y0, atol=1e-014, rtol=1e-014, events=event)
    print(f"Final value is {output.y[:,-1]}")
    print(f"Did we get an event? 1 if yes: {output.status}")
    print(f"The final integration time was {output.t_events}")
    norms = []
    for i in range(output1.y.shape[1]):
        norms.append(np.linalg.norm(output1.y[:3, i]) * MERCURY_RADIUS)

    print(f"The minimum altitude is {min(norms)-MERCURY_RADIUS}")
    plt.plot(output.y[0,:], output.y[1, :], '-k')
    plt.plot(output.y[0,0], output.y[1, 0], '-bx')
    last_list = np.zeros( (2, 2) )
    last_list[0, :] = output1.y[:2, -1]
    last_list[1, :] = np.array([0., 0.])
    plt.plot(output1.y[:, 0], output1.y[:, 1], '-g')   
    plt.plot(output1.y[0,-1], output1.y[1, -1], 'rx')
    plt.show()
    for i in range(len(h_list)):
        h_list[i] = h_list[i] * 1000.
    plt.plot((output.t/84600. * time_scale), h_list[:output.t.shape[0]], '-.k')
    plt.plot(output1.t / 86400. * time_scale, h_list[output.t.shape[0]:], '-k')

    plt.ylim( (-18, 0) ) 
    plt.ylabel('Kepler Energy [VU^2]')
    plt.xlabel('Time since epoch [days]')
    plt.show()

    plt.plot((output.t/84600. * time_scale), radius_list[:output.t.shape[0]], '-.k')
    plt.plot(output1.t / 86400. * time_scale, radius_list[output.t.shape[0]:], '-k')
    plt.ylim( (0, 200) )
    plt.show()

    output2 = scii.solve_ivp(force_function, (67.83 * 86400 / time_scale, -100 * 86400 / time_scale), output1.y[:,-1], atol=1e-014, rtol=1e-014, events=event)
    plt.plot(output.y[0, :], output.y[1, :], '-g')
    plt.plot(output1.y[0, :], output1.y[1, :], '-r')
    plt.plot(output2.y[0, :], output2.y[1, :], '-k')
    plt.show()

    return 0


if __name__ == "__main__":
    main()
