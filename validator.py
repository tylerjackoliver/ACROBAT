import numpy as np
import scipy.integrate as scii
import spiceypy as spice
import matplotlib.pyplot as plt
import copy

HOST = 'Mercury'
P_MUS = {'Venus':3.248599e5, 'Jupiter barycenter': 1.266865349e8, 'Saturn barycenter': 3.79311879e7, 'Mercury': 2.20329e4, 'Sun': 1.327124400189e11}
P_RADII = {'Mercury': 2439.7}
TIME_CONVERSION_COEFF = np.sqrt(P_RADII[HOST]  ** 3 / P_MUS[HOST])
EPOCH = 634798080.
MAX_RADIUS = 0.112e6 / P_RADII[HOST]
MAX_TIME = 8.0 * np.pi * ( (MAX_RADIUS * P_RADII[HOST]) ** (2./5.))
PREVIOUS_CONDITION_ONE = 1.0
TRAJECTORY_PLOT = 0


""" @brief Gets the vector from the HOST to the planet.
    @param[in] time Time (ephemeris seconds) at which to get the state
    @param[in] target Planet to target in the vector
    @returns NumPy ndarray containing the planet state
"""
def get_planet_vector(time: float, target: str):
    output, _ = spice.spkezr(target, time, "J2000", "NONE", HOST)
    return output


""" @brief Computes the right-hand-side of the ODE
    @param[in] t The current integration time
    @param[in] y The current state at time t
    @returns dy The derivative of y
"""
def force_function(t: float, y, negative=False, initial_condition, last_condition):

    # Get the derivative vector
    dy = np.zeros(y.shape)
    # Populate the first elements of dy with the last elements of y
    dy[:3] = y[3:]
    # Get the force on the s/c from the Host
    dy[3:] -= y[:3] / ( np.linalg.norm(y[:3]) ** 3 )
    # Get the force on the s/c d/t the Sun
    sun_vector = get_sun_vector(634798080 + t * TIME_CONVERSION_COEFF, 'Sun')
    sun_vector /= P_RADII[HOST]
    dy[3:] -= P_MUS['Sun']/P_MUS[HOST] * sun_vector[:3] / ( np.linalg.norm(sun_vector[:3]) ** 3 )
    # Difference vector
    diff_vec = y - sun_vector
    dy[3:] -= P_MUS['Sun']/P_RADII[HOST] * diff_vec[:3] / ( np.linalg.norm(diff_vec[:3]) ** 3)
    # Get force d/t additional planets
    for pname in ("Venus", "Jupiter barycenter", "Saturn barycenter"):
        p_vector = get_planet_vector(EPOCH + t * TIME_CONVERSION_COEFF, pname)
        p_vector /= P_RADII[HOST]
        diff_vec = y - p_vector
        dy[3:] -= P_MUS[pname] / P_MUS[HOST] * (p_vector[:3] / (np.linalg.norm(p_vector[:3]) ** 3) + diff_vec[:3] / (np.linalg.norm(diff_vec[:3]) ** 3) )
    
    return dy


""" @brief Determine whether any events have occured in the integration
    @param[in] t The current integration time
    @param[in] y The current integration state
    @returns Boolean flag corresponding to an event
"""
def find_event(t, y, negative, initial_condition, last_condition):
    radius = np.linalg.norm(y[:3])
    kepler_energy = np.linalg.norm(y[3:]) ** 2. / 2.0 - 1./np.linalg.norm(y[:3])
    
    crash = radius < 1.
    escape = radius > MAX_R
    
    if crash:
        return 1
    
    if not negative and escape:
        return 2
    
    angular_momentum = np.cross(initial_condition[:3], initial_condition[3:])
    condition_one = np.dot(y[:3], np.cross(angular_momentum, initial_condition[:3]))
    sign_change = PREVIOUS_CONDITION_ONE * condition_one < 0

    condition_two = np.dot(y[:3], initial_condition[:3])
    condition_three = np.dot(y[3:], y[:3]) * np.dot(previous_condition_one[:3], y[3:])

    if (sign_change and condition_two > 0 and condition_three > 0):
        return 3

    if (np.abs(t) >= MAX_TIME):
        return 4
    
    return 0


""" @brief Event controller for the ODE
    @param[in] t Current integration time
    @param[in] y Current integration state
    @param[in] negative Are we integrating in negative time
    @param[in] initial_condition The initial condition for the integration
    @param[in] last_condition The condition at previous plane crossing
    @param[in] previous_condition_one The value of the previous condition
"""
def has_an_event_occurred(t, y, negative, initial_condition, last_condition):
    event_code = find_event(t, y, negative, initial_condition, last_condition)
    if not event_code == 0:
        return 1
    else:
        return 0


# Make this function terminal
has_an_event_occurred.terminal = False


"""
    @brief Loads the ephemeris to be used
    @returns None
"""
def load_ephemeris() -> None:
    spice.furnsh('spice/naif0012.tls')
    spice.furnsh('spice/de430.bsp')


""" @brief Unload all the ephemeris
    @returns None
"""
def unload_ephemeris() -> None:
    spice.unload('spice/naif00012.tls')
    spice.unload('spice/de430.bsp')


""" @brief Load the trajectories to be computed.
    @returns The trajectories to be computed.
"""
def load_trajectories():
    return np.genfromtxt('validation_set', dtype=float, delimiter=',')


""" @brief Find the nutation angles
    @param[in] t The time at which to find the RADEC angles
    @returns A tuple of (alpha, delta)
"""
def RADEC(t):
    alpha0 = np.deg2rad(281.0097)
    delta0 = np.deg2rad(61.4143)

    secondsToCenturies = 1./(86400. * 36525)

    timeDelta = t * secondsToCenturies

    alphaPrec = np.deg2rad(-0.0328)
    deltaNut = np.deg2rad(-.0049)

    alpha = alpha0 + alphaPrec * timeDelta
    delta = delta0 + deltaNut * timeDelta 

    return alpha, delta

""" @brief Converts from a BME frame to an EME frame
    @param[in] The state in the BME frame
    @param[in] t Time in the BME frame
    @returns Returns the state in the EME frame
"""
def bme_to_eme(y, t):

    alpha, delta = RADEC(t)

    sina = np.sin(alpha)
    sind = np.sin(delta)

    cosa = np.cos(alpha)
    cosd = np.cos(delta)

    rot = np.zeros( (3, 3) )

    rot[0, 0] = -sina rot[0,1]= -cosa * sind rot[0,2]= cosa * cosd
    rot[1, 0] = cosa  rot[1,1]= -sina * sind rot[1,2]= cosd * sina
    rot[2, 0] = 0.0   rot[2,1]= cosd         rot[2,2]= sind

    return rot * y


def integrate_and_plot_trajectory(y0):
    pass
    return None

def main():
    load_ephemeris()
    trajectories_to_integrate = load_trajectories()

    #

    unload_ephemeris()
    return 0

