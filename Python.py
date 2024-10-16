import numpy as np
import matplotlib.pyplot as plt

# Function to classify and simulate the orbit
def orbit_simulation():
    # Input for current radial distance, radial velocity, and tangential velocity
    r_current = float(input('Enter the current radial distance from the center of Earth in km: ')) * 1e3  # Convert km to meters
    v_radial = float(input('Enter the radial velocity (velocity towards/away from Earth) in km/s: '))* 1e3 
    v_tangential = float(input('Enter the tangential velocity (velocity perpendicular to the radial direction) in km/s: '))* 1e3 

    # Earth's parameters
    mu_earth = 3.986e14  # Standard gravitational parameter of Earth (m^3/s^2)
    radius_earth = 6371e3  # Radius of the Earth in meters

    # Orbital energy equation to calculate semi-major axis (a)
    v_total = np.sqrt(v_radial**2 + v_tangential**2)  # Total velocity (magnitude)
    specific_orbital_energy = (v_total**2 / 2) - (mu_earth / r_current)
    semi_major_axis = -mu_earth / (2 * specific_orbital_energy)


    # Eccentricity calculation
    eccentricity = np.sqrt(1 + (2 * specific_orbital_energy * (r_current**2 * v_tangential**2)) / mu_earth**2)

    # Calculate perigee and apogee distances for elliptical orbits
    if eccentricity < 1:
        perigee = semi_major_axis * (1 - eccentricity)  # Closest point to Earth
        apogee = semi_major_axis * (1 + eccentricity)   # Farthest point from Earth
    else:
        perigee = None  # Perigee doesn't exist for hyperbolic orbits
        apogee = None  # Apogee doesn't exist for hyperbolic orbits

    # Velocities at perigee and apogee (m/s)
    if eccentricity < 1:
        v_perigee = np.sqrt(mu_earth * (2 / perigee - 1 / semi_major_axis))
        v_apogee = np.sqrt(mu_earth * (2 / apogee - 1 / semi_major_axis))
    else:
        v_perigee = None
        v_apogee = None

    # Orbital period for elliptical orbits
    if eccentricity < 1:
        orbital_period = 2 * np.pi * np.sqrt(semi_major_axis**3 / mu_earth)
    else:
        orbital_period = None  # No orbital period for hyperbolic orbits

    # Display calculated parameters
    print(f'Semi-major axis: {semi_major_axis / 1e3:.2f} km')
    print(f'Eccentricity: {eccentricity:.4f}')

    # Orbit shape classification based on eccentricity
    if eccentricity == 0:
        print('Orbit shape: Circle')
    elif eccentricity > 0 and eccentricity < 1:
        print('Orbit shape: Ellipse')
    elif eccentricity == 1:
        print('Orbit shape: Parabola')
    elif eccentricity > 1:
        print('Orbit shape: Hyperbola')

    if eccentricity < 1:
        print(f'Perigee: {perigee / 1e3:.2f} km')
        print(f'Apogee: {apogee / 1e3:.2f} km')
        print(f'Maximum velocity at perigee: {v_perigee:.2f} m/s')
        print(f'Minimum velocity at apogee: {v_apogee:.2f} m/s')
        print(f'Orbital period: {orbital_period:.2f} seconds')
    else:
        print('Perigee, apogee, and orbital period do not apply to hyperbolic orbits.')

    # Simulation and plotting for elliptical, parabolic, and hyperbolic orbits
    theta_steps = 1000  # Number of steps for the angular component

    if eccentricity < 1:  # Elliptical orbit
        theta = np.linspace(0, 2 * np.pi, theta_steps)  # True anomaly (angle around the orbit)
        r = (semi_major_axis * (1 - eccentricity**2)) / (1 + eccentricity * np.cos(theta))  # Orbit equation

        x = r * np.cos(theta)  # x-coordinates of orbit
        y = r * np.sin(theta)  # y-coordinates of orbit

        # Plot the orbit
        plt.figure()
        plt.plot(x / 1e3, y / 1e3, 'b-', linewidth=1.5)  # Plot in km
        plt.plot(0, 0, 'ro', markerfacecolor='r')  # Earth's center
        plt.title('Satellite Orbit (Elliptical)')
        plt.xlabel('x (km)')
        plt.ylabel('y (km)')
        plt.grid(True)
        plt.axis('equal')

    elif eccentricity == 1:  # Parabolic orbit
        theta = np.linspace(-np.pi/2, np.pi/2, theta_steps)  # Restricted angle for open trajectory
        r = (semi_major_axis * (1 - eccentricity**2)) / (1 + eccentricity * np.cos(theta))  # Orbit equation

        x = r * np.cos(theta)  # x-coordinates of orbit
        y = r * np.sin(theta)  # y-coordinates of orbit

        # Plot the orbit
        plt.figure()
        plt.plot(x / 1e3, y / 1e3, 'g-', linewidth=1.5)  # Plot in km
        plt.plot(0, 0, 'ro', markerfacecolor='r')  # Earth's center
        plt.title('Satellite Orbit (Parabolic)')
        plt.xlabel('x (km)')
        plt.ylabel('y (km)')
        plt.grid(True)
        plt.axis('equal')

    elif eccentricity > 1:  # Hyperbolic orbit
        theta = np.linspace(-np.pi/2, np.pi/2, theta_steps)  # Restricted angle for open trajectory
        r = (semi_major_axis * (1 - eccentricity**2)) / (1 + eccentricity * np.cos(theta))  # Orbit equation

        x = r * np.cos(theta)  # x-coordinates of orbit
        y = r * np.sin(theta)  # y-coordinates of orbit

        # Plot the orbit
        plt.figure()
        plt.plot(x / 1e3, y / 1e3, 'r-', linewidth=1.5)  # Plot in km
        plt.plot(0, 0, 'ro', markerfacecolor='r')  # Earth's center
        plt.title('Satellite Orbit (Hyperbolic)')
        plt.xlabel('x (km)')
        plt.ylabel('y (km)')
        plt.grid(True)
        plt.axis('equal')

    # Show all plots
    plt.show()

# Run the simulation
orbit_simulation()
