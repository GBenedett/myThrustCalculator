from math import sqrt, pi, degrees, sin, cos
import numpy as np
import matplotlib.pyplot as plt

# read data
def read_data_from_txt(file_name):
    with open(file_name, "r") as f:
        lines = f.readlines()

    settings = {}
    for line in lines:
        line = line.strip()
        if line.startswith("#"):
            continue
        elif "=" in line:
            key, value = line.split("=")
            settings[key.strip()] = float(value.strip())

    free_stream_velocity = settings["free_stream_velocity"]
    rotational_velocity = settings["rotational_velocity"]
    radius = settings["radius"]
    blades_number = settings["blades_number"]

    data = np.genfromtxt(file_name, skip_header=5, usecols=[0, 1, 2])
    r_R, c_R, beta = data.T
    return (
        free_stream_velocity,
        rotational_velocity,
        radius,
        blades_number,
        r_R,
        c_R,
        beta,
    )


(
    free_stream_velocity,
    rotational_velocity,
    radius,
    blades_number,
    r_R,
    c_R,
    beta,
) = read_data_from_txt("propeller1.txt")

density = 1

# preliminary calculation
chord = c_R * radius

r_hub = r_R[0] * radius
advanced_ratio = free_stream_velocity / (rotational_velocity * (2 * radius))

print(f"hub radius= {r_hub}")
print(f"advanced ratio= {advanced_ratio}")

MAXITER = 1000
a = 0.1
b = 0.01

a_vect = []
b_vect = []

for _ in range(1, MAXITER + 1):

    axial_velocity = free_stream_velocity * (1 + a)
    tangential_velocity = 2 * pi * r_R * rotational_velocity * (1 - b)

    relative_velocity = np.sqrt(axial_velocity**2 + tangential_velocity**2)

    phi = np.arctan2(axial_velocity, tangential_velocity)

    alpha = beta - np.degrees(phi)

    Cl = 2 * pi * np.radians(alpha)

    Cd = 0.008 - 0.003 * Cl + 0.01 * Cl**2

    # calculation blade element theory

    thrust_distribution = (
        0.5
        * density
        * relative_velocity**2
        * blades_number
        * chord
        * ((Cl * np.cos(phi)) - (Cd * np.sin(phi)))
        * (r_R * radius)
    )

    thrust_coefficient_distribution = (thrust_distribution) / (
        density * rotational_velocity**2 * (radius * 2) ** 4
    )
    total_thrust_coefficient = (r_R[1] - r_R[0]) * np.sum(
        thrust_coefficient_distribution
    )
    total_thrust = (r_R[1] - r_R[0]) * np.sum(thrust_distribution)

    torque_distribution = (
        0.5
        * density
        * relative_velocity**2
        * blades_number
        * chord
        * ((Cl * np.cos(phi)) + (Cd * np.sin(phi)))
        * (r_R * radius) ** 2
    )

    total_torque = (r_R[1] - r_R[0]) * np.sum(torque_distribution)

    axial_momentum = thrust_distribution / (
        4 * pi * (r_R * radius) * density * free_stream_velocity**2 * (1 + a)
    )
    angular_momentum = torque_distribution / (
        4
        * pi
        * (r_R * radius) ** 3
        * density
        * free_stream_velocity
        * (1 + a)
        * rotational_velocity
    )
    anew = 0.5 * (a + axial_momentum)
    bnew = 0.5 * (b + angular_momentum)

    if (np.abs(anew - a) < 1e-5).all() and (np.abs(bnew - b) < 1e-5).all():
        break

    a = anew
    b = bnew
a_vect.append(a)
b_vect.append(b)

print(f"total thrust= {total_thrust}")
print(f"total thrust coefficient= {total_thrust_coefficient}")

# Plot

plt.plot(r_R, thrust_distribution)
plt.xlabel("Adimensional Radius")
plt.ylabel("Local Thrust")
plt.title("Local Thrust vs Adimensional Radius")
plt.show()

plt.plot(r_R, thrust_coefficient_distribution)
plt.xlabel("Adimensional Radius")
plt.ylabel("Ct")
plt.title("Ct vs Adimensional Radius")
plt.show()

plt.plot(
    (rotational_velocity * r_R * radius / free_stream_velocity),
    np.transpose(a_vect),
    label="a",
)
plt.plot(
    (rotational_velocity * r_R * radius / free_stream_velocity),
    np.transpose(b_vect),
    label="b",
)
plt.xlabel("Chi")
plt.ylabel("a")
plt.title("a vs Chi")
plt.legend()
plt.show()
