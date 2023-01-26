from math import sqrt, pi, degrees, sin, cos
import numpy as np
import matplotlib.pyplot as plt

######### Function
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


def plot_data(
    r_R,
    thrust_distribution,
    thrust_coefficient_distribution,
    beta,
    alpha,
    phi,
    Cd,
    Cl,
    c_R,
):
    fig, axs = plt.subplots(2, 1, figsize=(10, 10))
    axs[0].plot(r_R, thrust_distribution)
    axs[0].set_xlabel("Adimensional Radius")
    axs[0].set_ylabel("Thrust")
    axs[0].set_title("Thrust vs Adimensional Radius")

    axs[1].plot(r_R, thrust_coefficient_distribution)
    axs[1].set_xlabel("Adimensional Radius")
    axs[1].set_ylabel("Thrust coefficient")
    axs[1].set_title("Thrust coefficient vs Adimensional Radius")
    plt.show()

    plt.plot(r_R, beta, label="beta")
    plt.plot(r_R, alpha, label="alpha")
    plt.plot(r_R, phi, label="phi")
    plt.xlabel("Adimensional Radius")
    plt.ylabel("Angles")
    plt.title("Beta vs Adimensional Radius")
    plt.legend()
    plt.show()

    plt.plot(r_R, Cd, label="$C_D$")
    plt.plot(r_R, Cl, label="$C_L$")
    plt.xlabel("Adimensional Radius")
    plt.ylabel("$C_D$ $C_L$")
    plt.title("Cd vs Adimensional Radius")
    plt.legend()
    plt.show()

    plt.plot(r_R, c_R)
    plt.xlabel("Adimensional Radius")
    plt.ylabel("Chord")
    plt.ylim((0, 0.8))
    plt.title("Chord vs Adimensional Radius")
    plt.legend()
    plt.show()


######### Main
(
    free_stream_velocity,
    rotational_velocity,
    radius,
    blades_number,
    r_R,
    c_R,
    beta,
) = read_data_from_txt("data.txt")

density = 1

# preliminary calculation
chord = c_R * radius
tangential_velocity = 2 * pi * r_R * rotational_velocity

relative_velocity = np.sqrt(free_stream_velocity**2 + tangential_velocity**2)

phi = np.arctan2(free_stream_velocity, tangential_velocity)

alpha = beta - np.degrees(phi)

Cl = 2 * pi * np.radians(alpha)

Cd = 0.008 - 0.003 * Cl + 0.01 * Cl**2

r_hub = r_R[0] * radius
advanced_ratio = free_stream_velocity / (rotational_velocity * (2 * radius))
print(f"hub radius= {r_hub}")
print(f"advanced ratio= {advanced_ratio}")

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
    density * rotational_velocity**2 * (radius * 2) ** 2
)
total_thrust_coefficient = (r_R[1] - r_R[0]) * np.sum(thrust_coefficient_distribution)
total_thrust = (r_R[1] - r_R[0]) * np.sum(thrust_distribution)

print(f"total thrust= {total_thrust}")
print(f"total thrust coefficient= {total_thrust_coefficient}")

plot_data(
    r_R,
    thrust_distribution,
    thrust_coefficient_distribution,
    beta,
    alpha,
    phi,
    Cd,
    Cl,
    c_R,
)
