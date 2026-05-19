import numpy as np

def P2R(abs: float, angle: float):
    angle = np.deg2rad(angle)  # Converter angulo de graus para radianos
    return abs*np.cos(angle) + 1j*abs*np.sin(angle)

x = P2R(np.inf, 0)
print(x)