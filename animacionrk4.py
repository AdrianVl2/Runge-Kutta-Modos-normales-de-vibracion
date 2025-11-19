import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Rectangle


# ---------------------------------------------------------
#   RK4 DIRECTO PARA DOS MASAS-RESORTES ACOPLADAS (2° orden)
# ---------------------------------------------------------
def rk4_segundo_orden_masas(k1, k2, k3, m1, m2,
                            t0, tf, h,
                            x10, v10, x20, v20,
                            mass_width=0.1):
    t = np.arange(t0, tf + h, h)
    N = len(t)

    x1 = np.zeros(N)
    v1 = np.zeros(N)
    x2 = np.zeros(N)
    v2 = np.zeros(N)

    x1[0] = x10
    v1[0] = v10
    x2[0] = x20
    v2[0] = v20

    # Aceleraciones
    def a1(tt, X1, V1, X2, V2):
        return (-(k1 + k2) * X1 + k2 * X2) / m1

    def a2(tt, X1, V1, X2, V2):
        return (k2 * X1 - (k2 + k3) * X2) / m2

    for n in range(N - 1):
        tn = t[n]

        # --- RK4 ---
        k11 = h * a1(tn, x1[n], v1[n], x2[n], v2[n])
        k21 = h * a2(tn, x1[n], v1[n], x2[n], v2[n])

        x1m = x1[n] + (h / 2) * v1[n] + (h ** 2 / 8) * k11
        v1m = v1[n] + k11 / 2
        x2m = x2[n] + (h / 2) * v2[n] + (h ** 2 / 8) * k21
        v2m = v2[n] + k21 / 2

        k12 = h * a1(tn + h / 2, x1m, v1m, x2m, v2m)
        k22 = h * a2(tn + h / 2, x1m, v1m, x2m, v2m)

        x1m = x1[n] + (h / 2) * (v1[n] + k11 / 2) + (h ** 2 / 8) * k12
        v1m = v1[n] + k12 / 2
        x2m = x2[n] + (h / 2) * (v2[n] + k21 / 2) + (h ** 2 / 8) * k22
        v2m = v2[n] + k22 / 2

        k13 = h * a1(tn + h / 2, x1m, v1m, x2m, v2m)
        k23 = h * a2(tn + h / 2, x1m, v1m, x2m, v2m)

        x1e = x1[n] + h * (v1[n] + k12 / 2) + (h ** 2 / 2) * k13
        v1e = v1[n] + k13
        x2e = x2[n] + h * (v2[n] + k22 / 2) + (h ** 2 / 2) * k23
        v2e = v2[n] + k23

        k14 = h * a1(tn + h, x1e, v1e, x2e, v2e)
        k24 = h * a2(tn + h, x1e, v1e, x2e, v2e)

        # actualizar velocidades
        v1[n + 1] = v1[n] + (k11 + 2 * k12 + 2 * k13 + k14) / 6
        v2[n + 1] = v2[n] + (k21 + 2 * k22 + 2 * k23 + k24) / 6

        # actualizar posiciones
        x1[n + 1] = x1[n] + h * (v1[n] + (h / 6) * (k11 + k12 + k13))
        x2[n + 1] = x2[n] + h * (v2[n] + (h / 6) * (k21 + k22 + k23))

    return t, x1, v1, x2, v2


# ---------------------------------------------------------
#        PARÁMETROS DEL SISTEMA
# ---------------------------------------------------------
k1 = 10;
k2 = 10;
k3 = 10;
m1 = 1;
m2 = 1
mass_width = 0.1

# RK4 numérico
t, x1, v1, x2, v2 = rk4_segundo_orden_masas(
    k1, k2, k3, m1, m2,
    0, 20, 0.0001,
    1, 0, 0, 0,  # Condiciones iniciales: x1_0=-1, v1_0=0, x2_0=0, v2_0=0
    mass_width
)

# ---------------------------------------------------------
#                ANIMACIÓN Y GRÁFICOS EN TIEMPO REAL
# ---------------------------------------------------------
skip = 100  # Reducimos un poco el skip para que los gráficos se actualicen más rápido
t_anim = t[::skip]
x1_anim = x1[::skip]
x2_anim = x2[::skip]

# Configuración de la figura y subplots
fig, (ax_anim, ax_plot) = plt.subplots(2, 1, figsize=(10, 7), gridspec_kw={'height_ratios': [1, 1]})

# --- Subplot de Animación ---
ax_anim.set_xlim(-2, 2)
ax_anim.set_ylim(-0.5, 0.5)
ax_anim.set_xlabel("Posición (m)")
ax_anim.set_title("Animación de masas acopladas con resortes")
ax_anim.set_yticks([])  # Eliminar ticks del eje Y para una vista más limpia

mass_height = 0.1
mass1_rect = Rectangle((0, -mass_height / 2), mass_width, mass_height, color='blue', label='Masa 1')
mass2_rect = Rectangle((0, -mass_height / 2), mass_width, mass_height, color='red', label='Masa 2')
ax_anim.add_patch(mass1_rect)
ax_anim.add_patch(mass2_rect)

spring1, = ax_anim.plot([], [], 'g', lw=2)
spring2, = ax_anim.plot([], [], 'g', lw=2)
spring3, = ax_anim.plot([], [], 'g', lw=2)

ax_anim.plot([-2, -2], [-0.3, 0.3], 'k', lw=3)  # Pared izquierda
ax_anim.plot([2, 2], [-0.3, 0.3], 'k', lw=3)  # Pared derecha

# --- Subplot de Gráficas en Tiempo Real ---
ax_plot.set_xlim(t_anim[0], t_anim[-1])
ax_plot.set_ylim(min(min(x1_anim), min(x2_anim)) * 1.1, max(max(x1_anim), max(x2_anim)) * 1.1)
ax_plot.set_xlabel("Tiempo (s)")
ax_plot.set_ylabel("Desplazamiento (m)")
ax_plot.set_title("Desplazamiento de las Masas en Tiempo Real")
ax_plot.grid(True)

line_x1, = ax_plot.plot([], [], 'b-', label='x1(t) (RK4)')
line_x2, = ax_plot.plot([], [], 'r-', label='x2(t) (RK4)')
ax_plot.legend()

# Líneas verticales para mostrar el tiempo actual en las gráficas
time_line_x1 = ax_plot.axvline(x=t_anim[0], color='blue', linestyle='--', lw=1, alpha=0.7)
time_line_x2 = ax_plot.axvline(x=t_anim[0], color='red', linestyle='--', lw=1, alpha=0.7)


def spring_line(x_start, x_end, y=0, turns=7, amplitude=0.08):  # Aumentar vueltas y amplitud para mejor visibilidad
    X = np.linspace(x_start, x_end, 2 * turns + 1)
    Y = np.zeros_like(X)
    Y[1::2] = amplitude
    Y[2::2] = -amplitude
    return X, Y + y  # Añadir un desplazamiento vertical 'y' para que el resorte esté centrado


def init():
    # Animación
    spring1.set_data([], [])
    spring2.set_data([], [])
    spring3.set_data([], [])
    mass1_rect.set_x(x1_anim[0] - mass_width / 2)
    mass2_rect.set_x(x2_anim[0] - mass_width / 2)

    # Gráficos en tiempo real
    line_x1.set_data([], [])
    line_x2.set_data([], [])
    time_line_x1.set_xdata(t_anim[0])
    time_line_x2.set_xdata(t_anim[0])

    return mass1_rect, mass2_rect, spring1, spring2, spring3, line_x1, line_x2, time_line_x1, time_line_x2


def update(frame):
    # --- Actualizar Animación ---
    current_x1 = x1_anim[frame]
    current_x2 = x2_anim[frame]
    current_t = t_anim[frame]

    mass1_rect.set_x(current_x1 - mass_width / 2)
    mass2_rect.set_x(current_x2 - mass_width / 2)

    x_s1, y_s1 = spring_line(-2, current_x1)
    x_s2, y_s2 = spring_line(current_x1 + mass_width / 2,
                             current_x2 - mass_width / 2)  # Resortes entre bordes de las masas
    x_s3, y_s3 = spring_line(current_x2 + mass_width / 2, 2)

    spring1.set_data(x_s1, y_s1)
    spring2.set_data(x_s2, y_s2)
    spring3.set_data(x_s3, y_s3)

    # --- Actualizar Gráficos en Tiempo Real ---
    line_x1.set_data(t_anim[:frame], x1_anim[:frame])
    line_x2.set_data(t_anim[:frame], x2_anim[:frame])

    time_line_x1.set_xdata(current_t)
    time_line_x2.set_xdata(current_t)

    # La función update debe devolver una tupla de todos los objetos que han sido modificados.
    return mass1_rect, mass2_rect, spring1, spring2, spring3, line_x1, line_x2, time_line_x1, time_line_x2


ani = FuncAnimation(fig, update, frames=len(t_anim),
                    init_func=init, blit=True, interval=20)

plt.tight_layout()  # Ajustar el diseño para evitar solapamientos
plt.show()


