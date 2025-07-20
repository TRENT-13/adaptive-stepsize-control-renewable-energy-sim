import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec


class RenewableEnergySystemAdaptive:
    def __init__(self):
        self.params = {
            'max_solar_capacity': 150.0,
            'max_wind_capacity': 120.0,
            'battery_capacity': 300.0,
            'grid_connection_limit': 100.0,
            'solar_efficiency': 0.20,
            'wind_efficiency': 0.38,
            'battery_charge_efficiency': 0.95,
            'battery_discharge_efficiency': 0.92,
            'temperature_sensitivity': 0.03,
            'wind_variability': 0.15,
        }

        self.generation_models = {
            'solar': lambda t, temp: (np.sin(2 * np.pi * t / 24) * 0.5 + 0.5) *
                                     (1 + self.params['temperature_sensitivity'] * temp),
            'wind': lambda t, wind_speed: (np.cos(2 * np.pi * t / 24) * 0.5 + 0.5) *
                                          (1 + self.params['wind_variability'] * wind_speed)
        }

        self.demand_model = {
            'residential': lambda t: 50 + 20 * np.sin(2 * np.pi * t / 24),
            'industrial': lambda t: 80 + 30 * np.cos(2 * np.pi * t / 24),
            'commercial': lambda t: 40 + 10 * np.sin(2 * np.pi * t / 12)
        }

    def system_dynamics(self, X, t, environmental_conditions):
        """System dynamics function"""
        solar_gen, wind_gen, battery_stored, grid_draw = X

        temp = environmental_conditions.get('temperature', 20)
        wind_speed = environmental_conditions.get('wind_speed', 5)

        solar_input = (self.params['max_solar_capacity'] *
                       self.params['solar_efficiency'] *
                       self.generation_models['solar'](t, temp))

        wind_input = (self.params['max_wind_capacity'] *
                      self.params['wind_efficiency'] *
                      self.generation_models['wind'](t, wind_speed))

        total_demand = (
                self.demand_model['residential'](t) +
                self.demand_model['industrial'](t) +
                self.demand_model['commercial'](t)
        )

        derivatives = np.array([
            solar_input * (1 - solar_gen / self.params['max_solar_capacity']) - solar_gen,
            wind_input * (1 - wind_gen / self.params['max_wind_capacity']) - wind_gen,
            (solar_gen + wind_gen - total_demand) * (
                0.9 if battery_stored < self.params['battery_capacity'] else 0
            ),
            min(max(total_demand - (solar_gen + wind_gen + battery_stored), 0),
                self.params['grid_connection_limit'])
        ], dtype=np.float64)

        return derivatives

    def simulate_with_adaptive_stepsize(self, T=24, initial_dt=0.1, tolerance=1e-6):
        """
        Simulate system using adaptive step size control
        """
        X0 = np.array([10.0, 10.0, 50.0, 0.0], dtype=np.float64)  # Initial state
        env_conditions = {'temperature': 20.0, 'wind_speed': 5.0}

        def system_wrapper(x, t):
            return self.system_dynamics(x, t, env_conditions)

        # Initialize storage for results
        t_points = [0.0]
        X_points = [X0]
        h_points = []  # Store step sizes for analysis

        t = 0.0
        X = X0
        h = initial_dt

        while t < T:
            # Prevent overshooting the end time
            if t + h > T:
                h = T - t

            # Take step with error estimate
            X_new, error = AdaptiveStepSizeControl.embedded_rk_step(
                system_wrapper, t, X, h
            )

            # Compute new step size
            h_new = AdaptiveStepSizeControl.step_size_control(error, h, order=4,
                                                              tolerance=tolerance)

            # Accept or reject step based on error
            if error <= tolerance:
                t += h
                X = X_new
                t_points.append(t)
                X_points.append(X)
                h_points.append(h)
                h = h_new
            else:
                h = h_new  # Reject step and try again with smaller step size

        return np.array(t_points), np.array(X_points), np.array(h_points)

    def visualize_adaptive_results(self, t, X, h):
        """
        Create visualization including step size adaptation
        """
        fig = plt.figure(figsize=(15, 15))
        gs = GridSpec(5, 1, figure=fig)

        variables = ['Solar Generation', 'Wind Generation',
                     'Battery Storage', 'Grid Draw']

        # Plot system variables
        for idx, var in enumerate(variables):
            ax = fig.add_subplot(gs[idx])
            ax.plot(t, X[:, idx], '-o', markersize=3, label=var)
            ax.set_title(f'{var} Evolution')
            ax.set_xlabel('Time (hours)')
            ax.set_ylabel('Value')
            ax.grid(True)

        # Plot step sizes
        ax_h = fig.add_subplot(gs[4])
        ax_h.semilogy(t[:-1], h, '-o', markersize=3, label='Step Size')
        ax_h.set_title('Adaptive Step Size Evolution')
        ax_h.set_xlabel('Time (hours)')
        ax_h.set_ylabel('Step Size (log scale)')
        ax_h.grid(True)

        plt.tight_layout()
        return fig

class AdaptiveStepSizeControl:
    @staticmethod
    def step_size_control(error, h, order, tolerance=1e-2):
        """
        Compute new step size based on error estimate
        Args:
            error: estimated error
            h: current step size
            order: order of the method
            tolerance: error tolerance
        """
        safety_factor = 0.9
        max_factor = 2.0
        min_factor = 0.1

        if error == 0:
            factor = max_factor
        else:
            factor = safety_factor * (tolerance / error) ** (1.0 / (order + 1))
            factor = min(max_factor, max(min_factor, factor))

        return h * factor

    @staticmethod
    def embedded_rk_step(f, t, y, h):
        """
        Embedded Runge-Kutta method (Fehlberg 4(5))
        """
        y = np.array(y, dtype=np.float64)

        a = np.array([
            [0, 0, 0, 0, 0],
            [1 / 4, 0, 0, 0, 0],
            [3 / 32, 9 / 32, 0, 0, 0],
            [1932 / 2197, -7200 / 2197, 7296 / 2197, 0, 0],
            [439 / 216, -8, 3680 / 513, -845 / 4104, 0]
        ], dtype=np.float64)

        b4 = np.array([25 / 216, 0, 1408 / 2565, 2197 / 4104, -1 / 5], dtype=np.float64)  # 4th order
        b5 = np.array([16 / 135, 0, 6656 / 12825, 28561 / 56430, -9 / 50, 2 / 55], dtype=np.float64)  # 5th order
        c = np.array([0, 1 / 4, 3 / 8, 12 / 13, 1], dtype=np.float64)

        # Compute stages
        k = np.zeros((6, len(y)), dtype=np.float64)
        k[0] = np.array(f(y, t), dtype=np.float64)

        for i in range(1, 5):
            yi = y.copy()
            for j in range(i):
                yi += h * a[i][j] * k[j]
            k[i] = np.array(f(yi, t + c[i] * h), dtype=np.float64)

        y4 = y + h * np.sum([b4[i] * k[i] for i in range(5)], axis=0)
        y5 = y + h * np.sum([b5[i] * k[i] for i in range(6)], axis=0)

        error = np.max(np.abs(y5 - y4))

        return y5, error



def run_adaptive_simulation():
    # Create system instance
    system = RenewableEnergySystemAdaptive()

    print("Running simulation with adaptive step size control...")
    t, X, h = system.simulate_with_adaptive_stepsize(T=24, initial_dt=0.1, tolerance=1e-2)

    print("\nGenerating visualizations...")
    system.visualize_adaptive_results(t, X, h)
    plt.show()

    print("\nSimulation Statistics:")
    print(f"Total time steps: {len(t)}")
    print(f"Average step size: {np.mean(h):.6f}")
    print(f"Minimum step size: {np.min(h):.6f}")
    print(f"Maximum step size: {np.max(h):.6f}")


if __name__ == "__main__":
    run_adaptive_simulation()