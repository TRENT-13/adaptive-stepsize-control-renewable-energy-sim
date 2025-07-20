import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

from AP8.method import AdaptiveStepSizeControl


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
