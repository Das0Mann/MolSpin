import numpy as np

def generate_reconstructed_hyperfine_fluctuations(initial_tensor, total_time, time_step, damping=0.1, temperature=300,
                                                   restoring_coeff=0.5, output_file="output.mst"):
    """
    Generates an MST file simulating symmetric thermal fluctuations of hyperfine tensors by updating the independent components.

    Parameters:
    - initial_tensor (3x3 numpy array): Initial hyperfine tensor.
    - total_time (float): Total simulation time in nanoseconds.
    - time_step (float): Time step in nanoseconds.
    - damping (float): Langevin damping coefficient (γ).
    - temperature (float): Temperature in Kelvin.
    - restoring_coeff (float): Restoring force coefficient (k_restore).
    - output_file (str): Name of the output MST file.
    """
    k_B = 1.380649e-23  # Boltzmann constant (J/K)
    D = k_B * temperature / damping  # Diffusion constant

    # Extract the six independent components
    A_xx, A_yy, A_zz = initial_tensor[0, 0], initial_tensor[1, 1], initial_tensor[2, 2]
    A_xy, A_xz, A_yz = initial_tensor[0, 1], initial_tensor[0, 2], initial_tensor[1, 2]

    # Time steps
    num_steps = int(total_time / time_step)
    times = np.arange(0, total_time, time_step)

    with open(output_file, "w") as f:
        for t in times:
            # Generate Gaussian noise
            noise_xx = np.random.normal(0, np.sqrt(2 * D * time_step))
            noise_yy = np.random.normal(0, np.sqrt(2 * D * time_step))
            noise_zz = np.random.normal(0, np.sqrt(2 * D * time_step))
            noise_xy = np.random.normal(0, np.sqrt(2 * D * time_step))
            noise_xz = np.random.normal(0, np.sqrt(2 * D * time_step))
            noise_yz = np.random.normal(0, np.sqrt(2 * D * time_step))

            # Update the independent components
            A_xx += -restoring_coeff * (A_xx - initial_tensor[0, 0]) * time_step + noise_xx
            A_yy += -restoring_coeff * (A_yy - initial_tensor[1, 1]) * time_step + noise_yy
            A_zz += -restoring_coeff * (A_zz - initial_tensor[2, 2]) * time_step + noise_zz
            A_xy += -restoring_coeff * (A_xy - initial_tensor[0, 1]) * time_step + noise_xy
            A_xz += -restoring_coeff * (A_xz - initial_tensor[0, 2]) * time_step + noise_xz
            A_yz += -restoring_coeff * (A_yz - initial_tensor[1, 2]) * time_step + noise_yz

            # Reconstruct the symmetric tensor
            hyperfine_tensor = np.array([
                [A_xx, A_xy, A_xz],
                [A_xy, A_yy, A_yz],
                [A_xz, A_yz, A_zz]
            ])

            # Flatten tensor and write to file
            tensor_flat = hyperfine_tensor.flatten()
            f.write(f"{t:.8f} " + " ".join(f"{x:.8e}" for x in tensor_flat) + "\n")

    print(f"MST file with reconstructed symmetric fluctuations generated: {output_file}")


# Example usage
initial_tensor = np.array([
    [-0.0003, -0.00004, -0.0001],
    [-0.00004, -0.00055, 0.00004],
    [-0.0001, 0.00004, -0.00038]
])

generate_reconstructed_hyperfine_fluctuations(
    initial_tensor=initial_tensor,
    total_time=10,       # Total time in nanoseconds
    time_step=0.01,       # Time step in nanoseconds
    damping=1e-8,          # Langevin damping coefficient
    temperature=300,      # Temperature in Kelvin
    restoring_coeff=0.05,  # Restoring force coefficient
    output_file="reconstructed_hyperfine.mst"
)
