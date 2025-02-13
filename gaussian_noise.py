import numpy as np
import matplotlib.pyplot as plt

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

    print("D = ", D) 
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

def plot_trj(file, T, N, ax, terms, pref=1, color="black"):

    dt = T/N
    print(dt)
    time = np.linspace(0, T, N)

    labels = ["1","2","3","4","5", "6","7","8","9"]

    if terms ==3:
        labels = ["field.x", "field.y", "field.z"]
    elif terms == 9: 
        labels = ["mat.xx", "mat.xy", "mat.xz", "mat.yx", "mat.yy", "mat.yz", "mat.zx", "mat.zy", "mat.zz"]

    stdevs= []

    for i in range(terms):
        
        # i=2
        d_file = open(file, 'r')
        next(d_file)
        sig = []
        for line in d_file:
            data = line.split()
            sig.append(np.double(data[i+1]))
        d_file.close() 
        sig=np.array(sig)*pref

        stdev = np.std(sig)
        stdevs.append(stdev)

        # fig, ax = plt.subplots(num = fig_no, figsize=(8,6))
        # plt.ylim(1.3, 2.2)
        # plt.xlim(0, 500)
        fontsize = 14
        linewidth = 2
        ax.tick_params(labelsize=fontsize * 1.1, length=10, width=linewidth)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_linewidth(linewidth)
        ax.spines['bottom'].set_linewidth(linewidth)

        ax.plot(time*1e9, sig[0:N], alpha =0.8, label=labels[i])
        ax.legend()

        ax.set_xlabel("Simulation Time / ns", fontsize=20)
        ax.set_ylabel("$T_{zz}(t)$ / mT", fontsize=20)

fig, ax = plt.subplots(num=1)
# plot_trj("Example/standard_examples/GaussianNoise.mst", 10, 5000, ax)
# plot_trj("Example/standard_examples/FieldBBNoise.mst", 10, 50, ax)
# plot_trj("Example/standard_examples/OUGeneral.mst", 5000, 5000, ax)
# plot_trj("Example/standard_examples/dist.txt", 4998, 4998, ax)
# plot_trj("Example/standard_examples/OUGeneral.mst", 500, 500, ax)
# plot_trj("zeeman_bb2.mst", 70, 700, ax, 3 )
plot_trj("hf1.mst", 5000, 5000, ax, terms=9)
fig.savefig("GaussianTensor.png", dpi=300, bbox_inches="tight")
  
# fig, ax = plt.subplots(num=2)
# plot_trj("reconstructed_hyperfine.mst", 10, 999, ax)
# fig.savefig("reconstructed_hyperfine.png", dpi=300, bbox_inches="tight")
def plot_result(file, ax, label=""):

    d_file = open(file, 'r')
    next(d_file)
    time = []
    sig = []
    for line in d_file:
        data = line.split()
        time.append(np.double(data[1]))
        sig.append(np.double(data[2]))

    ax.plot(np.array(time), np.array(sig), label=label)
    plt.legend()

fig, ax = plt.subplots(num=13)
# plot_result("dat_gaussian_test_mst.dat", ax, label="mst")
# plot_result("dat_gaussian_test.dat", ax, label="built-in")
# plot_result("example_timeevolution.dat", ax)
# plot_result("dat_BBfield_test.dat", ax, label="BB")
# plot_result("dat_BBfield_mst_test.dat", ax, label="BB mst")
# plot_result("dat_monochromatic_test.dat", ax, label="built-in")
# plot_result("dat_monochromatic_test_mst.dat", ax, label="built-in")
# plot_result("dat_tdMono_test.dat", ax, label="built-in")
# plot_result("dat_ougeneral_test.dat", ax, label="built-in")


fig.savefig("gaussian_Pst.png", dpi=300)

