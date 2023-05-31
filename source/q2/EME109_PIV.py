import numpy as np
import matplotlib.pyplot as plt

viscosity = 0.0010016           #N*s/m^2
density = 998.21                #kg/m^3

diameter_vals = np.array([0.009525, 0.009525, 0.012573, 0.012573])
velocity_vals = np.array([0.15, 0.5, 0.15, 0.5])

Re_vals = density * velocity_vals * diameter_vals / viscosity

print("Reynolds Numbers: ", np.round(Re_vals, 2))

St_vals = 0.2 * (1 - 20 / Re_vals)
print("Strouhal Numbers: ", np.round(St_vals, 4))

vortex_frequency = St_vals * velocity_vals / diameter_vals
print("Vortex Shedding Frequency: ", np.round(vortex_frequency, 3)) 


pixel_length_conv = 0.012573 / 62           #m/pixel

x_displacment_vals = pixel_length_conv * np.array([5, 11, 5.385165, 12])

time_change = 0.00625                       #s

exp_velocity_vals  = x_displacment_vals / time_change
print("Experimental Velocity (m/s): ", exp_velocity_vals)

exp_Re_vals = density * exp_velocity_vals * diameter_vals / viscosity

print("Experimental Reynolds Numbers: ", np.round(exp_Re_vals, 2))

exp_St_vals = 0.2 * (1 - 20 / exp_Re_vals)
print("Experimental Strouhal Numbers: ", np.round(exp_St_vals, 4))

exp_vortex_frequency = exp_St_vals * exp_velocity_vals / diameter_vals
print("Experimental Vortex Shedding Frequency: ", np.round(exp_vortex_frequency, 3)) 

visual_freq_vals = np.array([3, 8, 3, 5])

visual_St_vals = visual_freq_vals * diameter_vals / exp_velocity_vals
print("Visual Strouhal Numbers: ", np.round(visual_St_vals, 3))

plt.figure(figsize=(8,8))
plt.scatter(exp_Re_vals, exp_St_vals, label = "Experimental Velocity Values")
#plt.scatter(exp_Re_vals, visual_St_vals, label = "Experimental Frequency Values")
plt.scatter(Re_vals, St_vals, label = "Theoretical Data")
plt.legend(loc = 'best')
plt.xlabel("Reynolds Number")
plt.ylabel("Strouhal Number")
plt.title("Reynolds Number vs Strouhal Number")
plt.show()