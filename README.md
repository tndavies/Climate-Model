## Earth Climate Model

# Description

A simulation of Earth's surface temperature, using a 1D energy-balence model. Can be used to make temperature predications for various future greenhouse gas emission scenarios. e.g:


![gat_forecast](https://github.com/user-attachments/assets/8588cf63-fa68-418f-a638-8c2404ac6f46)


# Usage

- Import climate.py
- Define a simulation specification, as so
```
Sim_Specification mySpec(Duration=<# years to simulate>)
results = Simulate_Climate(mySpec)
```
You can then access the latitude grid used, temporal points, as well as the corressponding surface temperature distributions. 
(see Sim_Result dataclass in climate.py)

(plots.py demos some uses for the climate model)
