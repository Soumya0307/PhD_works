How to show:
**0. Setup**

Have two .xyz files for a molecule or reaction (Initial geometry and final geometry), then create an input file (.inp) for ORCA that commands to run the NEB calculation or whatever you require to calculate. Then run the ORCA inside the same folder where the .inp and .xyz files are to get the .out, .interp and other backup files. Command prompt: "Path where ORCA is stored" "input file name".inp -> "output file name".out 

**1. To get TS geometry**



If full experiment/run is allowed, then the files with ***'H3\_rxn\_NEB-TS\_converged.xyz'*** or ***'H3\_rxn\_NEB-CI\_converged.xyz'*** will be appeared. Those appeared only after full convergence.



If not, then look at NEB iteration table in the .out file (Here H3\_rxn.out). The table looks as follows





Starting iterations:



Optim.  Iteration  HEI  E(HEI)-E(0)  max(|Fp|)   RMS(Fp)    dS

Switch-on CI threshold               0.020000 

&nbsp;  LBFGS     0      5    0.160982    0.083566   0.017014  6.0889       

&nbsp;  LBFGS     1      5    0.160894    0.078382   0.016368  6.1029       

&nbsp;  LBFGS     2      5    0.160794    0.070905   0.015412  6.1411       

&nbsp;  LBFGS     3      5    0.160685    0.064380   0.013985  6.2028       

&nbsp;  LBFGS     4      5    0.160560    0.060227   0.013068  6.2511       

&nbsp;  LBFGS     5      5    0.160424    0.053172   0.011436  6.3282       

&nbsp;  LBFGS     6      5    0.160254    0.048278   0.010704  6.3771       

&nbsp;  LBFGS     7      5    0.160078    0.044139   0.009182  6.4553       

&nbsp;  LBFGS     8      5    0.159848    0.040571   0.008647  6.5038       

&nbsp;  LBFGS     9      5    0.159605    0.036194   0.007285  6.5818       

&nbsp;  LBFGS    10      5    0.159282    0.031784   0.006971  6.6311       

&nbsp;  LBFGS    11      5    0.158917    0.027326   0.005718  6.7094       

&nbsp;  LBFGS    12      5    0.158461    0.023026   0.005382  6.7629       

&nbsp;  LBFGS    13      5    0.157862    0.019436   0.004618  6.8324       



Image  5 will be converted to a climbing image in the next iteration (max(|Fp|) < 0.0200) 



Optim.  Iteration  CI   E(CI)-E(0)   max(|Fp|)   RMS(Fp)    dS     max(|FCI|)   RMS(FCI)

Convergence thresholds               0.020000   0.010000            0.002000    0.001000 

&nbsp;  LBFGS    14      5    0.157155    0.018582   0.005093  6.8726    0.009709    0.004383       

&nbsp;  LBFGS    15      5    0.156297    0.030009   0.005883  6.9268    0.010159    0.004651       

&nbsp;  LBFGS    16      5    0.155692    0.014072   0.004138  6.9885    0.010562    0.004845       

&nbsp;  LBFGS    17      5    0.154571    0.015548   0.004414  7.0439    0.011275    0.005201       

&nbsp;  LBFGS    18      5    0.153469    0.039243   0.006919  7.0333    0.011865    0.005548       

&nbsp;  LBFGS    19      5    0.151549    0.063733   0.010267  7.0730    0.013250    0.006155       

&nbsp;  LBFGS    20      5    0.151021    0.027083   0.005786  7.1199    0.013683    0.006322       

&nbsp;  LBFGS    21      5    0.149803    0.014434   0.004992  7.1914    0.014688    0.006709       

&nbsp;  LBFGS    22      5    0.147680    0.040882   0.007778  7.1966    0.016447    0.007382       

&nbsp;  LBFGS    23      5    0.144509    0.067063   0.011446  7.2581    0.019056    0.008381       

&nbsp;  LBFGS    24      5    0.143567    0.029293   0.007234  7.3170    0.019679    0.008658       

&nbsp;  LBFGS    25      4    0.144047    0.020070   0.006913  7.4120    0.050588    0.023056       

&nbsp;  LBFGS    26      4    0.144661    0.082744   0.013902  7.5357    0.052101    0.023671       

&nbsp;  LBFGS    27      4    0.148933    0.025511   0.007526  7.5062    0.049138    0.022037       

&nbsp;  LBFGS    28      4    0.155158    0.026809   0.008729  7.5328    0.045674    0.020016  



&nbsp;

The CI image shows 5. So, then, the file ***H3\_rxn\_im5.xyz*** will be the TS geometry.





**2. Animation on Avogadro:** Drag and paste ***H3\_rxn\_MEP\_ALL\_trj***



**3. Plot the MEP:** The ***.interp*** file





import re, numpy as np, matplotlib.pyplot as plt, os



path = "/mnt/data/H3\_rxn.interp"



\# Parse numeric blocks

with open(path, "r", encoding="utf-8", errors="replace") as f:

&nbsp;   lines = f.readlines()



float\_line = re.compile(r"^\\s\*\[-+0-9.]+(?:\[eEdD]\[-+]?\\d+)?\\s+\[-+0-9.]+(?:\[eEdD]\[-+]?\\d+)?")

blocks, cur = \[], \[]

for line in lines:

&nbsp;   if float\_line.match(line):

&nbsp;       cur.append(line.replace("D","E").replace("d","E"))

&nbsp;   else:

&nbsp;       if cur:

&nbsp;           blocks.append(cur); cur=\[]

if cur:

&nbsp;   blocks.append(cur)



data = np.loadtxt(blocks\[-1])

c1, c2, c3 = data\[:,0], data\[:,1], data\[:,2]  # inferred meaning below



\# Identify peak

imax = np.argmax(c3)

peak = (c1\[imax], c2\[imax], c3\[imax])

hartree\_to\_kcal = 627.509474



print("Last block columns ranges:")

print("  col1 min/max:", float(c1.min()), float(c1.max()))

print("  col2 min/max:", float(c2.min()), float(c2.max()))

print("  col3 min/max:", float(c3.min()), float(c3.max()))

print(f"Peak at: col1={peak\[0]:.4f}, col2={peak\[1]:.4f}, col3={peak\[2]:.6f} Eh = {peak\[2]\*hartree\_to\_kcal:.2f} kcal/mol")



\# Plot 1: normalized coordinate vs relative energy

plt.figure()

plt.plot(c1, c3, marker="o")

plt.xlabel("col1 (normalized path coordinate)")

plt.ylabel("col3 (relative energy, Eh)")

plt.title("MEP from H3\_rxn.interp (last block): Energy vs normalized coordinate")

plt.show()



\# Plot 2: path length vs relative energy (often the most MEP-like)

plt.figure()

plt.plot(c2, c3, marker="o")

plt.xlabel("col2 (path length coordinate)")

plt.ylabel("col3 (relative energy, Eh)")

plt.title("MEP from H3\_rxn.interp (last block): Energy vs path length")

plt.show()









