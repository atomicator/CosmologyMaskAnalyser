import numpy as np
import matplotlib.pyplot as plt
from toolkit import toolkit
import hmf
import copy

toolkit.plt_use_tex()


def load_richness_redshift(richness, redshift):
    if richness <= min_richness:
        return False
    else:
        sdss_richness.append(richness)
        sdss_redshift.append(redshift)
        return True

def test(*args):
    sdss_data.append(np.max(args))
    return True


def filter(richness):
    if richness > 10:
        sdss_data.append(richness)
        return True
    else:
        return False


def mass_to_richness(mass):
    a = 1
    return 30.0 * (mass / a) ** (3/4)

def richness_to_mass(richness):
    #return 10 ** 14.45 * (richness / 40) ** 1.29
    return (richness / 30) ** (4/3) * 3e14

def mass_to_act_y0(mass, z):
    matter_density = 0.22
    pivot_mass = 3e14 * 0.72 # assume h to be 0.72
    beta = 0.08
    return (4.95e-5) * (matter_density * (1 + z) ** 3 + (1 - matter_density)) * (mass / pivot_mass) ** (1 + beta)

bins = np.linspace(0, 10, 51)
min_richness = 0
act_cat = toolkit.StarCatalogue("../../data/DR5_cluster-catalog_v1.1.fits", hdu=1, table=True)
sdss_richness = []
sdss_redshift = []
act_cat.load_with_selection(load_richness_redshift, ("fixed_y_c", "fixed_SNR"), True)
sdss_richness = np.array(sdss_richness) * 1e-4
sdss_redshift = np.array(sdss_redshift)
act_noise = np.mean(sdss_richness / sdss_redshift)
act_snr = copy.copy(sdss_redshift)

sdss_height = np.zeros(len(bins) - 1)
for i in range(len(bins) - 1):
    sdss_height[i] = np.sum(np.float_(np.bitwise_and(act_snr > bins[i], act_snr < bins[i + 1]))) / (len(act_snr) / 100)
plt.stairs(sdss_height, bins, label="ACT")

print(act_noise)

#exit()


sdss_cat = toolkit.load_catalogue("sdss")
sdss_data = []
sdss_cat.load_with_selection(test, ("LAMBDA_CHISQ",), True)
sdss_data = np.array(sdss_data)
sdss_height = np.zeros(len(bins) - 1)
for i in range(len(bins) - 1):
    sdss_height[i] = np.sum(np.float_(np.bitwise_and(sdss_data > bins[i], sdss_data < bins[i + 1]))) / (len(sdss_data) / 100)
#plt.stairs(sdss_height, bins, label="SDSS")

#sdss_cat = toolkit.StarCatalogue("../../data/DR5_cluster-catalog_v1.1.fits", hdu=1, table=True)
#sdss_data = []
#sdss_cat.load_with_selection(test, ("RM_LAMBDA", "RMDESY3_LAMBDA_CHISQ", "CAMIRA_N_mem"), True)
#sdss_data = np.array(sdss_data)
#sdss_height = np.zeros(len(bins) - 1)
#for i in range(len(bins) - 1):
#    sdss_height[i] = np.sum(np.float_(np.bitwise_and(sdss_data > bins[i], sdss_data < bins[i + 1]))) / (len(sdss_data) / 100)
#plt.stairs(sdss_height, bins, label="ACT")
#print(np.sum(sdss_height))
#print(np.min(sdss_data), np.max(sdss_data))
#print(np.sum(sdss_data < 0))

min_richness = 0
sdss_cat = toolkit.load_catalogue("sdss")
sdss_richness = []
sdss_redshift = []
sdss_cat.load_with_selection(load_richness_redshift, ("LAMBDA_CHISQ", "Z"), True)
sdss_redshift = np.array(sdss_redshift)
sdss_richness = np.array(sdss_richness)
sdss_mass_200_m = richness_to_mass(sdss_richness)
data_test = hmf.halos.mass_definitions.MassDefinition()
mdef_1 = hmf.halos.mass_definitions.SOMean(overdensity=200)
mdef_2 = hmf.halos.mass_definitions.SOCritical(overdensity=500)
new_mass = []
for i in range(len(sdss_richness)):
    mnew, rnew, cnew = mdef_1.change_definition(m=sdss_mass_200_m[i], z=sdss_redshift[i], mdef=mdef_2)
    new_mass.append(mnew)
    print(i)
print(sdss_mass_200_m)
new_mass = np.array(new_mass)
y_0 = mass_to_act_y0(new_mass, sdss_redshift)
snr = y_0 / act_noise

print(snr)
print(np.min(snr), np.max(snr))

sdss_height = np.zeros(len(bins) - 1)
for i in range(len(bins) - 1):
    sdss_height[i] = np.sum(np.float_(np.bitwise_and(snr > bins[i], snr < bins[i + 1]))) / (len(snr) / 100)
plt.stairs(sdss_height, bins, label="SDSS, complete")

min_richness = 20
sdss_cat = toolkit.load_catalogue("sdss")
sdss_richness = []
sdss_redshift = []
sdss_cat.load_with_selection(load_richness_redshift, ("LAMBDA_CHISQ", "Z"), True)
sdss_redshift = np.array(sdss_redshift)
sdss_richness = np.array(sdss_richness)
sdss_mass_200_m = richness_to_mass(sdss_richness)
data_test = hmf.halos.mass_definitions.MassDefinition()
mdef_1 = hmf.halos.mass_definitions.SOMean(overdensity=200)
mdef_2 = hmf.halos.mass_definitions.SOCritical(overdensity=500)
new_mass = []
for i in range(len(sdss_richness)):
    mnew, rnew, cnew = mdef_1.change_definition(m=sdss_mass_200_m[i], z=sdss_redshift[i], mdef=mdef_2)
    new_mass.append(mnew)
    print(i)
print(sdss_mass_200_m)
new_mass = np.array(new_mass)
y_0 = mass_to_act_y0(new_mass, sdss_redshift)
snr = y_0 / act_noise

print(snr)
print(np.min(snr), np.max(snr))

sdss_height = np.zeros(len(bins) - 1)
for i in range(len(bins) - 1):
    sdss_height[i] = np.sum(np.float_(np.bitwise_and(snr > bins[i], snr < bins[i + 1]))) / (len(snr) / 100)
plt.stairs(sdss_height, bins, label=r"SDSS, $\lambda > 20$")

act_cat = toolkit.StarCatalogue("../../data/DR5_cluster-catalog_v1.1.fits", hdu=1, table=True)


plt.legend()
plt.title("The cluster count of the catalogues, as a function of ACT SNR")
#plt.yscale("log")
plt.ylabel("Percentage")
plt.xlabel("SNR")

"""
act_cat = toolkit.StarCatalogue("../../data/DR5_cluster-catalog_v1.1.fits", hdu=1, table=True)
sdss_richness = []
sdss_snr = []
sdss_redshift = []

def load_act_richness(flag, richness, snr, redshift):
    if flag == 0.0:
        return False
    else:
        sdss_richness.append(richness)
        sdss_snr.append(snr)
        sdss_redshift.append(redshift)


act_cat.load_with_selection(load_act_richness, ("RM", "RM_LAMBDA", "fixed_SNR", "redshift"), True)
sdss_snr = np.array(sdss_snr)
sdss_richness = np.array(sdss_richness)
sdss_redshift = np.array(sdss_redshift)
print(sdss_richness, sdss_snr)
sdss_mass_200_m = richness_to_mass(sdss_richness)
data_test = hmf.halos.mass_definitions.MassDefinition()
mdef_1 = hmf.halos.mass_definitions.SOMean(overdensity=200)
mdef_2 = hmf.halos.mass_definitions.SOMean(overdensity=500)
new_mass = []
for i in range(len(sdss_richness)):
    mnew, rnew, cnew = mdef_1.change_definition(m=sdss_mass_200_m[i], z=sdss_redshift[i], mdef=mdef_2)
    new_mass.append(mnew)
    print(i)
print(sdss_mass_200_m)
new_mass = np.array(new_mass)
y_0 = mass_to_act_y0(new_mass, sdss_redshift)
snr = y_0 / act_noise
plt.scatter(sdss_snr, sdss_richness, marker="+")
#plt.plot((0, 40), (0, 40), color="k")
#plt.xlim(0, 40)
#plt.ylim(0, 40)
plt.title("Matching clusters with both SNR and richness measurements")
plt.xlabel("ACT measured SNR")
#plt.ylabel("SDSS richness, mapped to ACT SNR")
plt.ylabel("SDSS richness")"""
plt.savefig("ACT_SNR.png", dpi=1000)
plt.show()
