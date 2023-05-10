import scipy.integrate as integrate
import scipy.signal as signal
import numpy as np
import matplotlib.pyplot as plt

# Inputs (Angles are in decimal degrees!)
# Graph Parameters (seconds)
start = 0
stop = 1000000
step = 1000
UScalar = 5     # Divides distances/speeds by 10^x. Keeps numbers computer-reasonable.
# Item Demand per Second
ips = 48
# Vessel Inputs
VCount = 60
dVSpeed = 1200
VCargo = 600
# Takeoff/Landing Adjustments
VDelay = 12
dDAdjust = 2500
# Inclination
dthetaA = 6.23
dthetaB = 6.88
# Longitude of Ascending Node
dlanA = 38.32
dlanB = 126.85
# Orbital Period
TA = 18238
TB = 27636
# Orbital Radius
dRA = 332520
dRB = 438680

# Conversion Values (d)
VSpeed = dVSpeed / np.power(10, UScalar)
RA = dRA / np.power(10, UScalar)
RB = dRB / np.power(10, UScalar)
DAdjust = dDAdjust / np.power(10, UScalar)
thetaA = np.radians(dthetaA + 90)
thetaB = np.radians(dthetaB + 90)
lanA = np.radians(dlanA + 90)
lanB = np.radians(dlanB + 90)
lambdaA = (2 * np.pi) / TA
lambdaB = (2 * np.pi) / TB
VFac = VCount * VCargo


# Calculates the difference in inclination between the two orbits being compared.
def findAngleDifference():
    # Calculate the normal vectors' cartesian values
    Nax = np.cos(lanA) * np.cos(thetaA)
    Nay = np.sin(lanA) * np.cos(thetaA)
    Naz = np.sin(thetaA)
    Nbx = np.cos(lanB) * np.cos(thetaB)
    Nby = np.sin(lanB) * np.cos(thetaB)
    Nbz = np.sin(thetaB)
    # Find the dot product of the vectors
    N_AdotB = (Nax * Nbx) + (Nay * Nby) + (Naz * Nbz)
    # Arc cosine of the dot product of these two vectors calculates the angle between them
    return np.arccos(N_AdotB)


thetaF = findAngleDifference()


# Convert the parameters into parametric equations modeling the orbits, and calculates the distance between
# them at a given time.
# DAdjust takes off some distance from each planet, and the distance is then converted into vessel travel time.
# VDelay adds some time to account for take off and landing, VFac then converts the time into items/s.
def itemSim(t):
    return (VFac / (4 * VDelay + 2 * ((np.sqrt(np.power(
        RB * np.sin(lambdaB * t) - RA * np.sin(lambdaA * t) * np.cos(thetaF), 2) + np.power(
        RB * np.cos(lambdaB * t) - RA * np.cos(lambdaA * t), 2) + np.power(
        0 - RA * np.sin(lambdaA * t) * np.sin(thetaF), 2)) - 2 * DAdjust) / VSpeed)))


# Build the time step array
T = np.arange(start, stop, step, np.int64)


# Fa and Fb integrates the item/s over time data into an item count.
# Fa subtracts the item/s that the factory requires, while Fb does not
def Fa(t):
    res = np.zeros_like(t)
    marker = 0
    for i, val in enumerate(t):
        if i > 0:
            res[i] = res[i-1] + integrate.quad(itemSim, marker, val)[0] - (ips * (val - marker))
            marker += step
        else:
            res[i] = integrate.quad(itemSim, marker, val)[0] - (ips * (val - marker))
    return res


def Fb(t):
    res = np.zeros_like(t)
    marker = 0
    for i, val in enumerate(t):
        if i > 0:
            res[i] = res[i-1] + integrate.quad(itemSim, marker, val)[0]
            marker += step
        else:
            res[i] = integrate.quad(itemSim, marker, val)[0]
    return res


# Simulate the storage while it is being drained and while it is not
simWithDrain = Fa(T)
simWithoutDrain = Fb(T)


# findPeakData finds the local minima and maxima, does a few basic operations, and returns the results.
# If there are no local extrema (Transportation always/never meets factory needs regardless of distance), returns Null
# Returns an array of:
# Maxima values array
# Minima values array
# Average gain per cycle
# Average cycle time as the average time between each cycle's minimum distance.
# (in ascending time) Avg of the local maxima - minima (loss per cycle while planets are distant, average needed buffer)
# Maximum of the values of maxima - minima (absolute minimum buffer size for 100% uptime)
# (in ascending time) Avg of local minima - maxima (gain per cycle while planets are close, buffer size+gain per cycle)
def findPeakData(Orbit, test=False):
    maximaIndices = signal.find_peaks(Orbit)[0]
    maxima = Orbit[maximaIndices]
    minima = Orbit[signal.find_peaks(-Orbit)[0]]
    if maxima.size == 0:
        if test:
            print("No extrema!")
        return None
    extIndDiff = maxima.size - minima.size
    maxToMin = np.subtract(minima[0:minima.size+(1*extIndDiff)], maxima[0:maxima.size-(1*extIndDiff)])
    minToMax = np.subtract(maxima[1:], minima[0:minima.size+(2*extIndDiff)])
    avgMaxToMin = np.average(maxToMin)
    avgMinToMax = np.average(minToMax)
    avgGainPerCycle = np.average(np.subtract(maxima[1:maxima.size], maxima[0:maxima.size-1]))
    avgCycleTime = np.average(np.subtract(maximaIndices[1:maximaIndices.size], maximaIndices[0:maximaIndices.size-1]))
    if test:
        print("Average gain per cycle: %.3f" % avgGainPerCycle)
        print("Average cycle time: %.3f" % (avgCycleTime*2))
        print("Average needed buffer per side: %.3f" % np.abs(avgMaxToMin))
        print("Maximum needed buffer per side: %d" % np.ceil(np.amax(np.abs(maxToMin))))
        print("Average buffer + gain: %.3f" % avgMinToMax)
    return [maxima, minima, avgGainPerCycle, avgCycleTime, avgMaxToMin, np.amax(maxToMin), avgMinToMax]


findPeakData(simWithDrain, True)
# slope, intercept = stat.linear_regression(OrbitWithNegative, T, proportional=True)
# print("The average increase in storage is: %f items/s" % slope)
plt.plot(T, simWithDrain)
plt.show()
