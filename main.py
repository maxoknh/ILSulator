import scipy.integrate as integrate
import scipy.signal as signal
import numpy as np
import matplotlib.pyplot as plt

# Only works with concentric orbits.

# Inputs

# Item Demand per Second
ips = 24
# Vessel Stats
VCount = 6
dVSpeed = 600
VCargo = 200

# Planet A:
# Orbit Radius (AU)
dRA = .105
# Orbital Period (seconds)
TA = 2048
# Orbit Inclination (degrees, minutes)
dInclinationA = 11
dInclinationAminutes = 47
# Longitude of Ascending Node (degrees, minutes)
dlanA = 83
dlanAminutes = 37

# Planet B:
# Orbit Radius (AU)
dRB = .065
# Orbital Period (seconds)
TB = 993
# Orbit Inclination (degrees, minutes)
dInclinationB = 7
dInclinationBminutes = 16
# Longitude of Ascending Node (degrees, minutes)
dlanB = 42
dlanBminutes = 41

# Simulation Time Parameters (seconds). Based on orbit with longest period.
# If step is too fine, (~5000 steps per revolution,) output seems to break.
revolutions = 5
stepsPerRevolution = 2000
start = 0
stop = np.maximum(TA,TB) * revolutions
step = stop/(revolutions * stepsPerRevolution)

# Takeoff/Landing Assumptions (seconds, meters)
VDelay = 12
dDAdjust = 2500

# Divides all numbers used by 10^x. Prevents integer overflow. May lower accuracy.
UScalar = 0






# Conversion Values (d)
VSpeed = dVSpeed / np.power(10, UScalar)
RA = dRA * 40000 / np.power(10, UScalar)
RB = dRB * 40000 / np.power(10, UScalar)
DAdjust = dDAdjust / np.power(10, UScalar)
thetaA = np.radians(dInclinationA + (dInclinationAminutes / 60) + 90)
thetaB = np.radians(dInclinationB + (dInclinationBminutes / 60) + 90)
lanA = np.radians(dlanA + (dlanAminutes / 60) + 90)
lanB = np.radians(dlanB + (dlanBminutes / 60) + 90)
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
    return (VFac / (4 * VDelay + 2 * (
        (np.sqrt(np.power(
        RB * np.sin(lambdaB * t) - 
        RA * np.sin(lambdaA * t) * np.cos(thetaF), 2) + 
        np.power(RB * np.cos(lambdaB * t) - RA * np.cos(lambdaA * t), 2) + 
        np.power(0 - RA * np.sin(lambdaA * t) * np.sin(thetaF), 2)) - 2 * DAdjust) / VSpeed)))


# Build the time step array
T = np.arange(start, stop, step, np.int64)


# Integrates the item/s data into an item count.
# One accounts for the item/s demand, while the other does not
def simWithDrain(t):
    res = np.zeros_like(t)
    for i in range(1, t.size, 1):
        res[i] = np.maximum(res[i-1] + integrate.quad(itemSim, t[i-1], t[i])[0] - (ips * (t[i] - t[i-1])), 0)
    return res


def simWithoutDrain(t):
    res = np.zeros_like(t)
    for i in range(1, t.size, 1):
        res[i] = res[i-1] + integrate.quad(itemSim, t[i-1], t[i])[0]
    return res


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
def findPeakData(data, t, verbo=True):
    maximaIndices = signal.find_peaks(data)[0]
    
    if maximaIndices.size == 0:
        if verbo: print("No extrema!")
        return None

    maxima = data[maximaIndices]
    maximaTimes = t[maximaIndices]    
    minima = data[signal.find_peaks(-data)[0]]
    comIndice = np.minimum(maxima.size, minima.size)
    maxToMin = np.subtract(maxima[0:comIndice], minima[0:comIndice])
    minToMax = np.subtract(minima[0:comIndice-1], maxima[1:comIndice])
    avgMaxToMin = np.average(maxToMin)
    avgMinToMax = np.average(minToMax)
    avgGainPerCycle = np.average(np.subtract(maxima[1:], maxima[0:maxima.size-1]))
    avgCycleTime = np.average(np.subtract(maximaTimes[1:], maximaTimes[0:maximaTimes.size-1]))
    
    if verbo:
        print("Average Peak to Peak Gain per Cycle: %.3f" % avgGainPerCycle)
        print("Maximum needed buffer per side/Maximum Peak to Trough Loss: %d" % np.ceil(np.amax(np.abs(maxToMin))))
        print("////////////////////////////////////////////////////////////////////")
        print("Average needed buffer per side: %.3f" % np.abs(avgMaxToMin))
        print("Average Trough to Peak Gain: %.3f" % np.abs(avgMinToMax))
        print("Average Peak Cycle time: %.3f" % (avgCycleTime))
        print("Number of maximums: %d" % (maxima.size))
        print("Number of minimums: %d" % (minima.size))

    return [maxima, minima, avgGainPerCycle, avgCycleTime, avgMaxToMin, np.amax(maxToMin), avgMinToMax]

finalData = simWithDrain(T)
findPeakData(finalData, T)
plt.plot(T, finalData)
plt.show()
