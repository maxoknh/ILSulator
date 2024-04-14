import scipy.integrate as integrate
import scipy.signal as signal
import numpy as np
import matplotlib.pyplot as plt

# Only works with concentric orbits.
# How to use:
# Choose your planets, copy their orbital stats into the area below
# Set your desired transfer rate between the planets
# Set your cargo vessel stats, which can be found on the right side of the upgrade tab under 'Interstellar Vessel Upgrades'
# Choose your number of cargo vessels you wish to use
# Repeat until graph does not intersect with zero.

# Inputs

# Item Demand per Second
ips = 48
# Vessel Count
VCount = 54
# Vessel Speed
dVSpeed = 1200
# Vessel Cargo Count
VCargo = 600
# Items already in buffer on receiving planet
itemStartingBuffer = 0

# Planet A:
# Orbit Radius (AU)
dRA = 8.313
# Orbital Period (seconds)
TA = 18238
# Orbit Inclination (degrees, minutes)
dInclinationA = 6
dInclinationAminutes = 13
# Longitude of Ascending Node (degrees, minutes)
dlanA = 38
dlanAminutes = 19

# Planet B:
# Orbit Radius (AU)
dRB = 10.967
# Orbital Period (seconds)
TB = 27636
# Orbit Inclination (degrees, minutes)
dInclinationB = 6
dInclinationBminutes = 53
# Longitude of Ascending Node (degrees, minutes)
dlanB = 126
dlanBminutes = 51

# Simulation Time Parameters (seconds). Based on orbit with longest period.
# Step count can greatly affect accuracy depending on orbit size!
# A few revolutions is usually enough to get a good enough idea, higher revolutions is better for planets with high inclinations and different LAN's.
# Recommended step count is at least 100, accuracy seems to be ~1-5% better at 1000 based on rough testing on a large system.
# It is not recommended to exceed 5000 steps per revolution, the program may hang for a long time.
revolutions = 30
stepsPerRevolution = 100
start = 0
stop = np.maximum(TA,TB) * revolutions
step = stop/(revolutions * stepsPerRevolution)

# Takeoff/Landing Assumptions (seconds, meters)
VDelay = 12
dDAdjust = 2500

# Divides all numbers used by 10^x. May help prevent errors. May lower accuracy.
# Not much purpose anymore, but it's here.
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
# them at a given time, and then converting that distance into an item/s count.
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
T = np.arange(start, stop, step, np.float64)


# Integrates the item/s data into an item count.
# One accounts for the item/s demand, while the other does not
def simWithDrain(t):
    res = np.zeros_like(t)
    res[0] = itemStartingBuffer
    for i in range(1, t.size, 1):
        res[i] = np.maximum(res[i-1] + integrate.quad(itemSim, t[i-1], t[i])[0] - (ips * (t[i] - t[i-1])), 0)
    return res


def simWithoutDrain(t):
    res = np.zeros_like(t)
    res[0] = itemStartingBuffer
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

    bufferDepletesFirstCycle = False
    for i in range (1, maximaIndices[1], 1):
        if data[i] == 0:
            bufferDepletesFirstCycle = True
            break

    bufferDepletes = False
    for i in range(maximaIndices[1], data.size, 1):
        if data[i] == 0:
            bufferDepletes = True
            break
    
    if verbo:
        print("Average Peak to Peak Gain per Cycle: %.3f" % avgGainPerCycle)
        print("Maximum needed buffer per side (Maximum Max-to-Min Loss): %d" % np.ceil(np.amax(np.abs(maxToMin))))
        if bufferDepletes and bufferDepletesFirstCycle:
            print("Vessel count cannot sustain this demand!")
        elif not bufferDepletes and bufferDepletesFirstCycle:
            print("Vessel count sustains demand. Starting buffer may be needed.")
        elif bufferDepletes and not bufferDepletesFirstCycle:
            print("Vessel count cannot sustain demand. Starting buffer will delay depletion by at least one cycle.")
        else:
            print("Vessel count always sustains demand.")
        print("Starting buffer needed (Max Peak - First Peak): %.3f" % (np.ceil(np.amax(np.abs(maxToMin))) - maxima[0]))
        print("/////////////////////////////////Other Statistics///////////////////////////////////")
        print("Vessel count: %d vessels" % VCount)
        print("First Peak Height: %.3f items" % maxima[0])
        print("Average needed buffer per side: %.3f items" % np.abs(avgMaxToMin))
        print("Average Min-to-Max Gain: %.3f items" % np.abs(avgMinToMax))
        print("Average Time Between Peaks: %.3f seconds" % (avgCycleTime))
        print("Number of maximums: %d" % (maxima.size))
        print("Number of minimums: %d" % (minima.size))

    return [maxima, minima, avgGainPerCycle, avgCycleTime, avgMaxToMin, np.amax(maxToMin), avgMinToMax]

finalData = simWithDrain(T)
findPeakData(finalData, T)
plt.plot(T, finalData)
plt.show()
