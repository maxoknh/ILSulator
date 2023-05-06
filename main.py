import scipy.integrate as integrate
import numpy as np
import matplotlib.pyplot as plt
import statistics as stat

# Inputs (Angles are in decimal degrees)
# Graph Parameters (seconds)
start = 0
stop = 1000000
step = 100
UScalar = 5
# Factory Demand Input
ips = 44
# Vessel Inputs
VCount = 100
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


# Conversion Values
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


# Find the normal of both orbits and find the angle between them
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


# Convert the orbits into parametric equations, and calculate the distance between them at a given time.
# DAdjust takes off some distance from each planet, and the distance is then converted into vessel travel time.
# VDelay adds some time to account for take off and landing, and then VFac converts this into
# VFac then converts the time into items/s.
def orbitSim(t):
    return (VFac / (4 * VDelay + 2 * ((np.sqrt(np.power(
            RB * np.sin(lambdaB * t) - RA * np.sin(lambdaA * t) * np.cos(thetaF), 2) + np.power(
            RB * np.cos(lambdaB * t) - RA * np.cos(lambdaA * t), 2) + np.power(
            0 - RA * np.sin(lambdaA * t) * np.sin(thetaF), 2)) - 2 * DAdjust) / VSpeed)))


# print("%f" % orbitSim(0))


# Build the time step array
T = np.arange(start, stop, step, np.int64)


# Integrates the orbital sim according to an array of input times
# If the step size reaches 500,000, the integral will be broken up into 500,000 step-sized chunks to preserve precision
def F(t):
    res = np.zeros_like(t)
    for i, val in enumerate(t):
        chunkSize = 500000
        chunk = val // chunkSize
        if chunk > 0:
            for v in range(chunk):
                a = chunkSize * v
                b = chunkSize * (v+1)
                y = integrate.quad(orbitSim, a, b)[0]
                res[i] = res[i]+y
            yf = integrate.quad(orbitSim, chunkSize*chunk, val)[0]
            res[i] = res[i]+yf
        else:
            y = integrate.quad(orbitSim, 0, val)[0]
            res[i] = y
        res[i] = res[i] - (ips * val)
    return res


def Fa(t):
    res = np.zeros_like(t)
    marker = 0
    for i, val in enumerate(t):
        if i > 0:
            res[i] = res[i-1] + integrate.quad(orbitSim, marker, val)[0]
            res[i] = res[i] - (ips * (val - marker))
            marker = marker + step
        else:
            res[i] = integrate.quad(orbitSim, marker, val)[0] - (ips * (val - marker))
    return res


def Fb(t):
    res = np.zeros_like(t)
    marker = 0
    for i, val in enumerate(t):
        if i > 0:
            res[i] = res[i-1] + integrate.quad(orbitSim, marker, val)[0]
            marker = marker + step
        else:
            res[i] = integrate.quad(orbitSim, marker, val)[0]
    return res



def test():
    out = F(T)
    for i, val in enumerate(T):
        print("%f, %f" % (val, out[i]))


OrbitWithNegative = Fa(T)
OrbitWithoutNegative = Fb(T)

def findPeaks(Orbit):
    mark = 0
    num = 0
    increasing = False
    res = np.zeros_like(Orbit)
    resInd = np.zeros_like(Orbit)
    for i, val in enumerate(Orbit):
        if val > mark:
            if not increasing:
                increasing = True
                res[num] = mark
                resInd[num] = i
                num += 1
            mark = val
        elif val < mark:
            if increasing:
                increasing = False
                res[num] = mark
                resInd[num] = i
                num += 1
            mark = val
    return np.trim_zeros(res), np.trim_zeros(resInd)


mixedPeaks = findPeaks(OrbitWithNegative)
def printPeaks():
    dbp = 0
    dap = 0
    avgBef = 0
    avgBefCount = 0
    avgAfter = 0
    avgAfterCount = 0
    maxBef = 0
    maxAft = 0
    for i, val in enumerate(mixedPeaks[0]):
        if i == mixedPeaks[0].size-1 or val > mixedPeaks[0][i + 1]:
            if i != 0:
                dbp = val - mixedPeaks[0][i - 1]
                avgBef += dbp
                avgBefCount += 1
                # print("Difference before peak: %d" % dbp)
                if dbp > maxBef:
                    maxBef = dbp
            if i < mixedPeaks[0].size-1:
                dap = val - mixedPeaks[0][i + 1]
                avgAfter += dap
                avgAfterCount += 1
                # print("Difference after peak: %d" % dap)
                if dap > maxAft:
                    maxAft = dap
    print("Average difference pre peak: %d" % (avgBef/avgBefCount))
    print("Average difference post peak: %d" % (avgAfter/avgAfterCount))
    print("Maximum difference pre peak: %d" % maxBef)
    print("Maximum difference post peak: %d" % maxAft)
    print("Average item gain per cycle: %d" % (avgBef/avgBefCount - avgAfter/avgAfterCount))


def splitPeaks():
    highs = np.zeros_like(mixedPeaks[0])
    highInd = np.zeros_like(mixedPeaks[0])
    lows = np.zeros_like(mixedPeaks[0])
    lowInd = np.zeros_like(mixedPeaks[0])
    for i, val in enumerate(mixedPeaks[0]):
        if i % 2 == 0:
            highs[i//2] = val
            highInd[i//2] = mixedPeaks[1][i]
        else:
            lows[i//2] = val
            lowInd[i//2] = mixedPeaks[1][i]
    return [(np.trim_zeros(highs), np.trim_zeros(highInd)), (np.trim_zeros(lows), np.trim_zeros(lowInd))]


peaks = splitPeaks()
highPeaks = peaks[0]
lowPeaks = peaks[1]
#print("Average time between peaks: %d" % np.average(np.append(highPeaks[1], lowPeaks[1])))
slope, intercept = stat.linear_regression(OrbitWithNegative, T, proportional=True)
print("The average increase in storage is: %f items/s" % slope)


#printPeaks()
# test()
plt.plot(T, Fa(T))
plt.show()
#plt.plot(findPeaks())
#plt.show()


