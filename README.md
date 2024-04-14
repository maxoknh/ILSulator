# ILSulator
ILSulator is a simple intrastellar item throughput calculator for Dyson Sphere Program's Interstellar Logistics Stations. 
ILSulator uses basic orbital parameters of two concentric planets and the stats of your cargo vessels to calculate the minimum number of vessels and smallest buffer size you need to sustain logistics from one planet to another.

ILSulator only works with two planets which orbit the same body. ILSulator cannot calculate warp.

ILSulator works by modeling the distance between two planets over time using their orbital parameters. Using this, the rate that items are moved from one planet to another can be calculated. 
The integral of this rate allows the user to determine the required buffer sizes (if needed) and vessel counts in order to sustain items on another planet. 
Note that this does not actually simulate the vessels flying back and forth, only the number of items that could be moved given any single slice of time accumulated over many steps. 
It will not be 100% accurate due to vessel travel time and the movement of planets creating a doppler effect, which I cannot be bothered to code. But it can produce a good approximation.

This project is currently usable but is unfinished and is not user friendly. It will require some knowledge of coding in order to work.
