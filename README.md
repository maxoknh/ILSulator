# ILSulator
ILSulator is a simple instrastellar item throughput calculator for Dyson Sphere Program's Interstellar Logistics Stations. ILSulator uses basic orbital parameters you can get from in DSP and the level of your vessels to calculate the minimum number of vessels and smallest buffer size you need to sustain logistics from one planet to another without warp.

ILSulator works by mathematically modeling the distance between two planets over time using their orbital parameters. Using this, the rate that items are moved from one planet to another can be calculated. The integral of this rate allows the user to determine the required buffer sizes (if needed) and vessel counts in order to maintain item throughput.

This project is currently usable but is unfinished and not user friendly, and requires some knowledge of coding in order to work.
