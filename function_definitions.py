import typing
import numpy
from math import log, exp


def correctionFactorCalc(arrangement: str, longitudinalCellNumber: int) -> float:
    """
    This checks the number of cells in the longitudinal direction.
    If the number is less than 16, a correction factor must be applied to the result.
    :param arrangement: str (aligned or staggered)
    :param longitudinalCellNumber: int (unitless)
    :return: float (unitless)
    """
    correctionFactorDict = {}
    if arrangement == "aligned":
        correctionFactorDict = {1: .7,
                                2: .8,
                                3: .86,
                                4: .9,
                                5: .92,
                                6: .935,
                                7: .95,
                                8: .96,
                                9: .965,
                                10: .97,
                                11: .975,
                                12: .98,
                                13: .98,
                                14: .983,
                                15: .986,
                                16: .99}
    elif arrangement == "staggered":
        correctionFactorDict = {1: .64,
                                2: .76,
                                3: .84,
                                4: .89,
                                5: .92,
                                6: .935,
                                7: .95,
                                8: .96,
                                9: .965,
                                10: .97,
                                11: .975,
                                12: .98,
                                13: .98,
                                14: .983,
                                15: .986,
                                16: .99}
    correctionFactor = correctionFactorDict.get(longitudinalCellNumber, 1)
    return correctionFactor


def findMaxReynolds(airDensity: float, cellDiameter: float, dynamicViscosity: float, transversePitch: float,
                    freestreamVelocity: float, diametricalPitch: float = 0) -> float:
    """
    This finds the maximum reynolds number in the system. If the system is staggered, diametrical pitch must be specified. Otherwise, it's default
    value is zero and is not used.
    :param airDensity: float (kg/m^3)
    :param cellDiameter: float (mm)
    :param dynamicViscosity: float (kg/ms)
    :param transversePitch: float (mm)
    :param freestreamVelocity: float (m/s)
    :param diametricalPitch: float (mm & optional)
    :return: float (unitless)
    """
    velocity1 = (transversePitch / (transversePitch - cellDiameter)) * freestreamVelocity
    velocity2 = (transversePitch / (2 * (diametricalPitch - cellDiameter))) * freestreamVelocity
    # this will return a negative number in the aligned case, which will me ignored by the max function used later in this function.
    velocityMax = max(velocity1, velocity2)
    maxReynolds = (airDensity * velocityMax * cellDiameter / 1000) / dynamicViscosity
    return maxReynolds


def prandtlNumberCalculation(surfaceTemperature: float, freestreamTemperature: float, filmTemperature: float = 0) \
        -> typing.Tuple[float, float, float]:
    """
    A prandtl number must be calculated for the air at both the cell temperature and the freestream temperature. These numbers
    are taken from this website:
    :param filmTemperature: float (celsius & optional)
    :param surfaceTemperature: float (celsius)
    :param freestreamTemperature: float (celsius)
    :return: Tuple [freestream, surface, film] (unitless)
    """
    celsiusDataPoints = [0.0,
                         6.9,
                         15.6,
                         26.9,
                         46.9,
                         66.9,
                         86.9,
                         106.9,
                         126.9]
    prandtlNumberDataPoints = [0.711,
                               0.710,
                               0.709,
                               0.707,
                               0.705,
                               0.703,
                               0.701,
                               0.700,
                               0.699]
    surfacePrandtl = numpy.interp(surfaceTemperature, celsiusDataPoints, prandtlNumberDataPoints)
    freestreamPrandtl = numpy.interp(freestreamTemperature, celsiusDataPoints, prandtlNumberDataPoints)
    filmPrandtl = 0

    if filmTemperature != 0:
        filmPrandtl = numpy.interp(filmTemperature, celsiusDataPoints, prandtlNumberDataPoints)
        return freestreamPrandtl, surfacePrandtl, filmPrandtl

    return freestreamPrandtl, surfacePrandtl, filmPrandtl


def constantCalculation(maxReynolds: float, transversePitch: float, longitudinalPitch: float, arrangement: str) -> \
        typing.Tuple[float, float]:
    """
    Calculates necessary constants using max reynolds number of the system.
    :param longitudinalPitch: float (mm)
    :param transversePitch: float (mm)
    :param arrangement: string (aligned or staggered)
    :param maxReynolds: float (unitless)
    :return: Tuple [constantOne, constantTwo] (unitless)
    """
    constantOne = 0
    constantTwo = 0
    if arrangement == "aligned":
        if 10 <= maxReynolds < 100:
            constantOne = .8
            constantTwo = .4

        elif 100 <= maxReynolds < 1000:
            constantOne = 0
            constantTwo = 0

        elif (1000 <= maxReynolds < 2 * (10 ** 5)) & (transversePitch / longitudinalPitch > .7):
            constantOne = .27
            constantTwo = .63

        elif (1000 <= maxReynolds < 2 * (10 ** 5)) & (transversePitch / longitudinalPitch < .7):
            print(
                "aligned tubes are inefficient in this geometry case. A staggered arrangement should be used instead.")

        elif 2 * (10 ** 5) <= maxReynolds < 2 * (10 ** 6):
            constantOne = .021
            constantTwo = .84

        elif maxReynolds > 2 * (10 ** 6):
            print("The Reynold's number for this case too large, try reducing flow rate")

        elif maxReynolds < 10:
            print("The Reynold's number for this case too small, try increasing flow rate")
    if arrangement == "staggered":
        if 10 <= maxReynolds < 100:
            constantOne = .9
            constantTwo = .4

        elif 100 <= maxReynolds < 1000:
            constantOne = 0
            constantTwo = 0

        elif (1000 <= maxReynolds < 2 * (10 ** 5)) & (transversePitch / longitudinalPitch < 2):
            constantOne = .35 * ((transversePitch / longitudinalPitch) ** .2)
            constantTwo = .6

        elif (1000 <= maxReynolds < 2 * (10 ** 5)) & (transversePitch / longitudinalPitch > 2):
            constantOne = .4
            constantTwo = .6

        elif 2 * (10 ** 5) <= maxReynolds < 2 * (10 ** 6):
            constantOne = .022
            constantTwo = .84

        elif maxReynolds > 2 * (10 ** 6):
            print("The Reynold's number for this case too large, try reducing flow rate")

        elif maxReynolds < 10:
            print("The Reynold's number for this case too small, try increasing flow rate")
    return constantOne, constantTwo


def calculateDynamicViscosity(surfaceTemp: float, freestreamTemp: float) -> float:
    """
    Calculates dynamic viscosity at the film temperature of the system. Data taken from Engineering Toolbox website.
    :param surfaceTemp: float (celsius)
    :param freestreamTemp: float (celsius)
    :return: float (Ns/m^2)
    """
    filmTemp = (surfaceTemp + freestreamTemp) / 2
    celsiusDataPoints = [0,
                         5,
                         10,
                         15,
                         20,
                         25,
                         30,
                         40,
                         50,
                         60,
                         80,
                         100,
                         125]
    dynamicViscosityDataPoints = [17.15,
                                  17.40,
                                  17.64,
                                  17.89,
                                  18.13,
                                  18.37,
                                  18.60,
                                  19.07,
                                  19.53,
                                  19.99,
                                  20.88,
                                  21.74,
                                  22.79]
    dynamicViscosity = numpy.interp(filmTemp, celsiusDataPoints, dynamicViscosityDataPoints) * (10**(-6))
    return dynamicViscosity


def nusseltNumberCalculation(constantOne: float, constantTwo: float, maxReynolds: float, surfacePrandtl: float,
                             freestreamPrandtl: float, surfaceTemperature: float, freestreamTemperature: float,
                             correctionFactor: float) -> float:
    """
    Calculation of the nusselt number for the given situation. If the constants given are both zero, the nusselt number
    for a cylinder is calculated.
    :param correctionFactor: float (unitless)
    :param freestreamTemperature: float (celsius)
    :param surfaceTemperature: float (celsius)
    :param constantOne: float (unitless)
    :param constantTwo: float (unitless)
    :param maxReynolds: float (unitless)
    :param surfacePrandtl: float (unitless)
    :param freestreamPrandtl: float (unitless)
    :return: float (unitless)
    """
    if (constantOne == 0) & (constantTwo == 0):
        filmPrandtl = (prandtlNumberCalculation(0, 0, (surfaceTemperature + freestreamTemperature) / 2))[2]
        nusselt = .683 * (maxReynolds ** .466) * (filmPrandtl ** 1 / 3)
        return nusselt * correctionFactor

    nusselt = constantOne * (maxReynolds ** constantTwo) * (freestreamPrandtl ** .36) * (
            (freestreamPrandtl / surfacePrandtl) ** .25)
    return nusselt * correctionFactor


def fluidThermalConductivityCalculation(surfaceTemperature: float, freestreamTemperature: float) -> float:
    """
    Calculates fluid thermal conductivity at the film temperature.
    :param surfaceTemperature: float (celsius)
    :param freestreamTemperature: float (celsius)
    :return: float (w/mk)
    """
    filmTemperature = (surfaceTemperature + freestreamTemperature) / 2
    celsiusDataPoints = [0,
                         5,
                         10,
                         15,
                         20,
                         25,
                         30,
                         40,
                         50,
                         60,
                         80,
                         100,
                         125]
    thermalConductivityDataPoints = [24.36,
                                     24.74,
                                     25.12,
                                     25.50,
                                     25.87,
                                     26.24,
                                     26.62,
                                     27.35,
                                     28.08,
                                     28.80,
                                     30.23,
                                     31.62,
                                     33.33]
    thermalConductivity = numpy.interp(filmTemperature, celsiusDataPoints, thermalConductivityDataPoints) / 1000
    return thermalConductivity


def calculateAverageConvectiveCoefficient(fluidThermalConductivity: float, cellDiameter: float,
                                          nusseltNumber: float) -> float:
    """
    Calculates average convective heat transfer coefficient.
    :param fluidThermalConductivity: float (w/mk)
    :param cellDiameter: float (mm)
    :param nusseltNumber: float (unitless)
    :return: float (w/m^2k)
    """
    averageConvectiveCoef = nusseltNumber * fluidThermalConductivity / (cellDiameter/1000)
    return averageConvectiveCoef


def calculateFluidDensity(surfaceTemperature: float, freestreamTemperature: float) -> float:
    """
    Calculates density of the fluid.
    :param surfaceTemperature: float (celsius)
    :param freestreamTemperature: float (celsius)
    :return: float (kg/m^3)
    """
    filmTemp = (surfaceTemperature + freestreamTemperature) / 2
    celsiusDataPoints = [0,
                         5,
                         10,
                         15,
                         20,
                         25,
                         30,
                         40,
                         50,
                         60,
                         80,
                         100,
                         125]
    densityDataPoints = [1.292,
                         1.268,
                         1.246,
                         1.225,
                         1.204,
                         1.184,
                         1.164,
                         1.127,
                         1.093,
                         1.060,
                         1.000,
                         0.9467,
                         0.8868]
    fluidDensity = numpy.interp(filmTemp, celsiusDataPoints, densityDataPoints)
    return fluidDensity


def calculateFluidSpecificHeat(surfaceTemperature: float, freestreamTemperature: float) -> float:
    """
    Calculate film temperature fluid specific heat.
    :param surfaceTemperature: float (celsius)
    :param freestreamTemperature: float (celsius)
    :return: float (J/KgK)
    """
    filmTemperature = (surfaceTemperature + freestreamTemperature) / 2
    celsiusDataPoints = [0.0,
                         6.9,
                         15.6,
                         26.9,
                         46.9,
                         66.9,
                         86.9,
                         107,
                         127]
    specificHeatDataPoints = [1.006,
                              1.006,
                              1.006,
                              1.006,
                              1.007,
                              1.009,
                              1.010,
                              1.012,
                              1.014]
    fluidSpecificHeat = numpy.interp(filmTemperature, celsiusDataPoints, specificHeatDataPoints)
    return fluidSpecificHeat * 1000


def calculateExitTemp(cellDiameter: float, cellNumber: float,
                      averageConvectiveCoef: float, filmTempDensity: float, freestreamVelocity: float,
                      cellNumberTransverse: float, transversePitch: float, filmTempSpecificHeat: float,
                      surfaceTemp: float, freestreamTemp: float) -> float:
    """
    Calculates the exit temperature of the air through the battery pack.
    :param freestreamTemp: float (celsius)
    :param surfaceTemp: float (celsius)
    :param cellDiameter: float (mm)
    :param cellNumber: float (unitless)
    :param averageConvectiveCoef: float (w/m^2K)
    :param filmTempDensity: float (kg/m^3)
    :param freestreamVelocity: float (m/s)
    :param cellNumberTransverse: float (unitless)
    :param transversePitch: float (mm)
    :param filmTempSpecificHeat: float (J/kgK)
    :return: float (celsius)
    """
    exitTemp = -((exp(
        (-3.14159 * cellDiameter / 1000 * cellNumber * averageConvectiveCoef) / (filmTempDensity * freestreamVelocity *
                                                                                 cellNumberTransverse * transversePitch
                                                                                 / 1000 * filmTempSpecificHeat)) * (surfaceTemp-freestreamTemp)) - surfaceTemp)
    return exitTemp


def calculateLogMeanTempDifference(surfaceTemperature: float, freestreamTemperature: float, exitTemperature: float) \
        -> float:
    """
    Calculates the log mean temperature difference of the system.
    :param surfaceTemperature: float (celsius)
    :param freestreamTemperature: float (celsius)
    :param exitTemperature: float (celsius)
    :return: float (celsius)
    """
    deltaOne = surfaceTemperature - freestreamTemperature
    deltaTwo = surfaceTemperature - exitTemperature
    logMeanTempDifference = (deltaOne - deltaTwo) / log(deltaOne / deltaTwo)
    return logMeanTempDifference


def calculateTotalHeatTransfer(cellNumber: float, averageConvectiveCoef: float, cellDiameter: float,
                               logMeanTempDif: float, cellLength: float) -> float:
    """
    Calculates the total heat transfer of the system.
    :param cellLength: float (m)
    :param cellNumber: float (unitless)
    :param averageConvectiveCoef: float (w/m^2K)
    :param cellDiameter: float (mm)
    :param logMeanTempDif: float (celsius)
    :return: float (W)
    """
    totalHeatTransfer = cellNumber * averageConvectiveCoef * 3.14159 * cellDiameter/1000 * logMeanTempDif * cellLength
    return totalHeatTransfer


