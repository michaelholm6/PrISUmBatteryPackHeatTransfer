from function_definitions import *

cellNumber = 180 #unitless
transversePitch = 20 #mm
longitudinalPitch = 18.536 #mm
cellDiameter = 18.5 #mm
cellLength = .32535 #m
numberLongitudinal = 45 #unitless
numberTransverse = 4 #unitless
freestreamTemp = 30 #Celsius
velocity = 2.5 #m/s
surfaceTemp = 60 #Celsius
arrangement = "aligned" #"staggered" or "aligned"
diametricalPitch = 0 #mm, if aligned, say this is zero


fluidDensity = calculateFluidDensity(surfaceTemp, freestreamTemp)
correctionFactor = correctionFactorCalc(arrangement, 10)
dynamicViscosity = calculateDynamicViscosity(surfaceTemp, freestreamTemp)
maxReynolds = findMaxReynolds(fluidDensity, cellDiameter, dynamicViscosity, transversePitch, velocity, diametricalPitch)
freestreamPrandtl, surfacePrandtl, filmPrandtl = prandtlNumberCalculation(surfaceTemp, freestreamTemp, 45)
constantOne, constantTwo = constantCalculation(maxReynolds, transversePitch, longitudinalPitch, arrangement)
nusseltNumber = nusseltNumberCalculation(constantOne, constantTwo, maxReynolds, surfacePrandtl, freestreamPrandtl, surfaceTemp, freestreamTemp, correctionFactor)
fluidThermalConductivity = fluidThermalConductivityCalculation(surfaceTemp, freestreamTemp)
hBar = calculateAverageConvectiveCoefficient(fluidThermalConductivity, cellDiameter, nusseltNumber)
specificHeat = calculateFluidSpecificHeat(surfaceTemp, freestreamTemp)
exitTemp = calculateExitTemp(cellDiameter, cellNumber, hBar, fluidDensity, velocity, numberTransverse, transversePitch, specificHeat, surfaceTemp, freestreamTemp)
logMeanTempDif = calculateLogMeanTempDifference(surfaceTemp, freestreamTemp, exitTemp)
totalHeatTransfer = calculateTotalHeatTransfer(cellNumber, hBar, cellDiameter, logMeanTempDif, cellLength)
print("Exit temperature of the air:" + str(exitTemp))
print("Total heat transfer:" + str(totalHeatTransfer))