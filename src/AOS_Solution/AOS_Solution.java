package AOS_Solution;

import Structs.*;

public class AOS_Solution {
    //Function to perform AquaCrop-OS solution for a single time step
    public static Object[] run(AOS_InitialiseStruct AOS_InitialiseStruct, ClockStruct AOS_ClockStruct, Crop Crop,
                               Soil Soil, Weather Weather, IrrMngtStruct IrrMngt, FieldMngtStruct FieldMngt,
                               GwStruct Groundwater, InitCondStruct InitCond, boolean GrowingSeason, CO2 CO2) {
        //Run simulations
        //1. Soil water balance
        Object[] a = AOS_SoilWaterBalance.run(AOS_ClockStruct, Crop, Soil, Weather, IrrMngt, FieldMngt, Groundwater, InitCond, GrowingSeason);
        InitCondStruct NewCond = (InitCondStruct) a[0];
        SoilWatOutStruct SoilWBOut = (SoilWatOutStruct) a[1];

        //2. Crop growth and yield formation
        NewCond = AOS_CropGrowthYieldForm.run(Crop, Soil, Weather, Groundwater, NewCond, SoilWBOut, GrowingSeason, CO2);

        //Update model outputs
        return AOS_UpdateOutputs.run(AOS_InitialiseStruct, AOS_ClockStruct, NewCond, SoilWBOut, IrrMngt, Weather.GDD, GrowingSeason);
    }
}