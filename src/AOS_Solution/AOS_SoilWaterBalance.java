package AOS_Solution;

import Structs.*;

public class AOS_SoilWaterBalance {
    //Function to execute AquaCrop-OS soil water balance module
    public static Object[] run(Crop Crop, Soil Soil, Weather Weather, IrrMngtStruct IrrMngt, FieldMngtStruct FieldMngt,
                               GwStruct Groundwater, InitCondStruct InitCond, boolean GrowingSeason) {
        //Unpack weather structure %%
        double P = Weather.Precip;
        double Et0 = Weather.RefET;
        int GDD = Weather.GDD;

        //Store initial conditions for updating
        InitCondStruct NewCond = InitCond;

        //Check for presence of groundwater table
        NewCond = AOS_CheckGroundwaterTable(Soil, Groundwater, NewCond);

        return new Object[]{Crop, Soil};
    }

    private static InitCondStruct AOS_CheckGroundwaterTable(Soil Soil, GwStruct Groundwater, InitCondStruct NewCond) {
        return null;
    }
}