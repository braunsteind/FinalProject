import AOS_Solution.AOS_Solution;
import Structs.*;

public class AOS_PerformTimeStep {
    public static void PerformTimeStep(AOS_Initialize aos_initialize) {
        //Setup parameters for time-step solution
        Object[] a = AOS_SetupSolution.run(aos_initialize.AOS_InitialiseStruct, aos_initialize.clockStruct);
        Crop Crop = (Crop) a[0];
        Soil Soil = (Soil) a[1];
        Weather Weather = (Weather) a[2];
        IrrMngtStruct IrrMngt = (IrrMngtStruct) a[3];
        FieldMngtStruct FieldMngt = (FieldMngtStruct) a[4];
        GwStruct Groundwater = (GwStruct) a[5];
        InitCondStruct InitCond = (InitCondStruct) a[6];
        boolean GrowingSeason = (boolean) a[7];
        CO2 CO2 = (CO2) a[8];

        //Get model solution for current time-step
        a = AOS_Solution.run(Crop, Soil, Weather, IrrMngt, FieldMngt, Groundwater, InitCond, GrowingSeason, CO2);
        InitCondStruct NewCond = (InitCondStruct) a[0];
        Outputs Outputs = (Structs.Outputs) a[1];

        //Update initial conditions and outputs
        AOS_InitialiseStruct AOS_InitialiseStruct = aos_initialize.AOS_InitialiseStruct;
        AOS_InitialiseStruct.InitialCondition = NewCond;
        AOS_InitialiseStruct.Outputs = Outputs;

        //Check model termination
        AOS_CheckModelTermination(aos_initialize);

        //Update time step
        AOS_UpdateTime();
    }

    private static void AOS_CheckModelTermination(AOS_Initialize aos_initialize) {
        ClockStruct AOS_ClockStruct = aos_initialize.clockStruct;
        AOS_InitialiseStruct AOS_InitialiseStruct = aos_initialize.AOS_InitialiseStruct;
        //Check if current time-step is the last
        int CurrentTime = aos_initialize.clockStruct.StepEndTime;
        if (CurrentTime < AOS_ClockStruct.SimulationEndDate) {
            AOS_ClockStruct.ModelTermination = false;
        } else if (CurrentTime >= AOS_ClockStruct.SimulationEndDate) {
            AOS_ClockStruct.ModelTermination = true;
        }

        //Check if at the end of last growing season
        //Allow model to exit early if crop has reached maturity or died, and in
        //the last simulated growing season

        if (AOS_InitialiseStruct.InitialCondition.HarvestFlag && AOS_ClockStruct.SeasonCounter == AOS_ClockStruct.nSeasons) {
            AOS_ClockStruct.ModelTermination = true;
        }
    }

    //Function to update current time in model
    private static void AOS_UpdateTime() {
        //TODO
    }
}