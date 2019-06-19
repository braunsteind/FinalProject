package AOS_Solution;

import Structs.*;


public class AOS_UpdateOutputs {

    //Function to store model outputs for current time step
    public static Object[] run(AOS_InitialiseStruct AOS_InitialiseStruct, ClockStruct AOS_ClockStruct, InitCondStruct InitCond, SoilWatOutStruct SoilWBOut, IrrMngtStruct IrrMngt,
                               double GDD, boolean GrowingSeason) {
        //Store initial conditions for updating
        InitCondStruct NewCond = InitCond;

        //Store daily outputs
        //Extract outputs
        Outputs Outputs = AOS_InitialiseStruct.Outputs;
        int row_day = AOS_ClockStruct.TimeStepCounter;
        int row_gs = AOS_ClockStruct.SeasonCounter;

        double IrrDay, IrrTot;
        if (GrowingSeason) {
            if (IrrMngt.IrrMethod == 4) {
                //Net irrigation
                IrrDay = SoilWBOut.IrrNet;
                IrrTot = NewCond.IrrNetCum;
            } else {
                //Irrigation
                IrrDay = SoilWBOut.Irr;
                IrrTot = NewCond.IrrCum;
            }
        } else {
            IrrDay = 0;
            IrrTot = 0;
        }

        //1. Store water contents
        Outputs.WaterContents[row_day - 1][3] = AOS_ClockStruct.TimeStepCounter;
        Outputs.WaterContents[row_day - 1][4] = GrowingSeason ? 1 : 0;
        for (int i = 5; i < Outputs.WaterContents[0].length; i++) {
            Outputs.WaterContents[row_day - 1][i] = NewCond.th[i - 5];
        }

        //2. Store water fluxes
        Outputs.WaterFluxes[row_day - 1][3] = AOS_ClockStruct.TimeStepCounter;
        Outputs.WaterFluxes[row_day - 1][4] = GrowingSeason ? 1 : 0;
        Outputs.WaterFluxes[row_day - 1][5] = SoilWBOut.Wr;
        Outputs.WaterFluxes[row_day - 1][6] = NewCond.zGW;
        Outputs.WaterFluxes[row_day - 1][7] = NewCond.SurfaceStorage;
        Outputs.WaterFluxes[row_day - 1][8] = IrrDay;
        Outputs.WaterFluxes[row_day - 1][9] = SoilWBOut.Infl;
        Outputs.WaterFluxes[row_day - 1][10] = SoilWBOut.Runoff;
        Outputs.WaterFluxes[row_day - 1][11] = SoilWBOut.DeepPerc;
        Outputs.WaterFluxes[row_day - 1][12] = SoilWBOut.CR;
        Outputs.WaterFluxes[row_day - 1][13] = SoilWBOut.GwIn;
        Outputs.WaterFluxes[row_day - 1][14] = SoilWBOut.Es;
        Outputs.WaterFluxes[row_day - 1][15] = SoilWBOut.EsPot;
        Outputs.WaterFluxes[row_day - 1][16] = SoilWBOut.Tr;
        Outputs.WaterFluxes[row_day - 1][17] = SoilWBOut.TrPot;

        //3. Store crop growth
        Outputs.CropGrowth[row_day - 1][3] = AOS_ClockStruct.TimeStepCounter;
        Outputs.CropGrowth[row_day - 1][4] = GrowingSeason ? 1 : 0;
        Outputs.CropGrowth[row_day - 1][5] = GDD;
        Outputs.CropGrowth[row_day - 1][6] = NewCond.GDDcum;
        Outputs.CropGrowth[row_day - 1][7] = NewCond.Zroot;
        Outputs.CropGrowth[row_day - 1][8] = NewCond.CC;
        Outputs.CropGrowth[row_day - 1][9] = NewCond.CC_NS;
        Outputs.CropGrowth[row_day - 1][10] = NewCond.B;
        Outputs.CropGrowth[row_day - 1][11] = NewCond.B_NS;
        Outputs.CropGrowth[row_day - 1][12] = NewCond.HI;
        Outputs.CropGrowth[row_day - 1][13] = NewCond.HIadj;
        Outputs.CropGrowth[row_day - 1][14] = NewCond.Y;

        //Store seasonal outputs
        //Check if it is the end of a season
        if (AOS_ClockStruct.SeasonCounter > 0) {
            if ((NewCond.CropMature || NewCond.CropDead || (AOS_ClockStruct.StepEndTime ==
                    AOS_ClockStruct.HarvestDate[AOS_ClockStruct.SeasonCounter - 1])) && !NewCond.HarvestFlag) {
                //Get planting and harvest dates
                int plant_sdate = 0;
                for (int i = 0; i < AOS_ClockStruct.TimeSpan.length; i++) {
                    if (AOS_ClockStruct.TimeSpan[i] == AOS_ClockStruct.PlantingDate[AOS_ClockStruct.SeasonCounter - 1]) {
                        plant_sdate = i + 1;
                        break;
                    }
                }
                //TODO fix to datestr
//                plant_cdate = datestr(AOS_ClockStruct.PlantingDate(AOS_ClockStruct.SeasonCounter),'dd/mm/yyyy');
                String plant_cdate = "01/05/2015";
//                harvest_cdate = datestr(AOS_ClockStruct.StepStartTime,'dd/mm/yyyy');
                String harvest_cdate = "08/09/2015";
                int harvest_sdate = AOS_ClockStruct.TimeStepCounter;

                //Store end of season outputs
                Outputs.FinalOutput = new String[8];
                Outputs.FinalOutput[0] = String.valueOf(AOS_ClockStruct.SeasonCounter);
                Outputs.FinalOutput[1] = AOS_InitialiseStruct.CropChoices[AOS_ClockStruct.SeasonCounter - 1].name;
                Outputs.FinalOutput[2] = plant_cdate;
                Outputs.FinalOutput[3] = String.valueOf(plant_sdate);
                Outputs.FinalOutput[4] = harvest_cdate;
                Outputs.FinalOutput[5] = String.valueOf(harvest_sdate);
                //TODO bad value - check why
                Outputs.FinalOutput[6] = String.valueOf(NewCond.Y);
                Outputs.FinalOutput[7] = String.valueOf(IrrTot);

                //Set harvest flag
                NewCond.HarvestFlag = true;
            }
        }

        return new Object[]{NewCond, Outputs};
    }
}