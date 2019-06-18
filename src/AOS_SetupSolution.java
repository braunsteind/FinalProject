import Structs.*;

public class AOS_SetupSolution {

    public static Object[] run(AOS_InitialiseStruct AOS_InitialiseStruct, ClockStruct AOS_ClockStruct) {
        //Extract initial conditions
        InitCondStruct InitCond = AOS_InitialiseStruct.InitialCondition;

        //Perform growing season check
        boolean GrowingSeason = CheckGrowingSeason(InitCond, AOS_ClockStruct);

        //Get model parameter structures
        Soil Soil = AOS_InitialiseStruct.Parameter.soil;
        GwStruct Groundwater = AOS_InitialiseStruct.Groundwater;

        IrrMngtStruct IrrMngt = new IrrMngtStruct();
        FieldMngtStruct FieldMngt = new FieldMngtStruct();
        Crop Crop;
        if (AOS_ClockStruct.SeasonCounter > 0) {
            Crop = AOS_InitialiseStruct.Parameter.crop[AOS_ClockStruct.SeasonCounter - 1];
            IrrMngt = AOS_InitialiseStruct.IrrigationManagement;

            if (GrowingSeason) {
                FieldMngt = AOS_InitialiseStruct.FieldManagement[AOS_ClockStruct.SeasonCounter];
            } else {
                //TODO check that Fallow value is 1
                int Fallow = 1;
                FieldMngt = AOS_InitialiseStruct.FieldManagement[Fallow];
            }
        } else {
            //Assign crop, irrigation management, and field management structures
            Crop = new Crop();
            Crop.Zmin = 0.3;
            Crop.Aer = 5;
            IrrMngt.IrrMethod = -99;
            //TODO check that Fallow value is 1
            int Fallow = 1;
            FieldMngt = AOS_InitialiseStruct.FieldManagement[Fallow];
        }

        CO2 CO2 = AOS_InitialiseStruct.Parameter.CO2;

        //Get weather and GDD data for current day
        Object[] a = ExtractWeatherData(AOS_InitialiseStruct, AOS_ClockStruct, Crop, InitCond, GrowingSeason);
        Weather Weather = (Structs.Weather) a[0];
        InitCondStruct NewCond = (InitCondStruct) a[1];

        //Check for germination
        NewCond = AOS_Germination(NewCond, Soil, Crop, Weather.GDD, GrowingSeason);

        return new Object[]{Crop, Soil, Weather, IrrMngt, FieldMngt, Groundwater, NewCond, GrowingSeason, CO2};
    }

    //Function to check if crop has germinated
    private static InitCondStruct AOS_Germination(InitCondStruct InitCond, Soil Soil, Crop Crop, double GDD, boolean GrowingSeason) {
        //Store initial conditions in new structure for updating
        InitCondStruct NewCond = InitCond;

        //Check for germination (if in growing season)
        if (GrowingSeason) {
            //Find compartments covered by top soil layer affecting germination
            int comp_sto = 0;
            for (int i = 0; i < Soil.comp.dzsum.length; i++) {
                if (Soil.comp.dzsum[i] >= Soil.zGerm) {
                    comp_sto = i;
                    break;
                }
            }
            //Calculate water content in top soil layer
            double Wr = 0;
            double WrFC = 0;
            double WrWP = 0;

            for (int ii = 0; ii < comp_sto; ii++) {
                //Get soil layer
                int layeri = Soil.comp.layer[ii];

                //Determine fraction of compartment covered by top soil layer
                double factor;
                if (Soil.comp.dzsum[ii] > Soil.zGerm) {
                    factor = 1 - ((Soil.comp.dzsum[ii] - Soil.zGerm) / Soil.comp.dz[ii]);
                } else {
                    factor = 1;
                }

                //Increment actual water storage (mm)
                Wr = Wr + (factor * 1000 * InitCond.th[ii] * Soil.comp.dz[ii]);
                //Increment water storage at field capacity (mm)
                WrFC = WrFC + (factor * 1000 * Soil.layer.th_fc[layeri] * Soil.comp.dz[ii]);
                //Increment water storage at permanent wilting point (mm)
                WrWP = WrWP + (factor * 1000 * Soil.layer.th_wp[layeri] * Soil.comp.dz[ii]);
            }

            //Limit actual water storage to not be less than zero
            if (Wr < 0) {
                Wr = 0;
            }
            //Calculate proportional water content
            double WcProp = 1 - ((WrFC - Wr) / (WrFC - WrWP));
            //Check if water content is above germination threshold
            if ((WcProp >= Crop.GermThr) && (!NewCond.Germination)) {
                //Crop has germinated
                NewCond.Germination = true;
                //If crop sown as seedling, turn on seedling protection
                if (Crop.PlantMethod == 1) {
                    NewCond.ProtectedSeed = true;
                } else {
                    //Crop is transplanted so no protection
                    NewCond.ProtectedSeed = false;
                }
            }

            //Increment delayed growth time counters if germination is yet to
            //occur, and also set seed protection to false if yet to germinate
            if (!NewCond.Germination) {
                NewCond.DelayedCDs = InitCond.DelayedCDs + 1;
                NewCond.DelayedGDDs = InitCond.DelayedGDDs + GDD;
                NewCond.ProtectedSeed = false;
            }
        } else {
            //Not in growing season so no germination calculation is performed.
            NewCond.Germination = false;
            NewCond.ProtectedSeed = false;
            NewCond.DelayedCDs = 0;
            NewCond.DelayedGDDs = 0;
        }

        return NewCond;
    }

    //Function to calculate number of growing degree days on current day
    private static Object[] GrowingDegreeDay(Crop Crop, InitCondStruct InitCond, double Tmax, double Tmin) {
        //Store initial conditions for updating
        InitCondStruct NewCond = InitCond;

        double GDD = 0;
        //Calculate GDDs
        if (Crop.GDDmethod == 1) {
            // Method 1
            double Tmean = (Tmax + Tmin) / 2;
            if (Tmean > Crop.Tupp) {
                Tmean = Crop.Tupp;
            } else if (Tmean < Crop.Tbase) {
                Tmean = Crop.Tbase;
            }
            GDD = Tmean - Crop.Tbase;

        } else if (Crop.GDDmethod == 2) {
            //Method 2
            if (Tmax > Crop.Tupp) {
                Tmax = Crop.Tupp;
            } else if (Tmax < Crop.Tbase) {
                Tmax = Crop.Tbase;
            }
            if (Tmin > Crop.Tupp) {
                Tmin = Crop.Tupp;
            } else if (Tmin < Crop.Tbase) {
                Tmin = Crop.Tbase;
            }
            double Tmean = (Tmax + Tmin) / 2;
            GDD = Tmean - Crop.Tbase;
        } else if (Crop.GDDmethod == 3) {
            //Method 3
            if (Tmax > Crop.Tupp) {
                Tmax = Crop.Tupp;
            } else if (Tmax < Crop.Tbase) {
                Tmax = Crop.Tbase;
            }
            if (Tmin > Crop.Tupp) {
                Tmin = Crop.Tupp;
            }
            double Tmean = (Tmax + Tmin) / 2;
            if (Tmean < Crop.Tbase) {
                Tmean = Crop.Tbase;
            }
            GDD = Tmean - Crop.Tbase;
        }

        // returns a struct of GDD & NewCond
        NewCond.GDDcum = InitCond.GDDcum + GDD;

        return new Object[]{GDD, NewCond};
    }

    //Function to extract weather data for current time step
    private static Object[] ExtractWeatherData(AOS_InitialiseStruct AOS_InitialiseStruct, ClockStruct AOS_ClockStruct, Crop Crop, InitCondStruct InitCond, boolean GrowingSeason) {
        //Store initial conditions for updating
        InitCondStruct NewCond = InitCond;
        //Extract weather dataset
        double[][] WeatherDB = AOS_InitialiseStruct.Weather;

        //Extract weather data for current time step
        //Get current date
        int Date = AOS_ClockStruct.StepStartTime;
        //Find row corresponding to the current date in dataset
        int Row = 0;
        for (int i = 0; i < WeatherDB.length; i++) {
            if (Date == WeatherDB[i][0]) {
                Row = i;
            }
        }

        //Get weather variables
        Weather Weather = new Weather();
        Weather.Tmin = WeatherDB[Row][1];
        Weather.Tmax = WeatherDB[Row][2];
        Weather.Precip = WeatherDB[Row][3];
        Weather.RefET = WeatherDB[Row][4];

        //Calculate growing degree days
        if (GrowingSeason) {
            //Calendar days after planting
            NewCond.DAP = NewCond.DAP + 1;
            //Growing degree days after planting
            Object[] a = GrowingDegreeDay(Crop, NewCond, Weather.Tmax, Weather.Tmin);
            Weather.GDD = (double) a[0];
            NewCond = (InitCondStruct) a[1];
        } else {
            //Calendar days after planting
            NewCond.DAP = 0;
            //Growing degree days after planting
            Weather.GDD = 0;
            NewCond.GDDcum = 0;
        }

        return new Object[]{Weather, NewCond};
    }

    //Function to check if current time-step is during growing season
    private static boolean CheckGrowingSeason(InitCondStruct InitCond, ClockStruct AOS_ClockStruct) {
        //Check if in growing season
        if (AOS_ClockStruct.SeasonCounter > 0) {
            int CurrentDate = AOS_ClockStruct.StepStartTime;
            double PlantingDate = AOS_ClockStruct.PlantingDate[AOS_ClockStruct.SeasonCounter - 1];
            double HarvestDate = AOS_ClockStruct.HarvestDate[AOS_ClockStruct.SeasonCounter - 1];

            return (CurrentDate >= PlantingDate) && (CurrentDate <= HarvestDate) && (!InitCond.CropMature) && (!InitCond.CropDead);
        }
        return false;
    }
}