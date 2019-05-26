import Structs.InitialiseStrcut;
import Structs.SetupSolutionStruct;

public class AOS_SetupSolution {


    //global variables - where should they be?
    //public ClockStruct cs;
    //public InitialiseStrcut aos_is;


    /**
     * Function to check if crop has germinated
     *
     * This function corresponds to the AOS_Germination.m file
     *
     *
     * @param InitCond
     * @param Soil
     * @param Crop
     * @param GDD
     * @param GrowingSeason
     * @return
     */
    private NewCond Germination(Object InitCond, Object Soil, Object Crop,Object GDD, boolean GrowingSeason) {

        //Store initial conditions in new structure for updating
        Object NewCond = InitCond;

        //Check for germination (if in growing season)
        if (GrowingSeason == true) {
            //Find compartments covered by top soil layer affecting germination
            comp_sto = find(Soil.Comp.dzsum>=Soil.zGerm,1,'first');

            //Calculate water content in top soil layer
            double Wr = 0;
            double WrFC = 0;
            double WrWP = 0;

            for (int ii = 1; ii < comp_sto; ii++) {
                //Get soil layer
                layeri = Soil.Comp.Layer(ii);

                //Determine fraction of compartment covered by top soil layer
                if (Soil.Comp.dzsum(ii) > Soil.zGerm) {
                    double factor = 1-((Soil.Comp.dzsum(ii)-Soil.zGerm)/Soil.Comp.dz(ii));
                } else {
                    double factor = 1;
                }

                //Increment actual water storage (mm)
                Wr = Wr+(factor*1000*InitCond.th(ii)*Soil.Comp.dz(ii));
                //Increment water storage at field capacity (mm)
                WrFC = WrFC+(factor*1000*Soil.Layer.th_fc(layeri)*Soil.Comp.dz(ii));
                //Increment water storage at permanent wilting point (mm)
                WrWP = WrWP+(factor*1000*Soil.Layer.th_wp(layeri)*Soil.Comp.dz(ii));
            }

            //Limit actual water storage to not be less than zero
            if (Wr < 0) {
                Wr = 0;
            }
            //Calculate proportional water content
            double WcProp = 1-((WrFC-Wr)/(WrFC-WrWP));
            //Check if water content is above germination threshold
            if ((WcProp >= Crop.GermThr) && (NewCond.Germination == false)) {
                //Crop has germinated
                NewCond.Germination = true;
                //If crop sown as seedling, turn on seedling protection
                if (Crop.PlantMethod == true) {
                    NewCond.ProtectedSeed = true;
                } else {
                    //Crop is transplanted so no protection
                    NewCond.ProtectedSeed = false;
                }
            }

            //Increment delayed growth time counters if germination is yet to
            //occur, and also set seed protection to false if yet to germinate
            if (NewCond.Germination == false) {
                NewCond.DelayedCDs = InitCond.DelayedCDs+1;
                NewCond.DelayedGDDs = InitCond.DelayedGDDs+GDD;
                NewCond.ProtectedSeed = false;
            }


        } else  {
            //Not in growing season so no germination calculation is performed.
            NewCond.Germination = false;
            NewCond.ProtectedSeed = false;
            NewCond.DelayedCDs = 0;
            NewCond.DelayedGDDs = 0;
        }
    }



    /**
     * Function to calculate number of growing degree days on current day
     *
     * This function corresponds to the AOS_GrowingDegreeDay.m file
     *
     * @param Crop
     * @param InitCond
     * @param Tmax
     * @param Tmin
     * @return
     */
    public AnotherStruct GrowingDegreeDay(Object Crop, Object InitCond, Object Tmax,Object Tmin) {
        //Store initial conditions for updating
        Object NewCond = InitCond;

        //Calculate GDDs
        if (Crop.GDDmethod == 1) {
            // Method 1
            Tmean = (Tmax+Tmin)/2;
            Tmean(Tmean>Crop.Tupp) = Crop.Tupp;
            Tmean(Tmean<Crop.Tbase) = Crop.Tbase;
            GDD = Tmean-Crop.Tbase;

        } else if (Crop.GDDmethod == 2) {
            //Method 2
            Tmax(Tmax>Crop.Tupp) = Crop.Tupp;
            Tmax(Tmax<Crop.Tbase) = Crop.Tbase;
            Tmin(Tmin>Crop.Tupp) = Crop.Tupp;
            Tmin(Tmin<Crop.Tbase) = Crop.Tbase;
            Tmean = (Tmax+Tmin)/2;
            GDD = Tmean-Crop.Tbase;
        } else if (Crop.GDDmethod == 3) {
            //Method 3
            Tmax(Tmax>Crop.Tupp) = Crop.Tupp;
            Tmax(Tmax<Crop.Tbase) = Crop.Tbase;
            Tmin(Tmin>Crop.Tupp) = Crop.Tupp;
            Tmean = (Tmax+Tmin)/2;
            Tmean(Tmean<Crop.Tbase) = Crop.Tbase;
            GDD = Tmean-Crop.Tbase;
        }

        // returns a struct of GDD & NewCond
        NewCond.GDDcum = InitCond.GDDcum+GDD;
    }


    /**
     * Function to extract weather data for current time step
     *
     * This function corresponds to AOS_ExtractWeatherData.m
     *
     * @param Crop
     * @param InitCond
     * @param GrowingSeason
     * @return
     */
    public SomeStruct ExtractWeatherData(Object Crop, Object InitCond, boolean GrowingSeason) {

        //Store initial conditions for updating
        Object NewCond = InitCond;
        //Extract weather dataset
        Object WeatherDB = AOS_InitialiseStruct.Weather;

        //Extract weather data for current time step

        //Get current date
        Object Date = AOS_ClockStruct.StepStartTime;

        //Find row corresponding to the current date in dataset
        Object Row = WeatherDB(:,1)==Date;

        //Get weather variables
        Weather.Tmin = WeatherDB(Row,2);
        Weather.Tmax = WeatherDB(Row,3);
        Weather.Precip = WeatherDB(Row,4);
        Weather.RefET = WeatherDB(Row,5);


        if (GrowingSeason == true) {
            //Calendar days after planting
            NewCond.DAP = NewCond.DAP+1;
            //Growing degree days after planting
            [Weather.GDD,NewCond] = this.GrowingDegreeDay(Crop, NewCond, Weather.Tmax, Weather.Tmin);
        } else {
            //Calendar days after planting
            NewCond.DAP = 0;
            //Growing degree days after planting
            Weather.GDD = 0;
            NewCond.GDDcum = 0;
        }

        //Should return a struct of Weather & NewCond
    }




    /**
     * function CheckGrowingSeason is the equvilant of
     * the AOS_CheckGrowingSeason() in matalb.
     *
     * if needed we can take it out from this class and make it
     * public so that other function could use it as well
     */
    public boolean CheckGrowingSeason() {

        Object InitCond = AOS_InitialiseStruct.InitialCondition;
        if (AOS_ClockStruct.SeasonCounter > 0) {

            Object CurrentDate = AOS_ClockStruct.StepStartTime;
            Object PlantingDate = AOS_ClockStruct.PlantingDate(AOS_ClockStruct.SeasonCounter);
            Object HarvestDate = AOS_ClockStruct.HarvestDate(AOS_ClockStruct.SeasonCounter);

            if ((CurrentDate >= PlantingDate) && (CurrentDate <= HarvestDate)&&(InitCond.CropMature == false) && (InitCond.CropDead == false)) {
                return true;
            } else {
                return false;
            }
        } else {
            return false;
        }
    }


    /**
     * Function to setup parameters for time-step solution
     *
     * This function represents the work flow of AOS_SetupSolution function
     *
     * @param s
     */
    public void run(SetupSolutionStruct s) {

        // should be Another type than Object
        Object InitCond = AOS_InitialiseStruct.InitialCondition;
        boolean GrowingSeason = this.CheckGrowingSeason();


        Object Soil = AOS_InitialiseStruct.Parameter.Soil;
        Object Groundwater = AOS_InitialiseStruct.Groundwater;

        if (AOS_ClockStruct.SeasonCounter > 0) {

            // what exactly does the dot (.) after Crop and IrrigationManagement say?
            Object Crop = AOS_InitialiseStruct.Parameter.Crop.(AOS_InitialiseStruct.CropChoices{AOS_ClockStruct.SeasonCounter});
            Object IrrMngt = AOS_InitialiseStruct.IrrigationManagement.(AOS_InitialiseStruct.CropChoices{AOS_ClockStruct.SeasonCounter});

            if (GrowingSeason == true) {

            } else  {
                Object FieldMngt = AOS_InitialiseStruct.FieldManagement.(AOS_InitialiseStruct.CropChoices{AOS_ClockStruct.SeasonCounter});
            }

        } else {
            Object Crop = struct('Zmin',0.3,'Aer',5);
            Object IrrMngt = struct('IrrMethod',-99);
            Object FieldMngt = AOS_InitialiseStruct.FieldManagement.Fallow;
        }

        Object CO2 = AOS_InitialiseStruct.Parameter.CO2;

        //Get weather and GDD data for current day
        [Weather,NewCond] = AOS_ExtractWeatherData(Crop,InitCond,GrowingSeason);

        //Check for germination
        NewCond = AOS_Germination(NewCond,Soil,Crop,Weather.GDD,GrowingSeason);

    }
}