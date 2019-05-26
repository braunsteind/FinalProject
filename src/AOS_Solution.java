import Structs.SetupSolutionStruct;
import Structs.SolutionReturnStruct;

public class AOS_Solution {

    private SetupSolutionStruct sss;
    //should be completed ********************

    public AOS_Solution(SetupSolutionStruct sss) {
        this.sss = sss;
    }






    /**
     * Function to get irrigation depth for current day
     *
     * this function corresponds to AOS_Irrigation.m
     *
     * @param GrowingSeason
     * @return
     */
    private SomeStruct Irrigation(InitCond,IrrMngt,Crop,Soil, boolean GrowingSeason,Rain,Runoff) {

        //Declare global variables
        //global AOS_ClockStruct - using the global variable

        //Store intial conditions for updating
        Object NewCond = InitCond;
        //Determine irrigation depth (mm/day) to be applied
        if (GrowingSeason == true) {
            //Calculate root zone water content and depletion
            [~,Dr,TAW,thRZ] = AOS_RootZoneWater(Soil,Crop,NewCond);   //// NEW FUNCTION FOR AOS_RooTZoneWater should be created later
            //Use root zone depletions and TAW only for triggering irrigation
            Dr = Dr.Rz;
            TAW = TAW.Rz;

            //Determine adjustment for inflows and outflows on current day
            if (thRZ.Act > thRZ.FC) {
                rootdepth = max(InitCond.Zroot,Crop.Zmin);
                AbvFc = (thRZ.Act-thRZ.FC)*1000*rootdepth;
            } else {
                AbvFc = 0;
            }
            WCadj = InitCond.Tpot+InitCond.Epot-Rain+Runoff-AbvFc;

            //Update growth stage if it is first day of a growing season
            if (NewCond.DAP == 1) {
                NewCond.GrowthStage = 1;
            }
            //Run irrigation depth calculation
            if (IrrMngt.IrrMethod == 0) {   //Rainfed - no irrigation
                double Irr = 0;
            } else if (IrrMngt.IrrMethod == 1) {    //Irrigation - soil moisture
                //Get soil moisture target for current growth stage
                SMT = IrrMngt.SMT(NewCond.GrowthStage);
                //Determine threshold to initiate irrigation
                IrrThr = (1-SMT/100)*TAW;
                //Adjust depletion for inflows and outflows today
                Dr = Dr+WCadj;
                if (Dr < 0) {
                    Dr = 0;
                }
                //Check if depletion exceeds threshold
                if (Dr > IrrThr) {
                    //Irrigation will occur
                    double IrrReq = Math.max(0,Dr);
                    //Adjust irrigation requirements for application efficiency
                    EffAdj = ((100-IrrMngt.AppEff)+100)/100;
                    IrrReq = IrrReq*EffAdj;
                    //Limit irrigation to maximum depth
                    double Irr = Math.min(IrrMngt.MaxIrr,IrrReq);
                } else {
                    //No irrigation
                    double Irr = 0;
                }

            } else if (IrrMngt.IrrMethod == 2) {    //Irrigation - fixed interval
                //Get number of days in growing season so far (subtract 1 so that
                //always irrigate first on day 1 of each growing season)
                int nDays = NewCond.DAP-1;
                //Adjust depletion for inflows and outflows today
                double Dr = Dr+WCadj;
                if (Dr < 0) {
                    Dr = 0;
                }
                if (rem(nDays,IrrMngt.IrrInterval) == 0) {
                    //Irrigation occurs
                    double IrrReq = Math.max(0,Dr);
                    //Adjust irrigation requirements for application efficiency
                    EffAdj = ((100-IrrMngt.AppEff)+100)/100;
                    IrrReq = IrrReq*EffAdj;
                    //Limit irrigation to maximum depth
                    double Irr = Math.min(IrrMngt.MaxIrr,IrrReq);
                } else {
                    //No irrigation
                    double Irr = 0;
                }
            } else if (IrrMngt.IrrMethod == 3) {   //Irrigation - pre-defined schedule
                //Get current date
                CurrentDate = AOS_ClockStruct.StepStartTime;
                //Find irrigation value corresponding to current date
                Irr = IrrMngt.IrrigationSch((IrrMngt.IrrigationSch(:,1)==CurrentDate),2);  //??????
            } else if (IrrMngt.IrrMethod == 4) {    //Irrigation - net irrigation
                //Net irrigation calculation performed after transpiration, so
                //irrigation is zero here
                double Irr = 0;
            }
            //Update cumulative irrigation counter for growing season
            NewCond.IrrCum = NewCond.IrrCum+Irr;
        } else if (GrowingSeason == false) {
            //No irrigation outside growing season
            double Irr = 0;
            NewCond.IrrCum = 0;
        }

        //here should be returned struct of NewCond & Irr
    }










    /**
     * Function to partition rainfall into surface runoff and infiltration
     * using the curve number approach
     *
     * This function corresponds to AOS_RainfallPartition.m
     *
     * @return
     */
    private SomeStruct RainfallPartition(P,Soil,FieldMngt,InitCond) {

        //Store initial conditions for updating
        NewCond = InitCond;

        //Calculate runoff
        if ((FieldMngt.SRinhb == 'N') && ((FieldMngt.Bunds == 'N') || (FieldMngt.zBund < 0.001))) {
            //Surface runoff is not inhibited and no soil bunds are on field
            //Reset submerged days
            NewCond.DaySubmerged = 0;
            //Adjust curve number for field management practices
            CN = Soil.CN*(1+(FieldMngt.CNadjPct/100));

            if (Soil.AdjCN == 1) {  //Adjust CN for antecedent moisture
                //Calculate upper and lowe curve number bounds
                CNbot = round(1.4*(exp(-14*log(10)))+(0.507*CN)-(0.00374*CN^2)+(0.0000867*CN^3));
                CNtop = round(5.6*(exp(-14*log(10)))+(2.33*CN)-(0.0209*CN^2)+(0.000076*CN^3));

                //Check which compartment cover depth of top soil used to adjust curve number
                comp_sto = find(Soil.Comp.dzsum>=Soil.zCN,1,'first');
                if (isempty(comp_sto)) {
                    comp_sto = Soil.nComp;
                }
                //Calculate weighting factors by compartment
                double xx = 0;
                wrel = zeros(1,comp_sto);
                for (int ii = 1; i < comp_sto; i++) {
                    if (Soil.Comp.dzsum(ii) > Soil.zCN) {
                        Soil.Comp.dzsum(ii) = Soil.zCN;
                    }
                    double wx = 1.016*(1-exp(-4.16*(Soil.Comp.dzsum(ii)/Soil.zCN)));
                    wrel(ii) = wx-xx;
                    if (wrel(ii) < 0) {
                        wrel(ii) = 0;
                    } else if (wrel(ii) > 1) {
                        wrel(ii) = 1;
                    }
                    xx = wx;
                }
                //Calculate relative wetness of top soil
                double wet_top = 0;
                for (int ii = 1; i < comp_sto; i++) {
                    layeri = Soil.Comp.Layer(ii);
                    th = max(Soil.Layer.th_wp(layeri),InitCond.th(ii));
                    wet_top = wet_top+(wrel(ii)*((th-Soil.Layer.th_wp(layeri))/(Soil.Layer.th_fc(layeri)-Soil.Layer.th_wp(layeri))));
                }
                //Calculate adjusted curve number
                if (wet_top > 1) {
                    wet_top = 1;
                } else if (wet_top < 0) {
                    wet_top = 0;
                }
                CN = round(CNbot+(CNtop-CNbot)*wet_top);
            }
            //Partition rainfall into runoff and infiltration (mm)
            S = (25400/CN)-254;
            term = P-((5/100)*S);
            if (term <= 0) {
                double Runoff = 0;
                Infl = P;
            } else {
                double Runoff = (term^2)/(P+(1-(5/100))*S);
                Infl = P-Runoff;
            }

        } else {
            //Bunds on field, therefore no surface runoff
            double Runoff = 0;
            Infl = P;
        }
        // HERE SHOULD BE A RETURNED STRUCT OF Runoff, Infl & NewCond
    }










    /**
     * Function to calculate pre-irrigation when in net irrigation mode
     *
     * This function corresponds to AOS_PreIrrigation.m
     * @param GrowingSeason
     * @return
     */
    private SomeStruct PreIrrigation(Soil,Crop,IrrMngt,InitCond,boolean GrowingSeason) {

        //Store initial conditions for updating
        Object NewCond = InitCond;
        //Calculate pre-irrigation needs
        if (GrowingSeason == true) {
            if ((IrrMngt.IrrMethod ~= 4) || (NewCond.DAP ~= 1)) {
                //No pre-irrigation as not in net irrigation mode or not on first day
                //of the growing season
                double PreIrr = 0;
            } else {
                //Determine compartments covered by the root zone
                rootdepth = max(InitCond.Zroot,Crop.Zmin);
                rootdepth = round((rootdepth*100))/100;
                comp_sto = find(Soil.Comp.dzsum>=rootdepth,1,'first');
                //Calculate pre-irrigation requirements
                double PreIrr = 0;
                for (ii = 1; i< comp_sto; i++) {
                    //Get soil layer
                    layeri = Soil.Comp.Layer(ii);
                    //Determine critical water content threshold
                    thCrit = Soil.Layer.th_wp(layeri)+((IrrMngt.NetIrrSMT/100)*(Soil.Layer.th_fc(layeri)-Soil.Layer.th_wp(layeri)));
                    //Check if pre-irrigation is required
                    if (NewCond.th(ii) < thCrit) {
                        double PreIrr = PreIrr+((thCrit-NewCond.th(ii))*1000*Soil.Comp.dz(ii));
                        NewCond.th(ii) = thCrit;
                    }
                }
            }

        } else {
            double PreIrr = 0
        }
        //here should be the returned struct
    }







    /**
     * Function to execute AquaCrop-OS soil water balance module
     *
     * This function corresponds to AOS_SoilWaterBalance.m
     * @return
     */
    private SoilWaterBalanceStruct SoilWaterBalance(Crop,SoilWeather,IrrMngt,FieldMngt,Groundwater,InitCond, boolean GrowingSeason) {

        //Unpack weather structure
        P = Weather.Precip;
        Et0 = Weather.RefET;
        GDD = Weather.GDD;

        //Store initial conditions for updating
        Object NewCond = InitCond;

        //Check for presence of groundwater table
        NewCond = AOS_CheckGroundwaterTable(Soil,Groundwater,NewCond);

        //Pre-irrigation
        [NewCond,PreIrr] = this.PreIrrigation(Soil,Crop,IrrMngt,NewCond,GrowingSeason);
        //Drainage
        [NewCond,DeepPerc,FluxOut] = HelperUtils.Drainage(Soil,NewCond);
        //Rainfall partitioning
        [Runoff,Infl,NewCond] = this.RainfallPartition(P,Soil,FieldMngt,NewCond);
        //Irrigation
        [NewCond,Irr] = this.Irrigation(NewCond,IrrMngt,Crop,Soil,GrowingSeason,P,Runoff);

    }


    /**
     *
     * @return
     */
    private NewCond CropGrowthYieldForm(Crop,Soil,Weather,Groundwater,NewCond,SoilWBOut,GrowingSeason,CO2) {




    }


    /**
     *
     *
     * @return
     */
    private OutputsStruct UpdateOutputs(NewCond,SoilWBOut,IrrMngt,Weather.GDD,GrowingSeason) {






    }





    //SolutuinReturnStruct should be constructed from:
    //NewCond & Outputs
    public SolutionReturnStruct run() {
        // ##### Here should be the content of AOS_Solution

        //Run simulations
        //1. Soil water balance
        [NewCond,SoilWBOut] = this.SoilWaterBalance(this.sss);

        //2. Crop growth and yield formation
        Object NewCond = AOS_CropGrowthYieldForm(Crop,Soil,Weather,Groundwater,NewCond,SoilWBOut,GrowingSeason,CO2);

        //Update model outputs
        SolutionReturnStruct srs = AOS_UpdateOutputs(NewCond,SoilWBOut,IrrMngt,Weather.GDD,GrowingSeason);

        //returning struct of NewCond & Outputs
    }
}
