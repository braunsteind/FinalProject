package AOS_Solution;

import Structs.*;

import static java.lang.Math.*;

public class AOS_SoilWaterBalance {
    //Function to execute AquaCrop-OS soil water balance module
    public static Object[] run(ClockStruct AOS_ClockStruct, Crop Crop, Soil Soil, Weather Weather, IrrMngtStruct IrrMngt, FieldMngtStruct FieldMngt,
                               GwStruct Groundwater, InitCondStruct InitCond, boolean GrowingSeason) {
        //Unpack weather structure
        double P = Weather.Precip;
        double Et0 = Weather.RefET;
        int GDD = Weather.GDD;

        //Store initial conditions for updating
        InitCondStruct NewCond = InitCond;

        //Check for presence of groundwater table
        //TODO
        NewCond = AOS_CheckGroundwaterTable(Soil, Groundwater, NewCond);

        //Pre-irrigation
        //TODO
        Object[] a = AOS_PreIrrigation(Soil, Crop, IrrMngt, NewCond, GrowingSeason);
        NewCond = (InitCondStruct) a[0];
        double PreIrr = (double) a[1];

        //Drainage
        //TODO
        a = AOS_Drainage(Soil, NewCond);
        NewCond = (InitCondStruct) a[0];
        double DeepPerc = (double) a[1];
        double[] FluxOut = (double[]) a[2];

        //Rainfall partitioning
        a = AOS_RainfallPartition(P, Soil, FieldMngt, NewCond);
        double Runoff = (double) a[0];
        double Infl = (double) a[1];
        NewCond = (InitCondStruct) a[2];

        //Irrigation
        a = AOS_Irrigation(AOS_ClockStruct, NewCond, IrrMngt, Crop, Soil, GrowingSeason, P, Runoff);

        return new Object[]{Crop, Soil};
    }

    //Function to check for presence of a groundwater table, and, if present,
    //to adjust compartment water contents and field capacities where necessary
    private static InitCondStruct AOS_CheckGroundwaterTable(Soil Soil, GwStruct Groundwater, InitCondStruct InitCond) {
        //Store initial conditions for updating %%
        InitCondStruct NewCond = InitCond;

        //Perform calculations (if variable water table is present)
        //TODO change the if
        //if ((Groundwater.WaterTable == 1) && (strcmp(Groundwater.Method,'Variable'))) {
        if ((Groundwater.WaterTable == 1)) {
            //TODO
        }
        return NewCond;
    }

    //Function to calculate pre-irrigation when in net irrigation mode
    private static Object[] AOS_PreIrrigation(Soil Soil, Crop Crop, IrrMngtStruct IrrMngt, InitCondStruct InitCond, boolean GrowingSeason) {
        //Store initial conditions for updating
        InitCondStruct NewCond = InitCond;
        double PreIrr = 0;

        //Calculate pre-irrigation needs
        if (GrowingSeason) {
            if ((IrrMngt.IrrMethod != 4) || (NewCond.DAP != 1)) {
                //No pre-irrigation as not in net irrigation mode or not on first day
                // of the growing season
                PreIrr = 0;
            } else {
                //TODO
                //Determine compartments covered by the root zone
//                double rootdepth = max(InitCond.Zroot, Crop.Zmin);
//                rootdepth = Math.round((rootdepth * 100)) / 100;
//                comp_sto = find(Soil.Comp.dzsum >= rootdepth, 1, 'first');
            }
        }

        return new Object[]{NewCond, PreIrr};
    }

    //Function to redistribute stored soil water
    private static Object[] AOS_Drainage(Soil Soil, InitCondStruct InitCond) {
        //Store initial conditions in new structure for updating
        InitCondStruct NewCond = InitCond;

        //Preallocate arrays
        double[] thnew = new double[Soil.nComp];
        double[] FluxOut = new double[Soil.nComp];

        //Initialise counters and states
        double drainsum = 0;
        double dthdt;

        //Calculate drainage and updated water contents
        for (int ii = 0; ii < Soil.nComp; ii++) {
            //Specify layer for compartment
            int layeri = Soil.comp.layer[ii];

            //Calculate drainage ability of compartment ii
            if (InitCond.th[ii] <= InitCond.th_fc_Adj[ii]) {
                dthdt = 0;
            } else if (InitCond.th[ii] >= Soil.layer.th_s[layeri]) {
                dthdt = Soil.layer.tau[layeri] * (Soil.layer.th_s[layeri] - Soil.layer.th_fc[layeri]);
                if ((InitCond.th[ii] - dthdt) < InitCond.th_fc_Adj[ii]) {
                    dthdt = InitCond.th[ii] - InitCond.th_fc_Adj[ii];
                }
            } else {
                dthdt = Soil.layer.tau[layeri] * (Soil.layer.th_s[layeri] - Soil.layer.th_fc[layeri]) * ((Math.exp(InitCond.th[ii] -
                        Soil.layer.th_fc[layeri]) - 1) / (Math.exp(Soil.layer.th_s[layeri] -
                        Soil.layer.th_fc[layeri]) - 1));
                if ((InitCond.th[ii] - dthdt) < InitCond.th_fc_Adj[ii]) {
                    dthdt = InitCond.th[ii] - InitCond.th_fc_Adj[ii];
                }
            }

            //Drainage from compartment ii (mm)
            double draincomp = dthdt * Soil.comp.dz[ii] * 1000;

            //Check drainage ability of compartment ii against cumulative drainage
            //from compartments above
            double excess = 0;
            double prethick = Soil.comp.dzsum[ii] - Soil.comp.dz[ii];
            double drainmax = dthdt * 1000 * prethick;
            boolean drainability = drainsum <= drainmax;

            //Drain compartment ii
            if (drainability) {
                //No storage needed. Update water content in compartment ii
                thnew[ii] = InitCond.th[ii] - dthdt;

                //Update cumulative drainage (mm)
                drainsum = drainsum + draincomp;

                //Restrict cumulative drainage to saturated hydraulic
                //conductivity and adjust excess drainage flow
                if (drainsum > Soil.layer.Ksat[layeri]) {
                    excess = excess + drainsum - Soil.layer.Ksat[layeri];
                    drainsum = Soil.layer.Ksat[layeri];
                }
            } else {
                //Storage is needed
                dthdt = drainsum / (1000 * prethick);

                //Calculate value of theta (thX) needed to provide a
                //drainage ability equal to cumulative drainage
                double thX;
                if (dthdt <= 0) {
                    thX = InitCond.th_fc_Adj[ii];
                } else if (Soil.layer.tau[layeri] > 0) {
                    double A = 1 + ((dthdt * (Math.exp(Soil.layer.th_s[layeri] - Soil.layer.th_fc[layeri]) - 1)) /
                            (Soil.layer.tau[layeri] * (Soil.layer.th_s[layeri] - Soil.layer.th_fc[layeri])));
                    thX = Soil.layer.th_fc[layeri] + Math.log(A);
                    if (thX < InitCond.th_fc_Adj[ii]) {
                        thX = InitCond.th_fc_Adj[ii];
                    }
                } else {
                    //TODO check if layer is actually layeri
                    //thX = Soil.layer.th_s[layer]+0.01;
                    thX = Soil.layer.th_s[layeri] + 0.01;
                }

                //Check thX against hydraulic properties of current soil layer
                if (thX <= Soil.layer.th_s[layeri]) {
                    //Increase compartment ii water content with cumulative
                    //drainage
                    thnew[ii] = InitCond.th[ii] + (drainsum / (1000 * Soil.comp.dz[ii]));
                    //TODO continue
                }
            }
        }


        double DeepPerc = drainsum;
        return new Object[]{NewCond, DeepPerc, FluxOut};
    }

    //Function to partition rainfall into surface runoff and infiltration using the curve number approach
    private static Object[] AOS_RainfallPartition(double P, Soil Soil, FieldMngtStruct FieldMngt, InitCondStruct InitCond) {
        //Store initial conditions for updating
        InitCondStruct NewCond = InitCond;
        double Runoff = 0, Infl = 0;

        //Calculate runoff
        if ((FieldMngt.SRinhb.compareTo("N") == 0) && ((FieldMngt.Bunds.compareTo("N") == 0) || (FieldMngt.zBund < 0.001))) {
            //Surface runoff is not inhibited and no soil bunds are on field
            //Reset submerged days
            NewCond.DaySubmerged = 0;
            //Adjust curve number for field management practices
            double CN = Soil.CN * (1 + (FieldMngt.CNadjPct / 100));
            if (Soil.AdjCN == 1) { //Adjust CN for antecedent moisture
                //Calculate upper and lowe curve number bounds
                double CNbot = round(1.4 * (Math.exp(-14 * Math.log(10))) + (0.507 * CN) -
                        Math.pow(0.00374 * CN, 2) + Math.pow(0.0000867 * CN, 3));
                double CNtop = round(5.6 * (Math.exp(-14 * Math.log(10))) + (2.33 * CN) -
                        Math.pow(0.0209 * CN, 2) + Math.pow(0.000076 * CN, 3));
                //Check which compartment cover depth of top soil used to adjust curve number
                int comp_sto = 0;
                for (int i = 0; i < Soil.comp.dzsum.length; i++) {
                    if (Soil.comp.dzsum[i] >= Soil.zCN) {
                        comp_sto = i + 1;
                        break;
                    }
                }
                if (comp_sto != 0) {
                    comp_sto = Soil.nComp;
                }
                //Calculate weighting factors by compartment
                double xx = 0;
                double[] wrel = new double[comp_sto];
                for (int ii = 0; ii < comp_sto; ii++) {
                    if (Soil.comp.dzsum[ii] > Soil.zCN) {
                        Soil.comp.dzsum[ii] = Soil.zCN;
                    }
                    double wx = 1.016 * (1 - Math.exp(-4.16 * (Soil.comp.dzsum[ii] / Soil.zCN)));
                    wrel[ii] = wx - xx;
                    if (wrel[ii] < 0) {
                        wrel[ii] = 0;
                    } else if (wrel[ii] > 1) {
                        wrel[ii] = 1;
                    }
                    xx = wx;
                }
                //Calculate relative wetness of top soil
                double wet_top = 0;
                for (int ii = 0; ii < comp_sto; ii++) {
                    int layeri = Soil.comp.layer[ii];
                    double th = Math.max(Soil.layer.th_wp[layeri], InitCond.th[ii]);
                    wet_top = wet_top + (wrel[ii] * ((th - Soil.layer.th_wp[layeri]) /
                            (Soil.layer.th_fc[layeri] - Soil.layer.th_wp[layeri])));
                }
                //Calculate adjusted curve number
                if (wet_top > 1)
                    wet_top = 1;
                else if (wet_top < 0)
                    wet_top = 0;
                CN = round(CNbot + (CNtop - CNbot) * wet_top);
            }
            //Partition rainfall into runoff and infiltration (mm)
            double S = (25400 / CN) - 254;
            double term = P - ((5 / 100) * S);
            if (term <= 0) {
                Runoff = 0;
                Infl = P;
            } else {
                Runoff = Math.pow(term, 2) / (P + 0.95 * S);
                Infl = P - Runoff;
            }
        } else {
            //Bunds on field, therefore no surface runoff
            Runoff = 0;
            Infl = P;
        }

        return new Object[]{Runoff, Infl, NewCond};
    }

    //Function to get irrigation depth for current day
    private static Object[] AOS_Irrigation(ClockStruct AOS_ClockStruct, InitCondStruct InitCond, IrrMngtStruct IrrMngt, Crop Crop, Soil Soil,
                                           boolean GrowingSeason, double Rain, double Runoff) {
        //Store intial conditions for updating
        InitCondStruct NewCond = InitCond;
        double Irr = 0;

        //Determine irrigation depth (mm/day) to be applied
        if (GrowingSeason) {
            //Calculate root zone water content and depletion
            Object[] a = AOS_RootZoneWater(Soil, Crop, NewCond);
            TAW DrReturn = (TAW) a[0];
            TAW TAWReturn = (TAW) a[1];
            thRZStruct thRZ = (thRZStruct) a[2];
            //Use root zone depletions and TAW only for triggering irrigation
            double Dr = DrReturn.Rz;
            double TAW = TAWReturn.Rz;

            //Determine adjustment for inflows and outflows on current day %
            double rootdepth, AbvFc;
            if (thRZ.Act > thRZ.FC) {
                rootdepth = max(InitCond.Zroot, Crop.Zmin);
                AbvFc = (thRZ.Act - thRZ.FC) * 1000 * rootdepth;
            } else {
                AbvFc = 0;
            }
            double WCadj = InitCond.Tpot + InitCond.Epot - Rain + Runoff - AbvFc;

            //Update growth stage if it is first day of a growing season
            if (NewCond.DAP == 1) {
                NewCond.GrowthStage = 1;
            }
            //Run irrigation depth calculation
            if (IrrMngt.IrrMethod == 0) { //Rainfed - no irrigation
                Irr = 0;
            } else if (IrrMngt.IrrMethod == 1) { //Irrigation - soil moisture
                //Get soil moisture target for current growth stage
                double SMT = IrrMngt.SMT[NewCond.GrowthStage];
                //Determine threshold to initiate irrigation
                double IrrThr = (1 - SMT / 100) * TAW;
                //Adjust depletion for inflows and outflows today
                Dr = Dr + WCadj;
                if (Dr < 0) {
                    Dr = 0;
                }
                //Check if depletion exceeds threshold
                if (Dr > IrrThr) {
                    //Irrigation will occur
                    double IrrReq = max(0, Dr);
                    //Adjust irrigation requirements for application efficiency
                    double EffAdj = ((100 - IrrMngt.AppEff) + 100) / 100;
                    IrrReq = IrrReq * EffAdj;
                    //Limit irrigation to maximum depth
                    Irr = min(IrrMngt.MaxIrr, IrrReq);
                } else {
                    //No irrigation
                    Irr = 0;
                }
            } else if (IrrMngt.IrrMethod == 2) { //Irrigation - fixed interval
                //Get number of days in growing season so far (subtract 1 so that
                //always irrigate first on day 1 of each growing season)
                int nDays = NewCond.DAP - 1;
                //Adjust depletion for inflows and outflows today
                Dr = Dr + WCadj;
                if (Dr < 0) {
                    Dr = 0;
                }
                if (nDays % IrrMngt.IrrInterval == 0) {
                    //Irrigation occurs
                    double IrrReq = max(0, Dr);
                    //Adjust irrigation requirements for application efficiency
                    double EffAdj = ((100 - IrrMngt.AppEff) + 100) / 100;
                    IrrReq = IrrReq * EffAdj;
                    //Limit irrigation to maximum depth
                    Irr = min(IrrMngt.MaxIrr, IrrReq);
                } else {
                    //No irrigation
                    Irr = 0;
                }
            } else if (IrrMngt.IrrMethod == 3) { //Irrigation - pre-defined schedule
                //Get current date
                int CurrentDate = AOS_ClockStruct.StepStartTime;
                //Find irrigation value corresponding to current date
                //TODO
//                Irr = IrrMngt.IrrigationSch((IrrMngt.IrrigationSch(:,1)==CurrentDate),2);
            } else if (IrrMngt.IrrMethod == 4) { //Irrigation - net irrigation
                //Net irrigation calculation performed after transpiration, so
                //irrigation is zero here
                Irr = 0;
            }
            //Update cumulative irrigation counter for growing season
            NewCond.IrrCum = NewCond.IrrCum + Irr;
        } else {
            //No irrigation outside growing season
            Irr = 0;
            NewCond.IrrCum = 0;
        }

        return new Object[]{NewCond, Irr};
    }

    //Function to calculate actual and total available water in the root zone at current time step
    private static Object[] AOS_RootZoneWater(Soil Soil, Crop Crop, InitCondStruct InitCond) {
        //Calculate root zone water content and available water
        //Compartments covered by the root zone
        double rootdepth = max(InitCond.Zroot, Crop.Zmin);
        rootdepth = round((rootdepth * 100)) / 100;
        int sum = 0;
        for (int i = 0; i < Soil.comp.dzsum.length; i++) {
            if (Soil.comp.dzsum[i] < rootdepth) {
                sum++;
            }
        }
        int comp_sto = min(sum + 1, Soil.nComp);
        //Initialise counters
        double WrAct = 0;
        double WrS = 0;
        double WrFC = 0;
        double WrWP = 0;
        double WrDry = 0;
        double WrAer = 0;

        double factor;
        for (int ii = 0; ii < comp_sto; ii++) {
            //Specify layer
            int layeri = Soil.comp.layer[ii];
            //Fraction of compartment covered by root zone
            if (Soil.comp.dzsum[ii] > rootdepth) {
                factor = 1 - ((Soil.comp.dzsum[ii] - rootdepth) / Soil.comp.dz[ii]);
            } else {
                factor = 1;
            }
            //Actual water storage in root zone (mm)
            WrAct = WrAct + (factor * 1000 * InitCond.th[ii] * Soil.comp.dz[ii]);
            //Water storage in root zone at saturation (mm)
            WrS = WrS + (factor * 1000 * Soil.layer.th_s[layeri] * Soil.comp.dz[ii]);
            //Water storage in root zone at field capacity (mm)
            WrFC = WrFC + (factor * 1000 * Soil.layer.th_fc[layeri] * Soil.comp.dz[ii]);
            // Water storage in root zone at permanent wilting point (mm)
            WrWP = WrWP + (factor * 1000 * Soil.layer.th_wp[layeri] * Soil.comp.dz[ii]);
            // Water storage in root zone at air dry (mm)
            WrDry = WrDry + (factor * 1000 * Soil.layer.th_dry[layeri] * Soil.comp.dz[ii]);
            // Water storage in root zone at aeration stress threshold (mm)
            WrAer = WrAer + (factor * 1000 * (Soil.layer.th_s[layeri] - (Crop.Aer / 100)) * Soil.comp.dz[ii]);
        }
        if (WrAct < 0) {
            WrAct = 0;
        }
        //Calculate total available water (m3/m3)
        TAW TAW = new TAW();
        TAW.Rz = max(WrFC - WrWP, 0);
        //Calculate soil water depletion (mm)
        TAW Dr = new TAW();
        Dr.Rz = min(WrFC - WrAct, TAW.Rz);

        thRZStruct thRZ = new thRZStruct();
        //Actual root zone water content (m3/m3)
        thRZ.Act = WrAct / (rootdepth * 1000);
        //Root zone water content at saturation (m3/m3)
        thRZ.S = WrS / (rootdepth * 1000);
        //Root zone water content at field capacity (m3/m3)
        thRZ.FC = WrFC / (rootdepth * 1000);
        //Root zone water content at permanent wilting point (m3/m3)
        thRZ.WP = WrWP / (rootdepth * 1000);
        //Root zone water content at air dry (m3/m3)
        thRZ.Dry = WrDry / (rootdepth * 1000);
        //Root zone water content at aeration stress threshold (m3/m3)
        thRZ.Aer = WrAer / (rootdepth * 1000);

        //Calculate top soil water content and available water
        if (rootdepth > Soil.zTop) {
            //Determine compartments covered by the top soil
            double ztopdepth = Soil.zTop;
            ztopdepth = round((ztopdepth * 100)) / 100;
            sum = 0;
            for (int i = 0; i < Soil.comp.dzsum.length; i++) {
                if (Soil.comp.dzsum[i] < rootdepth) {
                    sum++;
                }
            }
            comp_sto = sum + 1;
            //Initialise counters
            double WrAct_Zt = 0;
            double WrFC_Zt = 0;
            double WrWP_Zt = 0;
            //Calculate water storage in top soil
            for (int ii = 0; ii < comp_sto; ii++) {
                //Specify layer
                int layeri = Soil.comp.layer[ii];
                //Fraction of compartment covered by root zone
                if (Soil.comp.dzsum[ii] > ztopdepth) {
                    factor = 1 - ((Soil.comp.dzsum[ii] - ztopdepth) / Soil.comp.dz[ii]);
                } else {
                    factor = 1;
                }
                //Actual water storage in top soil (mm)
                WrAct_Zt = WrAct_Zt + (factor * 1000 * InitCond.th[ii] * Soil.comp.dz[ii]);
                //Water storage in top soil at field capacity (mm)
                WrFC_Zt = WrFC_Zt + (factor * 1000 * Soil.layer.th_fc[layeri] * Soil.comp.dz[ii]);
                // Water storage in top soil at permanent wilting point (mm)
                WrWP_Zt = WrWP_Zt + (factor * 1000 * Soil.layer.th_wp[layeri] * Soil.comp.dz[ii]);
            }
            //Ensure available water in top soil is not less than zero
            if (WrAct_Zt < 0) {
                WrAct_Zt = 0;
            }
            //Calculate total available water in top soil (m3/m3)
            TAW.Zt = max(WrFC_Zt - WrWP_Zt, 0);
            //Calculate depletion in top soil (mm)
            Dr.Zt = min(WrFC_Zt - WrAct_Zt, TAW.Zt);
        } else {
            // Set top soil depletions and TAW to root zone values
            Dr.Zt = Dr.Rz;
            TAW.Zt = TAW.Rz;
        }

        return new Object[]{Dr, TAW, thRZ};
    }
}