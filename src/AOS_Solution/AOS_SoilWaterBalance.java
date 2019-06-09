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
        NewCond = (InitCondStruct) a[0];
        double Irr = (double) a[1];

        //Infiltration
        a = AOS_Infiltration(Soil, NewCond, Infl, Irr, IrrMngt, FieldMngt, FluxOut, DeepPerc, Runoff, GrowingSeason);
        NewCond = (InitCondStruct) a[0];
        DeepPerc = (double) a[1];
        Runoff = (double) a[2];
        Infl = (double) a[3];
        FluxOut = (double[]) a[4];

        //Capillary rise
        a = AOS_CapillaryRise(Soil, Groundwater, NewCond, FluxOut);
        NewCond = (InitCondStruct) a[0];
        double CR = (double) a[1];

        //Soil evaporation
        a = AOS_SoilEvaporation(AOS_ClockStruct, Soil, IrrMngt, FieldMngt, NewCond, Et0, Infl, P, Irr);
        NewCond = (InitCondStruct) a[0];
        double Es = (double) a[1];
        double EsPot = (double) a[2];




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

    //Function to get irrigation depth for current day
    private static Object[] AOS_Infiltration(Soil Soil, InitCondStruct InitCond, double Infl, double Irr,
                                             IrrMngtStruct IrrMngt, FieldMngtStruct FieldMngt, double[] FluxOut,
                                             double DeepPerc0, double Runoff0, boolean GrowingSeason) {
        //Store initial conditions in new structure for updating
        InitCondStruct NewCond = InitCond;
        double[] thnew = NewCond.th;

        //Update infiltration rate for irrigation
        //Note: irrigation amount adjusted for specified application efficiency
        if (GrowingSeason) {
            Infl = Infl + (Irr * (IrrMngt.AppEff / 100));
        }

        //Determine surface storage (if bunds are present)
        double RunoffIni = 0, ToStore = 0;
        if (FieldMngt.Bunds.compareTo("Y") == 0) {
            //TODO
        } else if (FieldMngt.Bunds.compareTo("N") == 0) {
            //No bunds on field
            if (Infl > Soil.layer.Ksat[1]) {
                //Infiltration limited by saturated hydraulic conductivity of top soil layer
                ToStore = Soil.layer.Ksat[1];
                //Additional water runs off
                RunoffIni = Infl - Soil.layer.Ksat[1];
            } else {
                //All water infiltrates
                ToStore = Infl;
                RunoffIni = 0;
            }

            //Update surface storage
            NewCond.SurfaceStorage = 0;
            //Add any water remaining behind bunds to surface runoff (needed for
            //days when bunds are removed to maintain water balance)
            RunoffIni = RunoffIni + InitCond.SurfaceStorage;
        }

        //Initialise counters
        int ii = 0;
        double Runoff = 0;

        //Infiltrate incoming water
        double DeepPerc = 0;
        if (ToStore > 0) {
            //TODO
        } else {
            //No infiltration
            DeepPerc = 0;
            Runoff = 0;
        }

        //Update total runoff
        Runoff = Runoff + RunoffIni;

        //Update surface storage (if bunds are present)
        if (Runoff > RunoffIni) {
            if (FieldMngt.Bunds.compareTo("Y") == 0) {
                if (FieldMngt.zBund > 0.001) {
                    //Increase surface storage
                    NewCond.SurfaceStorage = NewCond.SurfaceStorage + (Runoff - RunoffIni);
                    //Limit surface storage to bund height
                    if (NewCond.SurfaceStorage > (FieldMngt.zBund * 1000)) {
                        //Additional water above top of bunds becomes runoff
                        Runoff = RunoffIni + (NewCond.SurfaceStorage - (FieldMngt.zBund * 1000));
                        //Set surface storage to bund height
                        NewCond.SurfaceStorage = FieldMngt.zBund * 1000;
                    } else {
                        //No additional overtopping of bunds
                        Runoff = RunoffIni;
                    }
                }
            }
        }
        //Store updated water contents
        NewCond.th = thnew;

        //Update deep percolation, surface runoff, and infiltration values
        DeepPerc = DeepPerc + DeepPerc0;
        Infl = Infl - Runoff;
        double RunoffTot = Runoff + Runoff0;

        return new Object[]{NewCond, DeepPerc, RunoffTot, Infl, FluxOut};
    }

    //Function to calculate capillary rise from a shallow groundwater table
    private static Object[] AOS_CapillaryRise(Soil Soil, GwStruct Groundwater, InitCondStruct InitCond, double[] FluxOut) {
        //Store initial conditions for updating
        InitCondStruct NewCond = InitCond;
        double CrTot = 0;

        //Get groundwater table elevation on current day %%
        double zGW = NewCond.zGW;

        //Calculate capillary rise
        if (Groundwater.WaterTable == 0) {//No water table present
            //Capillary rise is zero
            CrTot = 0;
        } else if (Groundwater.WaterTable == 1) { //Water table present
            //TODO
        }

        return new Object[]{NewCond, CrTot};
    }

    //Function to calculate daily soil evaporation in AOS
    private static Object[] AOS_SoilEvaporation(ClockStruct AOS_ClockStruct, Soil Soil, IrrMngtStruct IrrMngt,
                                                FieldMngtStruct FieldMngt, InitCondStruct InitCond, double Et0,
                                                double Infl, double Rain, double Irr) {
        double EsAct = 0, EsPot = 0;

        //Store initial conditions in new structure that will be updated
        InitCondStruct NewCond = InitCond;

        //Initialise Wevap structur
        WevapStruct Wevap = new WevapStruct();

        //Prepare stage 2 evaporation (REW gone)
        //Only do this if it is first day of simulation, or if it is first day of
        //growing season and not simulating off-season
        if (AOS_ClockStruct.TimeStepCounter == 1 || (NewCond.DAP == 1 && AOS_ClockStruct.OffSeason.compareTo("N") == 0)) {
            //Reset storage in surface soil layer to zero
            NewCond.Wsurf = 0;
            //Set evaporation depth to minimum
            NewCond.EvapZ = Soil.EvapZmin;
            //Trigger stage 2 evaporation
            NewCond.Stage2 = true;
            //Get relative water content for start of stage 2 evaporation
            Wevap = AOS_EvapLayerWaterContent(NewCond, Soil, Wevap);
            NewCond.Wstage2 = (Wevap.Act - (Wevap.Fc - Soil.REW)) / (Wevap.Sat - (Wevap.Fc - Soil.REW));
            NewCond.Wstage2 = round((100 * NewCond.Wstage2)) / 100;
            if (NewCond.Wstage2 < 0) {
                NewCond.Wstage2 = 0;
            }
        }

        //Prepare soil evaporation stage 1
        //Adjust water in surface evaporation layer for any infiltration
        if (Rain > 0 || (Irr > 0 && IrrMngt.IrrMethod != 4)) {
            //Only prepare stage one when rainfall occurs, or when irrigation is trigerred (not in net irrigation mode)
            if (Infl > 0) {
                //Update storage in surface evaporation layer for incoming infiltration
                NewCond.Wsurf = Infl;
                //Water stored in surface evaporation layer cannot exceed REW
                if (NewCond.Wsurf > Soil.REW) {
                    NewCond.Wsurf = Soil.REW;
                }
                //Reset variables
                NewCond.Wstage2 = 0;
                NewCond.EvapZ = Soil.EvapZmin;
                NewCond.Stage2 = false;
            }
        }

        //Calculate potential soil evaporation rate (mm/day)
        EsPot = NewCond.Ke * Et0;

        //Adjust potential soil evaporation for mulches and/or partial wetting
        //Mulches
        double EsPotMul = 0;
        if (NewCond.SurfaceStorage < 0.000001) {
            if (FieldMngt.Mulches.compareTo("N") == 0) {
                //No mulches present
                EsPotMul = EsPot;
            } else if (FieldMngt.Mulches.compareTo("Y") == 0) {
                //Mulches present
                EsPotMul = EsPot * (1 - FieldMngt.fMulch * (FieldMngt.MulchPct / 100));
            }
        } else {
            //Surface is flooded - no adjustment of potential soil evaporation for mulches
            EsPotMul = EsPot;
        }

        //Partial surface wetting by irrigation
        double EsPotIrr;
        if (Irr > 0 && IrrMngt.IrrMethod != 4) {
            //Only apply adjustment if irrigation occurs and not in net irrigation mode
            if (Rain > 1 || NewCond.SurfaceStorage > 0) {
                //No adjustment for partial wetting - assume surface is fully wet
                EsPotIrr = EsPot;
            } else {
                //Adjust for proprtion of surface area wetted by irrigation
                EsPotIrr = EsPot * (IrrMngt.WetSurf / 100);
            }
        } else {
            //No adjustment for partial surface wetting
            EsPotIrr = EsPot;
        }

        //Assign minimum value (mulches and partial wetting don't combine)
        EsPot = min(EsPotIrr, EsPotMul);

        //Surface evaporation
        //Initialise actual evaporation counter
        EsAct = 0;
        //Evaporate surface storage
        if (NewCond.SurfaceStorage > 0) {
            if (NewCond.SurfaceStorage > EsPot) {
                //All potential soil evaporation can be supplied by surface storage
                EsAct = EsPot;
                //Update surface storage
                NewCond.SurfaceStorage = NewCond.SurfaceStorage - EsAct;
            } else {
                //Surface storage is not sufficient to meet all potential soil evaporation
                EsAct = NewCond.SurfaceStorage;
                //Update surface storage, evaporation layer depth, stage
                NewCond.SurfaceStorage = 0;
                NewCond.Wsurf = Soil.REW;
                NewCond.Wstage2 = 0;
                NewCond.EvapZ = Soil.EvapZmin;
                NewCond.Stage2 = false;
            }
        }

        //Stage 1 evaporation
        //Determine total water to be extracted
        double ToExtract = EsPot - EsAct;
        //Determine total water to be extracted in stage one (limited by surface layer water storage)
        double ExtractPotStg1 = min(ToExtract, NewCond.Wsurf);
        //Extract water
        if (ExtractPotStg1 > 0) {
            //TODO
        }

        //Stage 2 evaporation
        //Extract water
        if (ToExtract > 0.0001) {
            //Start stage 2
            NewCond.Stage2 = true;
            //Get sub-daily evaporative demand
            double Edt = ToExtract / (double) AOS_ClockStruct.EvapTimeSteps;
            //Loop sub-daily steps
            for (int jj = 0; jj < AOS_ClockStruct.EvapTimeSteps; jj++) {
                //Get current water storage (mm)
                Wevap = AOS_EvapLayerWaterContent(NewCond, Soil, Wevap);
                //Get water storage (mm) at start of stage 2 evaporation
                double Wupper = NewCond.Wstage2 * (Wevap.Sat - (Wevap.Fc - Soil.REW)) + (Wevap.Fc - Soil.REW);
                //Get water storage (mm) when there is no evaporation
                double Wlower = Wevap.Dry;
                //Get relative depletion of evaporation storage in stage 2
                double Wrel = (Wevap.Act - Wlower) / (Wupper - Wlower);
                //Check if need to expand evaporation layer
                double Wcheck;
                if (Soil.EvapZmax > Soil.EvapZmin) {
                    Wcheck = Soil.fWrelExp * ((Soil.EvapZmax - NewCond.EvapZ) / (Soil.EvapZmax - Soil.EvapZmin));
                    while (Wrel < Wcheck && NewCond.EvapZ < Soil.EvapZmax) {
                        //Expand evaporation layer by 1 mm
                        NewCond.EvapZ = NewCond.EvapZ + 0.001;
                        //Update water storage (mm) in evaporation layer
                        Wevap = AOS_EvapLayerWaterContent(NewCond, Soil, Wevap);
                        Wupper = NewCond.Wstage2 * (Wevap.Sat - (Wevap.Fc - Soil.REW)) + (Wevap.Fc - Soil.REW);
                        Wlower = Wevap.Dry;
                        //Update relative depletion of evaporation storage
                        Wrel = (Wevap.Act - Wlower) / (Wupper - Wlower);
                        Wcheck = Soil.fWrelExp * ((Soil.EvapZmax - NewCond.EvapZ) / (Soil.EvapZmax - Soil.EvapZmin));
                    }
                }
                //Get stage 2 evaporation reduction coefficient
                double Kr = (exp(Soil.fevap * Wrel) - 1) / (exp(Soil.fevap) - 1);
                if (Kr > 1) {
                    Kr = 1;
                }
                //Get water to extract (mm)
                double ToExtractStg2 = Kr * Edt;

                //Extract water from compartments
                int sum = 0;
                for (int i = 0; i < Soil.comp.dzsum.length; i++) {
                    if (Soil.comp.dzsum[i] < NewCond.EvapZ) {
                        sum++;
                    }
                }
                int comp_sto = sum + 1;
                int comp = 0;
                double factor = 0;
                while (ToExtractStg2 > 0 && comp < comp_sto) {
                    //Increment compartment counter
                    comp = comp + 1;
                    //Specify layer number
                    int layeri = Soil.comp.layer[comp];
                    //Determine proportion of compartment in evaporation layer
                    if (Soil.comp.dzsum[comp] > NewCond.EvapZ) {
                        factor = 1 - ((Soil.comp.dzsum[comp] - NewCond.EvapZ) / Soil.comp.dz[comp]);
                    } else {
                        factor = 1;
                    }
                    //Water storage (mm) at air dry
                    double Wdry = 1000 * Soil.layer.th_dry[layeri] * Soil.comp.dz[comp];
                    //Available water (mm)
                    double W = 1000 * NewCond.th[comp] * Soil.comp.dz[comp];
                    //Water available in compartment for extraction (mm)
                    double AvW = (W - Wdry) * factor;
                    if (AvW >= ToExtractStg2) {
                        //Update actual evaporation
                        EsAct = EsAct + ToExtractStg2;
                        //Update depth of water in current compartment
                        W = W - ToExtractStg2;
                        //Update total water to be extracted
                        ToExtract = ToExtract - ToExtractStg2;
                        //Update water to be extracted from surface layer (stage 1)
                        ToExtractStg2 = 0;
                    } else {
                        //Update actual evaporation
                        EsAct = EsAct + AvW;
                        //Update depth of water in current compartment
                        W = W - AvW;
                        //Update water to be extracted from surface layer (stage 1)
                        ToExtractStg2 = ToExtractStg2 - AvW;
                        //Update total water to be extracted
                        ToExtract = ToExtract - AvW;
                    }
                    //Update water content
                    NewCond.th[comp] = W / (1000 * Soil.comp.dz[comp]);
                }
            }
        }

        //Store potential evaporation for irrigation calculations on next day
        NewCond.Epot = EsPot;

        return new Object[]{NewCond, EsAct, EsPot};
    }

    //Function to get water contents in the evaporation layer
    private static WevapStruct AOS_EvapLayerWaterContent(InitCondStruct InitCond, Soil Soil, WevapStruct Wevap) {
        //Determine actual water content (mm)
        //Find soil compartments covered by evaporation layer
        int sum = 0;
        for (int i = 0; i < Soil.comp.dzsum.length; i++) {
            if (Soil.comp.dzsum[i] < InitCond.EvapZ) {
                sum++;
            }
        }
        int comp_sto = sum + 1;

        //Initialise variables
        Wevap.Act = 0;
        Wevap.Sat = 0;
        Wevap.Fc = 0;
        Wevap.Wp = 0;
        Wevap.Dry = 0;

        for (int ii = 0; ii < comp_sto; ii++) {
            //Specify layer number
            int layeri = Soil.comp.layer[ii];
            //Determine fraction of soil compartment covered by evaporation layer
            double factor = 0;
            if (Soil.comp.dzsum[ii] > InitCond.EvapZ) {
                factor = 1 - ((Soil.comp.dzsum[ii] - InitCond.EvapZ) / Soil.comp.dz[ii]);
            } else {
                factor = 1;
            }
            //Actual water storage in evaporation layer (mm)
            Wevap.Act = Wevap.Act + (factor * 1000 * InitCond.th[ii] * Soil.comp.dz[ii]);
            //Water storage in evaporation layer at saturation (mm)
            Wevap.Sat = Wevap.Sat + (factor * 1000 * Soil.layer.th_s[layeri] * Soil.comp.dz[ii]);
            //Water storage in evaporation layer at field capacity (mm)
            Wevap.Fc = Wevap.Fc + (factor * 1000 * Soil.layer.th_fc[layeri] * Soil.comp.dz[ii]);
            //Water storage in evaporation layer at permanent wilting point (mm)
            Wevap.Wp = Wevap.Wp + (factor * 1000 * Soil.layer.th_wp[layeri] * Soil.comp.dz[ii]);
            //Water storage in evaporation layer at air dry (mm)
            Wevap.Dry = Wevap.Dry + (factor * 1000 * Soil.layer.th_dry[layeri] * Soil.comp.dz[ii]);
        }

        if (Wevap.Act < 0) {
            Wevap.Act = 0;
        }

        return Wevap;
    }
}