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

        //Pre-irrigation
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
                double CNbot = Math.round(1.4 * (Math.exp(-14 * Math.log(10))) + (0.507 * CN) -
                        Math.pow(0.00374 * CN, 2) + Math.pow(0.0000867 * CN, 3));
                double CNtop = Math.round(5.6 * (Math.exp(-14 * Math.log(10))) + (2.33 * CN) -
                        Math.pow(0.0209 * CN, 2) + Math.pow(0.000076 * CN, 3));
            }
        }

        return new Object[]{Runoff, Infl, NewCond};
    }
}