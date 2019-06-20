package AOS_Solution;

import Structs.*;

import static java.lang.Math.*;

class Common {
    //Function to calculate actual and total available water in the root zone at current time step
    static Object[] AOS_RootZoneWater(Soil Soil, Crop Crop, InitCondStruct InitCond) {
        //Calculate root zone water content and available water
        //Compartments covered by the root zone
        double rootdepth = max(InitCond.Zroot, Crop.Zmin);
        rootdepth = round(rootdepth * 100.0) / 100.0;
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
            int layeri = Soil.comp.layer[ii] - 1;
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
        Dr Dr = new Dr();
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
            ztopdepth = round((ztopdepth * 100)) / 100.0;
            sum = 0;
            for (int i = 0; i < Soil.comp.dzsum.length; i++) {
                if (Soil.comp.dzsum[i] < ztopdepth) {
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
                int layeri = Soil.comp.layer[ii] - 1;
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

        return new Object[]{WrAct, Dr, TAW, thRZ};
    }

    //Function to calculate water stress coefficients
    static KswStruct AOS_WaterStress(Crop Crop, InitCondStruct InitCond, double Dr, double TAW, double Et0, boolean beta) {
        //Calculate relative root zone water depletion for each stress type
        //Number of stress variables
        int nstress = Crop.p_up.length;

        //Store stress thresholds
        double[] p_up = Crop.p_up.clone();
        double[] p_lo = Crop.p_lo.clone();
        if (Crop.ETadj == 1) {
            //Adjust stress thresholds for Et0 on current day (don't do this for
            //pollination water stress coefficient)
            for (int ii = 0; ii < 3; ii++) {
                p_up[ii] = p_up[ii] + (0.04 * (5 - Et0)) * (log10(10 - 9 * p_up[ii]));
                p_lo[ii] = p_lo[ii] + (0.04 * (5 - Et0)) * (log10(10 - 9 * p_lo[ii]));
            }
        }
        //Adjust senescence threshold if early sensescence is triggered
        if (beta && InitCond.tEarlySen > 0) {
            p_up[2] = p_up[2] * (1 - Crop.beta / 100);
        }

        //Limit values
        for (int i = 0; i < p_up.length; i++) {
            if (p_up[i] < 0) {
                p_up[i] = 0;
            } else if (p_up[i] > 1) {
                p_up[i] = 1;
            }
        }
        for (int i = 0; i < p_lo.length; i++) {
            if (p_lo[i] < 0) {
                p_lo[i] = 0;
            } else if (p_lo[i] > 1) {
                p_lo[i] = 1;
            }
        }

        //Calculate relative depletion
        double[] Drel = new double[nstress];
        for (int ii = 0; ii < nstress; ii++) {
            if (Dr <= (p_up[ii] * TAW)) {
                //No water stress
                Drel[ii] = 0;
            } else if ((Dr > (p_up[ii] * TAW)) && (Dr < (p_lo[ii] * TAW))) {
                //Partial water stress
                Drel[ii] = 1 - ((p_lo[ii] - (Dr / TAW)) / (p_lo[ii] - p_up[ii]));
            } else if (Dr >= (p_lo[ii] * TAW)) {
                //Full water stress
                Drel[ii] = 1;
            }
        }

        //Calculate root zone water stress coefficients
        double[] Ks = new double[3];
        for (int ii = 0; ii < 3; ii++) {
            Ks[ii] = 1 - ((exp(Drel[ii] * Crop.fshape_w[ii]) - 1) / (exp(Crop.fshape_w[ii]) - 1));
        }

        KswStruct Ksw = new KswStruct();
        //Water stress coefficient for leaf expansion
        Ksw.Exp = Ks[0];
        //Water stress coefficient for stomatal closure
        Ksw.Sto = Ks[1];
        //Water stress coefficient for senescence
        Ksw.Sen = Ks[2];
        //Water stress coefficient for pollination failure
        Ksw.Pol = 1 - Drel[3];
        //Mean water stress coefficient for stomatal closure
        Ksw.StoLin = 1 - Drel[1];
        return Ksw;
    }
}
