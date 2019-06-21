package AOS_Solution;

import Structs.*;

public class AOS_CropGrowthYieldForm {


    public static Double nthRoot(Double base, Double n) {
        return Math.pow(Math.E, Math.log(base) / n);
    }


    private static Object find(Double[] arr, double val) {

        int result = 0;

        for (int i = 0; i < arr.length; i++) {
            if (arr[i] >= val) {
                // adapted to matlab arrays
                return i + 1;
            }
        }
        return result;
    }


    private static Object[] UpdateCCxCDC(double CCprev, double CDC, double CCx, double dt) {

        Object[] returnStruct = new Object[2];
        //Get adjusted CCx
        double CCXadj = CCprev / (1 - 0.05 * (Math.exp(dt * ((CDC * 3.33) / (CCx + 2.29))) - 1));
        //Get adjusted CDC
        double CDCadj = CDC * ((CCXadj + 2.29) / (CCx + 2.29));

        returnStruct[0] = CCXadj;
        returnStruct[1] = CDCadj;

        return returnStruct;
    }


    private static double CCRequiredTime(double CCprev, double CCo, double CCx, double CGC,
                                         double CDC, double dt, double tSum, String mode) {

        double CGCx, tReq = 1;
        //Get CGC and/or time (GDD or CD) required to reach CC on previous day
        if (mode.equals("CGC")) {
            if (CCprev <= (CCx / 2)) {
                CGCx = (Math.log(CCprev / CCo)) / (tSum - dt);
            } else {
                CGCx = (Math.log((0.25 * CCx * CCx / CCo) / (CCx - CCprev))) / (tSum - dt);
            }
            tReq = (tSum - dt) * (CGCx / CGC);
        } else if (mode.equals("CDC")) {
            tReq = (Math.log(1 + (1 - CCprev / CCx) / 0.05)) / (CDC / CCx);
        }
        return tReq;
    }


    private static double CCDevelopment(double CCo, double CCx, double CGC,
                                        double CDC, double dt, String mode, double CCx0) {
        //Initialise output
        double CC = 0;  // in matlab its CC = [];

        //Calculate new canopy cover
        if (mode.equals("Growth")) {
            //Calculate canopy growth
            //Exponential growth stage
            CC = CCo * Math.exp(CGC * dt);
            if (CC > (CCx / 2)) {
                //Exponential decay stage
                CC = CCx - 0.25 * (CCx / CCo) * CCx * Math.exp(-CGC * dt);
            }
            //Limit CC to CCx
            if (CC > CCx) {
                CC = CCx;
            }
        } else if (mode.equals("Decline")) {
            //Calculate canopy decline
            if (CCx < 0.001) {
                CC = 0;
            } else {
                CC = CCx * (1 - 0.05 * (Math.exp(dt * CDC * 3.33 * ((CCx + 2.29) / (CCx0 + 2.29)) / (CCx + 2.29)) - 1));
            }
        }
        return CC;
    }


    private static double AdjustCCx(double CCprev, double CCo, double CCx, double CGC,
                                    double CDC, double dt, double tSum, Crop crop) {
        //Get time required to reach CC on previous day
        double tCCtmp = CCRequiredTime(CCprev, CCo, CCx, CGC, CDC, dt, tSum, "CGC");

        double CCxAdj = 0;
        //Determine CCx adjusted
        if (tCCtmp > 0) {
            tCCtmp = tCCtmp + (crop.CanopyDevEnd - tSum) + dt;
            CCxAdj = CCDevelopment(CCo, CCx, CGC, CDC, tCCtmp, "Growth", crop.CCx);
        }
        return CCxAdj;
    }


    private static InitCondStruct HIadjPostAnthesis(InitCondStruct InitCond, Crop crop, KswStruct ksw) {

        double dCor, DayCor, tmax2;
        //Store initial conditions in a structure for updating
        InitCondStruct NewCond = InitCond;

        //Calculate harvest index adjustment
        //1. Adjustment for leaf expansion
        double tmax1 = crop.CanopyDevEndCD - crop.HIstartCD;
        double DAP = NewCond.DAP - InitCond.DelayedCDs;
        if ((DAP <= (crop.CanopyDevEndCD + 1)) && (tmax1 > 0) && (NewCond.Fpre > 0.99) && (NewCond.CC > 0.001) && (crop.a_HI > 0)) {
            dCor = (1 + (1 - ksw.Exp) / crop.a_HI);
            NewCond.sCor1 = InitCond.sCor1 + (dCor / tmax1);
            DayCor = DAP - 1 - crop.HIstartCD;
            NewCond.fpost_upp = (tmax1 / DayCor) * NewCond.sCor1;
        }
        //2. Adjustment for stomatal closure
        tmax2 = crop.YldFormCD;
        DAP = NewCond.DAP - InitCond.DelayedCDs;

        if ((DAP <= (crop.HIendCD + 1)) && (tmax2 > 0) && (NewCond.Fpre > 0.99) && (NewCond.CC > 0.001) && (crop.b_HI > 0)) {
            dCor = (Math.exp(0.1 * Math.log(ksw.Sto))) * (1 - (1 - ksw.Sto) / crop.b_HI);
            NewCond.sCor2 = InitCond.sCor2 + (dCor / tmax2);
            DayCor = DAP - 1 - crop.HIstartCD;
            NewCond.fpost_dwn = (tmax2 / DayCor) * NewCond.sCor2;
        }
        //Determine total multiplier
        if ((tmax1 == 0) && (tmax2 == 0)) {
            NewCond.Fpost = 1;
        } else {
            if (tmax2 == 0) {
                NewCond.Fpost = NewCond.fpost_upp;
            } else {
                if (tmax1 == 0) {
                    NewCond.Fpost = NewCond.fpost_dwn;
                } else if (tmax1 <= tmax2) {
                    NewCond.Fpost = NewCond.fpost_dwn * (((tmax1 * NewCond.fpost_upp) + (tmax2 - tmax1)) / tmax2);
                } else {
                    NewCond.Fpost = NewCond.fpost_upp * (((tmax2 * NewCond.fpost_dwn) + (tmax1 - tmax2)) / tmax1);
                }
            }
        }
        return NewCond;
    }


    private static InitCondStruct HIadjPollination(InitCondStruct InitCond, Crop crop, KswStruct ksw,
                                                   Kst kst, double HIt) {

        double FracFlow = 1, t1, F1 = 1, t2, F2, t2Pct, t1Pct, F, dFpol, Ks;
        //Store initial conditions for updating
        InitCondStruct NewCond = InitCond;
        //Caclulate harvest index adjustment for pollination
        //Get fractional flowering
        if (HIt == 0) {
            //No flowering yet
            FracFlow = 0;
        } else if (HIt > 0) {
            //Fractional flowering on previous day
            t1 = HIt - 1;
            if (t1 == 0) {
                F1 = 0;
            } else {
                t1Pct = 100 * (t1 / crop.FloweringCD);
                if (t1Pct > 100) {
                    t1Pct = 100;
                }
            }
            if (F1 < 0) {
                F1 = 0;
            }
            //Fractional flowering on current day
            t2 = HIt;
            if (t2 == 0) {
                F2 = 0;
            } else {
                t2Pct = 100 * (t2 / crop.FloweringCD);
                if (t2Pct > 100) {
                    t2Pct = 100;
                }
                F2 = 0.00558 * Math.exp(0.63 * Math.log(t2Pct)) - (0.000969 * t2Pct) - 0.00383;
            }
            if (F2 < 0) {
                F2 = 0;
            }
            //Weight values
            if (Math.abs(F1 - F2) < 0.0000001) {
                F = 0;
            } else {
                F = 100 * ((F1 + F2) / 2) / crop.FloweringCD;
            }
            FracFlow = F;
        }
        //Calculate pollination adjustment for current day
        if (InitCond.CC < crop.CCmin) {
            //No pollination can occur as canopy cover is smaller than minimum threshold
            dFpol = 0;
        } else {
            Ks = Math.min(ksw.Pol, kst.PolC);
            Ks = Math.min(Ks, kst.PolH);
            dFpol = Ks * FracFlow * (1 + (crop.exc / 100));
        }

        //Calculate pollination adjustment to date
        NewCond.Fpol = InitCond.Fpol + dFpol;
        if (NewCond.Fpol > 1) {
            //Crop has fully pollinated
            NewCond.Fpol = 1;
        }
        return NewCond;
    }


    public static InitCondStruct HIadjPreAnthesis(InitCondStruct InitCond, Crop crop) {
        //Store initial conditions in structure for updating
        InitCondStruct NewCond = InitCond;

        //Calculate adjustment
        //Get parameters
        double Br = InitCond.B / InitCond.B_NS;
        double Br_range = Math.log(crop.dHI_pre) / 5.62;
        double Br_upp = 1;
        double Br_low = 1 - Br_range;
        double Br_top = Br_upp - (Br_range / 3);

        //Get biomass ratios
        double ratio_low = (Br - Br_low) / (Br_top - Br_low);
        double ratio_upp = (Br - Br_top) / (Br_upp - Br_top);

        //Calculate adjustment factor
        if ((Br >= Br_low) && (Br < Br_top)) {
            NewCond.Fpre = 1 + (((1 + Math.sin((1.5 - ratio_low) * Math.PI)) / 2) * (crop.dHI_pre / 100));
        } else if ((Br > Br_top) && (Br <= Br_upp)) {
            NewCond.Fpre = 1 + (((1 + Math.sin((0.5 + ratio_upp) * Math.PI)) / 2) * (crop.dHI_pre / 100));
        } else {
            NewCond.Fpre = 1;
        }
        if (NewCond.CC <= 0.01) {
            //No green canopy cover left at start of flowering so no harvestable crop will develop
            NewCond.Fpre = 0;
        }
        return NewCond;
    }


    private static Kst TemperatureStress(Crop crop, double Tmax, double Tmin) {

        //Calculate temperature stress coefficients affecting crop pollination
        //Get parameters for logistic curve
        double KsPol_up = 1;
        double KsPol_lo = 0.001;
        double Trel;

        Kst kst = new Kst();
        //Calculate effects of heat stress on pollination
        if (crop.PolHeatStress == 0) {
            //No heat stress effects on pollination
            kst.PolH = 1;
        } else if (crop.PolHeatStress == 1) {
            //Pollination affected by heat stress
            if (Tmax <= crop.Tmax_lo) {
                kst.PolH = 1;
            } else if (Tmax >= crop.Tmax_up) {
                kst.PolH = 0;
            } else {
                Trel = (Tmax - crop.Tmax_lo) / (crop.Tmax_up - crop.Tmax_lo);
                kst.PolH = (KsPol_up * KsPol_lo) / (KsPol_lo + (KsPol_up - KsPol_lo) * Math.exp(-crop.fshape_b * (1 - Trel)));
            }
        }
        //Calculate effects of cold stress on pollination
        if (crop.PolColdStress == 0) {
            //No cold stress effects on pollination
            kst.PolC = 1;
        } else if (crop.PolColdStress == 1) {
            //Pollination affected by cold stress
            if (Tmin >= crop.Tmin_up) {
                kst.PolC = 1;
            } else if (Tmin <= crop.Tmin_lo) {
                kst.PolC = 0;
            } else {
                Trel = (crop.Tmin_up - Tmin) / (crop.Tmin_up - crop.Tmin_lo);
                kst.PolC = (KsPol_up * KsPol_lo) / (KsPol_lo + (KsPol_up - KsPol_lo) * Math.exp(-crop.fshape_b * (1 - Trel)));
            }
        }
        return kst;

    }


    private static InitCondStruct CanopyCover(Crop crop, Soil soil, InitCondStruct InitCond,
                                              double GDD, double Et0, CO2 co2, boolean GrowingSeason) {

        double tmp_tCC, dtCC = 0, tCCadj = 0, CDCadj = 0, DCadj, CCsen, CCXadj, DAPadj, tAdj = 0, KeMax, KeMin, mult = 1, CCxActAdj;
        //Store initial conditions in a new structure for updating
        InitCondStruct NewCond = InitCond;
        NewCond.CCprev = InitCond.CC;
        double dr, taw;
        //Calculate canopy development (if in growing season)
        if (GrowingSeason) {
            //Calculate root zone water content
            //[~,Dr,TAW,~] = AOS_RootZoneWater(Soil,Crop,NewCond);
            Object[] rootZoneWaterRetObj = Common.AOS_RootZoneWater(soil, crop, NewCond);
            //Check whether to use root zone or top soil depletions for calculating water stress
            Dr DrTemp = (Dr) rootZoneWaterRetObj[1];
            TAW TAWtemp = (TAW) rootZoneWaterRetObj[2];
            if ((DrTemp.Rz / TAWtemp.Rz) <= (DrTemp.Zt / TAWtemp.Zt)) {
                //Root zone is wetter than top soil, so use root zone value
                dr = DrTemp.Rz;
                taw = TAWtemp.Rz;
            } else {
                //Top soil is wetter than root zone, so use top soil values
                dr = DrTemp.Zt;
                taw = TAWtemp.Zt;
            }
            //Determine if water stress is occurring
            boolean beta = true;
            KswStruct ksw = Common.AOS_WaterStress(crop, NewCond, dr, taw, Et0, beta);
            //Get canopy cover growth time

            if (crop.CalendarType == 1) {
                dtCC = 1;
                tCCadj = NewCond.DAP - NewCond.DelayedCDs;
            } else if (crop.CalendarType == 2) {
                dtCC = GDD;
                tCCadj = NewCond.GDDcum - NewCond.DelayedGDDs;
            }
            //Canopy development (potential)
            if ((tCCadj < crop.Emergence) || (Math.round(tCCadj) > crop.Maturity)) {
                //No canopy development before emergence/germination or after maturity
                NewCond.CC_NS = 0;
            } else if (tCCadj < crop.CanopyDevEnd) {
                //Canopy growth can occur
                if (InitCond.CC_NS <= crop.CC0) {
                    //Very small initial CC
                    NewCond.CC_NS = crop.CC0 * Math.exp(crop.CGC * dtCC);
                } else {
                    //Canopy growing
                    tmp_tCC = tCCadj - crop.Emergence;
                    NewCond.CC_NS = CCDevelopment(crop.CC0, 0.98 * crop.CCx, crop.CGC, crop.CDC, tmp_tCC, "Growth", crop.CCx);
                }
                //Update maximum canopy cover size in growing season
                NewCond.CCxAct_NS = NewCond.CC_NS;
            } else if (tCCadj > crop.CanopyDevEnd) {
                //No more canopy growth is possible or canopy in decline
                //Set CCx for calculation of withered canopy effects
                NewCond.CCxW_NS = NewCond.CCxAct_NS;
                if (tCCadj < crop.Senescence) {
                    //Mid-season stage - no canopy growth
                    NewCond.CC_NS = InitCond.CC_NS;
                    //Update maximum canopy cover size in growing season
                    NewCond.CCxAct_NS = NewCond.CC_NS;
                } else {
                    //Late-season stage - canopy decline
                    tmp_tCC = tCCadj - crop.Senescence;
                    NewCond.CC_NS = CCDevelopment(crop.CC0, NewCond.CCxAct_NS, crop.CGC, crop.CDC, tmp_tCC, "Decline", crop.CCx);
                }
            }
            //Canopy development (actual)
            if ((tCCadj < crop.Emergence) || (Math.round(tCCadj) > crop.Maturity)) {
                //No canopy development before emergence/germination or after maturity
                NewCond.CC = 0;
            } else if (tCCadj < crop.CanopyDevEnd) {
                //Canopy growth can occur
                if (InitCond.CC <= NewCond.CC0adj || ((InitCond.ProtectedSeed) && (InitCond.CC <= (1.25 * NewCond.CC0adj)))) {
                    //Very small initial CC or seedling in protected phase of
                    //growth. In this case, assume no leaf water expansion stress
                    if (InitCond.ProtectedSeed) {
                        tmp_tCC = tCCadj - crop.Emergence;
                        NewCond.CC = CCDevelopment(crop.CC0, crop.CCx, crop.CGC, crop.CDC, tmp_tCC, "Growth", crop.CCx);
                        //Check if seed protection should be turned off
                        if (NewCond.CC > (1.25 * NewCond.CC0adj)) {
                            //Turn off seed protection - lead expansion stress can occur on future time steps.
                            NewCond.ProtectedSeed = false;
                        }
                    } else {
                        NewCond.CC = NewCond.CC0adj * Math.exp(crop.CGC * dtCC);
                    }
                } else {
                    //Canopy growing
                    if ((InitCond.CC < (0.9799 * crop.CCx))) {
                        //Adjust canopy growth coefficient for leaf expansion water stress effects
                        double CGCadj = crop.CGC * ksw.Exp;
                        if (CGCadj > 0) {
                            //Adjust CCx for change in CGC
                            CCXadj = AdjustCCx(InitCond.CC, NewCond.CC0adj, crop.CCx, CGCadj, crop.CDC, dtCC, tCCadj, crop);
                            if (CCXadj < 0) {
                                NewCond.CC = InitCond.CC;
                            } else if (Math.abs(InitCond.CC - (0.9799 * crop.CCx)) < 0.001) {
                                //Approaching maximum canopy cover size
                                tmp_tCC = tCCadj - crop.Emergence;
                                NewCond.CC = CCDevelopment(crop.CC0, crop.CCx, crop.CGC, crop.CDC, tmp_tCC, "Growth", crop.CCx);
                            } else {
                                //Determine time required to reach CC on previous, day, given CGCAdj value
                                double tReq = CCRequiredTime(InitCond.CC, NewCond.CC0adj, CCXadj, CGCadj, crop.CDC, dtCC, tCCadj, "CGC");
                                if (tReq > 0) {
                                    //Calclate GDD's for canopy growth
                                    tmp_tCC = tReq + dtCC;
                                    //Determine new canopy size
                                    NewCond.CC = CCDevelopment(NewCond.CC0adj, CCXadj, CGCadj, crop.CDC, tmp_tCC, "Growth", crop.CCx);
                                } else {
                                    //No canopy growth
                                    NewCond.CC = InitCond.CC;
                                }
                            }
                        } else {
                            //No canopy growth
                            NewCond.CC = InitCond.CC;
                            //Update CC0
                            if (NewCond.CC > NewCond.CC0adj) {  //ended with ; in matlab
                                NewCond.CC0adj = crop.CC0;
                            } else {
                                NewCond.CC0adj = NewCond.CC;
                            }
                        }
                    } else {
                        //Canopy approaching maximum size
                        tmp_tCC = tCCadj - crop.Emergence;
                        NewCond.CC = CCDevelopment(crop.CC0, crop.CCx, crop.CGC, crop.CDC, tmp_tCC, "Growth", crop.CCx);
                        NewCond.CC0adj = crop.CC0;
                    }
                }
                if (NewCond.CC > InitCond.CCxAct) {
                    //Update actual maximum canopy cover size during growing season
                    NewCond.CCxAct = NewCond.CC;
                }
            } else if (tCCadj > crop.CanopyDevEnd) {
                //No more canopy growth is possible or canopy is in decline
                if (tCCadj < crop.Senescence) {
                    //Mid-season stage - no canopy growth
                    NewCond.CC = InitCond.CC;
                    if (NewCond.CC > InitCond.CCxAct) {
                        //Update actual maximum canopy cover size during growing season
                        NewCond.CCxAct = NewCond.CC;
                    }
                } else {
                    //Late-season stage - canopy decline
                    //Adjust canopy decline coefficient for difference between actual
                    //and potential CCx
                    CDCadj = crop.CDC * ((NewCond.CCxAct + 2.29) / (crop.CCx + 2.29));
                    //Determine new canopy size
                    tmp_tCC = tCCadj - crop.Senescence;
                    NewCond.CC = CCDevelopment(NewCond.CC0adj, NewCond.CCxAct, crop.CGC, CDCadj, tmp_tCC, "Decline", crop.CCx);
                }
                //Check for crop growth termination
                if ((NewCond.CC < 0.001) && (!InitCond.CropDead)) {
                    //Crop has died
                    NewCond.CC = 0;
                    NewCond.CropDead = true;
                }
            }
            //Canopy senescence due to water stress (actual)
            if (tCCadj >= crop.Emergence) {
                if ((tCCadj < crop.Senescence) || (InitCond.tEarlySen > 0)) {
                    //Check for early canopy senescence  due to severe water stress.
                    if ((ksw.Sen < 1) && (!InitCond.ProtectedSeed)) {
                        //Early canopy senescence
                        NewCond.PrematSenes = true;
                        if (InitCond.tEarlySen == 0) {
                            //No prior early senescence
                            NewCond.CCxEarlySen = InitCond.CC;
                        }
                        //Increment early senescence GDD counter
                        NewCond.tEarlySen = InitCond.tEarlySen + dtCC;
                        //Adjust canopy decline coefficient for water stress
                        beta = false;
                        ksw = Common.AOS_WaterStress(crop, NewCond, dr, taw, Et0, beta);
                        if (ksw.Sen > 0.99999) {
                            DCadj = 0.0001;
                        } else {
                            CDCadj = (1 - Math.pow(ksw.Sen, 8)) * crop.CDC;
                        }
                        //Get new canpy cover size after senescence
                        if (NewCond.CCxEarlySen < 0.001) {
                            CCsen = 0;
                        } else {
                            //Get time required to reach CC at end of previous day, given CDCadj
                            double tReq = (Math.log(1 + (1 - InitCond.CC / NewCond.CCxEarlySen) / 0.05)) / ((CDCadj * 3.33) / (NewCond.CCxEarlySen + 2.29));
                            //Calculate GDD's for canopy decline
                            tmp_tCC = tReq + dtCC;
                            //Determine new canopy size
                            CCsen = NewCond.CCxEarlySen * (1 - 0.05 * (Math.exp(tmp_tCC * ((CDCadj * 3.33) / (NewCond.CCxEarlySen + 2.29))) - 1));
                            if (CCsen < 0) {
                                CCsen = 0;
                            }
                        }
                        //Update canopy cover size
                        if (tCCadj < crop.Senescence) {
                            //Limit CC to CCx
                            if (CCsen > crop.CCx) {
                                CCsen = crop.CCx;
                            }
                            //CC cannot be greater than value on previous day
                            NewCond.CC = CCsen;
                            if (NewCond.CC > InitCond.CC) {
                                NewCond.CC = InitCond.CC;
                            }
                            //Update maximum canopy cover size during growing season
                            NewCond.CCxAct = NewCond.CC;
                            //Update CC0 if current CC is less than initial canopy cover size at planting
                            if (NewCond.CC < crop.CC0) {
                                NewCond.CC0adj = NewCond.CC;
                            } else {
                                NewCond.CC0adj = crop.CC0;
                            }
                        } else {
                            //Update CC to account for canopy cover senescence due to water stress
                            if (CCsen < NewCond.CC) {
                                NewCond.CC = CCsen;
                            }
                        }
                        //Check for crop growth termination
                        if ((NewCond.CC < 0.001) && (!InitCond.CropDead)) {
                            //Crop has died
                            NewCond.CC = 0;
                            NewCond.CropDead = true;
                        }
                    } else {
                        //No water stress
                        NewCond.PrematSenes = false;
                        if ((tCCadj > crop.Senescence) && (InitCond.tEarlySen > 0)) {
                            //Rewatering of canopy in late season Get new values for CCx and CDC
                            tmp_tCC = tCCadj - dtCC - crop.Senescence;
                            //[CCXadj,CDCadj] = UpdateCCxCDC(InitCond.CC, crop.CDC, crop.CCx, tmp_tCC);
                            Object[] updateCCxCDCstruct = UpdateCCxCDC(InitCond.CC, crop.CDC, crop.CCx, tmp_tCC);
                            CCXadj = (double) updateCCxCDCstruct[0];
                            CDCadj = (double) updateCCxCDCstruct[2];
                            //Get new CC value for end of current day
                            tmp_tCC = tCCadj - crop.Senescence;
                            NewCond.CC = CCDevelopment(NewCond.CC0adj, CCXadj, crop.CGC, CDCadj, tmp_tCC, "Decline", crop.CCx);
                            //Check for crop growth termination
                            if ((NewCond.CC < 0.001) && (!InitCond.CropDead)) {
                                NewCond.CC = 0;
                                NewCond.CropDead = true;
                            }
                        }
                        //Reset early senescence counter
                        NewCond.tEarlySen = 0;
                    }
                    //Adjust CCx for effects of withered canopy
                    if (NewCond.CC > InitCond.CCxW) {
                        NewCond.CCxW = NewCond.CC;
                    }
                }
            }
            //Calculate canopy size adjusted for micro-advective effects
            //Check to ensure potential CC is not slightly lower than actual
            if (NewCond.CC_NS < NewCond.CC) {
                NewCond.CC_NS = NewCond.CC;
                if (tCCadj < crop.CanopyDevEnd) {
                    NewCond.CCxAct_NS = NewCond.CC_NS;
                }
            }
            //Actual (with water stress)
            NewCond.CCadj = (1.72 * NewCond.CC) - Math.pow(NewCond.CC, 2) + (0.3 * Math.pow(NewCond.CC, 3));
            //Potential (without water stress)
            NewCond.CCadj_NS = (1.72 * NewCond.CC_NS) - Math.pow(NewCond.CC_NS, 2) + (0.3 * Math.pow(NewCond.CC_NS, 3));

            //Update parameters for transpiration calculation
            //1. Kcb with no prior water stress
            //Update ageing days counter
            DAPadj = NewCond.DAP - NewCond.DelayedCDs;
            if (DAPadj > crop.MaxCanopyCD) {
                NewCond.AgeDays_NS = DAPadj - crop.MaxCanopyCD;
            }
            //Update crop coefficient for ageing of canopy
            if (NewCond.AgeDays_NS > 5) {
                NewCond.Kcb_NS = crop.Kcb - ((NewCond.AgeDays_NS - 5) * (crop.fage / 100)) * NewCond.CCxW_NS;
            } else {
                NewCond.Kcb_NS = crop.Kcb;
            }
            //Update crop coefficient for CO2 concentration
            if (co2.CurrentConc > co2.RefConc) {
                NewCond.Kcb_NS = NewCond.Kcb_NS * (1 - 0.05 * ((co2.CurrentConc - co2.RefConc) / (550 - co2.RefConc)));
            }
            //Correct crop coefficient for dying green canopy effects
            if (NewCond.CC_NS < NewCond.CCxW_NS) {
                if ((NewCond.CCxW_NS > 0.001) && (NewCond.CC_NS > 0.001)) {
                    NewCond.Kcb_NS = NewCond.Kcb_NS * Math.pow((NewCond.CC_NS / NewCond.CCxW_NS), crop.a_Tr);
                }
            }
            //Adjust for current potential canopy size
            NewCond.Kcb_NS = NewCond.Kcb_NS * (NewCond.CCadj_NS);

            //2. Kcb with potential prior water stress and/or delayed development
            //Update ageing days counter
            DAPadj = NewCond.DAP - NewCond.DelayedCDs - 1;
            if (DAPadj > crop.MaxCanopyCD) {
                NewCond.AgeDays = DAPadj - crop.MaxCanopyCD;
            }
            //Update crop coefficient for ageing of canopy
            if (NewCond.AgeDays > 5) {
                NewCond.Kcb = crop.Kcb - ((NewCond.AgeDays - 5) * (crop.fage / 100)) * NewCond.CCxW;
            } else {
                NewCond.Kcb = crop.Kcb;
            }
            //Update crop coefficient for CO2 concentration
            if (co2.CurrentConc > co2.RefConc) {
                NewCond.Kcb = NewCond.Kcb * (1 - 0.05 * ((co2.CurrentConc - co2.RefConc) / (550 - co2.RefConc)));
            }
            //Correct crop coefficient for dying green canopy effects
            if (NewCond.CC < NewCond.CCxW) {
                if ((NewCond.CCxW > 0.001) && (NewCond.CC > 0.001)) {
                    NewCond.Kcb = NewCond.Kcb * Math.pow((NewCond.CC / NewCond.CCxW), crop.a_Tr);
                }
            }
            //Adjust for current canopy size
            NewCond.Kcb = NewCond.Kcb * (NewCond.CCadj);

            //Update parameters for soil evaporation coefficient
            //Adjust time for any delayed development
            if (crop.CalendarType == 1) {
                tAdj = NewCond.DAP - NewCond.DelayedCDs;
            } else if (crop.CalendarType == 2) {
                tAdj = NewCond.GDDcum - NewCond.DelayedGDDs;
            }
            //Calculate maximum potential soil evaporation coefficient
            KeMax = soil.Kex * (1 - NewCond.CCxW * (soil.fwcc / 100));
            //Calculate actual soil evaporation coefficient (given current canopy cover size)
            NewCond.Ke = soil.Kex * (1 - NewCond.CCadj);
            //Adjust soil evaporation coefficient for effects of withered canopy
            if ((tAdj > crop.Senescence) && (NewCond.CCxAct > 0)) {
                if (NewCond.CC > (NewCond.CCxAct / 2)) {
                    if (NewCond.CC > NewCond.CCxAct) {
                        mult = 0;
                    } else {
                        mult = (NewCond.CCxAct - NewCond.CC) / (NewCond.CCxAct / 2);
                    }
                } else {
                    mult = 1;
                }
                NewCond.Ke = NewCond.Ke * (1 - NewCond.CCxAct * (soil.fwcc / 100) * mult);
                CCxActAdj = (1.72 * NewCond.CCxAct) - Math.pow(NewCond.CCxAct, 2) + 0.3 * Math.pow(NewCond.CCxAct, 3);
                KeMin = soil.Kex * (1 - CCxActAdj);
                if (KeMin < 0) {
                    KeMin = 0;
                }
                if (NewCond.Ke < KeMin) {
                    NewCond.Ke = KeMin;
                } else if (NewCond.Ke > KeMax) {
                    NewCond.Ke = KeMax;
                }
            }
            if (NewCond.PrematSenes) {
                if (NewCond.Ke > KeMax) {
                    NewCond.Ke = KeMax;
                }
            }
        } else {
            //No canopy outside growing season - set various values to zero
            NewCond.CC = 0;
            NewCond.CCadj = 0;
            NewCond.CC_NS = 0;
            NewCond.CCadj_NS = 0;
            NewCond.CCxW = 0;
            NewCond.CCxAct = 0;
            NewCond.CCxW_NS = 0;
            NewCond.CCxAct_NS = 0;
            NewCond.Kcb = 0;
            NewCond.Kcb_NS = 0;
            NewCond.Ke = soil.Kex;
        }
        return NewCond;
    }


    private static InitCondStruct HarvestIndex(Soil soil, Crop crop, InitCondStruct InitCond,
                                               double Et0, double Tmax, double Tmin, boolean GrowingSeason) {
        //Store initial conditions for updating
        InitCondStruct NewCond = InitCond;
        //Calculate harvest index build up (if in growing season)
        double HIadj = 0, HImax, taw, dr;
        if (GrowingSeason) {
            //Calculate root zone water content

            //[~,Dr,TAW,~] = Common.AOS_RootZoneWater(soil, crop, NewCond);           // Object []
            Object[] rootZoneWaterRetObj = Common.AOS_RootZoneWater(soil, crop, NewCond);
            //Check whether to use root zone or top soil depletions for calculating water stress
            Dr Dr = (Dr) rootZoneWaterRetObj[1];
            TAW TAW = (TAW) rootZoneWaterRetObj[2];
            if ((Dr.Rz / TAW.Rz) <= (Dr.Zt / TAW.Zt)) {
                //Root zone is wetter than top soil, so use root zone value
                dr = Dr.Rz;
                taw = TAW.Rz;
            } else {
                //Top soil is wetter than root zone, so use top soil values
                dr = Dr.Zt;
                taw = TAW.Zt;
            }

            //Calculate water stress
            boolean beta = true;
            KswStruct ksw = Common.AOS_WaterStress(crop, NewCond, dr, taw, Et0, beta);

            //Calculate temperature stress
            Kst kst = TemperatureStress(crop, Tmax, Tmin);

            //Get reference harvest index on current day
            double HIi = NewCond.HIref;
            //Get time for harvest index build-up
            double HIt = NewCond.DAP - NewCond.DelayedCDs - crop.HIstartCD - 1;
            //Calculate harvest index
            if ((NewCond.YieldForm) && (HIt >= 0)) {
                //Root/tuber or fruit/grain crops
                if ((crop.CropType == 2) || (crop.CropType == 3)) {
                    //Detemine adjustment for water stress before anthesis
                    if (!InitCond.PreAdj) {
                        InitCond.PreAdj = true;

                        NewCond = HIadjPreAnthesis(NewCond, crop);

                    }

                    //Determine adjustment for crop pollination failure
                    if (crop.CropType == 3) {   //Adjustment only for fruit/grain crops
                        if ((HIt > 0) && (HIt <= crop.FloweringCD)) {
                            NewCond = HIadjPollination(InitCond, crop, ksw, kst, HIt);
                        }
                        HImax = NewCond.Fpol * crop.HI0;
                    } else {
                        //No pollination adjustment for root/tuber crops
                        HImax = crop.HI0;
                    }

                    //Determine adjustments for post-anthesis water stress
                    if (HIt > 0) {
                        NewCond = HIadjPostAnthesis(NewCond, crop, ksw);
                    }
                    //Limit HI to maximum allowable increase due to pre- and
                    //post-anthesis water stress combinations
                    double HImult = NewCond.Fpre * NewCond.Fpost;
                    if (HImult > 1 + (crop.dHI0 / 100)) {
                        HImult = 1 + (crop.dHI0 / 100);
                    }
                    //Determine harvest index on current day, adjusted for stress effects
                    if (HImax >= HIi) {
                        HIadj = HImult * HIi;
                    } else {
                        HIadj = HImult * HImax;
                    }
                } else if (crop.CropType == 1) {
                    //Leafy vegetable crops - no adjustment, harvest index equal to reference value for current day
                    HIadj = HIi;
                }
            } else {
                //No build-up of harvest index if outside yield formation period
                HIi = InitCond.HI;
                HIadj = InitCond.HIadj;
            }
            //Store final values for current time step
            NewCond.HI = HIi;
            NewCond.HIadj = HIadj;
        } else {
            //No harvestable crop outside of a growing season
            NewCond.HI = 0;
            NewCond.HIadj = 0;
        }
        return NewCond;
    }


    private static InitCondStruct BiomassAccumulation(Crop crop, InitCondStruct InitCond, double Tr,
                                                      double TrPot, double Et0, boolean GrowingSeason) {
        //Store initial conditions in a new structure for updating
        InitCondStruct NewCond = InitCond;
        //Calculate biomass accumulation (if in growing season)
        if (GrowingSeason) {
            double WPadj = 0;
            //Get time for harvest index build-up
            int HIt = NewCond.DAP - NewCond.DelayedCDs - crop.HIstartCD - 1;

            if (((crop.CropType == 2) || (crop.CropType == 3)) && (NewCond.HIref > 0)) {
                //Adjust WP for reproductive stage
                double fswitch;
                if (crop.Determinant == 1) {
                    fswitch = NewCond.PctLagPhase / 100;
                } else {
                    if (HIt < (crop.YldFormCD / 3)) {
                        fswitch = HIt / (crop.YldFormCD / 3);
                    } else {
                        fswitch = 1;
                    }
                }
                WPadj = crop.WP * (1 - (1 - crop.WPy / 100) * fswitch);
            } else {
                WPadj = crop.WP;
            }
            //Adjust WP for CO2 effects
            WPadj = WPadj * crop.fCO2;
            //Calculate biomass accumulation on current day No water stress
            double dB_NS = WPadj * (TrPot / Et0);
            //With water stress
            Double dB = new Double(WPadj * (Tr / Et0));
            if (dB.isNaN()) {
                dB = new Double(0);
            }
            //Update biomass accumulation
            NewCond.B = NewCond.B + dB;
            NewCond.B_NS = NewCond.B_NS + dB_NS;
        } else {
            //No biomass accumulation outside of growing season
            NewCond.B = 0;
            NewCond.B_NS = 0;
        }
        return NewCond;
    }


    /**
     * Function to calculate reference (no adjustment for stress effects)
     * harvest index on current day
     *
     * @param InitCond
     * @param crop
     * @param GrowingSeason
     * @return
     */
    private static InitCondStruct HIrefCurrentDay(InitCondStruct InitCond, Crop crop, boolean GrowingSeason) {

        //Store initial conditions for updating
        InitCondStruct NewCond = InitCond;
        //Calculate reference harvest index (if in growing season)
        if (GrowingSeason) {
            //Check if in yield formation period
            double tAdj = NewCond.DAP - NewCond.DelayedCDs;
            NewCond.YieldForm = false;
            if (tAdj > crop.HIstartCD) {
                NewCond.YieldForm = true;
            }

            //Get time for harvest index calculation
            double HIt = NewCond.DAP - NewCond.DelayedCDs - crop.HIstartCD - 1;

            if (HIt <= 0) {
                //Yet to reach time for HI build-up
                NewCond.HIref = 0;
                NewCond.PctLagPhase = 0;
            } else {
                if (NewCond.CCprev <= (crop.CCmin * crop.CCx)) {
                    //HI cannot develop further as canopy cover is too small
                    NewCond.HIref = InitCond.HIref;
                } else {
                    //Check crop type
                    if ((crop.CropType == 1) || (crop.CropType == 2)) { // in matlab this line ended with ; why?
                        //If crop type is leafy vegetable or root/tuber, then proceed with
                        //logistic growth (i.e. no linear switch)
                        NewCond.PctLagPhase = 100; // No lag phase
                        //Calculate reference harvest index for current day
                        NewCond.HIref = (crop.HIini * crop.HI0) / (crop.HIini + (crop.HI0 - crop.HIini) * Math.exp(-crop.HIGC * HIt));
                        //Harvest index apprAOShing maximum limit
                        if (NewCond.HIref >= (0.9799 * crop.HI0)) {
                            NewCond.HIref = crop.HI0;
                        }
                    } else if (crop.CropType == 3) {
                        //If crop type is fruit/grain producing, check for linear switch
                        if (HIt < crop.tLinSwitch) {
                            //Not yet reached linear switch point, therefore proceed with logistic build-up
                            NewCond.PctLagPhase = 100 * (HIt / crop.tLinSwitch);
                            //Calculate reference harvest index for current day (logistic build-up)
                            NewCond.HIref = (crop.HIini * crop.HI0) / (crop.HIini + (crop.HI0 - crop.HIini) * Math.exp(-crop.HIGC * HIt));
                        } else {
                            //Linear switch point has been reached
                            NewCond.PctLagPhase = 100;
                            //Calculate reference harvest index for current day (logistic portion)
                            NewCond.HIref = (crop.HIini * crop.HI0) / (crop.HIini + (crop.HI0 - crop.HIini) * Math.exp(-crop.HIGC * crop.tLinSwitch));
                            //Calculate reference harvest index for current day (total - logistic portion + linear portion)
                            NewCond.HIref = NewCond.HIref + (crop.dHILinear * (HIt - crop.tLinSwitch));
                        }
                    }
                    //Limit HIref and round off computed value
                    if (NewCond.HIref > crop.HI0) {
                        NewCond.HIref = crop.HI0;
                    } else if (NewCond.HIref <= (crop.HIini + 0.004)) {
                        NewCond.HIref = 0;
                    } else if (((crop.HI0 - NewCond.HIref) < 0.004)) {
                        NewCond.HIref = crop.HI0;
                    }
                }
            }
        } else {
            //Reference harvest index is zero outside of growing season
            NewCond.HIref = 0;
        }
        return NewCond;
    }


    /**
     * Function to calculate number of growing degree days on current day
     *
     * @param crop
     * @param InitCond
     * @param GrowingSeason
     * @return
     */
    private static InitCondStruct GrowthStage(Crop crop, InitCondStruct InitCond
            , boolean GrowingSeason) {

        //Store initial conditions in new structure for updating
        InitCondStruct NewCond = InitCond;
        //Get growth stage (if in growing season
        if (GrowingSeason) {
            double tAdj = 0;
            //Adjust time for any delayed growth
            if (crop.CalendarType == 1) {
                tAdj = NewCond.DAP - NewCond.DelayedCDs;
            } else if (crop.CalendarType == 2) {
                tAdj = NewCond.GDDcum - NewCond.DelayedGDDs;
            }
            //Update growth stage
            if (tAdj <= crop.Canopy10Pct) {
                NewCond.GrowthStage = 1;
            } else if (tAdj <= crop.MaxCanopy) {
                NewCond.GrowthStage = 2;
            } else if (tAdj <= crop.Senescence) {
                NewCond.GrowthStage = 3;
            } else if (tAdj > crop.Senescence) {
                NewCond.GrowthStage = 4;
            }
        } else {
            //Not in growing season so growth stage is set to dummy value
            NewCond.GrowthStage = 0;
        }
        return NewCond;
    }


    /**
     * Function to calculate root zone expansio
     *
     * @param crop
     * @param soil
     * @param GrowingSeason
     * @return
     */
    private static InitCondStruct RootDevelopment(Crop crop, Soil soil, GwStruct Groundwater, InitCondStruct InitCond,
                                                  double GDD, boolean GrowingSeason) {
        //Store initial conditions for updating
        InitCondStruct NewCond = InitCond;

        //Calculate root expansion (if in growing season)
        if (GrowingSeason) {
            //If today is first day of season, root depth is equal to minimum depth
            if (NewCond.DAP == 1) {
                InitCond.Zroot = crop.Zmin;
            }
            //Adjust time for any delayed development
            double tAdj = 0;    // initialized with 0 to prevent errors
            if (crop.CalendarType == 1) {
                tAdj = NewCond.DAP - NewCond.DelayedCDs;
            } else if (crop.CalendarType == 2) {
                tAdj = NewCond.GDDcum - NewCond.DelayedGDDs;
            }
            //Calculate root expansion
            double Zini = crop.Zmin * (crop.PctZmin / 100);
            long t0 = Math.round((crop.Emergence / 2));
            double tmax = crop.MaxRooting;
            double tOld = 0;    // initialized with 0 to prevent errors
            if (crop.CalendarType == 1) {
                tOld = tAdj - 1;
            } else if (crop.CalendarType == 2) {
                tOld = tAdj - GDD;
            }
            //Potential root depth on previous day
            double ZrOld = 0;
            double X;
            if (tOld >= tmax) {
                ZrOld = crop.Zmax;
            } else if (tOld <= t0) {
                ZrOld = Zini;
            } else {
                X = (tOld - t0) / (tmax - t0);
                ZrOld = Zini + (crop.Zmax - Zini) * nthRoot(X, crop.fshape_r);
            }
            if (ZrOld < crop.Zmin) {
                ZrOld = crop.Zmin;
            }
            //Potential root depth on current day
            double Zr;
            if (tAdj >= tmax) {
                Zr = crop.Zmax;
            } else if (tAdj <= t0) {
                Zr = Zini;
            } else {
                X = (tAdj - t0) / (tmax - t0);
                Zr = Zini + (crop.Zmax - Zini) * nthRoot(X, crop.fshape_r);
            }
            if (Zr < crop.Zmin) {
                Zr = crop.Zmin;
            }
            //Store Zr as potential value
            double ZrPot = Zr;
            //Determine rate of change
            double dZr = Zr - ZrOld;

            //Adjust expansion rate for presence of restrictive soil horizons
            if (Zr > crop.Zmin) {
                int layeri = 1 - 1;
                double Zsoil = soil.layer.dz[layeri];
                while ((Zsoil <= crop.Zmin) && (layeri < soil.nLayer)) {
                    layeri = (layeri + 1) - 1;
                    Zsoil = Zsoil + soil.layer.dz[layeri];
                }
                double ZrAdj = crop.Zmin;
                double ZrRemain = Zr - crop.Zmin;
                double deltaZ = Zsoil - crop.Zmin;
                double ZrTest, ZrOUT = 0;
                boolean EndProf = false;
                while (!EndProf) {
                    ZrTest = ZrAdj + (ZrRemain * (soil.layer.Penetrability[layeri] / 100));
                    if (((layeri + 1) == soil.nLayer) || (soil.layer.Penetrability[layeri] == 0) || (ZrTest <= Zsoil)) {
                        ZrOUT = ZrTest;
                        EndProf = true;
                    } else {
                        ZrAdj = Zsoil;
                        ZrRemain = ZrRemain - (deltaZ / (soil.layer.Penetrability[layeri] / 100.0));
                        layeri = (layeri + 1) - 1;
                        Zsoil = Zsoil + soil.layer.dz[layeri];
                        deltaZ = soil.layer.dz[layeri];
                    }
                }
                //Correct Zr and dZr for effects of restrictive horizons
                Zr = ZrOUT;
                dZr = Zr - ZrOld;
            }
            //Adjust rate of expansion for any stomatal water stress
            if (NewCond.TrRatio < 0.9999) {
                if (crop.fshape_ex >= 0) {
                    dZr = dZr * NewCond.TrRatio;
                } else {
                    double fAdj = (Math.exp(NewCond.TrRatio * crop.fshape_ex) - 1) / (Math.exp(crop.fshape_ex) - 1);
                    dZr = dZr * fAdj;
                }
            }
            //Adjust rate of root expansion for dry soil at expansion front
            if (dZr > 0.001) {
                //Define water stress threshold for inhibition of root expansion
                double pZexp = crop.p_up[1] + ((1 - crop.p_up[1]) / 2.0);
                //Define potential new root depth
                double ZiTmp = InitCond.Zroot + dZr;
                //Find compartment that root zone will expand in to
                int compi = (int) find(soil.comp.dzsum, Math.round(ZiTmp)); //#########
                //Get TAW in compartment
                int layeri = soil.comp.layer[compi - 1] - 1;
                double TAWcompi = (soil.layer.th_fc[layeri] - soil.layer.th_wp[layeri]);
                //Define stress threshold
                double thThr = soil.layer.th_fc[layeri] - (pZexp * TAWcompi);
                //Check for stress conditions
                if (NewCond.th[compi - 1] < thThr) {
                    //Root expansion limited by water content at expansion front
                    if (NewCond.th[compi - 1] <= soil.layer.th_wp[layeri]) {
                        //Expansion fully inhibited
                        dZr = 0;
                    } else {
                        //Expansion partially inhibited
                        double Wrel = (soil.layer.th_fc[layeri] - NewCond.th[compi - 1]) / TAWcompi;
                        double Drel = 1 - ((1 - Wrel) / (1 - pZexp));
                        double Ks = 1 - ((Math.exp(Drel * crop.fshape_w[1]) - 1) / (Math.exp(crop.fshape_w[1]) - 1));
                        dZr = dZr * Ks;
                    }
                }
            }
            //Adjust for early senescence
            if ((NewCond.CC <= 0) && (NewCond.CC_NS > 0.5)) {
                dZr = 0;
            }
            //Adjust root expansion for failure to germinate (roots cannot expand
            //if crop has not germinated)
            if (!InitCond.Germination) {
                dZr = 0;
            }
            //Get new rooting depth
            NewCond.Zroot = InitCond.Zroot + dZr;

            //Adjust root density if deepening is restricted due to dry subsoil
            //and/or restrictive layers
            if (NewCond.Zroot < ZrPot) {
                NewCond.rCor = (2 * (ZrPot / NewCond.Zroot) * ((crop.SxTop + crop.SxBot) / 2) - crop.SxTop) / crop.SxBot;
                if (NewCond.Tpot > 0) {
                    NewCond.rCor = NewCond.rCor * NewCond.TrRatio;
                    if (NewCond.rCor < 1) {
                        NewCond.rCor = 1;
                    }
                }
            } else {
                NewCond.rCor = 1;
            }
            //Limit rooting depth if groundwater table is present (roots cannot develop below the water table)
            if ((Groundwater.WaterTable == 1) && (NewCond.zGW > 0)) {
                if (NewCond.Zroot > NewCond.zGW) {
                    NewCond.Zroot = NewCond.zGW;
                    if (NewCond.Zroot < crop.Zmin) {
                        NewCond.Zroot = crop.Zmin;
                    }
                }
            }
        } else {
            // No root system outside of the growing season
            NewCond.Zroot = 0;
        }
        return NewCond;
    }


    public static InitCondStruct run(Crop Crop, Soil Soil, Weather Weather, GwStruct Groundwater,
                                     InitCondStruct InitCond, SoilWatOutStruct SoilWatOut, boolean GrowingSeason, CO2 CO2) {
        //Unpack weather structure
        double Et0 = Weather.RefET;
        double GDD = Weather.GDD;
        double Tmin = Weather.Tmin;
        double Tmax = Weather.Tmax;

        //Unpack soil water balance output structure
        double Tr = SoilWatOut.Tr;
        double TrPot = SoilWatOut.TrPot_NS;


        //Store initial conditions for updating %%
        InitCondStruct NewCond = InitCond; //###########################
        //Root development
        NewCond = RootDevelopment(Crop, Soil, Groundwater, NewCond, GDD, GrowingSeason);

        //Update growth stage
        NewCond = GrowthStage(Crop, NewCond, GrowingSeason);

        //Canopy cover development
        NewCond = CanopyCover(Crop, Soil, NewCond, GDD, Et0, CO2, GrowingSeason);

        //Reference harvest index
        NewCond = HIrefCurrentDay(NewCond, Crop, GrowingSeason);

        //Biomass accumulation
        NewCond = BiomassAccumulation(Crop, NewCond, Tr, TrPot, Et0, GrowingSeason);

        //Harvest index
        NewCond = HarvestIndex(Soil, Crop, NewCond, Et0, Tmax, Tmin, GrowingSeason);

        //Crop yield
        if (GrowingSeason) {
            //Calculate crop yield (tonne/ha)
            NewCond.Y = (NewCond.B / 100) * NewCond.HIadj;
            //Check if crop has reached maturity
            if (((Crop.CalendarType == 1) && (NewCond.DAP >= Crop.Maturity)) ||
                    ((Crop.CalendarType == 2) && (NewCond.GDDcum >= Crop.Maturity))) {
                //Crop has reached maturity
                NewCond.CropMature = true;
            }
        } else {
            //Crop yield is zero outside of growing season
            NewCond.Y = 0;
        }
        return NewCond;
    }
}