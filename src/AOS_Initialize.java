import Structs.*;

import java.io.*;
import java.sql.Date;
import java.text.DecimalFormat;
import java.time.LocalDate;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

public class AOS_Initialize {

    public static final int MILLISECONDS_IN_DAY = 86400030;
    public static final int DATE_ADD = 719530;
    public ClockStruct clockStruct;
    public AOS_InitialiseStruct AOS_InitialiseStruct;

    public AOS_Initialize() {
        FileLocation fileLocation = AOS_ReadFileLocations();
        clockStruct = AOS_ReadClockParameters(fileLocation);
        double[][] WeatherStruct = AOS_ReadWeatherInputs(fileLocation, clockStruct);
        ParamStruct paramStruct = AOS_ReadModelParameters(fileLocation, clockStruct);
        IrrMngtStruct irrMngtStruct = AOS_ReadIrrigationManagement(fileLocation, paramStruct);
        FieldMngtStruct[] fieldMngtStruct = AOS_ReadFieldManagement(fileLocation, paramStruct);
        GwStruct gwStruct = AOS_ReadGroundwaterTable(fileLocation);
        paramStruct = AOS_ComputeVariables(paramStruct, WeatherStruct, clockStruct, gwStruct, fileLocation);

        //Define initial conditions
        InitCondStruct initCondStruct = AOS_ReadModelInitialConditions(paramStruct, gwStruct, fieldMngtStruct, fileLocation);

        //Pack output structure
        AOS_InitialiseStruct = new AOS_InitialiseStruct();
        AOS_InitialiseStruct.Parameter = paramStruct;
        AOS_InitialiseStruct.IrrigationManagement = irrMngtStruct;
        AOS_InitialiseStruct.FieldManagement = fieldMngtStruct;
        AOS_InitialiseStruct.Groundwater = gwStruct;
        AOS_InitialiseStruct.InitialCondition = initCondStruct;
        AOS_InitialiseStruct.CropChoices = paramStruct.crop;
        AOS_InitialiseStruct.Weather = WeatherStruct;
        AOS_InitialiseStruct.FileLocation = fileLocation;

        //Setup output files
        //Define output file location
        String FileLoc = fileLocation.output;
        //Setup blank matrices to store outputs
        AOS_InitialiseStruct.Outputs.WaterContents = new int[clockStruct.TimeSpan.length][5 + paramStruct.soil.nComp];
        for (int i = 0; i < clockStruct.TimeSpan.length; i++) {
            for (int j = 3; j < 5 + paramStruct.soil.nComp; j++) {
                if (j == 4) {
                    continue;
                }
                AOS_InitialiseStruct.Outputs.WaterContents[i][j] = -999;
            }
        }
        AOS_InitialiseStruct.Outputs.WaterFluxes = new int[clockStruct.TimeSpan.length][18];
        for (int i = 0; i < clockStruct.TimeSpan.length; i++) {
            for (int j = 3; j < 18; j++) {
                if (j == 4) {
                    continue;
                }
                AOS_InitialiseStruct.Outputs.WaterFluxes[i][j] = -999;
            }
        }
        AOS_InitialiseStruct.Outputs.CropGrowth = new int[clockStruct.TimeSpan.length][15];
        for (int i = 0; i < clockStruct.TimeSpan.length; i++) {
            for (int j = 3; j < 15; j++) {
                if (j == 4) {
                    continue;
                }
                AOS_InitialiseStruct.Outputs.CropGrowth[i][j] = -999;
            }
        }
        //TODO check it at the end
        AOS_InitialiseStruct.Outputs.FinalOutput = new String[clockStruct.nSeasons];

        //Store dates in daily matrices
        String s = clockStruct.SimulationStartTime;
        String e = clockStruct.SimulationEndTime;
        LocalDate start = LocalDate.parse(s);
        LocalDate end = LocalDate.parse(e);
        int i = 0;
        while (!start.isAfter(end)) {
            String[] split = start.toString().split("-");
            for (int j = 0; j < 3; j++) {
                AOS_InitialiseStruct.Outputs.WaterContents[i][j] = Integer.parseInt(split[j]);
                AOS_InitialiseStruct.Outputs.WaterFluxes[i][j] = Integer.parseInt(split[j]);
                AOS_InitialiseStruct.Outputs.CropGrowth[i][j] = Integer.parseInt(split[j]);
            }
            start = start.plusDays(1);
            i++;
        }

        if (AOS_InitialiseStruct.FileLocation.writeDaily.compareTo("Y") == 0) {
            //Water contents (daily)
            String[] names = new String[paramStruct.soil.nComp];
            DecimalFormat df = new DecimalFormat("#.###");
            for (int ii = 0; ii < paramStruct.soil.nComp; ii++) {
                double z = paramStruct.soil.comp.dzsum[ii] - (paramStruct.soil.comp.dz[ii] / 2);
                z = Double.parseDouble(df.format(z));
                names[ii] = String.valueOf(z).concat("m");
            }

            String path = FileLoc.concat(fileLocation.outputFilename + "_WaterContents.csv");
            String[] cHeader = new String[5 + names.length];
            cHeader[0] = "Year";
            cHeader[1] = ",Month";
            cHeader[2] = ",Day";
            cHeader[3] = ",SimDay";
            cHeader[4] = ",Season";
            for (i = 0; i < names.length; i++) {
                cHeader[i + 5] = "," + names[i];
            }
            write(path, cHeader);

            //Hydrological fluxes (daily)
            path = FileLoc.concat(fileLocation.outputFilename + "_WaterFluxes.csv");
            write(path, new String[]{"Year,Month,Day,SimDay,Season,wRZ,zGW,wSurf,Irr,Infl,RO,DP,CR,GWin,Es,EsX,Tr,TrX"});

            //Crop growth (daily)
            path = FileLoc.concat(fileLocation.outputFilename + "CropGrowth.csv");
            write(path, new String[]{"Year,Month,Day,SimDay,Season,GDD,TotGDD,Zr,CC,CCPot,Bio,BioPot,HI,HIadj,Yield"});
        }

        //Final output (at end of each growing season)
        String path = FileLoc.concat(fileLocation.outputFilename + "_FinalOutput.csv");
        write(path, new String[]{"Season,CropType,PlantDate,PlantSimDate,HarvestDate,HarvestSimDate,Yield,TotIrr"});
    }

    private void write(String path, String[] text) {
        File file = new File(path);
        file.getParentFile().mkdirs();
        BufferedWriter writer;
        try {
            writer = new BufferedWriter(new FileWriter(path));
            for (String s : text) {
                writer.write(s);
            }
            writer.close();
        } catch (IOException e1) {
            e1.printStackTrace();
        }
    }

    //Function to set up initial model conditions
    public InitCondStruct AOS_ReadModelInitialConditions(ParamStruct paramStruct, GwStruct
            gwStruct, FieldMngtStruct[] fieldMngtStruct, FileLocation fileLocation) {
        //Define initial conditions
        InitCondStruct InitCondStruct = new InitCondStruct();

        //Counters
        InitCondStruct.AgeDays = 0;
        InitCondStruct.AgeDays_NS = 0;
        InitCondStruct.AerDays = 0;
        InitCondStruct.IrrCum = 0;
        InitCondStruct.DelayedGDDs = 0;
        InitCondStruct.DelayedCDs = 0;
        InitCondStruct.PctLagPhase = 0;
        InitCondStruct.tEarlySen = 0;
        InitCondStruct.GDDcum = 0;
        InitCondStruct.DaySubmerged = 0;
        InitCondStruct.IrrNetCum = 0;
        InitCondStruct.DAP = 0;
        InitCondStruct.Epot = 0;
        InitCondStruct.Tpot = 0;

        //States
        InitCondStruct.PreAdj = false;
        InitCondStruct.CropMature = false;
        InitCondStruct.CropDead = false;
        InitCondStruct.Germination = false;
        InitCondStruct.PrematSenes = false;
        InitCondStruct.HarvestFlag = false;

        //Harvest index
        InitCondStruct.Stage = 1;
        InitCondStruct.Fpre = 1;
        InitCondStruct.Fpost = 1;
        InitCondStruct.fpost_dwn = 1;
        InitCondStruct.fpost_upp = 1;
        InitCondStruct.HIcor_Asum = 0;
        InitCondStruct.HIcor_Bsum = 0;
        InitCondStruct.Fpol = 0;
        InitCondStruct.sCor1 = 0;
        InitCondStruct.sCor2 = 0;

        //Growth stage
        InitCondStruct.GrowthStage = 0;

        //Aeration stress (compartment level)
        InitCondStruct.AerDaysComp = new double[paramStruct.soil.nComp];

        //Transpiration
        InitCondStruct.TrRatio = 1;

        //Soil evaporation
        InitCondStruct.Ke = paramStruct.soil.Kex;

        //Crop growth
        InitCondStruct.CC = 0;
        InitCondStruct.CCadj = 0;
        InitCondStruct.CC_NS = 0;
        InitCondStruct.CCadj_NS = 0;
        InitCondStruct.Zroot = 0;
        InitCondStruct.B = 0;
        InitCondStruct.B_NS = 0;
        InitCondStruct.HI = 0;
        InitCondStruct.HIadj = 0;
        InitCondStruct.CCxAct = 0;
        InitCondStruct.CCxAct_NS = 0;
        InitCondStruct.CCxW = 0;
        InitCondStruct.CCxW_NS = 0;
        InitCondStruct.CCxEarlySen = 0;
        InitCondStruct.CCprev = 0;
        InitCondStruct.ProtectedSeed = false;
        InitCondStruct.Kcb = 0;
        InitCondStruct.Kcb_NS = 0;
        InitCondStruct.rCor = 1;

        if (clockStruct.SeasonCounter == 0) {
            InitCondStruct.Zroot = 0;
            InitCondStruct.CC0adj = 0;
        } else if (clockStruct.SeasonCounter == 1) {
            InitCondStruct.Zroot = paramStruct.crop[clockStruct.SeasonCounter - 1].Zmin;
            InitCondStruct.CC0adj = paramStruct.crop[clockStruct.SeasonCounter - 1].CC0;
        }

        //Initial surface storage between any soil bunds
        if (clockStruct.SeasonCounter == 0) {
            //First day of simulation is in fallow period
            //TODO
        } else if (clockStruct.SeasonCounter == 1) {
            //First day of simulation is in first growing season
            //Get relevant field management structure parameters
            FieldMngtStruct FieldMngtTmp = fieldMngtStruct[clockStruct.SeasonCounter];
            if (FieldMngtTmp.Bunds.compareTo("Y") == 0 && FieldMngtTmp.zBund > 0.001) {
                //Get initial storage between surface bunds
                InitCondStruct.SurfaceStorage = FieldMngtTmp.BundWater;
                if (InitCondStruct.SurfaceStorage > FieldMngtTmp.zBund) {
                    InitCondStruct.SurfaceStorage = FieldMngtTmp.zBund;
                }
            } else {
                //No surface bunds
                InitCondStruct.SurfaceStorage = 0;
            }
        }

        //Check for presence of groundwater table
        if (gwStruct.WaterTable == 0) { //No water table present
            //Set initial groundwater level to dummy value
            InitCondStruct.zGW = -999;
            InitCondStruct.WTinSoil = false;
            //Set adjusted field capacity to default field capacity
            InitCondStruct.th_fc_Adj = paramStruct.soil.comp.th_fc;
        } else if (gwStruct.WaterTable == 1) { //Water table is present
            //Set initial groundwater level
            //TODO
        }

        //Define initial water contents
        //Read input file
        String Location = fileLocation.input;
        String filename = Location.concat("\\" + fileLocation.initialWCFilename);
        try {
            //check the file exists
            new FileReader(filename);
        } catch (FileNotFoundException e) {
            //Can't find text file defining soil initial water content
            System.out.println(e.getMessage());
            //Throw error message
            return null;
        }

        //Load data
        List<String> dataArray = new LinkedList<>();
        try {
            BufferedReader br = new BufferedReader(new FileReader(new File(filename)));
            String st;
            while ((st = br.readLine()) != null) {
                dataArray.add(st);
            }
        } catch (IOException e) {
            System.out.println(e.getMessage());
        }


        String TypeStr = dataArray.get(1).split(",")[1];
        String MethodStr = dataArray.get(2).split(",")[1];
        //TODO check for input with more than one value
        String[] Data_Pts = dataArray.get(5).split(",");

        //Extract data
        double[] Locs = new double[1];
        Locs[0] = Double.valueOf(Data_Pts[0]);
        double[] Vals = new double[(int) Locs[0]];
        if (MethodStr.compareTo("Depth") == 0) {
            for (int i = 0; i < Locs.length; i++) {
                Locs[i] *= 100;
            }
            Locs = round(Locs);
            for (int i = 0; i < Locs.length; i++) {
                Locs[i] /= 100;
            }
        }

        //Define soil compartment depths and layers
        Integer[] SoilLayers = paramStruct.soil.comp.layer;
        Double[] SoilDepths = cumsum(paramStruct.soil.comp.dz);
        for (int i = 0; i < SoilDepths.length; i++) {
            SoilDepths[i] *= 100;
        }
        SoilDepths = round(SoilDepths);
        for (int i = 0; i < SoilDepths.length; i++) {
            SoilDepths[i] /= 100;
        }


        //Assign data
        if (TypeStr.compareTo("Num") == 0) {
            //Values are defined as numbers (m3/m3) so no calculation required
            //TODO
        } else if (TypeStr.compareTo("Pct") == 0) {
            //TODO
        } else if (TypeStr.compareTo("Prop") == 0) {
            //Values are specified as soil hydraulic properties (SAT, FC, or WP).
            //Extract and assign value for each soil layer
            String ValsTmp = Data_Pts[1];
            //TODO change to ValsTmp array
            for (int ii = 0; ii < 1; ii++) {
                int s = ValsTmp.length();
                if (MethodStr.compareTo("Depth") == 0) {
                    //Find layer at specified depth
                    //TODO
                    if (Locs[ii] < SoilDepths[SoilDepths.length - 1]) {
                    }
                } else if (MethodStr.compareTo("Layer") == 0) {
                    //Calculate moisture content at specified layer
                    double LayTmp = Locs[ii];
                    if (ValsTmp.compareTo("SAT") == 0) {
                        //Calculate moisture content at specified layer
                        Vals[ii] = paramStruct.soil.layer.th_s[(int) LayTmp - 1];
                    } else if (ValsTmp.compareTo("FC") == 0) {
                        Vals[ii] = paramStruct.soil.layer.th_fc[(int) LayTmp - 1];
                    } else if (ValsTmp.compareTo("WP") == 0) {
                        Vals[ii] = paramStruct.soil.layer.th_wp[(int) LayTmp - 1];
                    }
                }
            }
        }

        //Interpolate values to all soil compartments
        double[] thini = new double[paramStruct.soil.nComp];
        if (MethodStr.compareTo("Layer") == 0) {
            for (int ii = 0; ii < Vals.length; ii++) {
                for (int j = 0; j < paramStruct.soil.comp.layer.length; j++) {
                    if (paramStruct.soil.comp.layer[j] == Locs[ii]) {
                        thini[j] = Vals[ii];
                    }
                }
            }
            InitCondStruct.th = thini;
        } else if (MethodStr.compareTo("Depth") == 0) {
            //Add zero point
            if (Locs[0] > 0) {
                //TODO
            }
        }

        //If groundwater table is present and calculating water contents based on field capacity,
        // then reset value to account for possible changes in field capacity caused by capillary rise effects
        if (gwStruct.WaterTable == 1) {
            //TODO
        }

        //If groundwater table is present in soil profile then set all water
        // contents below the water table to saturation
        if (InitCondStruct.WTinSoil) {
            //Find compartment mid-points
            Double[] zBot = cumsum(paramStruct.soil.comp.dz);
            //TODO
        }

        InitCondStruct.thini = InitCondStruct.th.clone();

        return InitCondStruct;
    }

    /**
     * Function to compute additional variables needed to run AOS
     */
    public ParamStruct AOS_ComputeVariables(ParamStruct paramStruct, double[][] WeatherStruct, ClockStruct
            clockStruct, GwStruct gwStruct, FileLocation fileLocation) {
        //Compute water contents and saturated hydraulic conductivity
        if (paramStruct.soil.CalcSHP == 0) {
            //Read soil texture file
            String fileName = fileLocation.input.concat("\\" + fileLocation.soilHydrologyFilename);

            try {
                //check the file exists
                new FileReader(fileName);
            } catch (FileNotFoundException e) {
                //Can't find text file defining soil hydraulic properties
                System.out.println(e.getMessage());
                return paramStruct;
            }

            //Load data
            List<String> dataArray = new LinkedList<>();
            try {
                BufferedReader br = new BufferedReader(new FileReader(new File(fileName)));
                String st;
                while ((st = br.readLine()) != null) {
                    dataArray.add(st);
                }
            } catch (IOException e) {
                System.out.println(e.getMessage());
            }

            String[] split = dataArray.get(1).split(",");
            //Assign data
            paramStruct.soil.layer.dz = new double[paramStruct.soil.nComp];
            paramStruct.soil.layer.dz[0] = Double.parseDouble(split[1]);
            paramStruct.soil.layer.th_s = new double[paramStruct.soil.nComp];
            paramStruct.soil.layer.th_s[0] = Double.parseDouble(split[2]);
            paramStruct.soil.layer.th_fc = new double[paramStruct.soil.nComp];
            paramStruct.soil.layer.th_fc[0] = Double.parseDouble(split[3]);
            paramStruct.soil.layer.th_wp = new double[paramStruct.soil.nComp];
            paramStruct.soil.layer.th_wp[0] = Double.parseDouble(split[4]);
            paramStruct.soil.layer.Ksat = new double[paramStruct.soil.nComp];
            paramStruct.soil.layer.Ksat[0] = Double.parseDouble(split[5]);
            paramStruct.soil.layer.Penetrability = new double[paramStruct.soil.nComp];
            paramStruct.soil.layer.Penetrability[0] = Double.parseDouble(split[6]);
            //Calculate additional variables
            paramStruct.soil.layer.th_dry = new double[paramStruct.soil.nComp];
            paramStruct.soil.layer.th_dry[0] = paramStruct.soil.layer.th_wp[0] / 2;
        }

        paramStruct.soil.comp.th_fc = new double[paramStruct.soil.nComp];
        //Assign field capacity values to each soil compartment
        for (int i = 0; i < paramStruct.soil.nComp; i++) {
            int layeri = paramStruct.soil.comp.layer[i] - 1;
            paramStruct.soil.comp.th_fc[i] = paramStruct.soil.layer.th_fc[layeri];
        }

        //Calculate capillary rise parameters for all soil layers
        //Only do calculation if water table is present. Calculations use equations
        //described in Raes et al. (2012)
        if (gwStruct.WaterTable == 1) {
            //TODO
        }

        //Calculate drainage characteristic (tau)
        //Calculations use equation given by Raes et al. 2012
        paramStruct.soil.layer.tau = new double[paramStruct.soil.nLayer];
        for (int i = 0; i < paramStruct.soil.nLayer; i++) {
            paramStruct.soil.layer.tau[i] = 0.0866 * (Math.pow(paramStruct.soil.layer.Ksat[i], 0.35));
            paramStruct.soil.layer.tau[i] = Math.round((100 * paramStruct.soil.layer.tau[i])) / 100.0;
            if (paramStruct.soil.layer.tau[i] > 1) {
                paramStruct.soil.layer.tau[i] = 1;
            } else if (paramStruct.soil.layer.tau[i] < 0) {
                paramStruct.soil.layer.tau[i] = 0;
            }
        }


        //Calculate readily evaporable water in surface layer %%
        if (paramStruct.soil.AdjREW == 0) {
            paramStruct.soil.REW = Math.round((1000 * (paramStruct.soil.layer.th_fc[0] - paramStruct.soil.layer.th_dry[0]) * paramStruct.soil.EvapZsurf));
        }

        //Calculate additional parameters for all crop types in mix
        int nCrops = paramStruct.crop.length;
        for (int i = 0; i < nCrops; i++) {
            //Fractional canopy cover size at emergence
            paramStruct.crop[i].CC0 = Math.round(10000 * (paramStruct.crop[i].PlantPop *
                    paramStruct.crop[i].SeedSize) * Math.pow(10, -8)) / 10000.0;
            //Root extraction terms
            double SxTopQ = paramStruct.crop[i].SxTopQ;
            double SxBotQ = paramStruct.crop[i].SxBotQ;
            double S1 = paramStruct.crop[i].SxTopQ;
            double S2 = paramStruct.crop[i].SxBotQ;

            double SxTop, SxBot, xx, SS1, SS2;
            if (S1 == S2) {
                SxTop = S1;
                SxBot = S2;
            } else {
                if (SxTopQ < SxBotQ) {
                    S1 = SxBotQ;
                    S2 = SxTopQ;
                }
                xx = 3 * (S2 / (S1 - S2));
                if (xx < 0.5) {
                    SS1 = (4 / 3.5) * S1;
                    SS2 = 0;
                } else {
                    SS1 = (xx + 3.5) * (S1 / (xx + 3));
                    SS2 = (xx - 0.5) * (S2 / xx);
                }
                if (SxTopQ > SxBotQ) {
                    SxTop = SS1;
                    SxBot = SS2;
                } else {
                    SxTop = SS2;
                    SxBot = SS1;
                }
            }
            paramStruct.crop[i].SxTop = SxTop;
            paramStruct.crop[i].SxBot = SxBot;

            //Water stress thresholds
            paramStruct.crop[i].p_up = new double[]{paramStruct.crop[i].p_up1, paramStruct.crop[i].p_up2,
                    paramStruct.crop[i].p_up3, paramStruct.crop[i].p_up4};
            paramStruct.crop[i].p_lo = new double[]{paramStruct.crop[i].p_lo1, paramStruct.crop[i].p_lo2,
                    paramStruct.crop[i].p_lo3, paramStruct.crop[i].p_lo4};
            paramStruct.crop[i].fshape_w = new double[]{paramStruct.crop[i].fshape_w1, paramStruct.crop[i].fshape_w2,
                    paramStruct.crop[i].fshape_w3, paramStruct.crop[i].fshape_w4};

            //Flowering function
            if (paramStruct.crop[i].CropType == 3) {
                //TODO
            }

            //Crop calendar
            //TODO
            paramStruct.crop[i] = AOS_ComputeCropCalendar(paramStruct.crop[i], WeatherStruct);

            //Harvest index growth coefficient
            //TODO

            //Days to linear HI switch point
            if (paramStruct.crop[i].CropType == 3) {
                //Determine linear switch point and HIGC rate for fruit/grain crops
                double tLin = 0, HIGClin = 0;
                //TODO
                //[tLin,HIGClin] = AOS_CalculateHILinear(paramStruct.crop[i]);
                paramStruct.crop[i].tLinSwitch = tLin;
                paramStruct.crop[i].dHILinear = HIGClin;
            } else {
                //No linear switch for leafy vegetable or root/tiber crops
                paramStruct.crop[i].tLinSwitch = 0;
                paramStruct.crop[i].dHILinear = 0;
            }
        }

        //Calculate WP adjustment factor for elevation in CO2 concentration
        //Load CO2 data
        String fileName = fileLocation.input.concat("\\" + fileLocation.CO2Filename);

        try {
            //check the file exists
            new FileReader(fileName);
        } catch (FileNotFoundException e) {
            //Can't find text file defining CO2 concentrations
            System.out.println(e.getMessage());
            return paramStruct;
        }

        //Load data
        List<String> dataArray = new LinkedList<>();
        try {
            BufferedReader br = new BufferedReader(new FileReader(new File(fileName)));
            String st;
            while ((st = br.readLine()) != null) {
                dataArray.add(st);
            }
        } catch (IOException e) {
            System.out.println(e.getMessage());
        }

        int[] Yrs = new int[dataArray.size() - 1];
        double[] CO2 = new double[dataArray.size() - 1];
        for (int i = 1; i < dataArray.size(); i++) {
            String[] split = dataArray.get(i).split(",");
            //Years
            Yrs[i - 1] = Integer.parseInt(split[0]);
            //CO2 concentrations (ppm)
            CO2[i - 1] = Double.parseDouble(split[1]);
        }

        List<Integer> YrsVec = new LinkedList<>();
        for (int i = Integer.parseInt(clockStruct.SimulationStartTime.split("-")[0]); i <= Integer.parseInt(clockStruct.SimulationEndTime.split("-")[0]); i++) {
            YrsVec.add(i);
        }

        //Store data
        paramStruct.CO2.Data = new Data[YrsVec.size()];
        for (int i = 0; i < YrsVec.size(); i++) {
            double CO2conc = 0.0;
            //TODO interp1
//            CO2conc = interp1(Yrs, CO2, YrsVec);
            paramStruct.CO2.Data[i] = new Data(YrsVec.get(i), CO2conc);
        }

        //Define reference CO2 concentration
        paramStruct.CO2.RefConc = 369.41;

        //Get CO2 concentration for first year
        paramStruct.CO2.CurrentConc = paramStruct.CO2.Data[0].value;

        //Get CO2 weighting factor for first year
        double CO2ref = paramStruct.CO2.RefConc;
        double CO2conc = paramStruct.CO2.CurrentConc;

        double fw, ftype;
        if (CO2conc <= CO2ref) {
            fw = 0;
        } else {
            if (CO2conc >= 550) {
                fw = 1;
            } else {
                fw = 1 - ((550 - CO2conc) / (550 - CO2ref));
            }
        }

        //Determine adjustment for each crop in first year of simulation
        for (int i = 0; i < nCrops; i++) {
            //Determine initial adjustment
            double fCO2 = (CO2conc / CO2ref) / (1 + (CO2conc - CO2ref) * ((1 - fw) *
                    paramStruct.crop[i].bsted + fw * ((paramStruct.crop[i].bsted *
                    paramStruct.crop[i].fsink) + (paramStruct.crop[i].bface *
                    (1 - paramStruct.crop[i].fsink)))));
            //Consider crop type
            if (paramStruct.crop[i].WP >= 40) {
                //No correction for C4 crops
                ftype = 0;
            } else if (paramStruct.crop[i].WP <= 20) {
                //Full correction for C3 crops
                ftype = 1;
            } else {
                ftype = (40 - paramStruct.crop[i].WP) / (40 - 20);
            }
            //Total adjustment
            paramStruct.crop[i].fCO2 = 1 + ftype * (fCO2 - 1);
        }

        return paramStruct;
    }

    //Function to compute additional parameters needed to define crop phenological calendar
    private Crop AOS_ComputeCropCalendar(Crop Crop, double[][] Weather) {
        //Define crop calendar mode
        double Mode = Crop.CalendarType;

        //Calculate variables
        if (Mode == 1) { //Growth in calendar days
            //Time from sowing to end of vegatative growth period
            if (Crop.Determinant == 1) {
                Crop.CanopyDevEnd = Math.round(Crop.HIstart + (Crop.Flowering / 2));
            } else {
                Crop.CanopyDevEnd = Crop.Senescence;
            }

            //Time from sowing to 10% canopy cover (non-stressed conditions)
            Crop.Canopy10Pct = Math.round(Crop.Emergence + (Math.log(0.1 / Crop.CC0) / Crop.CGC));

            //Time from sowing to maximum canopy cover (non-stressed conditions)
            //TODO
        } else if (Mode == 2) { //Growth in growing degree days
            //Time from sowing to end of vegatative growth period
            if (Crop.Determinant == 1) {
                Crop.CanopyDevEnd = Math.round(Crop.HIstart + (Crop.Flowering / 2));
            } else {
                Crop.CanopyDevEnd = Crop.Senescence;
            }

            //Time from sowing to 10% canopy cover (non-stressed conditions)
            Crop.Canopy10Pct = Math.round(Crop.Emergence + (Math.log(0.1 / Crop.CC0) / Crop.CGC));

            //Time from sowing to maximum canopy cover (non-stressed conditions)
            Crop.MaxCanopy = Math.round(Crop.Emergence + (Math.log((0.25 * Crop.CCx * Crop.CCx / Crop.CC0) /
                    (Crop.CCx - (0.98 * Crop.CCx))) / Crop.CGC));

            //Time from sowing to end of yield formation
            Crop.HIend = Crop.HIstart + Crop.YldForm;

            //Time from sowing to end of flowering (if fruit/grain crop)
            if (Crop.CropType == 3) {
                Crop.FloweringEnd = Crop.HIstart + Crop.Flowering;
            }

            //Additional variables
            //Extract weather data for first growing season that crop is planted
            //TODO check it
//            idx = find(strcmp(CropChoices,CropName)==1,1);
            double tSta = clockStruct.PlantingDate[0];
            double tSto = clockStruct.HarvestDate[0];
            double[] Dates = new double[Weather.length];
            for (int i = 0; i < Dates.length; i++) {
                Dates[i] = Weather[i][0];
            }
            int StaRow = 0, StoRow = 0;
            for (int i = 0; i < Dates.length; i++) {
                if (Dates[i] == tSta) {
                    StaRow = i;
                }
                if (Dates[i] == tSto) {
                    StoRow = i;
                }
            }
            double[] Tmin = new double[(StoRow - StaRow) + 1];
            double[] Tmax = new double[(StoRow - StaRow) + 1];
            for (int i = 0; i <= StoRow - StaRow; i++) {
                Tmin[i] = Weather[i][1];
                Tmax[i] = Weather[i][2];
            }
            //Calculate GDD's
            double[] GDD = new double[0];
            if (Crop.GDDmethod == 1) {
                //TODO
            } else if (Crop.GDDmethod == 2) {
                for (int i = 0; i < Tmax.length; i++) {
                    if (Tmax[i] > Crop.Tupp) {
                        Tmax[i] = Crop.Tupp;
                    }
                    if (Tmax[i] < Crop.Tbase) {
                        Tmax[i] = Crop.Tbase;
                    }
                    if (Tmin[i] > Crop.Tupp) {
                        Tmin[i] = Crop.Tupp;
                    }
                    if (Tmin[i] < Crop.Tbase) {
                        Tmin[i] = Crop.Tbase;
                    }
                }
                double[] Tmean = new double[Tmax.length];
                GDD = new double[Tmean.length];
                for (int i = 0; i < Tmean.length; i++) {
                    Tmean[i] = (Tmax[i] + Tmin[i]) / 2;
                    GDD[i] = Tmean[i] - Crop.Tbase;
                }

            } else if (Crop.GDDmethod == 3) {
                //TODO
            }
            double[] GDDcum = cumsum(GDD);

            //Find calendar days for some variables
            boolean set1 = false, set2 = false, set3 = false, set4 = false;
            for (int i = 0; i < GDDcum.length; i++) {
                //1. Calendar days from sowing to maximum canopy cover
                if (GDDcum[i] > Crop.MaxCanopy && !set1) {
                    Crop.MaxCanopyCD = i;
                    set1 = true;
                }
                //1. Calendar days from sowing to end of vegetative growth
                if (GDDcum[i] > Crop.CanopyDevEnd && !set2) {
                    Crop.CanopyDevEndCD = i;
                    set2 = true;
                }
                //2. Calendar days from sowing to start of yield formation
                if (GDDcum[i] > Crop.HIstart && !set3) {
                    Crop.HIstartCD = i;
                    set3 = true;
                }
                //3. Calendar days from sowing to end of yield formation
                if (GDDcum[i] > Crop.HIend && !set4) {
                    Crop.HIendCD = i;
                    set4 = true;
                }
            }
            //4. Duration of yield formation in calendar days
            Crop.YldFormCD = Crop.HIendCD - Crop.HIstartCD;

            if (Crop.CropType == 3) {
                //1. Calendar days from sowing to end of flowering
                int FloweringEnd = 0;
                for (int i = 0; i < GDDcum.length; i++) {
                    if (GDDcum[i] > Crop.FloweringEnd) {
                        FloweringEnd = i;
                        break;
                    }
                }
                //2. Duration of flowering in calendar days
                Crop.FloweringCD = FloweringEnd - Crop.HIstartCD;
            }
        }
        return Crop;
    }


    /**
     * Function to read input file and initialise groundwater table parameters
     */
    public GwStruct AOS_ReadGroundwaterTable(FileLocation fileLocation) {
        //Read input file location
        String location = fileLocation.input;
        //Define empty structure
        GwStruct gwStruct = new GwStruct();

        //Read groundwater table input file
        //Open file
        String fileName = location.concat("\\" + fileLocation.groundwaterFilename);
        try {
            //check the file exists
            new FileReader(fileName);
        } catch (FileNotFoundException e) {
            //Can't find text file defining irrigation management
            System.out.println(e.getMessage());
            return gwStruct;
        }

        //Load data
        List<String> dataArray = new LinkedList<>();
        try {
            BufferedReader br = new BufferedReader(new FileReader(new File(fileName)));
            String st;
            while ((st = br.readLine()) != null) {
                dataArray.add(st);
            }
        } catch (IOException e) {
            System.out.println(e.getMessage());
        }

        String WT = dataArray.get(1).split(",")[1];
        String Method = dataArray.get(2).split(",")[1];

        if (WT.compareTo("N") == 0) {
            //No water table present (don't read the rest of the input file)
            gwStruct.WaterTable = 0;
        } else if (WT.compareTo("Y") == 0) {
            //TODO
            gwStruct.WaterTable = 1;
        }
        return gwStruct;
    }

    /**
     * Function to read input files and initialise field management parameters
     */
    public FieldMngtStruct[] AOS_ReadFieldManagement(FileLocation fileLocation, ParamStruct paramStruct) {
        //Get input file location
        String location = fileLocation.input;
        //Read field management parameter input files (growing seasons) %%
        //Check for number of crop types

        //TODO change
        //Crops = fieldnames(ParamStruct.Crop);
        int nCrops = 1;

        //Create blank structure
        FieldMngtStruct[] fieldMngtStruct = new FieldMngtStruct[nCrops + 1];
        int i;
        for (i = 0; i < nCrops; i++) {
            //Open file
            String fileName = location.concat("\\" + paramStruct.crop[i].FieldMngtFile);
            try {
                //check the file exists
                new FileReader(fileName);
            } catch (FileNotFoundException e) {
                //Can't find text file defining irrigation management
                System.out.println(e.getMessage());
                return fieldMngtStruct;
            }

            //Load data
            List<String> dataArray = new LinkedList<>();
            try {
                BufferedReader br = new BufferedReader(new FileReader(new File(fileName)));
                String st;
                while ((st = br.readLine()) != null) {
                    dataArray.add(st);
                }
            } catch (IOException e) {
                System.out.println(e.getMessage());
            }

            for (int j = 0; j < dataArray.size(); j++) {
                String temp = dataArray.get(j).split(",")[1];
                dataArray.set(j, temp);
            }
            dataArray.remove(0);

            fieldMngtStruct[i] = new FieldMngtStruct();
            fieldMngtStruct[i].Mulches = dataArray.get(0);
            fieldMngtStruct[i].Bunds = dataArray.get(1);
            fieldMngtStruct[i].CNadj = dataArray.get(2);
            fieldMngtStruct[i].SRinhb = dataArray.get(3);
            fieldMngtStruct[i].MulchPct = Double.parseDouble(dataArray.get(4));
            fieldMngtStruct[i].fMulch = Double.parseDouble(dataArray.get(5));
            fieldMngtStruct[i].zBund = Double.parseDouble(dataArray.get(6));
            fieldMngtStruct[i].BundWater = Double.parseDouble(dataArray.get(7));
            fieldMngtStruct[i].CNadjPct = Double.parseDouble(dataArray.get(8));
        }

        //Read field management practice input file (fallow periods)
        //Open file
        String fileName = location.concat("\\" + fileLocation.fieldMngtFallowFilename);
        try {
            //check the file exists
            new FileReader(fileName);
        } catch (FileNotFoundException e) {
            //Can't find text file defining soil parameters
            System.out.println(e.getMessage());
            return fieldMngtStruct;
        }

        //Load data
        List<String> dataArray = new LinkedList<>();
        try {
            BufferedReader br = new BufferedReader(new FileReader(new File(fileName)));
            String st;
            while ((st = br.readLine()) != null) {
                dataArray.add(st);
            }
        } catch (IOException e) {
            System.out.println(e.getMessage());
        }

        for (int j = 0; j < dataArray.size(); j++) {
            String temp = dataArray.get(j).split(",")[1];
            dataArray.set(j, temp);
        }
        dataArray.remove(0);

        fieldMngtStruct[i] = new FieldMngtStruct();
        fieldMngtStruct[i].Mulches = dataArray.get(0);
        fieldMngtStruct[i].Bunds = dataArray.get(1);
        fieldMngtStruct[i].CNadj = dataArray.get(2);
        fieldMngtStruct[i].SRinhb = dataArray.get(3);
        fieldMngtStruct[i].MulchPct = Double.parseDouble(dataArray.get(4));
        fieldMngtStruct[i].fMulch = Double.parseDouble(dataArray.get(5));
        fieldMngtStruct[i].zBund = Double.parseDouble(dataArray.get(6));
        fieldMngtStruct[i].BundWater = Double.parseDouble(dataArray.get(7));
        fieldMngtStruct[i].CNadjPct = Double.parseDouble(dataArray.get(8));

        return fieldMngtStruct;
    }

    /**
     * Function to read and process input weather time-series
     */
    private double[][] AOS_ReadWeatherInputs(FileLocation fileLocation, ClockStruct clockStruct) {
        double[][] weatherDB = new double[(clockStruct.SimulationEndDate - clockStruct.SimulationStartDate) + 1][5];

        //Read input file location
        String location = fileLocation.input;
        //Read weather data inputs
        //Open file
        String fileName = location.concat("\\" + fileLocation.weatherFilename);

        try {
            //check the file exists
            new FileReader(fileName);
        } catch (FileNotFoundException e) {
            //Can't find text file defining locations of input and output folders.
            System.out.println(e.getMessage());
            return weatherDB;
        }

        //Load data
        List<String> dataArray = new LinkedList<>();
        try {
            BufferedReader br = new BufferedReader(new FileReader(new File(fileName)));
            String st;
            while ((st = br.readLine()) != null) {
                dataArray.add(st);
            }
        } catch (IOException e) {
            System.out.println(e.getMessage());
        }

        //Extract data
        //TODO find out how to calculate StartRow and EndRow
        int StartRow = 11809, EndRow = 11991;
        for (int i = 0; i <= clockStruct.SimulationEndDate - clockStruct.SimulationStartDate; i++) {
            weatherDB[i][0] = clockStruct.SimulationStartDate + i;
            weatherDB[i][1] = Double.parseDouble(dataArray.get(StartRow + i).split(",")[3]);
            weatherDB[i][2] = Double.parseDouble(dataArray.get(StartRow + i).split(",")[4]);
            weatherDB[i][3] = Double.parseDouble(dataArray.get(StartRow + i).split(",")[5]);
            weatherDB[i][4] = Double.parseDouble(dataArray.get(StartRow + i).split(",")[6]);
        }

        return weatherDB;
    }

    /**
     * Function to read input and output file locations.
     */
    private FileLocation AOS_ReadFileLocations() {
        FileLocation fileLocation = new FileLocation();
        //Read AOS file location input file
        String fileName = "FileLocations.csv";

        try {
            //check the file exists
            new FileReader(fileName);
        } catch (FileNotFoundException e) {
            //Can't find text file defining locations of input and output folders.
            System.out.println(e.getMessage());
            return fileLocation;
        }

        //Load data
        List<String> dataArray = new LinkedList<>();
        try {
            BufferedReader br = new BufferedReader(new FileReader(new File(fileName)));
            String st;
            while ((st = br.readLine()) != null) {
                dataArray.add(st);
            }
        } catch (IOException e) {
            System.out.println(e.getMessage());
        }

        fileLocation.input = dataArray.get(1).split(",")[1];
        fileLocation.output = dataArray.get(2).split(",")[1] + "\\";

        fileName = fileLocation.input.concat("\\FileSetup.csv");
        try {
            //check the file exists
            new FileReader(fileName);
        } catch (FileNotFoundException e) {
            //Can't find text file defining locations of input and output folders.
            System.out.println(e.getMessage());
            return fileLocation;
        }

        dataArray.clear();
        try {
            BufferedReader br = new BufferedReader(new FileReader(new File(fileName)));
            String st;
            while ((st = br.readLine()) != null) {
                dataArray.add(st);
            }
        } catch (IOException e) {
            System.out.println(e.getMessage());
        }

        fileLocation.clockFilename = dataArray.get(1).split(",")[1];
        fileLocation.weatherFilename = dataArray.get(2).split(",")[1];
        fileLocation.cropFilename = dataArray.get(3).split(",")[1];
        fileLocation.soilFilename = dataArray.get(4).split(",")[1];
        fileLocation.fieldMngtFallowFilename = dataArray.get(5).split(",")[1];
        fileLocation.initialWCFilename = dataArray.get(6).split(",")[1];
        fileLocation.groundwaterFilename = dataArray.get(7).split(",")[1];
        fileLocation.CO2Filename = dataArray.get(8).split(",")[1];
        fileLocation.outputFilename = dataArray.get(9).split(",")[1] + "\\";
        fileLocation.writeDaily = dataArray.get(10).split(",")[1];

        return fileLocation;
    }

    /**
     * Function to read input files and initialise model clock parameters.
     */
    private ClockStruct AOS_ReadClockParameters(FileLocation fileLocation) {
        ClockStruct clockStruct = new ClockStruct();

        //Read input file location
        String location = fileLocation.input;
        String fileName = location.concat("\\" + fileLocation.clockFilename);
        //Read clock parameter input file
        try {
            //check the file exists
            new FileReader(fileName);
        } catch (FileNotFoundException e) {
            //Can't find text file defining locations of input and output folders.
            System.out.println(e.getMessage());
            return clockStruct;
        }

        //Load data
        List<String> dataArray = new LinkedList<>();
        try {
            BufferedReader br = new BufferedReader(new FileReader(new File(fileName)));
            String st;
            while ((st = br.readLine()) != null) {
                dataArray.add(st);
            }
        } catch (IOException e) {
            System.out.println(e.getMessage());
        }

        //Create and assign variables
        clockStruct.SimulationStartTime = dataArray.get(1).split(",")[1];
        clockStruct.SimulationEndTime = dataArray.get(2).split(",")[1];
        clockStruct.OffSeason = dataArray.get(3).split(",")[1];

        //Define clock parameters
        //Initialise time step counter
        clockStruct.TimeStepCounter = 1;
        //Initialise model termination condition
        clockStruct.ModelTermination = false;
        // Simulation start time as serial date number
        clockStruct.SimulationStartDate = (int) ((Date.valueOf(clockStruct.SimulationStartTime).getTime() / MILLISECONDS_IN_DAY) + DATE_ADD);
        //Simulation end time as serial date number
        clockStruct.SimulationEndDate = (int) ((Date.valueOf(clockStruct.SimulationEndTime).getTime() / MILLISECONDS_IN_DAY) + DATE_ADD);
        //Time step (years)
        clockStruct.TimeStep = 1;
        //Total numbers of time steps (days)
        clockStruct.nSteps = clockStruct.SimulationEndDate - clockStruct.SimulationStartDate;

        //Time spans
        int[] TimeSpan = new int[clockStruct.nSteps + 1];
        TimeSpan[0] = clockStruct.SimulationStartDate;
        TimeSpan[clockStruct.nSteps - 1] = clockStruct.SimulationEndDate;
        for (int ss = 1; ss < clockStruct.nSteps; ss++) {
            TimeSpan[ss] = TimeSpan[ss - 1] + 1;
        }
        clockStruct.TimeSpan = TimeSpan;
        //Time at start of current time step
        clockStruct.StepStartTime = clockStruct.TimeSpan[clockStruct.TimeStepCounter - 1];
        //Time at end of current time step
        clockStruct.StepEndTime = clockStruct.TimeSpan[clockStruct.TimeStepCounter];
        //Number of time-steps (per day) for soil evaporation calculation
        clockStruct.EvapTimeSteps = 20;

        return clockStruct;
    }

    private ParamStruct AOS_ReadModelParameters(FileLocation fileLocation, ClockStruct clockStruct) {
        ParamStruct paramStruct = new ParamStruct();

        //Read input file location
        String location = fileLocation.input;
        //Read soil parameter input file
        //Open file
        String fileName = location.concat("\\" + fileLocation.soilFilename);

        try {
            //check the file exists
            new FileReader(fileName);
        } catch (FileNotFoundException e) {
            //Can't find text file defining locations of input and output folders.
            System.out.println(e.getMessage());
            return paramStruct;
        }

        //Load data
        List<String> dataArray = new LinkedList<>();
        try {
            BufferedReader br = new BufferedReader(new FileReader(new File(fileName)));
            String st;
            while ((st = br.readLine()) != null) {
                dataArray.add(st);
            }
        } catch (IOException e) {
            System.out.println(e.getMessage());
        }

        //Create assign string variables
        fileLocation.soilProfileFilename = dataArray.get(1).split(",")[1];
        fileLocation.soilTextureFilename = dataArray.get(2).split(",")[1];
        fileLocation.soilHydrologyFilename = dataArray.get(3).split(",")[1];
        paramStruct.soil.CalcSHP = Double.parseDouble(dataArray.get(4).split(",")[1]);
        paramStruct.soil.Zsoil = Double.parseDouble(dataArray.get(5).split(",")[1]);
        paramStruct.soil.nComp = Integer.parseInt(dataArray.get(6).split(",")[1]);
        paramStruct.soil.nLayer = Integer.parseInt(dataArray.get(7).split(",")[1]);
        paramStruct.soil.AdjREW = Double.parseDouble(dataArray.get(8).split(",")[1]);
        paramStruct.soil.REW = Double.parseDouble(dataArray.get(9).split(",")[1]);
        paramStruct.soil.CN = Double.parseDouble(dataArray.get(10).split(",")[1]);
        paramStruct.soil.zRes = Double.parseDouble(dataArray.get(11).split(",")[1]);

        //Assign default program properties (should not be changed without expert knowledge)
        paramStruct.soil.EvapZsurf = 0.04;  //Thickness of soil surface skin evaporation layer (m)
        paramStruct.soil.EvapZmin = 0.15;   //Minimum thickness of full soil surface evaporation layer (m)
        paramStruct.soil.EvapZmax = 0.30;   //Maximum thickness of full soil surface evaporation layer (m)
        paramStruct.soil.Kex = 1.1;         //Maximum soil evaporation coefficient
        paramStruct.soil.fevap = 4.0;         //Shape factor describing reduction in soil evaporation in stage 2.
        paramStruct.soil.fWrelExp = 0.4;    //Proportional value of Wrel at which soil evaporation layer expands
        paramStruct.soil.fwcc = 50.0;         //Maximum coefficient for soil evaporation reduction due to sheltering effect of withered canopy
        paramStruct.soil.zCN = 0.3;         //Thickness of soil surface (m) used to calculate water content to adjust curve number
        paramStruct.soil.zGerm = 0.3;       //Thickness of soil surface (m) used to calculate water content for germination
        paramStruct.soil.AdjCN = 1.0;         //Adjust curve number for antecedent moisture content (0: No, 1: Yes)
        paramStruct.soil.fshape_cr = 16.0;    //Capillary rise shape factor
        paramStruct.soil.zTop = 0.1;        //Thickness of soil surface layer for water stress comparisons (m)

        //Read soil profile input file
        //Open file
        fileName = location.concat("\\" + fileLocation.soilProfileFilename);

        try {
            //check the file exists
            new FileReader(fileName);
        } catch (FileNotFoundException e) {
            //Can't find text file defining locations of input and output folders.
            System.out.println(e.getMessage());
            return paramStruct;
        }

        //Load data
        dataArray = new LinkedList<>();
        try {
            BufferedReader br = new BufferedReader(new FileReader(new File(fileName)));
            String st;
            while ((st = br.readLine()) != null) {
                dataArray.add(st);
            }
        } catch (IOException e) {
            System.out.println(e.getMessage());
        }

        //Create vector of soil compartments sizes and associated layers
        ArrayList<Double> dz = new ArrayList<>();
        for (int i = 1; i < dataArray.size(); i++) {
            dz.add(Double.parseDouble(dataArray.get(i).split(",")[1]));
        }
        paramStruct.soil.comp.dz = dz.toArray(new Double[dz.size()]);
        Double[] dzsum = div(100, round(dot(100.0, cumsum(paramStruct.soil.comp.dz))));
        paramStruct.soil.comp.dzsum = dzsum;
        ArrayList<Integer> Layer = new ArrayList<>();
        for (int i = 1; i < dataArray.size(); i++) {
            Layer.add(Integer.parseInt(dataArray.get(i).split(",")[2]));
        }
        paramStruct.soil.comp.layer = Layer.toArray(new Integer[Layer.size()]);

        //Read crop mix input file
        //Open file
        fileName = location.concat("\\" + fileLocation.cropFilename);
        try {
            //check the file exists
            new FileReader(fileName);
        } catch (FileNotFoundException e) {
            //Can't find text file defining locations of input and output folders.
            System.out.println(e.getMessage());
            return paramStruct;
        }

        //Load data
        dataArray = new LinkedList<>();
        try {
            BufferedReader br = new BufferedReader(new FileReader(new File(fileName)));
            String st;
            while ((st = br.readLine()) != null) {
                dataArray.add(st);
            }
        } catch (IOException e) {
            System.out.println(e.getMessage());
        }

        //Number of crops
        int nCrops = Integer.parseInt(dataArray.get(1).split(",")[1]);
        //Crop rotation filename
        String Rotation = dataArray.get(2).split(",")[1];
        //Crop rotation filename
        String RotationFilename = dataArray.get(3).split(",")[1];
        //Crop information (type and filename)
        String[] CropInfo = dataArray.get(6).split(",");

        //Read crop parameter input files
        //Create blank structure

        //create crop parameter structure
        paramStruct.crop = new Crop[nCrops];

        //Loop crop types
        for (int i = 0; i < nCrops; i++) {
            //Open file
            fileName = location.concat("\\" + CropInfo[1]);
            try {
                //check the file exists
                new FileReader(fileName);
            } catch (FileNotFoundException e) {
                //Can't find text file defining locations of input and output folders.
                System.out.println(e.getMessage());
                return paramStruct;
            }

            //Load data
            dataArray = new LinkedList<>();
            try {
                BufferedReader br = new BufferedReader(new FileReader(new File(fileName)));
                String st;
                while ((st = br.readLine()) != null) {
                    dataArray.add(st);
                }
            } catch (IOException e) {
                System.out.println(e.getMessage());
            }

            String CropType = dataArray.get(1).split(",")[1];
            String PlantMethod = dataArray.get(2).split(",")[1];
            String CalendarType = dataArray.get(3).split(",")[1];
            String SwitchGDD = dataArray.get(4).split(",")[1];
            String PlantingDateStr = dataArray.get(5).split(",")[1];
            String HarvestDateStr = dataArray.get(6).split(",")[1];
            ArrayList<String> temp = new ArrayList<>();
            for (int j = 7; j < dataArray.size(); j++) {
                temp.add(dataArray.get(j).split(",")[1]);
            }
            String[] DataArray = temp.toArray(new String[temp.size()]);

            //Create crop parameter structure
            paramStruct.crop[i] = new Crop();
            paramStruct.crop[i].Emergence = Double.parseDouble(DataArray[0]);
            paramStruct.crop[i].MaxRooting = Double.parseDouble(DataArray[1]);
            paramStruct.crop[i].Senescence = Double.parseDouble(DataArray[2]);
            paramStruct.crop[i].Maturity = Double.parseDouble(DataArray[3]);
            paramStruct.crop[i].HIstart = Double.parseDouble(DataArray[4]);
            paramStruct.crop[i].Flowering = Double.parseDouble(DataArray[5]);
            paramStruct.crop[i].YldForm = Double.parseDouble(DataArray[6]);
            paramStruct.crop[i].GDDmethod = Double.parseDouble(DataArray[7]);
            paramStruct.crop[i].Tbase = Double.parseDouble(DataArray[8]);
            paramStruct.crop[i].Tupp = Double.parseDouble(DataArray[9]);
            paramStruct.crop[i].PolHeatStress = Double.parseDouble(DataArray[10]);
            paramStruct.crop[i].Tmax_up = Double.parseDouble(DataArray[11]);
            paramStruct.crop[i].Tmax_lo = Double.parseDouble(DataArray[12]);
            paramStruct.crop[i].PolColdStress = Double.parseDouble(DataArray[13]);
            paramStruct.crop[i].Tmin_up = Double.parseDouble(DataArray[14]);
            paramStruct.crop[i].Tmin_lo = Double.parseDouble(DataArray[15]);
            paramStruct.crop[i].TrColdStress = Double.parseDouble(DataArray[16]);
            paramStruct.crop[i].GDD_up = Double.parseDouble(DataArray[17]);
            paramStruct.crop[i].GDD_lo = Double.parseDouble(DataArray[18]);
            paramStruct.crop[i].Zmin = Double.parseDouble(DataArray[19]);
            paramStruct.crop[i].Zmax = Double.parseDouble(DataArray[20]);
            paramStruct.crop[i].fshape_r = Double.parseDouble(DataArray[21]);
            paramStruct.crop[i].SxTopQ = Double.parseDouble(DataArray[22]);
            paramStruct.crop[i].SxBotQ = Double.parseDouble(DataArray[23]);
            paramStruct.crop[i].SeedSize = Double.parseDouble(DataArray[24]);
            paramStruct.crop[i].PlantPop = Double.parseDouble(DataArray[25]);
            paramStruct.crop[i].CCx = Double.parseDouble(DataArray[26]);
            paramStruct.crop[i].CDC = Double.parseDouble(DataArray[27]);
            paramStruct.crop[i].CGC = Double.parseDouble(DataArray[28]);
            paramStruct.crop[i].Kcb = Double.parseDouble(DataArray[29]);
            paramStruct.crop[i].fage = Double.parseDouble(DataArray[30]);
            paramStruct.crop[i].WP = Double.parseDouble(DataArray[31]);
            paramStruct.crop[i].WPy = Double.parseDouble(DataArray[32]);
            paramStruct.crop[i].fsink = Double.parseDouble(DataArray[33]);
            paramStruct.crop[i].HI0 = Double.parseDouble(DataArray[34]);
            paramStruct.crop[i].dHI_pre = Double.parseDouble(DataArray[35]);
            paramStruct.crop[i].a_HI = Double.parseDouble(DataArray[36]);
            paramStruct.crop[i].b_HI = Double.parseDouble(DataArray[37]);
            paramStruct.crop[i].dHI0 = Double.parseDouble(DataArray[38]);
            paramStruct.crop[i].Determinant = Double.parseDouble(DataArray[39]);
            paramStruct.crop[i].exc = Double.parseDouble(DataArray[40]);
            paramStruct.crop[i].p_up1 = Double.parseDouble(DataArray[41]);
            paramStruct.crop[i].p_up2 = Double.parseDouble(DataArray[42]);
            paramStruct.crop[i].p_up3 = Double.parseDouble(DataArray[43]);
            paramStruct.crop[i].p_up4 = Double.parseDouble(DataArray[44]);
            paramStruct.crop[i].p_lo1 = Double.parseDouble(DataArray[45]);
            paramStruct.crop[i].p_lo2 = Double.parseDouble(DataArray[46]);
            paramStruct.crop[i].p_lo3 = Double.parseDouble(DataArray[47]);
            paramStruct.crop[i].p_lo4 = Double.parseDouble(DataArray[48]);
            paramStruct.crop[i].fshape_w1 = Double.parseDouble(DataArray[49]);
            paramStruct.crop[i].fshape_w2 = Double.parseDouble(DataArray[50]);
            paramStruct.crop[i].fshape_w3 = Double.parseDouble(DataArray[51]);
            paramStruct.crop[i].fshape_w4 = Double.parseDouble(DataArray[52]);

            //Add additional parameters
            paramStruct.crop[i].CropType = Double.parseDouble(CropType);
            paramStruct.crop[i].PlantMethod = Double.parseDouble(PlantMethod);
            paramStruct.crop[i].CalendarType = Double.parseDouble(CalendarType);
            paramStruct.crop[i].SwitchGDD = Double.parseDouble(SwitchGDD);
            paramStruct.crop[i].PlantingDate = PlantingDateStr;
            paramStruct.crop[i].HarvestDate = HarvestDateStr;
            //Add irrigation management information
            paramStruct.crop[i].IrrigationFile = CropInfo[2];
            paramStruct.crop[i].FieldMngtFile = CropInfo[3];
            //Assign default program properties (should not be changed without expert knowledge)
            paramStruct.crop[i].fshape_b = 13.8135; //Shape factor describing the reduction in biomass production for insufficient growing degree days
            paramStruct.crop[i].PctZmin = 70; //Initial percentage of minimum effective rooting depth
            paramStruct.crop[i].fshape_ex = -6; //Shape factor describing the effects of water stress on root expansion
            paramStruct.crop[i].ETadj = 1; //Adjustment to water stress thresholds depending on daily ET0 (0 = No, 1 = Yes)
            paramStruct.crop[i].Aer = 5; //Vol (%) below saturation at which stress begins to occur due to deficient aeration
            paramStruct.crop[i].LagAer = 3; //Number of days lag before aeration stress affects crop growth
            paramStruct.crop[i].beta = 12; //Reduction (%) to p_lo3 when early canopy senescence is triggered
            paramStruct.crop[i].a_Tr = 1; //Exponent parameter for adjustment of Kcx once senescence is triggered
            paramStruct.crop[i].GermThr = 0.2; //Proportion of total water storage needed for crop to germinate
            paramStruct.crop[i].CCmin = 0.05; //Minimum canopy size below which yield formation cannot occur
            paramStruct.crop[i].MaxFlowPct = 100 / 3; //Proportion of total flowering time (%) at which peak flowering occurs
            paramStruct.crop[i].HIini = 0.01; //Initial harvest index
            paramStruct.crop[i].bsted = 0.000138; //WP co2 adjustment parameter given by Steduto et al. 2007
            paramStruct.crop[i].bface = 0.001165; //WP co2 adjustment parameter given by FACE experiments
            paramStruct.crop[i].fsink /= 100; //Convert from %
        }

        double[] PlantDates = new double[1];
        double[] HarvestDates = new double[1];
        //Find planting and harvest dates
        if (nCrops > 1 || Rotation.compareTo("Y") == 0) {
            //Crop rotation occurs during the simulation period
            //Open rotation time-series file
            fileName = location.concat("\\" + fileLocation.cropFilename);
            try {
                //check the file exists
                new FileReader(fileName);
            } catch (FileNotFoundException e) {
                //Can't find text file defining locations of input and output folders.
                System.out.println(e.getMessage());
                return paramStruct;
            }

            //Load data
            dataArray = new LinkedList<>();
            try {
                BufferedReader br = new BufferedReader(new FileReader(new File(fileName)));
                String st;
                while ((st = br.readLine()) != null) {
                    dataArray.add(st);
                }
            } catch (IOException e) {
                System.out.println(e.getMessage());
            }
            //Extract data
            //TODO
//                PlantDates = datenum(DataArray{1,1},'dd/mm/yyyy');
//                String HarvestDates = datenum(DataArray{1,2},'dd/mm/yyyy');
//                int CropChoices = DataArray{1,3};
        } else if (nCrops == 1) {
            //Only one crop type considered during simulation - i.e. no rotations
            //either within or between yars
            //Get start and end years for full simulation
            String[] split = clockStruct.SimulationStartTime.split("-");
            String[] SimStaDate = {split[0], split[1], split[2]};
            split = clockStruct.SimulationEndTime.split("-");
            String[] SimEndDate = {split[0], split[1], split[2]};
            //Get temporary crop structure
            Crop CropTemp = paramStruct.crop[0];
            //Does growing season extend across multiple calendar years
            int[] PlantDatesSplit = {Integer.parseInt(CropTemp.PlantingDate.split("/")[0]), Integer.parseInt(CropTemp.PlantingDate.split("/")[1])};
            int[] HarvestDateSplit = {Integer.parseInt(CropTemp.HarvestDate.split("/")[0]), Integer.parseInt(CropTemp.HarvestDate.split("/")[1])};
            String[] YrsPlant = new String[1];
            String[] YrsHarvest = new String[1];
            if (PlantDatesSplit[1] < HarvestDateSplit[1] || (PlantDatesSplit[1] == HarvestDateSplit[1] && PlantDatesSplit[0] < HarvestDateSplit[0])) {
                ArrayList<String> YrsPlantTemp = new ArrayList<>();
                for (int i = Integer.parseInt(SimStaDate[0]); i <= Integer.parseInt(SimEndDate[0]); i++) {
                    YrsPlantTemp.add(String.valueOf(i));
                }
                YrsPlant = YrsPlantTemp.toArray(new String[YrsPlantTemp.size()]);
                YrsHarvest = YrsPlant.clone();
            } else {
                //TODO
//                YrsPlant = SimStaDate(1):SimEndDate(1) - 1;
//                YrsHarvest = SimStaDate(1) + 1:SimEndDate(1);
            }
//            Correct for partial first growing season (may occur when simulating
//            off-season soil water balance)
            //datenum(strcat(CropTemp.PlantingDate, '/', num2str(YrsPlant(1))), 'dd/mm/yyyy') < clockStruct.SimulationStartDate
            if (false) {
                //TODO
//                YrsPlant = YrsPlant(2:end);
//                YrsHarvest = YrsHarvest(2:end);
            }
            //Define blank variables
            PlantDates = new double[YrsPlant.length];
            HarvestDates = new double[YrsHarvest.length];
            String[] CropChoices = new String[YrsPlant.length];
            //Determine planting and harvest dates
            for (int ii = 0; ii < YrsPlant.length; ii++) {
                PlantDates[ii] = clockStruct.SimulationStartDate;
                //TODO check why -30
                HarvestDates[ii] = clockStruct.SimulationEndDate - 30;
                //TODO check it
                CropChoices[ii] = CropInfo[0];
            }
        }
        //Update clock parameters
        //Store planting and harvest dates
        clockStruct.PlantingDate = PlantDates;
        clockStruct.HarvestDate = HarvestDates;
        clockStruct.nSeasons = PlantDates.length;
//        //Initialise growing season counter
        if (clockStruct.StepStartTime == clockStruct.PlantingDate[0]) {
            clockStruct.SeasonCounter = 1;
        } else {
            clockStruct.SeasonCounter = 0;
        }
        return paramStruct;
    }

    public IrrMngtStruct AOS_ReadIrrigationManagement(FileLocation fileLocation, ParamStruct paramStruct) {
        //Read AOS input file location %%
        String location = fileLocation.input;

        //Read irrigation management input files %%
        //Check for number of crop types
        int nCrops = paramStruct.crop.length;
        //Create blank structure
        IrrMngtStruct irrMngtStruct = new IrrMngtStruct();

        for (int ii = 0; ii < nCrops; ii++) {
            //Open file
            String fileName = location.concat("\\" + paramStruct.crop[ii].IrrigationFile);
            try {
                //check the file exists
                new FileReader(fileName);
            } catch (FileNotFoundException e) {
                //Can't find text file defining locations of input and output folders.
                System.out.println(e.getMessage());
                return irrMngtStruct;
            }

            //Load data
            List<String> dataArray = new LinkedList<>();
            try {
                BufferedReader br = new BufferedReader(new FileReader(new File(fileName)));
                String st;
                while ((st = br.readLine()) != null) {
                    dataArray.add(st);
                }
            } catch (IOException e) {
                System.out.println(e.getMessage());
            }

            //Create and assign numeric variables
            irrMngtStruct.IrrMethod = Integer.parseInt(dataArray.get(2).split(",")[1]);
            irrMngtStruct.IrrInterval = Integer.parseInt(dataArray.get(3).split(",")[1]);
            irrMngtStruct.SMT1 = Integer.parseInt(dataArray.get(4).split(",")[1]);
            irrMngtStruct.SMT2 = Integer.parseInt(dataArray.get(5).split(",")[1]);
            irrMngtStruct.SMT3 = Integer.parseInt(dataArray.get(6).split(",")[1]);
            irrMngtStruct.SMT4 = Integer.parseInt(dataArray.get(7).split(",")[1]);
            irrMngtStruct.MaxIrr = Integer.parseInt(dataArray.get(8).split(",")[1]);
            irrMngtStruct.AppEff = Integer.parseInt(dataArray.get(9).split(",")[1]);
            irrMngtStruct.NetIrrSMT = Integer.parseInt(dataArray.get(10).split(",")[1]);
            irrMngtStruct.WetSurf = Integer.parseInt(dataArray.get(11).split(",")[1]);

            //Consolidate soil moisture targets in to one variable
            irrMngtStruct.SMT = new double[]{irrMngtStruct.SMT1, irrMngtStruct.SMT2, irrMngtStruct.SMT3, irrMngtStruct.SMT4};

            if (irrMngtStruct.IrrMethod == 3) {
                //TODO
            }
        }
        return irrMngtStruct;
    }

    //array operators
    private Double[] cumsum(Double[] in) {
        Double[] out = new Double[in.length];
        Double total = 0.0;
        for (int i = 0; i < in.length; i++) {
            total += in[i];
            out[i] = total;
        }
        return out;
    }

    private double[] cumsum(double[] in) {
        double[] out = new double[in.length];
        double total = 0.0;
        for (int i = 0; i < in.length; i++) {
            total += in[i];
            out[i] = total;
        }
        return out;
    }

    private Double[] dot(Double num, Double[] array) {
        Double[] newArray = new Double[array.length];
        for (int i = 0; i < array.length; i++) {
            newArray[i] = num * array[i];
        }
        return newArray;
    }

    private Double[] round(Double[] array) {
        Double[] newArray = new Double[array.length];
        for (int i = 0; i < array.length; i++) {
            Long temp = Math.round(array[i]);
            newArray[i] = temp.doubleValue();
        }
        return newArray;
    }

    private double[] round(double[] array) {
        double[] newArray = new double[array.length];
        for (int i = 0; i < array.length; i++) {
            Long temp = Math.round(array[i]);
            newArray[i] = temp.doubleValue();
        }
        return newArray;
    }

    private Double[] div(double num, Double[] array) {
        Double[] newArray = new Double[array.length];
        for (int i = 0; i < array.length; i++) {
            newArray[i] = array[i] / num;
        }
        return newArray;
    }
}