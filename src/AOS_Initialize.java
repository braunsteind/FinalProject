import Structs.ClockStruct;
import Structs.Crop;
import Structs.FileLocation;
import Structs.ParamStruct;

import java.io.*;
import java.sql.Date;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

public class AOS_Initialize {

    public static final int MILLISECONDS_IN_DAY = 86400000;
    public static final int DATE_ADD = 719530;

    public AOS_Initialize() {
        FileLocation fileLocation = AOS_ReadFileLocations();
        ClockStruct clockStruct = AOS_ReadClockParameters(fileLocation);
        double[][] WeatherStruct = AOS_ReadWeatherInputs(fileLocation, clockStruct);
        AOS_ReadModelParameters(fileLocation, clockStruct);
    }

    /**
     * Function to read and process input weather time-series
     */
    private double[][] AOS_ReadWeatherInputs(FileLocation fileLocation, ClockStruct clockStruct) {
        double[][] weatherDB = new double[clockStruct.nSteps][5];

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
        fileLocation.output = dataArray.get(2).split(",")[1];

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
        fileLocation.cO2Filename = dataArray.get(8).split(",")[1];
        fileLocation.outputFilename = dataArray.get(9).split(",")[1];
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

    private void AOS_ReadModelParameters(FileLocation fileLocation, ClockStruct clockStruct) {
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
            return;
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
        paramStruct.Soil.put("CalcSHP", Double.parseDouble(dataArray.get(4).split(",")[1]));
        paramStruct.Soil.put("Zsoil", Double.parseDouble(dataArray.get(5).split(",")[1]));
        paramStruct.Soil.put("nComp", Double.parseDouble(dataArray.get(6).split(",")[1]));
        paramStruct.Soil.put("nLayer", Double.parseDouble(dataArray.get(7).split(",")[1]));
        paramStruct.Soil.put("AdjREW", Double.parseDouble(dataArray.get(8).split(",")[1]));
        paramStruct.Soil.put("REW", Double.parseDouble(dataArray.get(9).split(",")[1]));
        paramStruct.Soil.put("CN", Double.parseDouble(dataArray.get(10).split(",")[1]));
        paramStruct.Soil.put("zRes", Double.parseDouble(dataArray.get(11).split(",")[1]));

        //Assign default program properties (should not be changed without expert knowledge)
        paramStruct.Soil.put("EvapZsurf", 0.04);  //Thickness of soil surface skin evaporation layer (m)
        paramStruct.Soil.put("EvapZmin", 0.15);   //Minimum thickness of full soil surface evaporation layer (m)
        paramStruct.Soil.put("EvapZmax", 0.30);   //Maximum thickness of full soil surface evaporation layer (m)
        paramStruct.Soil.put("Kex", 1.1);         //Maximum soil evaporation coefficient
        paramStruct.Soil.put("fevap", 4.0);         //Shape factor describing reduction in soil evaporation in stage 2.
        paramStruct.Soil.put("fWrelExp", 0.4);    //Proportional value of Wrel at which soil evaporation layer expands
        paramStruct.Soil.put("fwcc", 50.0);         //Maximum coefficient for soil evaporation reduction due to sheltering effect of withered canopy
        paramStruct.Soil.put("zCN", 0.3);         //Thickness of soil surface (m) used to calculate water content to adjust curve number
        paramStruct.Soil.put("zGerm", 0.3);       //Thickness of soil surface (m) used to calculate water content for germination
        paramStruct.Soil.put("AdjCN", 1.0);         //Adjust curve number for antecedent moisture content (0: No, 1: Yes)
        paramStruct.Soil.put("fshape_cr", 16.0);    //Capillary rise shape factor
        paramStruct.Soil.put("zTop", 0.1);        //Thickness of soil surface layer for water stress comparisons (m)

        //Read soil profile input file
        //Open file
        fileName = location.concat("\\" + fileLocation.soilProfileFilename);

        try {
            //check the file exists
            new FileReader(fileName);
        } catch (FileNotFoundException e) {
            //Can't find text file defining locations of input and output folders.
            System.out.println(e.getMessage());
            return;
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
        paramStruct.Comp.put("dz", dz.toArray(new Double[dz.size()]));
        Double[] dzsum = div(100, round(dot(100.0, cumsum(paramStruct.Comp.get("dz")))));
        paramStruct.Comp.put("dzsum", dzsum);
        ArrayList<Double> Layer = new ArrayList<>();
        for (int i = 1; i < dataArray.size(); i++) {
            Layer.add(Double.parseDouble(dataArray.get(i).split(",")[2]));
        }
        paramStruct.Comp.put("Layer", Layer.toArray(new Double[Layer.size()]));

        //Read crop mix input file
        //Open file
        fileName = location.concat("\\" + fileLocation.cropFilename);
        try {
            //check the file exists
            new FileReader(fileName);
        } catch (FileNotFoundException e) {
            //Can't find text file defining locations of input and output folders.
            System.out.println(e.getMessage());
            return;
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

        //declare
        int[] PlantDates = new int[1];
        int[] HarvestDates = new int[1];

        //Read crop parameter input files
        //Create blank structure

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
                return;
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
            paramStruct.crop = new Crop();
            paramStruct.crop.Emergence = Double.parseDouble(DataArray[0]);
            paramStruct.crop.MaxRooting = Double.parseDouble(DataArray[1]);
            paramStruct.crop.Senescence = Double.parseDouble(DataArray[2]);
            paramStruct.crop.Maturity = Double.parseDouble(DataArray[3]);
            paramStruct.crop.HIstart = Double.parseDouble(DataArray[4]);
            paramStruct.crop.Flowering = Double.parseDouble(DataArray[5]);
            paramStruct.crop.YldForm = Double.parseDouble(DataArray[6]);
            paramStruct.crop.GDDmethod = Double.parseDouble(DataArray[7]);
            paramStruct.crop.Tbase = Double.parseDouble(DataArray[8]);
            paramStruct.crop.Tupp = Double.parseDouble(DataArray[9]);
            paramStruct.crop.PolHeatStress = Double.parseDouble(DataArray[10]);
            paramStruct.crop.Tmax_up = Double.parseDouble(DataArray[11]);
            paramStruct.crop.Tmax_lo = Double.parseDouble(DataArray[12]);
            paramStruct.crop.PolColdStress = Double.parseDouble(DataArray[13]);
            paramStruct.crop.Tmin_up = Double.parseDouble(DataArray[14]);
            paramStruct.crop.Tmin_lo = Double.parseDouble(DataArray[15]);
            paramStruct.crop.TrColdStress = Double.parseDouble(DataArray[16]);
            paramStruct.crop.GDD_up = Double.parseDouble(DataArray[17]);
            paramStruct.crop.GDD_lo = Double.parseDouble(DataArray[18]);
            paramStruct.crop.Zmin = Double.parseDouble(DataArray[19]);
            paramStruct.crop.Zmax = Double.parseDouble(DataArray[20]);
            paramStruct.crop.fshape_r = Double.parseDouble(DataArray[21]);
            paramStruct.crop.SxTopQ = Double.parseDouble(DataArray[22]);
            paramStruct.crop.SxBotQ = Double.parseDouble(DataArray[23]);
            paramStruct.crop.SeedSize = Double.parseDouble(DataArray[24]);
            paramStruct.crop.PlantPop = Double.parseDouble(DataArray[25]);
            paramStruct.crop.CCx = Double.parseDouble(DataArray[26]);
            paramStruct.crop.CDC = Double.parseDouble(DataArray[27]);
            paramStruct.crop.CGC = Double.parseDouble(DataArray[28]);
            paramStruct.crop.Kcb = Double.parseDouble(DataArray[29]);
            paramStruct.crop.fage = Double.parseDouble(DataArray[30]);
            paramStruct.crop.WP = Double.parseDouble(DataArray[31]);
            paramStruct.crop.WPy = Double.parseDouble(DataArray[32]);
            paramStruct.crop.fsink = Double.parseDouble(DataArray[33]);
            paramStruct.crop.HI0 = Double.parseDouble(DataArray[34]);
            paramStruct.crop.dHI_pre = Double.parseDouble(DataArray[35]);
            paramStruct.crop.a_HI = Double.parseDouble(DataArray[36]);
            paramStruct.crop.b_HI = Double.parseDouble(DataArray[37]);
            paramStruct.crop.dHI0 = Double.parseDouble(DataArray[38]);
            paramStruct.crop.Determinant = Double.parseDouble(DataArray[39]);
            paramStruct.crop.exc = Double.parseDouble(DataArray[40]);
            paramStruct.crop.p_up1 = Double.parseDouble(DataArray[41]);
            paramStruct.crop.p_up2 = Double.parseDouble(DataArray[42]);
            paramStruct.crop.p_up3 = Double.parseDouble(DataArray[43]);
            paramStruct.crop.p_up4 = Double.parseDouble(DataArray[44]);
            paramStruct.crop.p_lo1 = Double.parseDouble(DataArray[45]);
            paramStruct.crop.p_lo2 = Double.parseDouble(DataArray[46]);
            paramStruct.crop.p_lo3 = Double.parseDouble(DataArray[47]);
            paramStruct.crop.p_lo4 = Double.parseDouble(DataArray[48]);
            paramStruct.crop.fshape_w1 = Double.parseDouble(DataArray[49]);
            paramStruct.crop.fshape_w2 = Double.parseDouble(DataArray[50]);
            paramStruct.crop.fshape_w3 = Double.parseDouble(DataArray[51]);
            paramStruct.crop.fshape_w4 = Double.parseDouble(DataArray[52]);

            //Add additional parameters
            paramStruct.crop.CropType = Double.parseDouble(CropType);
            paramStruct.crop.PlantMethod = Double.parseDouble(PlantMethod);
            paramStruct.crop.CalendarType = Double.parseDouble(CalendarType);
            paramStruct.crop.SwitchGDD = Double.parseDouble(SwitchGDD);
            paramStruct.crop.PlantingDate = PlantingDateStr;
            paramStruct.crop.HarvestDate = HarvestDateStr;
            //Add irrigation management information
            paramStruct.crop.IrrigationFile = CropInfo[2];
            paramStruct.crop.FieldMngtFile = CropInfo[3];
            //Assign default program properties (should not be changed without expert knowledge)
            paramStruct.crop.fshape_b = 13.8135; //Shape factor describing the reduction in biomass production for insufficient growing degree days
            paramStruct.crop.PctZmin = 70; //Initial percentage of minimum effective rooting depth
            paramStruct.crop.fshape_ex = -6; //Shape factor describing the effects of water stress on root expansion
            paramStruct.crop.ETadj = 1; //Adjustment to water stress thresholds depending on daily ET0 (0 = No, 1 = Yes)
            paramStruct.crop.Aer = 5; //Vol (%) below saturation at which stress begins to occur due to deficient aeration
            paramStruct.crop.LagAer = 3; //Number of days lag before aeration stress affects crop growth
            paramStruct.crop.beta = 12; //Reduction (%) to p_lo3 when early canopy senescence is triggered
            paramStruct.crop.a_Tr = 1; //Exponent parameter for adjustment of Kcx once senescence is triggered
            paramStruct.crop.GermThr = 0.2; //Proportion of total water storage needed for crop to germinate
            paramStruct.crop.CCmin = 0.05; //Minimum canopy size below which yield formation cannot occur
            paramStruct.crop.MaxFlowPct = 100 / 3; //Proportion of total flowering time (%) at which peak flowering occurs
            paramStruct.crop.HIini = 0.01; //Initial harvest index
            paramStruct.crop.bsted = 0.000138; //WP co2 adjustment parameter given by Steduto et al. 2007
            paramStruct.crop.bface = 0.001165; //WP co2 adjustment parameter given by FACE experiments
            paramStruct.crop.fsink /= 100; //Convert from %

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
                    return;
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
                PlantDates = datenum(DataArray {
                    1, 1
                },'dd/mm/yyyy');
//                String HarvestDates = datenum(DataArray{1,2},'dd/mm/yyyy');
//                int CropChoices = DataArray{1,3};
            } else if (nCrops == 1) {
                //Only one crop type considered during simulation - i.e. no rotations
                //either within or between yars
                //Get start and end years for full simulation
//                String SimStaDate = datevec(clockStruct.SimulationStartDate);
//                String SimEndDate = datevec(clockStruct.SimulationEndDate);
                //Get temporary crop structure
                Crop CropTemp = paramStruct.crop;
                //Does growing season extend across multiple calendar years
                if (datenum(CropTemp.PlantingDate, 'dd/mm') < datenum(CropTemp.HarvestDate, 'dd/mm')) {
//                    YrsPlant = SimStaDate(1):SimEndDate(1);
//                    YrsHarvest = YrsPlant;
                } else {
//                    YrsPlant = SimStaDate(1):SimEndDate(1) - 1;
//                    YrsHarvest = SimStaDate(1) + 1:SimEndDate(1);
                }
                //orrect for partial first growing season (may occur when simulating
                //off-season soil water balance)
                if (datenum(strcat(CropTemp.PlantingDate, '/', num2str(YrsPlant(1))), 'dd/mm/yyyy') < AOS_ClockStruct.SimulationStartDate) {
//                    YrsPlant = YrsPlant(2:end);
//                    YrsHarvest = YrsHarvest(2:end);
                }
                //Define blank variables
                PlantDates = zeros(length(YrsPlant), 1);
//                HarvestDates = zeros(length(YrsHarvest), 1);
                Crop CropChoices;
                //Determine planting and harvest dates
                for (int ii = 0; ii < YrsPlant; ii++) {
//                    PlantDates(ii) = datenum(strcat(CropTemp.PlantingDate, '/', num2str(YrsPlant(ii))), 'dd/mm/yyyy');
//                    HarvestDates(ii) = datenum(strcat(CropTemp.HarvestDate, '/', num2str(YrsHarvest(ii))), 'dd/mm/yyyy');
                    //TODO check it
                    CropChoices = paramStruct.crop;
                }
            }
        }
        //Update clock parameters %%
        //Store planting and harvest dates
        clockStruct.PlantingDate = PlantDates;
        clockStruct.HarvestDate = HarvestDates;
        clockStruct.nSeasons = PlantDates.length;
        //Initialise growing season counter
        if (clockStruct.StepStartTime == clockStruct.PlantingDate[0]) {
            clockStruct.SeasonCounter = 1;
        } else {
            clockStruct.SeasonCounter = 0;
        }
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

    private Double[] div(double num, Double[] array) {
        Double[] newArray = new Double[array.length];
        for (int i = 0; i < array.length; i++) {
            newArray[i] = array[i] / num;
        }
        return newArray;
    }
}