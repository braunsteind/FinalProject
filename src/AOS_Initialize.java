import Structs.ClockStruct;
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
        AOS_ReadModelParameters(fileLocation);
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

    private void AOS_ReadModelParameters(FileLocation fileLocation) {
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
            //TODO add the rest of the file to DataArray
            String[] DataArray = dataArray.get(7).split(",");
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