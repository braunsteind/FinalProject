import java.io.*;
import java.sql.Date;
import java.util.LinkedList;
import java.util.List;

public class AOS_Initialize {

    public static final int MILLISECONDS_IN_DAY = 86400000;

    public AOS_Initialize() {
        FileLocation fileLocation = AOS_ReadFileLocations();
        ClockStruct clockStruct = AOS_ReadClockParameters(fileLocation);
        int x = 5;
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
        clockStruct.SimulationStartDate = Date.valueOf(clockStruct.SimulationStartTime);
        //Simulation end time as serial date number
        clockStruct.SimulationEndDate = Date.valueOf(clockStruct.SimulationEndTime);
        //Time step (years)
        clockStruct.TimeStep = 1;
        //Total numbers of time steps (days)
        clockStruct.nSteps = (int) ((clockStruct.SimulationEndDate.getTime() - clockStruct.SimulationStartDate.getTime()) / MILLISECONDS_IN_DAY);

        //Time spans
        int[] TimeSpan = new int[clockStruct.nSteps];
//        TimeSpan[0] = clockStruct.SimulationStartDate;
//        TimeSpan[clockStruct.nSteps] = clockStruct.SimulationEndDate;
        TimeSpan[0] = 0;
        TimeSpan[clockStruct.nSteps - 1] = clockStruct.nSteps;
        for (int ss = 1; ss < clockStruct.nSteps; ss++) {
            TimeSpan[ss] = TimeSpan[ss - 1] + 1;
        }
        clockStruct.TimeSpan = TimeSpan;
        //Time at start of current time step
        clockStruct.StepStartTime = clockStruct.TimeSpan[clockStruct.TimeStepCounter];
        //Time at end of current time step
        clockStruct.StepEndTime = clockStruct.TimeSpan[clockStruct.TimeStepCounter + 1];
        //Number of time-steps (per day) for soil evaporation calculation
        clockStruct.EvapTimeSteps = 20;

        return clockStruct;
    }
}