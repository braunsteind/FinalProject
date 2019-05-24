import Structs.AOS_InitialiseStruct;
import Structs.ClockStruct;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class AOS_Finish {

    //Write outputs
    public static void finish(AOS_InitialiseStruct AOS_InitialiseStruct, ClockStruct AOS_ClockStruct) {
        //Define output file location and name
        String FileLoc = AOS_InitialiseStruct.FileLocation.output;
        String FileName = AOS_InitialiseStruct.FileLocation.outputFilename;

        //Write outputs
        if (AOS_ClockStruct.ModelTermination) {
            if (AOS_InitialiseStruct.FileLocation.writeDaily.compareTo("Y") == 0) {
                //TODO
                //Water contents
//                dlmwrite(FileLoc.concat( FileName+"_WaterContents.csv")), AOS_InitialiseStruct.Outputs.WaterContents, '-append');
                //Water fluxes

                //Crop growth
            }

            //Final output
            String[] FinalOut = AOS_InitialiseStruct.Outputs.FinalOutput;
            String path = FileLoc.concat(FileName + "_FinalOutput.csv");
            write(path, FinalOut);
        }
    }

    private static void write(String path, String[] text) {
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
}