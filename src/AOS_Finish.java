import Structs.AOS_InitialiseStruct;
import Structs.ClockStruct;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;

public class AOS_Finish {


    private static String[][] parseToString(double[][] arr) {
        String[][] data = new String[arr.length][arr[0].length];
        for (int i = 0; i < arr.length; i++) {
            for (int j = 0; j < arr[0].length; j++) {
                String str = Double.toString(arr[i][j]);
                data[i][j] = str;
            }
        }
        return data;
    }


    //Write outputs
    public static void finish(AOS_InitialiseStruct AOS_InitialiseStruct, ClockStruct AOS_ClockStruct) {
        //Define output file location and name
        String FileLoc = AOS_InitialiseStruct.FileLocation.output;
        String FileName = AOS_InitialiseStruct.FileLocation.outputFilename;

        //Write outputs
        if (AOS_ClockStruct.ModelTermination) {
            if (AOS_InitialiseStruct.FileLocation.writeDaily.compareTo("Y") == 0) {
                try {
                    //Water contents
                    String[][] data = parseToString(AOS_InitialiseStruct.Outputs.WaterContents);
                    StringBuilder text = new StringBuilder();
                    text.append("\n");
                    for (int i = 0; i < AOS_InitialiseStruct.Outputs.WaterContents.length; i++) {
                        for (int j = 0; j < AOS_InitialiseStruct.Outputs.WaterContents[i].length; j++) {
                            text.append(data[i][j]);
                            text.append(",");
                        }
                        text.append("\n");
                    }
                    Files.write(Paths.get(FileLoc.concat(FileName.concat("_WaterContents.csv"))), text.toString().getBytes(), StandardOpenOption.APPEND);

                    //Water fluxes
                    data = parseToString(AOS_InitialiseStruct.Outputs.WaterFluxes);
                    text = new StringBuilder();
                    text.append("\n");
                    for (int i = 0; i < AOS_InitialiseStruct.Outputs.WaterFluxes.length; i++) {
                        for (int j = 0; j < AOS_InitialiseStruct.Outputs.WaterFluxes[i].length; j++) {
                            text.append(data[i][j]);
                            text.append(",");
                        }
                        text.append("\n");
                    }
                    Files.write(Paths.get(FileLoc.concat(FileName.concat("_WaterFluxes.csv"))), text.toString().getBytes(), StandardOpenOption.APPEND);

                    //Crop growth
                    data = parseToString(AOS_InitialiseStruct.Outputs.CropGrowth);
                    text = new StringBuilder();
                    text.append("\n");
                    for (int i = 0; i < AOS_InitialiseStruct.Outputs.CropGrowth.length; i++) {
                        for (int j = 0; j < AOS_InitialiseStruct.Outputs.CropGrowth[i].length; j++) {
                            text.append(data[i][j]);
                            text.append(",");
                        }
                        text.append("\n");
                    }
                    Files.write(Paths.get(FileLoc.concat(FileName.concat("_CropGrowth.csv"))), text.toString().getBytes(), StandardOpenOption.APPEND);

                } catch (Exception e) {
                    System.out.println(e.getMessage());
                }
            }

            try {
                //Final output
                StringBuilder text = new StringBuilder();
                text.append("\n");
                for (int i = 0; i < AOS_InitialiseStruct.Outputs.FinalOutput.length; i++) {
                    text.append(AOS_InitialiseStruct.Outputs.FinalOutput[i]);
                    text.append(",");
                }
                Files.write(Paths.get(FileLoc.concat(FileName.concat("_FinalOutput.csv"))), text.toString().getBytes(), StandardOpenOption.APPEND);
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }
}