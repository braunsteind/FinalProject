import Structs.AOS_InitialiseStruct;
import Structs.ClockStruct;

import java.io.*;

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
                    //TODO
                    ;
                    //Water contents
                    //                dlmwrite(FileLoc.concat( FileName+"_WaterContents.csv")), AOS_InitialiseStruct.Outputs.WaterContents, '-append');


                    FileWriter csvWriter = new FileWriter("C:\\Users\\ShlomiZ\\Desktop\\ProjectGmar\\FinalProject\\Output\\Sample\\_WaterContents.csv");
                    csvWriter.append("Year,Month,Day,SimDay,Season,GDD,TotGDD,Zr,CC,CCPot,Bio,BioPot,HI,HIadj,Yield");
                    csvWriter.append("\n");
                    int i, j;
                    String[][] data = parseToString(AOS_InitialiseStruct.Outputs.WaterContents);
                    for (i = 0; i < AOS_InitialiseStruct.Outputs.WaterContents.length; i++) {
                        for (j = 0; j < AOS_InitialiseStruct.Outputs.WaterContents[i].length; j++) {
                            csvWriter.append(data[i][j]);
                            csvWriter.append(",");
                        }
                        csvWriter.append("\n");
                    }
                    csvWriter.flush();
                    csvWriter.close();


                    String _WaterFluxesPath = "C:\\Users\\ShlomiZ\\Desktop\\ProjectGmar\\FinalProject\\Output\\Sample\\_WaterFluxes.csv";
                    csvWriter = new FileWriter(_WaterFluxesPath);
                    csvWriter.append("Year,Month,Day,SimDay,Season,wRZ,zGW,wSurf,Irr,Infl,RO,DP,CR,GWin,Es,EsX,Tr,TrX");
                    csvWriter.append("\n");
                    //Water fluxes
                    data = parseToString(AOS_InitialiseStruct.Outputs.WaterFluxes);
                    for (i = 0; i < AOS_InitialiseStruct.Outputs.WaterFluxes.length; i++) {
                        for (j = 0; j < AOS_InitialiseStruct.Outputs.WaterFluxes[i].length; j++) {
                            csvWriter.append(data[i][j]);
                            csvWriter.append(",");
                        }
                        csvWriter.append("\n");
                    }
                    csvWriter.flush();
                    csvWriter.close();


                    String CropGrowthPath = "C:\\Users\\ShlomiZ\\Desktop\\ProjectGmar\\FinalProject\\Output\\Sample\\CropGrowth.csv";
                    csvWriter = new FileWriter(CropGrowthPath);
                    //Crop growth
                    csvWriter.append("Year,Month,Day,SimDay,Season,GDD,TotGDD,Zr,CC,CCPot,Bio,BioPot,HI,HIadj,Yield");
                    csvWriter.append("\n");
                    data = parseToString(AOS_InitialiseStruct.Outputs.CropGrowth);
                    for (i = 0; i < AOS_InitialiseStruct.Outputs.CropGrowth.length; i++) {
                        for (j = 0; j < AOS_InitialiseStruct.Outputs.CropGrowth[i].length; j++) {
                            csvWriter.append(data[i][j]);
                            csvWriter.append(",");
                        }
                        csvWriter.append("\n");
                    }
                    csvWriter.flush();
                    csvWriter.close();

                } catch (Exception e) {
                    System.out.println(e);
                }
            }

            //Final output
            String[] FinalOut = AOS_InitialiseStruct.Outputs.FinalOutput;
            String path = FileLoc.concat(FileName + "_FinalOutput.csv");
            //write(path, FinalOut);

            try {
                FileWriter csvWriter = new FileWriter(path);
                csvWriter.append("Season,CropType,PlantDate,PlantSimDate,HarvestDate,HarvestSimDate,Yield,TotIrr");
                csvWriter.append("\n");
                csvWriter.append(AOS_InitialiseStruct.Outputs.FinalOutput[0]);
                csvWriter.append(",");
                csvWriter.append(AOS_InitialiseStruct.Outputs.FinalOutput[1]);
                csvWriter.append(",");
                csvWriter.append(AOS_InitialiseStruct.Outputs.FinalOutput[2]);
                csvWriter.append(",");
                csvWriter.append(AOS_InitialiseStruct.Outputs.FinalOutput[3]);
                csvWriter.append(",");
                csvWriter.append(AOS_InitialiseStruct.Outputs.FinalOutput[4]);
                csvWriter.append(",");
                csvWriter.append(AOS_InitialiseStruct.Outputs.FinalOutput[5]);
                csvWriter.append(",");
                csvWriter.append(AOS_InitialiseStruct.Outputs.FinalOutput[6]);
                csvWriter.append(",");
                csvWriter.append(AOS_InitialiseStruct.Outputs.FinalOutput[7]);

                csvWriter.flush();
                csvWriter.close();
            } catch (Exception e) {
                System.out.println(e);
            }
        }
    }





    /*
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
    } */
}