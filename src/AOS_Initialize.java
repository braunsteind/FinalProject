import java.io.*;
import java.util.LinkedList;
import java.util.List;

public class AOS_Initialize {
    private String fileLocation;

    public AOS_Initialize() {
        fileLocation = AOS_ReadFileLocations();
    }

    /**
     * Function to read input and output file locations.
     */
    private String AOS_ReadFileLocations() {
        //Read AOS file location input file
        String fileName = "FileLocations.csv";

        FileReader fr;
        try {
            fr = new FileReader(fileName);
        } catch (FileNotFoundException e) {
            //Can't find text file defining locations of input and output folders.
            System.out.println(e.getMessage());
            return "";
        }

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

        return fileName;
    }
}