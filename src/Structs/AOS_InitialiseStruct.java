package Structs;

public class AOS_InitialiseStruct {
    public ParamStruct Parameter;
    public IrrMngtStruct IrrigationManagement;
    public FieldMngtStruct[] FieldManagement;
    public GwStruct Groundwater;
    //TODO change
    public InitCondStruct InitialCondition = new InitCondStruct();
    public Crop[] CropChoices;
    public double[][] Weather;
    public FileLocation FileLocation;
    public Outputs Outputs = new Outputs();
}