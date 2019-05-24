package Structs;

public class AOS_InitialiseStruct {
    public ParamStruct Parameter;
    public IrrMngtStruct IrrigationManagement;
    public FieldMngtStruct[] FieldManagement;
    public GwStruct Groundwater;
    public InitCondStruct InitialCondition;
    public Crop[] CropChoices;
    public double[][] Weather;
    public FileLocation FileLocation;
    public Outputs Outputs = new Outputs();
}