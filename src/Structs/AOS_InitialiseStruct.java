public class InitialiseStrcut {

    public ParameterStruct Parameter;
    public IrrigManagementStrcut IrrigationManagement;
    public FieldManagementStruct FieldManagement;
    public GroundwaterStruct Groundwater;
    public InitialConditionStruct InitialCondition;
    public CropChoices cc; // in matlab CropChoices is 1x1 cell, need to be checked
    public double[][] Weather; //was 183X5
    public FileLocationStruct FileLocation;
    public OutputsStruct Outputs;
}