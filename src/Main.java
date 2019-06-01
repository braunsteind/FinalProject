public class Main {
    public static void main(String[] args) {
        //Run model
        //Initialise simulation
        AOS_Initialize aos_initialize = new AOS_Initialize();

        //Perform single time-step (day)
        while (!aos_initialize.clockStruct.ModelTermination) {
            //TODO
            AOS_PerformTimeStep.PerformTimeStep(aos_initialize);
        }

        //Finish simulation
        AOS_Finish.finish(aos_initialize.AOS_InitialiseStruct, aos_initialize.clockStruct);
    }
}