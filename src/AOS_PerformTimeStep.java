import Structs.InitialiseStrcut;
import Structs.SetupSolutionStruct;
import Structs.SolutionReturnStruct;

public class AOS_PerformTimeStep {

    // should be global - need to be checked
    public InitialiseStrcut aos_is;

    public void PerformTimeStep() {

        // original code - lines 8-9
        //Setup parameters for time-step solution
        SetupSolutionStruct setSolStruct = new SetupSolutionStruct();
        AOS_SetupSolution setSol = new AOS_SetupSolution();
        setSol.run(setSolStruct);

        // original code - lines 12-13
        //Get model solution for current time-step
        AOS_Solution sol = new AOS_Solution(setSolStruct);
        SolutionReturnStruct srs = sol.run();

        //editing the global variable - maybe from AOS_Initialize
        //Update initial conditions and outputs
        aos_is.InitialCondition = srs.NewCond;
        aos_is.Outputs = srs.Outputs;

        //Check model termination
        AOS_CheckModelTermination();
        //Update time step
        AOS_UpdateTime();

    }
}


