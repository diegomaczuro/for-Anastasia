#ifndef TESTMONODOMAIN3DEXAMPLETUTORIAL_HPP_
#define TESTMONODOMAIN3DEXAMPLETUTORIAL_HPP_


#include <cxxtest/TestSuite.h>
#include "MonodomainProblem.hpp"
#include "LuoRudy1991BackwardEuler.hpp"
#include "SimpleStimulus.hpp"
#include "TetrahedralMesh.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "BidomainProblem.hpp"


class CellFactory_measurement_of_speed : public AbstractCardiacCellFactory<1>
{
private:
    boost::shared_ptr<SimpleStimulus> mpStimulus;
public:
    CellFactory_measurement_of_speed()
    : AbstractCardiacCellFactory<1>(),
      mpStimulus(new SimpleStimulus(-100000, 2, 0 )) //nA/cm3, ms, ms 
    {
    }
    
    AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<1>* pNode)
    {
        double x = pNode->rGetLocation()[0];
        if (x < 1){
            return new CellLuoRudy1991FromCellMLBackwardEuler(mpSolver, mpStimulus);            
        }
        else{
            return new CellLuoRudy1991FromCellMLBackwardEuler(mpSolver, mpZeroStimulus); 
        }
        
        
    }

};



class TestClass_measurement_of_speed : public CxxTest::TestSuite
{
public:
    void Test_measurement_of_speed() throw (Exception)
    {
        DistributedTetrahedralMesh<1, 1> mesh;
        double h = 0.025; //размерноcть шага по сетке в см
        mesh.ConstructRegularSlabMesh(h, 21, 0.0, 0.0);
        HeartConfig::Instance()->SetOutputUsingOriginalNodeOrdering(true);
        HeartConfig::Instance()->SetSimulationDuration(500);//ms
        HeartConfig::Instance()->SetOutputDirectory("measurement_of_speed_output");
        HeartConfig::Instance()->SetOutputFilenamePrefix("results");
        HeartConfig::Instance()->SetVisualizeWithVtk(true);
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.6, 0.6, 0.6));//mS/cm
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(0.8, 0.8, 0.8));//mS/cm
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.05, 1.0);//ms
        HeartConfig::Instance()->SetCapacitance(0.7); // uF/cm^2
        

        CellFactory_measurement_of_speed cell_factory;
        BidomainProblem<1> bidomain_problem( &cell_factory );
        bidomain_problem.SetMesh(&mesh);
        bidomain_problem.SetWriteInfo();
        bidomain_problem.Initialise();
        bidomain_problem.Solve();
    }
};

#endif /*TESTMONODOMAIN3DEXAMPLETUTORIAL_HPP_*/