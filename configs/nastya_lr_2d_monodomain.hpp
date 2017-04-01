#include <cxxtest/TestSuite.h>
#include "MonodomainProblem.hpp"
#include "GenericMeshReader.hpp"
#include "SimpleStimulus.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "FixedModifier.hpp"
#include "LuoRudy1991BackwardEuler.hpp"
#include "BidomainProblem.hpp"

// read csv
#include <fstream>
#include <string>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

/** TEMPORARY FOR DEBUGGING */
///\todo #2739
#include <unistd.h>
#include <sys/resource.h>
#include "Debug.hpp"
#include <mpi.h>


//Write
#include "VtkMeshWriter.hpp"

const int MAGIC_NUMBER = 361;//TODO Тут написать сколько строк в data.csv

class FreiburgHeartCellFactory : public AbstractCardiacCellFactory<2> // <3> here
{
private:
    boost::numeric::ublas::matrix<double> m;

public:
    FreiburgHeartCellFactory()
            : AbstractCardiacCellFactory<2>()
    {
        m = boost::numeric::ublas::matrix<double>(MAGIC_NUMBER, 9);

        std::ifstream infile("heart/test/nastya2/chaste/MONODOMAIN_2D/data.csv");
        double idx, x, y, z, d_Ks, amplitude, duration, time, model;
        int line_num = 0;
        while (infile >> idx >> x >> y >> z >> d_Ks >> amplitude >> duration >> time >> model)
        {
            m(line_num, 0) = idx;
            m(line_num, 1) = x;
            m(line_num, 2) = y;
            m(line_num, 3) = z;
            m(line_num, 4) = d_Ks;
            m(line_num, 5) = amplitude;
            m(line_num, 6) = duration;
            m(line_num, 7) = time;
            m(line_num, 8) = model;
            line_num++;
        }

        PRINT_VARIABLE( m(5, 0));
        PRINT_VARIABLE( m(5, 1));
        PRINT_VARIABLE( m(5, 2));
        PRINT_VARIABLE( m(5, 3));
        PRINT_VARIABLE( m(5, 4));
        PRINT_VARIABLE( m(5, 5));
        PRINT_VARIABLE( m(5, 6));
        PRINT_VARIABLE( m(5, 7));
        PRINT_VARIABLE( m(5, 8));

    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<2>* pNode)
    {
        double x = pNode->rGetLocation()[0];
        double y = pNode->rGetLocation()[1];
        //double z = pNode->rGetLocation()[2];
        unsigned int node_index = pNode->GetIndex();

        AbstractCardiacCell* cell = NULL;

        for(int idx=0; idx < MAGIC_NUMBER; idx++){
            if ((x - m(idx,1))*(x - m(idx,1)) + (y - m(idx,2))*(y - m(idx,2))  < 0.1*0.1) { //TODO Тут менять ход волокон
                boost::shared_ptr <SimpleStimulus> mpStimulus(new SimpleStimulus(-20000.0, 2.0, 0.0)); // mpStimulus(new SimpleStimulus(m(idx,5), m(idx,6), m(idx,7)));
                //TODO: Memory leak
            	cell = new CellLuoRudy1991FromCellMLBackwardEuler(mpSolver, mpStimulus);
                break;
            }
        }

        if (cell == NULL)
        {
            cell = new CellLuoRudy1991FromCellMLBackwardEuler(mpSolver, mpZeroStimulus);
        }
        PRINT_VARIABLE(node_index);
        return cell;
    }
};

class TestClass_nastya_tnnp06 : public CxxTest::TestSuite
{
public:
    double GetMemoryUsage()
    {
        struct rusage rusage;
        getrusage( RUSAGE_SELF, &rusage );

        return (double)(rusage.ru_maxrss)/(1024);// Convert KB to MB
    }

    void Test_nastya_tnnp06() throw(Exception)
    {
        int rank, size;
        MPI_Comm_rank (PETSC_COMM_WORLD, &rank);	/* get current process id */
        MPI_Comm_size (PETSC_COMM_WORLD, &size);	/* get number of processes */
        PRINT_VARIABLE(rank);
        PRINT_VARIABLE(size);


        HeartConfig::Instance()->SetMeshFileName("heart/test/nastya2/chaste/MONODOMAIN_2D/data", cp::media_type::Orthotropic);


        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(3.4, 0.6));//, 0.6));     
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(1.2, 0.8));//, 0.8));//mS/cm

        HeartConfig::Instance()->SetSimulationDuration(15);//500ms
        HeartConfig::Instance()->SetOutputDirectory("nastya_lr_2d_monodomain");
        HeartConfig::Instance()->SetOutputFilenamePrefix("results");
        HeartConfig::Instance()->SetVisualizeWithVtk(false);
        HeartConfig::Instance()->SetOutputUsingOriginalNodeOrdering(false);

        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.05, 1.);
        HeartConfig::Instance()->SetCapacitance(0.7); // uF/cm^2

        FreiburgHeartCellFactory cell_factory;
        MonodomainProblem<2> monodomain_problem( &cell_factory );
        TRACE("monodomain_problem(.) complete!");

        monodomain_problem.SetWriteInfo();
        TRACE("SetWriteInfo() complete!");
        monodomain_problem.Initialise();
        TRACE("Initialise() complete!");
        monodomain_problem.Solve();


        AbstractTetrahedralMesh<2,2>* p_mesh = &(monodomain_problem.rGetMesh());
        VtkMeshWriter<2,2> vtk_writer("./nastya_lr_2d_monodomain/", "init_mesh_1", false);
        vtk_writer.WriteFilesUsingMesh(*p_mesh);

        VtkMeshWriter<2,2> vtk_writer2("./nastya_lr_2d_monodomain/", "init_mesh_2", false);
        std::string original_file = p_mesh->GetMeshFileBaseName();
        std::auto_ptr<AbstractMeshReader<2, 2> > p_original_mesh_reader = GenericMeshReader<2, 2>(original_file);
        vtk_writer2.WriteFilesUsingMeshReader(*p_original_mesh_reader);



    }
};
