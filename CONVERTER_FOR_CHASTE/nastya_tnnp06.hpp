#include <cxxtest/TestSuite.h>
#include "MonodomainProblem.hpp"
#include "ten_tusscher_model_2006_IK1Ko_epi_units.hpp"
#include "GenericMeshReader.hpp"
#include "SimpleStimulus.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "FixedModifier.hpp"
#include "TenTusscher2006EpiBackwardEuler.hpp"

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

const int MAGIC_NUMBER = 31779+1;//TODO Тут написать сколько строк в data.csv

class FreiburgHeartCellFactory : public AbstractCardiacCellFactory<3> // <3> here
{
private:
    boost::numeric::ublas::matrix<double> m;

public:
    FreiburgHeartCellFactory()
            : AbstractCardiacCellFactory<3>()
    {
        m = boost::numeric::ublas::matrix<double>(MAGIC_NUMBER, 9);

        std::ifstream infile("heart/test/nastya/chaste/data.csv");
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

    AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<3>* pNode)
    {
        double x = pNode->rGetLocation()[0];
        double y = pNode->rGetLocation()[1];
        double z = pNode->rGetLocation()[2];
        unsigned int node_index = pNode->GetIndex();

        AbstractCardiacCell* cell = NULL;

        for(int idx=0; idx < MAGIC_NUMBER; idx++){
            if ((x - m(idx,1))*(x - m(idx,1)) + (y - m(idx,2))*(y - m(idx,2)) + (z - m(idx,3))*(z - m(idx,3)) < 0.15*0.15) { //TODO Тут менять ход волокон
                boost::shared_ptr <SimpleStimulus> mpStimulus(new SimpleStimulus(-50000, 2.0, 0.0));
                //TODO: Memory leak
            	cell = new CML_tentusscher_panfilov_2006_epi_cell_pe_lut_be(mpSolver, mpStimulus);
            	cell->SetParameter("ScaleFactorGks", m(idx,4)/0.392);
                break;
            }
        }

        if (cell == NULL)
        {
            cell = new CML_tentusscher_panfilov_2006_epi_cell_pe_lut_be(mpSolver, mpZeroStimulus);
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


        HeartConfig::Instance()->SetMeshFileName("heart/test/nastya/chaste/data", cp::media_type::Axisymmetric);


        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(8.57, 3.87, 3.87));

        HeartConfig::Instance()->SetSimulationDuration(500);
        HeartConfig::Instance()->SetOutputDirectory("nastya_tnnp06");
        HeartConfig::Instance()->SetOutputFilenamePrefix("results");
        HeartConfig::Instance()->SetVisualizeWithVtk(false);
        HeartConfig::Instance()->SetOutputUsingOriginalNodeOrdering(false);

        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.05, 1.);

        FreiburgHeartCellFactory cell_factory;
        MonodomainProblem<3> monodomain_problem( &cell_factory );
        TRACE("monodomain_problem(.) complete!");
        //int number_of_cells = cell_factory.GetNumberOfCells();
        //PRINT_VARIABLE(10);

        monodomain_problem.SetWriteInfo();
        TRACE("SetWriteInfo() complete!");
        monodomain_problem.Initialise();
        TRACE("Initialise() complete!");
        monodomain_problem.Solve();


        /*
         * // Convert simulation data to VTK format
            Hdf5ToVtkConverter<ELEMENT_DIM,SPACE_DIM> converter(test_output,
                    HeartConfig::Instance()->GetOutputFilenamePrefix(),
                    mpMesh,
                    false,
                    HeartConfig::Instance()->GetOutputUsingOriginalNodeOrdering());
            std::string subdirectory_name = converter.GetSubdirectory();
            HeartConfig::Instance()->Write(false, subdirectory_name);
         *
         * */

        AbstractTetrahedralMesh<3,3>* p_mesh = &(monodomain_problem.rGetMesh());
        VtkMeshWriter<3,3> vtk_writer("./nastya_tnnp06/", "init_mesh_1", false);
        vtk_writer.WriteFilesUsingMesh(*p_mesh);

        VtkMeshWriter<3,3> vtk_writer2("./nastya_tnnp06/", "init_mesh_2", false);
        std::string original_file = p_mesh->GetMeshFileBaseName();
        std::auto_ptr<AbstractMeshReader<3, 3> > p_original_mesh_reader = GenericMeshReader<3, 3>(original_file);
        vtk_writer2.WriteFilesUsingMeshReader(*p_original_mesh_reader);


        //AbstractTetrahedralMesh<3,3>* p_mesh = &(monodomain_problem.rGetMesh());

        /** todo #2739
        double before_init = GetMemoryUsage();
        monodomain_problem.Initialise();
        double after_init = GetMemoryUsage();
        PRINT_VARIABLE(after_init - before_init);
        */

        //for (AbstractTetrahedralMesh<3,3>::NodeIterator iter = p_mesh->GetNodeIteratorBegin();
        //     iter != p_mesh->GetNodeIteratorEnd(); ++iter)
        //{
        //    TS_ASSERT_EQUALS(HeartRegionCode::IsRegionBath( iter->GetRegion() ), false);
        //}


    }
};
