// Copyright (c) 2002-2014, Boyce Griffith
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of New York University nor the names of its
//      contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

// Config files
//#include <IBAMR_config.h>
//#include <IBTK_config.h>
#include <SAMRAI_config.h>
//#include <cmath>   <<--- may not need to include this!

// Headers for basic PETSc functions
#include <petscsys.h>

// Headers for basic SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// Headers for basic libMesh objects
#include <libmesh/boundary_info.h>
#include <libmesh/equation_systems.h>
#include <libmesh/exodusII_io.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IBFEMethod.h>
#include <ibamr/INSCollocatedHierarchyIntegrator.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>
#include <ibamr/app_namespaces.h>
#include <ibtk/AppInitializer.h>
#include <ibtk/libmesh_utilities.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

//--------------------------------------------------------------------
//
// Elasticity model data by AP Hoover and LA Miller
//
//--------------------------------------------------------------------
namespace ModelData
{

    // Problem parameters.
    static const double mu = 10.0; //Elastic modulus
    static double kappa = 5.0e6;   // body force spring constant
    static double J_maximum=1;
    static double J_minimum=1;

    /**********************************************************************
     *
     * TARGET FORCE FUNCTION FOR MESH 1
     *
     ***********************************************************************/
    void
    target1_force_function(
        VectorValue<double>& F,
        const TensorValue<double>& /*FF*/,
        const libMesh::Point& X,
        const libMesh::Point& s,
        Elem* const /*elem*/,
        const std::vector<const std::vector<double>*>& /*var_data*/,
        const std::vector<const std::vector<VectorValue<double> >*>& /*grad_var_data*/,
        double time,
        void* /*ctx*/)
    {
       libMesh::Point s_dump;
       //double kappa;

        // s_dump: the target point position 
        // s: is the reference configuration. 
        // X: current position 

        // Use modular arithmetic to get effective time in [0,2]
        //double tTilde = time % 2;    <<--- Modulo operator only does ints in C++
        double tTilde = std::fmod(time,4);

        // Max Speed and Time-step
        double dt = 1.25e-4;
        double LTube = 0.225;     //Accounting for side bars; amount of tube to travel
        double maxSpeed = 0.3375; // MaxSpeed = 3/2*LTube 

        if(tTilde<1)      
        {
            /*------------------------------------------------------
            * SPECIFIES NEXT GEOMETRIC STATE BASED OFF PREVIOUS
            /*------------------------------------------------------*/
            //s_dump(0)=X(0)+ 4.0*(maxSpeed)*(tTilde)*(1.0-tTilde)*dt; // How to move geometry in x-Direction
            //s_dump(1)=X(1);                          // How to move geoemtry in y-Direction

            /*--------------------------------------------------------------
            * SPECIFIES NEXT GEOMETRIC STATE BASED OFF ORIGINAL STATE
            *       -> Distance from starting = INT_0^t v(t) dt
            *       -> v(t) = 4*Smax*(t-t^2)
            *       -> d(t) = 4*Smax*( 0.5*t^2 - 1/3*t^3 )
            /*--------------------------------------------------------------*/
            s_dump(0)=s(0)+ 4.0*(maxSpeed)*( 0.5*pow(tTilde,2) - 1.0/3.0*pow(tTilde,3) ); // How to move geometry in x-Direction
            s_dump(1)=s(1);     // How to move geoemtry in y-Direction
            s_dump(2)=s(2);     // How to move geometry in z-Direction
        }

        else
        {
            //s_dump(1)=X(1);    // Let's geometry just 'flow' from previous momentum
            //s_dump(0)=X(0);    // Let's geometry just 'flow' from previous momentum	      
                    
            s_dump(0)=s(0)+LTube;    // Specify where geometry should be (x)
            s_dump(1)=s(1);          // Specify where geometry should be (y)
            s_dump(2)=s(2);          // Specify where geometry should be (z)
        }
    

            F = kappa*(s_dump-X);  // target force

            return;

        }
    
        /*********************************************************************
        *
        * TARGET FORCE FUNCTION FOR MESH 2
        *
        ********************************************************************/
        void
        target2_force_function(
            VectorValue<double>& F,
            const TensorValue<double>& /*FF*/,
            const libMesh::Point& X,
            const libMesh::Point& s,
            Elem* const /*elem*/,
            const std::vector<const std::vector<double>*>& /*var_data*/,
            const std::vector<const std::vector<VectorValue<double> >*>& /*grad_var_data*/,
            double time,
            void* /*ctx*/)
        {
            libMesh::Point s_dump;
            //double kappa;

            // s_dump: the target point position 
            // s: is the reference configuration. 
            // X: current position 

            // Use modular arithmetic to get effective time in [0,2]
            //double tTilde = time % 2;    <<--- Modulo operator only does ints in C++
            double tTilde = std::fmod(time,4);

            // Max Speed and Time-step
            double dt = 1.25e-4;
            double LTube = 0.225;     //Accounting for side bars; amount of tube to travel
            double maxSpeed = 0.3375; // MaxSpeed = 3/2*LTube 

            if(tTilde<1)      
            {
                /*------------------------------------------------------
                * SPECIFIES NEXT GEOMETRIC STATE BASED OFF PREVIOUS
                /*------------------------------------------------------*/
                //s_dump(0)=X(0)+ 4.0*(maxSpeed)*(tTilde)*(1.0-tTilde)*dt; // How to move geometry in x-Direction
                //s_dump(1)=X(1);                          // How to move geoemtry in y-Direction

                /*--------------------------------------------------------------
                * SPECIFIES NEXT GEOMETRIC STATE BASED OFF ORIGINAL STATE
                *       -> Distance from starting = INT_0^t v(t) dt
                *       -> v(t) = 4*Smax*(t-t^2)
                *       -> d(t) = 4*Smax*( 0.5*t^2 - 1/3*t^3 )
                /*--------------------------------------------------------------*/
                s_dump(0)=s(0) - 4.0*(maxSpeed)*( 0.5*pow(tTilde,2) - 1.0/3.0*pow(tTilde,3) ); // How to move geometry in x-Direction
                s_dump(1)=s(1);     // How to move geoemtry in y-Direction
                s_dump(2)=s(2);     // How to move geometry in z-Direction
            }

            else
            {
                //s_dump(1)=X(1);    // Let's geometry just 'flow' from previous momentum
                //s_dump(0)=X(0);    // Let's geometry just 'flow' from previous momentum	      
                        
                s_dump(0)=s(0)-LTube;    // Specify where geometry should be (x)
                s_dump(1)=s(1);          // Specify where geometry should be (y)
                s_dump(2)=s(2);          // Specify where geometry should be (z)
            }
        

            F = kappa*(s_dump-X);  // target force

            return;
        }

        /*--------------------------------------------------------------------
        * TARGET FORCE FUNCTION 3: KEEP GEOMETRY STATIONARY!
        *----------------------------------------------------------------------*/
        void
        target3_force_function(
        VectorValue<double>& F,
            const TensorValue<double>& /*FF*/,
            const libMesh::Point& X,
            const libMesh::Point& s,
            Elem* const /*elem*/,
            const std::vector<const std::vector<double>*>& /*var_data*/,
            const std::vector<const std::vector<VectorValue<double> >*>& /*grad_var_data*/,
            double time,
        void* /*ctx*/)
        {

            // If trying to prescribe movement of geometry; defining points to use
            //libMesh::Point s_dump;

            // Target Point Force based on "Linear Spring"  
            F = kappa*(s-X);
            
            return;
        
        }

        /*********************************************************************
        *
        * DEVIATORIC STRESS FUNCTION
        *
        ********************************************************************/    
        void
        PK1_dev_stress_function(
            TensorValue<double>& PP,
            const TensorValue<double>& FF,
            const libMesh::Point& /*X*/,
            const libMesh::Point& /*s*/,
            Elem* const /*elem*/,
            const std::vector<const std::vector<double>*>& /*var_data*/,
            const std::vector<const std::vector<VectorValue<double> >*>& /*grad_var_data*/,
            double /*time*/,
            void* /*ctx*/)
        {
            
            //double mu;
            
            PP = (mu)*FF;
            
            return;
            
        }// PK1_dev_stress_function

        /*********************************************************************
        *
        * DILATATIONAL STRESS FUNCTION 
        *
        ********************************************************************/        
        void
        PK1_dil_stress_function(
            TensorValue<double>& PP,
            const TensorValue<double>& FF,
            const libMesh::Point& /*X*/,
            const libMesh::Point& /*s*/,
            Elem* const /*elem*/,
            const std::vector<const std::vector<double>*>& /*var_data*/,
            const std::vector<const std::vector<VectorValue<double> >*>& /*grad_var_data*/,
            double /*time*/,
            void* /*ctx*/)
        {
            double mu;
            
            // From Wing Example -> Dilational Stress
            //PP = (-1.0*(mu)+beta*log(FF.det()))*tensor_inverse_transpose(FF,NDIM);
                
            // Original Dilational Stress
            PP = -mu*tensor_inverse_transpose(FF,NDIM);

            return;
            
        }// PK1_dil_stress_function

} // Ends ModelData

using namespace ModelData;

// Function prototypes
static ofstream forcex_stream, forcey_stream, forcez_stream, velx_stream, vely_stream, velz_stream, workx_stream, worky_stream, workz_stream, U_L1_norm_stream, U_L2_norm_stream, U_max_norm_stream;
void postprocess_data(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                      Pointer<INSHierarchyIntegrator> navier_stokes_integrator,
                      Mesh& mesh,
                      EquationSystems* equation_systems,
                      const int iteration_num,
                      const double loop_time,
                      const string& data_dump_dirname);

/*******************************************************************************
 * For each run, the input filename and restart information (if needed) must   *
 * be given on the command line.  For non-restarted case, command line is:     *
 *                                                                             *
 *    executable <input file name>                                             *
 *                                                                             *
 * For restarted run, command line is:                                         *
 *                                                                             *
 *    executable <input file name> <restart directory> <restart number>        *
 *                                                                             *
 *******************************************************************************/
int
main(
    int argc,
    char* argv[])
{
    // Initialize libMesh, PETSc, MPI, and SAMRAI.
    LibMeshInit init(argc, argv);
    SAMRAI_MPI::setCommunicator(PETSC_COMM_WORLD);
    SAMRAI_MPI::setCallAbortInSerialInsteadOfExit();
    SAMRAIManager::startup();

    {// cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "IB.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Get various standard options set in the input file.
        const bool dump_viz_data = app_initializer->dumpVizData();
        const int viz_dump_interval = app_initializer->getVizDumpInterval();
        const bool uses_visit = dump_viz_data && app_initializer->getVisItDataWriter();
        const bool uses_exodus = dump_viz_data && !app_initializer->getExodusIIFilename().empty();
        
        /*-------------------------------------------------------------------------------
         *
         * MAKE SURE IS CORRECT FOR READING IN APPROPRIATE NUMBER OF MESHES (NAB)
        *
        *-------------------------------------------------------------------------------*/
        const string mesh1_exodus_filename = app_initializer->getExodusIIFilename("mesh1");
        const string mesh2_exodus_filename = app_initializer->getExodusIIFilename("mesh2");
        const string mesh3_exodus_filename = app_initializer->getExodusIIFilename("mesh3");
        const string mesh4_exodus_filename = app_initializer->getExodusIIFilename("mesh4");


        const bool dump_restart_data = app_initializer->dumpRestartData();
        const int restart_dump_interval = app_initializer->getRestartDumpInterval();
        const string restart_dump_dirname = app_initializer->getRestartDumpDirectory();
        const string restart_read_dirname = app_initializer->getRestartReadDirectory();
        const int restart_restore_num = app_initializer->getRestartRestoreNumber();

        const bool dump_postproc_data = app_initializer->dumpPostProcessingData();
        
        const int postproc_data_dump_interval = app_initializer->getPostProcessingDataDumpInterval();
        const string postproc_data_dump_dirname = app_initializer->getPostProcessingDataDumpDirectory();
        if (dump_postproc_data && (postproc_data_dump_interval > 0) && !postproc_data_dump_dirname.empty())
        {
            Utilities::recursiveMkdir(postproc_data_dump_dirname);
        }


        const bool dump_timer_data = app_initializer->dumpTimerData();
        const int timer_dump_interval = app_initializer->getTimerDumpInterval();
		
    
        const double dx = input_db->getDouble("DX");
        const double ds = input_db->getDouble("MFAC")*dx;
        //const string mesh_filename = input_db->getStringWithDefault("MESH_FILENAME","Leaf.e");
        //string elem_type = input_db->getString("ELEM_TYPE");

        /*-------------------------------------------------------------------------------
        *
        * READING IN MULTIPLE MESHES!!!!
        *       --> Note that boundary condition data must be registered with each FE
        *              system before calling IBFEMethod::initializeFEData().
        *
        *-------------------------------------------------------------------------------*/
        //Mesh mesh1(NDIM);
        //Mesh mesh2(NDIM);
        Mesh mesh1(init.comm(), NDIM);
        Mesh mesh2(init.comm(), NDIM);	
        Mesh mesh3(init.comm(), NDIM);
        Mesh mesh4(init.comm(), NDIM);
        mesh1.read("PlugLeft.e");
        mesh1.prepare_for_use();
        mesh2.read("PlugRight.e");
        mesh2.prepare_for_use();
        mesh3.read("PistonLeft.e");
        mesh3.prepare_for_use();
        mesh4.read("PistonRight.e");
        mesh4.prepare_for_use();
        
        /*---------------------------------
         * ORIGINAL (deprecated)!!!
         ---------------------------------*/
        //vector<Mesh*> meshes(2);
        //meshes[0] = &mesh1;
        //meshes[1] = &mesh2;

        /*-------------------------------------------------------------------------------
         *
         * NEEDED TO UPDATE FOR MULTIPLES MESHES FOR v0.13.0 (June 2023) by NAB
         * (fixed with AP Hoover's assistance from IBAMR Google Group discussion, see:
         *  https://groups.google.com/g/ibamr-users/c/ZrLq6NiB8NE/m/-GIrJHhRAgAJ  )
         *
         *-------------------------------------------------------------------------------*/
        int num_meshes = 4;
        vector<MeshBase*> meshes(num_meshes);
        meshes[0] = &mesh1;
        meshes[1] = &mesh2;
        meshes[2] = &mesh3;
        meshes[3] = &mesh4;

        for (int part = 0; part < num_meshes; ++part)
            {
              MeshBase* mesh = meshes[part];
              for (auto& elem : mesh->element_ptr_range())
              {
                elem->subdomain_id() = part;
              }
         }


        
        Pointer<INSHierarchyIntegrator> navier_stokes_integrator;
        const string solver_type = app_initializer->getComponentDatabase("Main")->getString("solver_type");
        if (solver_type == "STAGGERED")
        {
            navier_stokes_integrator = new INSStaggeredHierarchyIntegrator(
                "INSStaggeredHierarchyIntegrator", app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));
        }
        else if (solver_type == "COLLOCATED")
        {
            navier_stokes_integrator = new INSCollocatedHierarchyIntegrator(
                "INSCollocatedHierarchyIntegrator", app_initializer->getComponentDatabase("INSCollocatedHierarchyIntegrator"));
        }
        else
        {
            TBOX_ERROR("Unsupported solver type: " << solver_type << "\n" <<
                       "Valid options are: COLLOCATED, STAGGERED");
        }
		Pointer<IBFEMethod> ib_method_ops =
            new IBFEMethod("IBFEMethod",
                           app_initializer->getComponentDatabase("IBFEMethod"),
                           meshes,
                           app_initializer->getComponentDatabase("GriddingAlgorithm")->getInteger("max_levels"),
                           /*register_for_restart*/ true,
                           restart_read_dirname,
                           restart_restore_num);
        Pointer<IBHierarchyIntegrator> time_integrator = new IBExplicitHierarchyIntegrator(
            "IBHierarchyIntegrator", app_initializer->getComponentDatabase("IBHierarchyIntegrator"), ib_method_ops, navier_stokes_integrator);
        Pointer<CartesianGridGeometry<NDIM> > grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM> > patch_hierarchy = new PatchHierarchy<NDIM>(
            "PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM> > error_detector = new StandardTagAndInitialize<NDIM>(
            "StandardTagAndInitialize", time_integrator, app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM> > box_generator = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM> > load_balancer = new LoadBalancer<NDIM>(
            "LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM> > gridding_algorithm = new GriddingAlgorithm<NDIM>(
            "GriddingAlgorithm", app_initializer->getComponentDatabase("GriddingAlgorithm"), error_detector, box_generator, load_balancer);

        /*-------------------------------------------------------------------------------
         * CONFIGURE THE IBFE SOLVER 
         *-------------------------------------------------------------------------------*/
	    ib_method_ops->initializeFEEquationSystems();
		//
        IBFEMethod::LagBodyForceFcnData target1_force_data(target1_force_function);
        IBFEMethod::LagBodyForceFcnData target2_force_data(target2_force_function);  
        IBFEMethod::LagBodyForceFcnData target3_force_data(target3_force_function);  
		IBFEMethod::PK1StressFcnData PK1_dev_stress_data(PK1_dev_stress_function);
        IBFEMethod::PK1StressFcnData PK1_dil_stress_data(PK1_dil_stress_function);

       /*-------------------------------------------------------------------------------
        *
        * SPECIFY FORCE FUNCTIONS ON THE DIFFERENT FINITE ELEMENT MESHES (NAB)
        *    -> Target Force 1: Movements Left Wall Through Channel
        *    -> Target Force 2: Keeps channel walls rigid
        *
        *-------------------------------------------------------------------------------*/
        // TARGET FORCES
        ib_method_ops->registerLagBodyForceFunction(target1_force_data, 0); // Mesh 1 (plug left)
        ib_method_ops->registerLagBodyForceFunction(target2_force_data, 1); // Mesh 2 (plug right)
        ib_method_ops->registerLagBodyForceFunction(target3_force_data, 2); // Mesh 3 (rigid, stationary piece left)
        ib_method_ops->registerLagBodyForceFunction(target3_force_data, 3); // Mesh 4 (rigid, stationary piece right)
        // STANDARD PK1 DEVIATORIC STRESSES
        ib_method_ops->registerPK1StressFunction(PK1_dev_stress_data, 0); // Mesh 1
        ib_method_ops->registerPK1StressFunction(PK1_dev_stress_data, 1); // Mesh 2
        ib_method_ops->registerPK1StressFunction(PK1_dev_stress_data, 2); // Mesh 3
        ib_method_ops->registerPK1StressFunction(PK1_dev_stress_data, 3); // Mesh 4
        // STANDARD PK1 DILATION STRESSES
        ib_method_ops->registerPK1StressFunction(PK1_dil_stress_data, 0); // Mesh 1
        ib_method_ops->registerPK1StressFunction(PK1_dil_stress_data, 1); // Mesh 2
        ib_method_ops->registerPK1StressFunction(PK1_dil_stress_data, 2); // Mesh 3	
        ib_method_ops->registerPK1StressFunction(PK1_dil_stress_data, 3); // Mesh 4
		
       /*-------------------------------------------------------------------------------
        *
        * CONSTRUCT FEM SYSTEM OF EQUATIONS FOR EACH MESH
        *
        *-------------------------------------------------------------------------------*/
		ib_method_ops->initializeFEEquationSystems();
        EquationSystems* mesh1_equation_systems = ib_method_ops->getFEDataManager(0)->getEquationSystems();
        EquationSystems* mesh2_equation_systems = ib_method_ops->getFEDataManager(1)->getEquationSystems();
		EquationSystems* mesh3_equation_systems = ib_method_ops->getFEDataManager(2)->getEquationSystems();
        EquationSystems* mesh4_equation_systems = ib_method_ops->getFEDataManager(3)->getEquationSystems();

        PK1_dev_stress_data.quad_order = Utility::string_to_enum<libMeshEnums::Order>(input_db->getStringWithDefault("PK1_DEV_QUAD_ORDER","SEVENTH"));
        PK1_dil_stress_data.quad_order = Utility::string_to_enum<libMeshEnums::Order>(input_db->getStringWithDefault("PK1_DIL_QUAD_ORDER","SEVENTH"));
 

        // Create Eulerian initial condition specification objects.
        if (input_db->keyExists("VelocityInitialConditions"))
        {
            Pointer<CartGridFunction> u_init = new muParserCartGridFunction(
                "u_init", app_initializer->getComponentDatabase("VelocityInitialConditions"), grid_geometry);
            navier_stokes_integrator->registerVelocityInitialConditions(u_init);
        }

        if (input_db->keyExists("PressureInitialConditions"))
        {
            Pointer<CartGridFunction> p_init = new muParserCartGridFunction(
                "p_init", app_initializer->getComponentDatabase("PressureInitialConditions"), grid_geometry);
            navier_stokes_integrator->registerPressureInitialConditions(p_init);
        }

        // Create Eulerian boundary condition specification objects (when necessary).
        const IntVector<NDIM>& periodic_shift = grid_geometry->getPeriodicShift();
        vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM);
        if (periodic_shift.min() > 0)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                u_bc_coefs[d] = NULL;
            }
        }
        else
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                ostringstream bc_coefs_name_stream;
                bc_coefs_name_stream << "u_bc_coefs_" << d;
                const string bc_coefs_name = bc_coefs_name_stream.str();

                ostringstream bc_coefs_db_name_stream;
                bc_coefs_db_name_stream << "VelocityBcCoefs_" << d;
                const string bc_coefs_db_name = bc_coefs_db_name_stream.str();

                u_bc_coefs[d] = new muParserRobinBcCoefs(
                    bc_coefs_name, app_initializer->getComponentDatabase(bc_coefs_db_name), grid_geometry);
            }
            navier_stokes_integrator->registerPhysicalBoundaryConditions(u_bc_coefs);
        }

        // Create Eulerian body force function specification objects.
        if (input_db->keyExists("ForcingFunction"))
        {
            Pointer<CartGridFunction> f_fcn = new muParserCartGridFunction(
                "f_fcn", app_initializer->getComponentDatabase("ForcingFunction"), grid_geometry);
            time_integrator->registerBodyForceFunction(f_fcn);
        }

        // Set up visualization plot file writers.
        Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();
        if (uses_visit)
        {
            time_integrator->registerVisItDataWriter(visit_data_writer);
        }
        
	   /*-------------------------------------------------------------------------------
        *
        * SPECIFY CORRECT EXODUS NORMICLATURE FOR EACH MESH (NAB)
        *
        *-------------------------------------------------------------------------------*/
        // From circa 2016 UNC version --> autoibamr (Dec 2023) version!
        //AutoPtr<ExodusII_IO> mesh1_exodus_io(uses_exodus ? new ExodusII_IO(mesh1) : NULL);
        //AutoPtr<ExodusII_IO> mesh2_exodus_io(uses_exodus ? new ExodusII_IO(mesh2) : NULL);
        //std::unique_ptr<ExodusII_IO> mesh1_exodus_io(uses_exodus ? new ExodusII_IO(mesh1) : NULL);
        //std::unique_ptr<ExodusII_IO> mesh2_exodus_io(uses_exodus ? new ExodusII_IO(mesh2) : NULL);
        std::unique_ptr<ExodusII_IO> mesh1_exodus_io(uses_exodus ? new ExodusII_IO(*meshes[0]) : NULL);
        std::unique_ptr<ExodusII_IO> mesh2_exodus_io(uses_exodus ? new ExodusII_IO(*meshes[1]) : NULL);
        std::unique_ptr<ExodusII_IO> mesh3_exodus_io(uses_exodus ? new ExodusII_IO(*meshes[2]) : NULL);
        std::unique_ptr<ExodusII_IO> mesh4_exodus_io(uses_exodus ? new ExodusII_IO(*meshes[3]) : NULL);
    


        // Initialize hierarchy configuration and data on all patches.
        ib_method_ops->initializeFEData();
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        // Deallocate initialization objects.
        app_initializer.setNull();

        // Print the input database contents to the log file.
        plog << "Input database:\n";
        input_db->printClassData(plog);

        // Write out initial visualization data.
        int iteration_num = time_integrator->getIntegratorStep();
        double loop_time = time_integrator->getIntegratorTime();
        if (dump_viz_data)
        {
            pout << "\n\nWriting visualization files...\n\n";
            if (uses_visit)
            {
                time_integrator->setupPlotData();
                visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
            }
            if (uses_exodus) {
               /*-------------------------------------------------------------------------------
                *
                * FOR EXODUS OUTPUTS --> MAKE SURE TO SPECIFY EACH MESH (NAB)
                *
                *-------------------------------------------------------------------------------*/
                mesh1_exodus_io->write_timestep(mesh1_exodus_filename, *mesh1_equation_systems, iteration_num/viz_dump_interval+1, loop_time);
                mesh2_exodus_io->write_timestep(mesh2_exodus_filename, *mesh2_equation_systems, iteration_num/viz_dump_interval+1, loop_time);
                mesh3_exodus_io->write_timestep(mesh3_exodus_filename, *mesh3_equation_systems, iteration_num/viz_dump_interval+1, loop_time);
                mesh4_exodus_io->write_timestep(mesh4_exodus_filename, *mesh4_equation_systems, iteration_num/viz_dump_interval+1, loop_time);
            }
        }
        if (SAMRAI_MPI::getRank() == 0)
        {
            forcex_stream.open("F_x.curve", ios_base::out | ios_base::trunc);
            forcey_stream.open("F_y.curve", ios_base::out | ios_base::trunc);
            forcez_stream.open("F_z.curve", ios_base::out | ios_base::trunc);
	        velx_stream.open("U_x.curve", ios_base::out | ios_base::trunc);
            vely_stream.open("U_y.curve", ios_base::out | ios_base::trunc);
            velz_stream.open("U_z.curve", ios_base::out | ios_base::trunc);
	        workx_stream.open("W_x.curve", ios_base::out | ios_base::trunc);
            worky_stream.open("W_y.curve", ios_base::out | ios_base::trunc);
            workz_stream.open("W_z.curve", ios_base::out | ios_base::trunc);
	        U_L1_norm_stream.open("U_L1.curve", ios_base::out | ios_base::trunc);
	        U_L2_norm_stream.open("U_L2.curve", ios_base::out | ios_base::trunc);
	        U_max_norm_stream.open("U_max.curve", ios_base::out | ios_base::trunc);
        }
	
        // Main time step loop.
        static const double PI = acos(-1);
        double loop_time_end = time_integrator->getEndTime();
        double dt = 0.0;
        double max_J = 1;
        double min_J = 1;

        while (!MathUtilities<double>::equalEps(loop_time,loop_time_end) &&
               time_integrator->stepsRemaining())
        {
            iteration_num = time_integrator->getIntegratorStep();
            loop_time = time_integrator->getIntegratorTime();
	    
            pout <<                                                    "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "At beginning of timestep # " <<  iteration_num << "\n";
            pout << "Simulation time is " << loop_time              << "\n";

            dt = time_integrator->getMaximumTimeStepSize();
	    
	    time_integrator->advanceHierarchy(dt);
	    max_J=SAMRAI_MPI::maxReduction(J_maximum);
	    min_J=SAMRAI_MPI::minReduction(J_minimum);
	    pout << "Max Jacobian " <<  max_J << "\n";
	    pout << "Min Jacobian " <<  min_J << "\n";
            
            loop_time += dt;

            pout <<                                                    "\n";
            pout << "At end       of timestep # " <<  iteration_num << "\n";
            pout << "Simulation time is " << loop_time              << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout <<                                                    "\n";

            // At specified intervals, write visualization and restart files,
            // print out timer data, and store hierarchy data for post
            // processing.
            iteration_num += 1;
            const bool last_step = !time_integrator->stepsRemaining();
            if (dump_viz_data && (iteration_num%viz_dump_interval == 0 || last_step))
            {
                pout << "\nWriting visualization files...\n\n";
                if (uses_visit)
                {
                    time_integrator->setupPlotData();
                    visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
                }
                if (uses_exodus){
                   /*-------------------------------------------------------------------------------
                    *
                    * FOR EXODUS OUTPUTS --> MAKE SURE TO SPECIFY EACH MESH (NAB)
                    *
                    *-------------------------------------------------------------------------------*/
                    mesh1_exodus_io->write_timestep(mesh1_exodus_filename, *mesh1_equation_systems, iteration_num/viz_dump_interval+1, loop_time);
                    mesh2_exodus_io->write_timestep(mesh2_exodus_filename, *mesh2_equation_systems, iteration_num/viz_dump_interval+1, loop_time);
                    mesh3_exodus_io->write_timestep(mesh3_exodus_filename, *mesh3_equation_systems, iteration_num/viz_dump_interval+1, loop_time);
                    mesh4_exodus_io->write_timestep(mesh4_exodus_filename, *mesh4_equation_systems, iteration_num/viz_dump_interval+1, loop_time);
                }
            }
            if (dump_restart_data && (iteration_num%restart_dump_interval == 0 || last_step))
            {
                pout << "\nWriting restart files...\n\n";
                RestartManager::getManager()->writeRestartFile(restart_dump_dirname, iteration_num);
                ib_method_ops->writeFEDataToRestartFile(restart_dump_dirname, iteration_num);
            }
            if (dump_timer_data && (iteration_num%timer_dump_interval == 0 || last_step))
            {
                pout << "\nWriting timer data...\n\n";
                TimerManager::getManager()->print(plog);
            }
            if (dump_postproc_data && (iteration_num % postproc_data_dump_interval == 0 || last_step))
            {
                postprocess_data(patch_hierarchy,
                                 navier_stokes_integrator,
                                 mesh1,
                                 mesh1_equation_systems,
                                 iteration_num,
                                 loop_time,
                                 postproc_data_dump_dirname);
            }
        }

        // Close the logging streams.
        if (SAMRAI_MPI::getRank() == 0)
        {
            forcex_stream.close();
            forcey_stream.close();
            forcez_stream.close();
            velx_stream.close();
            vely_stream.close();
            velz_stream.close();
            workx_stream.close();
            worky_stream.close();
            workz_stream.close();
            U_L1_norm_stream.close();
            U_L2_norm_stream.close();
            U_max_norm_stream.close();
        }

        // Cleanup Eulerian boundary condition specification objects (when
        // necessary).
        for (unsigned int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];

    }// cleanup dynamically allocated objects prior to shutdown

    SAMRAIManager::shutdown();
    return 0;
}// main


void postprocess_data(Pointer<PatchHierarchy<NDIM> > /*patch_hierarchy*/,
                      Pointer<INSHierarchyIntegrator> /*navier_stokes_integrator*/,
                      Mesh& mesh,
                      EquationSystems* equation_systems,
                      const int /*iteration_num*/,
                      const double loop_time,
                      const string& /*data_dump_dirname*/)
{
    const unsigned int dim = mesh.mesh_dimension();
    {
        double F_integral[NDIM];

        for (unsigned int d = 0; d < NDIM; ++d) F_integral[d] = 0.0;

        System& F_system = equation_systems->get_system<System>(IBFEMethod::FORCE_SYSTEM_NAME);
//         System& U_system = equation_systems->get_system<System>(IBFEMethod::VELOCITY_SYSTEM_NAME);
//         System& X_system = equation_systems->get_system<System>(IBFEMethod::COORDS_SYSTEM_NAME);
//         System& dX_system = equation_systems->get_system<System>(IBFEMethod::COORD_MAPPING_SYSTEM_NAME);
        
        NumericVector<double>* F_vec = F_system.solution.get();
        NumericVector<double>* F_ghost_vec = F_system.current_local_solution.get();
        F_vec->localize(*F_ghost_vec);

        
        DofMap& F_dof_map = F_system.get_dof_map();
        std::vector<std::vector<unsigned int> > F_dof_indices(NDIM);
        unique_ptr<FEBase> fe(FEBase::build(dim, F_dof_map.variable_type(0)));
        unique_ptr<QBase> qrule = QBase::build(QGAUSS, dim, FIFTH);
        fe->attach_quadrature_rule(qrule.get());
        const std::vector<std::vector<double> >& phi = fe->get_phi();
        const std::vector<double>& JxW = fe->get_JxW();
        boost::multi_array<double, 2> F_node;
        
        


        
        const MeshBase::const_element_iterator el_begin = mesh.active_local_elements_begin();
        const MeshBase::const_element_iterator el_end = mesh.active_local_elements_end();
        for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
        {
            Elem* const elem = *el_it;
            fe->reinit(elem);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                F_dof_map.dof_indices(elem, F_dof_indices[d], d);
            }
            const int n_qp = qrule->n_points();
            const int n_basis = F_dof_indices[0].size();
            get_values_for_interpolation(F_node, *F_ghost_vec, F_dof_indices);
            for (int qp = 0; qp < n_qp; ++qp)
            {
                for (int k = 0; k < n_basis; ++k)
                {
                    for (int d = 0; d < NDIM; ++d)
                    {
                        F_integral[d] += F_node[k][d] * phi[k][qp] * JxW[qp];                        
                        
                    }
                    
                }
            }
        }

        SAMRAI_MPI::sumReduction(F_integral, NDIM);        
        
        if (SAMRAI_MPI::getRank() == 0)
        {
            forcex_stream << loop_time << " " << F_integral[0] << endl;
            forcey_stream << loop_time << " " << F_integral[1] << endl;
            forcez_stream << loop_time << " " << F_integral[2] << endl;
            
        }
    }
    
    {

    }
    return;
} // postprocess_data
