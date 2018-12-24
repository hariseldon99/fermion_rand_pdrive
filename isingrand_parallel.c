#include "isingrand_parallel.h"
#include "params.h"

/* Display program usage, and exit.
 */
void
display_usage (int argc, char *argv[])
{
  printf
    ("\nUsage: %s [--option1 value1] [--option2 value2] [-flag1 value1] [-flag2 value2] ...\n",
     argv[0]);
  puts ("Integrates a periodically driven disordered Ising model.");
  printf
    ("Current hardcoded lattice size is %d. To change, edit file 'params.h' and recompile\n",
     N);
  printf
    ("Current hardcoded # of random instances per MPI process is %d. To change, edit file 'params.h' and recompile\n",
     NRAND);
#if defined(FIELD_DISORDER_OFF)
  puts ("\nField disorder is OFF!\n");
#else
  puts ("\nField disorder is ON!\n");
#endif
#if defined(INITCONDS_GHZ)
  puts
    ("Current hardcoded initial state is the GHZ state. To change, edit file 'params.h' and recompile\n");
#else
  puts
    ("Current hardcoded initial state is the ordered ground state (all spins up). To change, edit file 'params.h' and recompile\n");
#endif
#if defined(PERIODIC_BC)
  puts ("Periodic Boundary Conditions");
#else
  puts ("Open Boundary Conditions");
#endif
  printf
    ("Current hardcoded integrator method is given at runtime. To change, edit file 'integrator.h' and recompile\n\n");
  puts (options);
  printf ("\n");
  puts (envv);
  printf ("\nNOTE:\n");
  printf
    (" For environment variable options, please consult the GNU Scientific Library Reference Manual.\n");
  exit (EXIT_FAILURE);
}

//Function for loading static array of GArrays in memory
void
load_arrays (GArray * garrays[NRAND])
{
  gint j;
  // put the load arrays stuff here.
  for (j = 0; j < NRAND; j++)
    {
      garrays[j] = g_array_new (FALSE, FALSE, sizeof (gdouble));
    }
}

/* Progress bar function for verbose mode.
   Process has done x out of n rounds,
   and we want a bar of width w and resolution r */
static inline void
loadbar (int x, int n, int r, int w)
{
  // Only update r times.
  if (x % (n / r) != 0)
    return;

  // Calculuate the ratio of complete-to-incomplete.
  float ratio = x / (float) n;
  int c = ratio * w;

  // Show the percentage complete.
  printf ("%3d%% [", (int) (ratio * 100));

  // Show the load bar.
  for (x = 0; x < c; x++)
    printf ("=");

  for (x = c; x < w; x++)
    printf (" ");

  // ANSI Control codes to go back to the
  // previous line and clear it.
  printf ("]\n\033[F\033[J");
}

#undef __FUNCT__
#define __FUNCT__ "main"
int
main (int argc, char *argv[])
{
  int mpi_numprocs, mpi_rank;
  time_t begin, end;
  paramspace_pt param;
  double random_no;
  //Create output block
  double rand_arrh[NRAND][N], rand_arrj[NRAND][N];

  double data_mag[2];
  GArray *magrand_out[NRAND];

  //Initiate MPI environment. This is for one process per node
  MPI_Init (&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &mpi_numprocs);
  MPI_Comm_rank (MPI_COMM_WORLD, &mpi_rank);
  //Now, each of these processes runs the code below.
  //Create array of NRAND files. One for each instance of randomness
  FILE *outfiles_mag[NRAND];

  load_arrays (magrand_out);

  char appender[100];
  char tempstring1[100];

  /************************************  Input Block  ************************************/

  //Initiate defaults
  int opt = 0;
  int noentry = 1;
  int longIndex = 0;
  globalArgs.gamma_a = 0.0;
  globalArgs.omega = 0.0;
  globalArgs.alpha = 0.0;
  globalArgs.finaltime = 1.0;
  globalArgs.tsteps = 100;
  globalArgs.outFileName_Mag = "mags.dat";
  globalArgs.check_unitarity = 0;
  globalArgs.verbosity = 0;
  globalArgs.inputFiles = NULL;
  globalArgs.numInputFiles = 0;

  //Process the arguments with getopt_long(), then populate globalArgs. 
  opt = getopt_long (argc, argv, optString, longOpts, &longIndex);
  while (opt != -1)
    {
      switch (opt)
	{

	case 'a':
	  globalArgs.gamma_a = atof (optarg);
	  param.gamma_a = globalArgs.gamma_a;
	  noentry = 0;
	  break;

	case 'f':
	  globalArgs.omega = atof (optarg);
	  param.omega = globalArgs.omega;
	  noentry = 0;
	  break;

	case 'p':
	  globalArgs.alpha = atof (optarg);
	  param.alpha = globalArgs.alpha;
	  noentry = 0;
	  break;

	case 't':
	  globalArgs.finaltime = atof (optarg);
	  noentry = 0;
	  break;

	case 'n':
	  globalArgs.tsteps = atoi (optarg);
	  noentry = 0;
	  break;

	case 'm':
	  globalArgs.outFileName_Mag = optarg;
	  noentry = 0;
	  break;

	case 'u':
	  globalArgs.check_unitarity++;
	  noentry = 0;

	case 'v':
	  globalArgs.verbosity++;
	  noentry = 0;
	  break;

	default:
	  display_usage (argc, argv);
	  exit (EXIT_FAILURE);
	}

      opt = getopt_long (argc, argv, optString, longOpts, &longIndex);

    }
  if ((noentry != 0) && (mpi_rank == 0))
    {
      display_usage (argc, argv);
      exit (EXIT_FAILURE);
    }

  globalArgs.inputFiles = argv + optind;
  globalArgs.numInputFiles = argc - optind;


  /**********************************  End Input Block  **********************************/

  (void) time (&begin);
  int i, j;

  double *y = malloc (DIM * sizeof (double));
  double *y_init = malloc (DIM * sizeof (double));	//This is the full state

  int status;
  double t_dyn, t_inc = globalArgs.finaltime / (globalArgs.tsteps - 1);

  int x = 0, res = BAR_RES, w = BAR_RES;	//loadbar point, resolution and width

/*******************************  Random Numbers Block  *******************************/
  gsl_rng *r;
  const gsl_rng_type *rtype;
  gsl_rng_env_setup ();

  rtype = gsl_rng_default;
  r = gsl_rng_alloc (rtype);

  long seed =
    (long int) cantor_pair (gsl_rng_default_seed, mpi_rank) + 100000;

  gsl_rng_set (r, seed);

  if ((globalArgs.verbosity != 0) && (mpi_rank == 0))
    {
      printf ("\n");
      puts ("Random Number Generator Info:");
      printf ("=============================================\n");
      printf ("Generator type \t=\t%s\n", gsl_rng_name (r));
      printf ("Seed value \t=\t%lu\n", gsl_rng_default_seed);
      printf ("First value \t=\t%lf\n", gsl_rng_uniform (r));
      printf ("=============================================\n");
    }

  int randcount;
  //Initialize NRAND * N random numbers in the range (-1,1) to use as disorder in lattice
  for (i = 0; i < NRAND; i++)
    for (j = 0; j < N; j++)
      {
	random_no = gsl_rng_uniform (r);	//This generates a rand number in the interval [0,1) with uniform probability
	random_no = 2.0 * random_no - 1.0;	//Scale it to the interval (-1,1)
	rand_arrh[i][j] = random_no;
#if defined(FIELD_DISORDER_OFF)
	rand_arrh[i][j] = 0.0;
#endif
	random_no = gsl_rng_uniform (r);	//This generates a random number in the interval [0,1) with uniform probability
	random_no = 2.0 * random_no - 1.0;	//Scale it to the interval (-1,1)
	rand_arrj[i][j] = random_no;
      }

  /******************************  End Random Numbers Block  *****************************/

  //This loop initializes N instances of random arrays and uses them in N integrations
  {
    for (randcount = 0; randcount < NRAND; randcount++)
      {
	for (i = 0; i < N; i++)
	  {
	    //Copy from random array of h and j for each instance of randcount
	    {
	      param.hrand[i] = param.alpha * rand_arrh[randcount][i];	//Scale random number by the perturbation
	      param.jrand[i] = param.alpha * rand_arrj[randcount][i];
	    }
	  }
   /******************************  Initial Conditions Block  *****************************/
#if defined(INITCONDS_GHZ)
	//This sets the initial conditions to a random but normalized state
	set_initconds_ghz (y_init);
#else
	//This sets the initial conditions to all ising spins up ie the ordered ground state
	set_initconds_allup (y_init);
#endif
  /******************************  End Initial Conditions Block  *****************************/

  /******************************  Set Initial Conditions Block  *****************************/
	for (i = 0; i < DIM; i++)
	  {
	    y[i] = y_init[i];
	  }

  /******************************  End Set Initial Conditions Block  *****************************/

  /********************************    Integration  Block    ********************************/
	if ((globalArgs.verbosity != 0) && (mpi_rank == 0))
	  {
	    printf ("\n");
	    printf ("Now beginning integration with magnetization = %lf\n",
		    magnetization (y));
	    printf ("=============================================\n");
	  }
	x = 0;			//x counts the number of executions of integrator  
	status = GSL_SUCCESS;
	for (t_dyn = 0.0; t_dyn < globalArgs.finaltime; t_dyn = t_dyn + t_inc)
	  {
	    {
	      if (status != GSL_SUCCESS)
		printf ("\n Error. Integrator failed at t=%lf", t_dyn);
	      else
		{
		  data_mag[0] = t_dyn;
		  data_mag[1] = magnetization (y);

		  if (globalArgs.check_unitarity != 0)
		    printf
		      ("\n Nonunitary parameter at t = %f is = %f",
		       t_dyn, check_unitarity (y));

		  {
		    g_array_append_vals (magrand_out[randcount], &data_mag,
					 2);
		  }

		}
	      {
		if ((globalArgs.verbosity != 0) && (mpi_rank == 0))
		  loadbar (x, globalArgs.tsteps, res, w);
	      }
	    }
	    status = integrate (y, t_dyn, t_dyn + t_inc, &param);
	    x++;		//x counts the number of executions of integrator  
	  }

	if ((globalArgs.verbosity != 0) && (mpi_rank == 0))
	  {
	    printf ("\n");
	    printf ("=============================================\n");
	    printf ("End integration with magnetization = %lf",
		    magnetization (y));
	  }
  /********************************  End Integration  Block  ********************************/
      }
  }
  (void) time (&end);

  /************************************  Output Block  ************************************/
  {

    //Output magrand_out, entrand_out, negrand_out. magrand_out[i] is the time, mag data of 
    //the ith random number iteration etc
    //First, open NRAND files with NRAND filenames
    for (i = 0; i < NRAND; i++)
      {
	sprintf (appender, ".%d_seed_%ld", i, seed);
	strcpy (tempstring1, globalArgs.outFileName_Mag);
	outfiles_mag[i] =
	  fopen (strcat (globalArgs.outFileName_Mag, appender), "w");
	strcpy (globalArgs.outFileName_Mag, tempstring1);
      }
    //Then, write output to each filenames
    for (i = 0; i < NRAND; i++)
      {
	for (j = 0; j < magrand_out[i]->len; j = j + 2)
	  {
	    fprintf (outfiles_mag[i], "%lf %lf\n",
		     g_array_index (magrand_out[i], gdouble, j),
		     g_array_index (magrand_out[i], gdouble, j + 1));
	  }
      }
    printf ("\n");
    //Print entered values
    printf ("\n");
    puts ("Runtime Parameters:");
    printf ("=============================================\n");
    printf ("Drive amplitude \t=\t%lf\n", param.gamma_a);
    printf ("Drive frequency \t=\t%lf\n", param.omega);
    printf ("Perturbation amplitude\t=\t%lf\n", param.alpha);
    printf ("Final time\t\t=\t%lf\n", globalArgs.finaltime);
    printf ("Time steps\t\t=\t%d\n", globalArgs.tsteps);
    printf ("# of rands\t\t=\t%d\n", NRAND * mpi_numprocs);
    printf ("Lattice size\t\t=\t%d\n", N);
    printf ("=============================================\n");
    puts ("Final Output:");
    printf ("=============================================\n");
    printf ("Total Runtime (root) \t=\t%lf\n", difftime (end, begin));
    printf ("Number of processes \t=\t%d\n", mpi_numprocs);
    printf ("=============================================\n");
  }
  /**********************************  End Output Block  **********************************/

  //fclose (globalArgs.outFile);
  gsl_rng_free (r);

  for (i = 0; i < NRAND; i++)
    {
      fclose (outfiles_mag[i]);
      g_array_free (magrand_out[i], TRUE);
    }
  free (y);
  free (y_init);

  MPI_Finalize ();
  return EXIT_SUCCESS;
}
