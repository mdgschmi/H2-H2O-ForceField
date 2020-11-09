/* C routines for H2H2OFF.py */


#include "MMTK/universe.h"
#include "MMTK/forcefield.h"
#include "MMTK/forcefield_private.h"
#include <math.h>
#include <stdio.h>

/* This function does the actual energy (and gradient) calculation.
   Everything else is just bookkeeping. */
static void
h2h2o_evaluator(PyFFEnergyTermObject *self,
		PyFFEvaluatorObject *eval,
		energy_spec *input,
		energy_data *energy)
/* The four parameters are pointers to structures that are
   defined in MMTK/forcefield.h.
   PyFFEnergyTermObject: All data relevant to this particular
   energy term.
   PyFFEvaluatorObject:  Data referring to the global energy
   evaluation process, e.g. parallelization
   options. Not used here.
   energy_spec:          Input parameters for this routine, i.e.
   atom positions and parallelization parameters.
   energy_data:          Storage for the results (energy terms,
   gradients, second derivatives).
*/
{
  int Nwater=20;
  double e=0.0;
  double totalgradx=0.0;
  double totalgrady=0.0;
  double totalgradz=0.0;
  int radindx,theindx,chiindx,arayindx,sign;
  int chiindxstore;
  int ihighr,ihight,ihighc;
  double vadj,dvdradj,dvdtadj,dvdcadj;

  double rmin=3.00;
  double deltar=0.046;
  vector3 *coordinates = (vector3 *)input->coordinates->data;
  vector3 *g;
  int i,j,n;
  int atom_index = (int)self->param[0];  /* atom index */
  //int atom_index2 = (int)self->param[1];  /* atom index */
  //int atom_index3 = (int)self->param[2];  /* atom index */
  
  //double x1 = coordinates[atom_index][0]*10.;
  //double y1 = coordinates[atom_index][1]*10.;
  //double z1 = coordinates[atom_index][2]*10.;
  
  PyArrayObject *pot_array = (PyArrayObject *)self->data[0];
  double *pot = (double *)pot_array->data;  /* atomic charges */

  PyArrayObject *dvdr_array = (PyArrayObject *)self->data[1];
  double *dvdr = (double *)dvdr_array->data;  /* atomic charges */

  PyArrayObject *dvdt_array = (PyArrayObject *)self->data[2];
  double *dvdt = (double *)dvdt_array->data;  /* atomic charges */

  PyArrayObject *dvdc_array = (PyArrayObject *)self->data[3];
  double *dvdc = (double *)dvdc_array->data;  /* atomic charges */

  PyArrayObject *com_array = (PyArrayObject *)self->data[4];
  double *r_com = (double *)com_array->data;

  PyArrayObject *rotmatrixarray = (PyArrayObject *)self->data[5];
  double *rotmat=(double *)rotmatrixarray->data;
  
  /* energy_terms is an array because each routine could compute
     several terms that should logically be kept apart. For example,
     a single routine calculates Lennard-Jones and electrostatic interactions
     in a single iteration over the nonbonded list. The separation of
     terms is only done for the benefit of user code (universe.energyTerms())
     returns the value of each term separately), the total energy is
     always the sum of all terms. Here we have only one energy term. */
  
  double rotpH2[3];

  double Rvec[3];
  double r;
  double chi;
  double theta;
  
  double drdx[3];
  double dthetadx[3];
  double dchidx[3];
  
  double rotgr[3];
  double gr[3];

  for (n=0;n<Nwater;n++)
    {
      
      for (i=0;i<3;i++)
	{
	  Rvec[i]=coordinates[atom_index][i]-r_com[i+3*n];
	}
      
      // Construct Rotation Matrix
      double RotMatrix[3][3]=
	{
	  {rotmat[0+9*n],rotmat[1+9*n],rotmat[2+9*n]},
	  {rotmat[3+9*n],rotmat[4+9*n],rotmat[5+9*n]},
	  {rotmat[6+9*n],rotmat[7+9*n],rotmat[8+9*n]}
	};
      
      // Rotation Matrix is Orthogonal by Construction
      // So Inverse is the Transpose
      double RotInverse[3][3]=
	{
	  {rotmat[0+9*n], rotmat[3+9*n], rotmat[6+9*n]},
	  {rotmat[1+9*n], rotmat[4+9*n], rotmat[7+9*n]},
	  {rotmat[2+9*n], rotmat[5+9*n], rotmat[8+9*n]}
	};
      
      // Rotate Water and para-hydrogen
      // from Space Fixed Frame to Internal Water Frame
      for (i=0;i<3;i++)
	{
	  rotpH2[i]=0.0;
	  for (j=0;j<3;j++)
	    {
	      rotpH2[i]+=RotMatrix[i][j]*Rvec[j];
	    }
	}
      
      // Calculate r
      r=sqrt(pow(rotpH2[0],2.)+pow(rotpH2[1],2.)+pow(rotpH2[2],2));
      
      // Calculate drdx for each component
      drdx[0]=rotpH2[0]/r;
      drdx[1]=rotpH2[1]/r;
      drdx[2]=rotpH2[2]/r;
      
      // Calculate theta  
      theta=acos(rotpH2[2]/r)*180./(M_PI);
      
      if (theta < 0) {
	theta=theta+360.;
      }
      
      
      // Calculate dthetadx for each component
      dthetadx[0]=rotpH2[0]*rotpH2[2]/(sqrt(pow(rotpH2[0],2.)+pow(rotpH2[1],2.)));
      dthetadx[0]=dthetadx[0]/(r*r);
      dthetadx[1]=rotpH2[1]*rotpH2[2]/(sqrt(pow(rotpH2[0],2.)+pow(rotpH2[1],2.)));
      dthetadx[1]=dthetadx[1]/(r*r);
      dthetadx[2]=-1.0*(sqrt(pow(rotpH2[0],2.)+pow(rotpH2[1],2.)))/(r*r);
      
      
      dthetadx[0]=dthetadx[0]*180./(M_PI);
      dthetadx[1]=dthetadx[1]*180./(M_PI);
      dthetadx[2]=dthetadx[2]*180./(M_PI);
      
      // Calculate chi  
      chi=atan(rotpH2[1]/rotpH2[0])*180./(M_PI);
      
      if (chi < 0){
	chi=chi+360.;
      }
      
      // Calculate dchidx for each component
      dchidx[0]=-1.0*rotpH2[1]/(pow(rotpH2[0],2.)+pow(rotpH2[1],2.));
      dchidx[1]=rotpH2[0]/(pow(rotpH2[0],2.)+pow(rotpH2[1],2.));
      dchidx[2]=0.0;
      
      dchidx[0]=dchidx[0]*180./(M_PI);
      dchidx[1]=dchidx[1]*180./(M_PI);
      
      //  printf("R : %e   Theta : %e   Chi : %e \n", r, theta, chi);
      
      radindx=(int)round((r/.0529177249-rmin)/deltar);
      theindx=(int)round(theta);
      chiindx=(int)round(chi);
      chiindxstore=chiindx;
      
      if (radindx<0)
	radindx=0;
      if (radindx>500)
	radindx=500;
      if (theindx<0)
	theindx=0;
      if (theindx>180)
	theindx=180;
      if (chiindx<0)
	chiindx=0;
      if(chiindx>360)
	chiindx=360;
      
      sign=1;
      
      if(chiindx>90)
	{
	  if(chiindx<180)
	    {
	      chiindx=180-chiindx;
	      sign=-1;
	    }
	  else
	    if(chiindx>270)
	      {
		chiindx=360-chiindx;
		sign=-1;
	      }
	    else
	      chiindx=chiindx-180;
	}
      
      arayindx=chiindx+91*theindx+91*181*radindx;
      ihighr=arayindx+181*91;
      ihight=arayindx+91;
      ihighc=arayindx+1;
      
      vadj=pot[arayindx];
      dvdradj=dvdr[arayindx];
      dvdtadj=dvdt[arayindx];
      dvdcadj=dvdc[arayindx];
      
      if (ihighr < 8251971)
	{
	  vadj=vadj+((pot[ihighr]-pot[arayindx])/deltar)*
	    (r/.0529177249-(radindx*deltar+rmin));
	  vadj=vadj+((pot[ihight]-pot[arayindx])/(1))*(theta-theindx);
	  vadj=vadj+((pot[ihighc]-pot[arayindx])/(1))*(chi-chiindxstore);
	  
	  dvdradj=dvdradj+((dvdr[ihighr]-dvdr[arayindx])/(deltar))*
	    (r/.0529177249-(radindx*deltar+rmin));
	  dvdradj=dvdradj+((dvdr[ihight]-dvdr[arayindx])/(1))*(theta-theindx);
	  dvdradj=dvdradj+((dvdr[ihighc]-dvdr[arayindx])/(1))*(chi-chiindxstore);
	  
	  dvdtadj=dvdtadj+((dvdt[ihighr]-dvdt[arayindx])/(deltar))*
	    (r/.0529177249-(radindx*deltar+rmin));
	  dvdtadj=dvdtadj+((dvdt[ihight]-dvdt[arayindx])/(1))*(theta-theindx);
	  dvdtadj=dvdtadj+((dvdt[ihighc]-dvdt[arayindx])/(1))*(chi-chiindxstore);
	  
	  dvdcadj=dvdcadj+((dvdc[ihighr]-dvdc[arayindx])/(deltar))*
	    (r/.0529177249-(radindx*deltar+rmin));
	  dvdcadj=dvdcadj+((dvdc[ihight]-dvdc[arayindx])/(1))*(theta-theindx);
	  dvdcadj=dvdcadj+((dvdc[ihighc]-dvdc[arayindx])/(1))*(chi-chiindxstore);
	}
      
      e=e+vadj;
      
      // now calculate the gradients!
      for (i=0; i < 3; i ++)
	{ 
	  rotgr[i]=(dvdradj*drdx[i] + dvdtadj*dthetadx[i]+dvdcadj*dchidx[i]*sign);
	}
      
      
      // Rotate the gradients to match the original coordinates
      for (i=0;i<3;i++)
	{
	  gr[i]=0.0;
	  for (j=0;j<3;j++)
	    {
	      gr[i]+=RotInverse[i][j]*rotgr[j];
	    }
	}
      totalgradx=totalgradx+gr[0];
      totalgrady=totalgrady+gr[1];
      totalgradz=totalgradz+gr[2];
      
      
    }
  
  energy->energy_terms[self->index] = e;//(double)beads;                                                                                                                                                                                           
  
  // If only the energy is asked for, stop here.                                                                                                                                                                                                   
  if (energy->gradients == NULL)
    return;
  
  
  // Add the gradient contribution to the global gradient array.
  //   It would be a serious error to use '=' instead of '+=' here,
  //   in that case all previously calculated forces would be erased.
  g = (vector3 *)((PyArrayObject*)energy->gradients)->data;
  
  g[atom_index][0]+=totalgradx;
  g[atom_index][1]+=totalgrady;
  g[atom_index][2]+=totalgradz;

}

/* A utility function that allocates memory for a copy of a string */
static char *
allocstring(char *string)
{
  char *memory = (char *)malloc(strlen(string)+1);
  if (memory != NULL)
    strcpy(memory, string);
  return memory;
}

/* The next function is meant to be called from Python. It creates the
   energy term object at the C level and stores all the parameters in
   there in a form that is convient to access for the C routine above.
   This is the routine that is imported into and called by the Python
   module, H2H2OFF.py. */
static PyObject *
H2H2OTerm(PyObject *dummy, PyObject *args)
{
  PyFFEnergyTermObject *self;
  PyArrayObject *pot;
  PyArrayObject *dvdr;
  PyArrayObject *dvdt;
  PyArrayObject *dvdc;
  PyArrayObject *r_com;
  PyArrayObject *rotmat;
  int atom_index;

  /* Create a new energy term object and return if the creation fails. */
  self = PyFFEnergyTerm_New();
  if (self == NULL)
    return NULL;

  /* Convert the parameters to C data types.*/
  if (!PyArg_ParseTuple(args, "O!O!O!O!O!O!O!i",
			&PyUniverseSpec_Type, &self->universe_spec,
			&PyArray_Type, &pot,
			&PyArray_Type, &dvdr,
			&PyArray_Type, &dvdt,
			&PyArray_Type, &dvdc,
			&PyArray_Type, &r_com,
			&PyArray_Type, &rotmat,
			&atom_index))
    return NULL;
  /* We keep a reference to the universe_spec in the newly created
     energy term object, so we have to increase the reference count. */
  Py_INCREF(self->universe_spec);
  /* A pointer to the evaluation routine. */
  self->eval_func = h2h2o_evaluator;
  /* The name of the energy term object. */
  self->evaluator_name = "h2h2o";
  /* The names of the individual energy terms - just one here. */
  self->term_names[0] = allocstring("h2h2o");
  if (self->term_names[0] == NULL)
    return PyErr_NoMemory();
  self->nterms = 1;
  /* self->param is a storage area for parameters. Note that there
     are only 40 slots (double) there, if you need more space, you can use
     self->data, an array for up to 40 Python object pointers. */
  self->param[0] = (double) atom_index;

  /* self->data is the other storage area for parameters. There are
     40 Python object slots there */
  self->data[0] = (PyObject *)pot;
  self->data[1] = (PyObject *)dvdr;
  self->data[2] = (PyObject *)dvdt;
  self->data[3] = (PyObject *)dvdc;
  self->data[4] = (PyObject *)r_com;
  self->data[5] = (PyObject *)rotmat;

  Py_INCREF(pot);
  Py_INCREF(dvdr);
  Py_INCREF(dvdt);
  Py_INCREF(dvdc);
  Py_INCREF(r_com);
  Py_INCREF(rotmat);

  /* Return the energy term object. */
  return (PyObject *)self;
}

/* This is a list of all Python-callable functions defined in this
   module. Each list entry consists of the name of the function object
   in the module, the C routine that implements it, and a "1" signalling
   new-style parameter passing conventions (only veterans care about the
   alternatives). The list is terminated by a NULL entry. */
static PyMethodDef functions[] = {
  {"H2H2OTerm", H2H2OTerm, 1},
  {NULL, NULL}		/* sentinel */
};


/* The initialization function for the module. This is the only function
   that must be publicly visible, everything else should be declared
   static to prevent name clashes with other modules. The name of this
   function must be "init" followed by the module name. */
DL_EXPORT(void)
initMMTK_h2_h2o(void)
{
  PyObject *m;

  /* Create the module and add the functions. */
  m = Py_InitModule("MMTK_h2_h2o", functions);

  /* Import the array module. */
#ifdef import_array
  import_array();
#endif

  /* Import MMTK modules. */
  import_MMTK_universe();
  import_MMTK_forcefield();

  /* Check for errors. */
  if (PyErr_Occurred())
    Py_FatalError("can't initialize module MMTK_h2_h2o");
}
