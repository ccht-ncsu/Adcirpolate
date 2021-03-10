# What is adcirpolate ?
*adcirpolate* is a module for interpolating hotstart files between different ADCIRC meshes.
Hence, we can switch between different meshes to avoid unnecessary computation on meshes with
high resolution in areas where no important phenomenon is going on. If you do not know about
ADCIRC, or its hotstart files, you are in the wrong place !

- [Prerequisites](##Prerequisites)
- [CMake and Fortran version requirements](##CMake-and-Fortran-version-requirements)
- [Running the code](##Running-the-code)

## Prerequisites
In order to run this code, you need MPI and ESMF libraries. The program is to be built with the same MPI library that is used to build ESMF. You also need to set the following environment variables, and run cmake:

    export FC=mpif90
    export ESMF_CONFIG_FILE=/path/to/esmf.mk
    cmake /path/to/CMakeLists.txt
    make all

## CMake and Fortran version requirements
The current version of Adcirpolate has been tested on multiple HPC clusters with different versions of CMake and GNU and Intel Fortran compilers. Understandably, some tunings to the CMakeLists.txt might be required for clusters with less up-to-date software. In case you are having trouble compiling the code on your machine, please submit an issue with output log from CMake or Fortran compiler.

## Running the code
Adcirpolate should be executed in a directory with two folders named `coarse` and `fine`. Typically, we transfer the results from hotstart file of the coarse mesh to a hotstart file of the fine mesh. The two folders contain the partitioned ADCIRC meshes. In the `coarse` directory, the hotstart file of the coarse mesh exists. After running, the program creates a hotstart file for the fine mesh.

## Compiling and Running ADCIRC with ADCIRPOLATE

**Step 1: Compile Necessary Source Code - ADCIRC**

Load file ADCIRC source code into workspace and unzip. Go to the work folder, this is where we will perform the compiling of the ADCIRC source code. Type:

	make clean 
	make clobber

Then, go into the cmplrflags.mk file and make sure the correct machine is set and there are debug lines (if needed line 116):

	DP            :=  -DREAL8 -DLINUX -DCSCA -DCMPI -DDEBUG_WARN_ELEV

I also had to switch the NetCDF location of files as follows (lines 179-181). This may be specific by machine: 

     ifeq ($(MACHINENAME),henry2)
        NETCDFHOME     :=/usr/local/apps/netcdf-centos7/4.6.1-intel2017/
        HDF5HOME       :=/usr/local/apps/hdf-centos7/hdf5/v1.10.2/
        FLIBS          := ${FLIBS} -I${NETCDFHOME}/include -L${NETCDFHOME}/lib -lnetcdff -lnetcdf -L${HDF5HOME}/lib -lhdf5_hl -lhdf5
     endif

Then we compile:

	cmake/3.16.3
	intel/2017.1.132
	intel_mpi/2017
	PrgEnv-intel/2017.1.132
	hdf5/1.10.2-intel2017
	netcdf/4.6.1-intel2017

	make adcprep MACHINENAME=henry2 NETCDF=enable NETCDF4=enable NETCDF4_COMPRESSION=enable
	make padcirc MACHINENAME=henry2 NETCDF=enable NETCDF4=enable NETCDF4_COMPRESSION=enable

You can also made it so all of the above modules are automatically loaded when opening the connection by adding it to a ./tcshrc file or similar. Now the ADCIRC source code will have the adcprep and padcirc executables. 

The next step is to load esmf. This takes quite some time. Gather the esmf files, create an esmf folder and then load the submission script and the zip file in this folder. Then submit the submission script and wait. 

Then, when this is complete, the next step is to install adcirpolate. This must be done AFTER esmf is installed because we need the location of a specific file for the compilation. This is similar to the esmf compiling, start by creating a folder adcirpolate and move the adcirpolate.zip and the submission script into the folder. THEN you need to go into the submission script and change the following line to match where the associated esmf file is located (all one line):

	setenv ESMF_CONFIG_FILE /usr/local/usrapps/jcdietri/esmf/installl_debug/lib/libO/Linux.intel.64.intelmpi.default/esmf.mk

**Step 2: Gather Run Input Files - Tidal Spinup**

For a tidal spinup file we just need tidal constituents for the time span of the combined simulations. Tidal spinup runs are typically 15 days in length. Files needed:

Mesh file - fort.13  
Nodal attributes file - fort.14  
Main control file - fort.15  

For Hurricane Matthew in particular (and possibly others) there is an HSOFS offset file for the vertical datum. 

HSOFS_Matthew_Offsetsurface  

Then the executable files and the associated submission script files. 

adcprep  
adcprep.csh  
padcirc  
padcirc.csh  

Once all of these files have been collected, check them over and make sure they are set up specifically for your storm and time span and that you have the correct mesh. Som components of note 

fort.15  

Some of the parameters to note:
IHOT 0 - since no hotstart  
RNDAY 15.0 for tidal spinup  
DT 1.0 for HSOFS, 0.5 for FEMA-HR  
NOUTE 3 for netCDF files  
NSPOOL/V/E/GV/GE - 3600 means every hour water elevations written  
NHSINC - 86400 means hotstart written once a day, every 24 hours  

Cores

Also, mind the number of cores within both submissions scripts as well as the writer cores in padcirc.csh. The number of cores in padcirc should match the number generated PE files from adcprep and number of writer cores specified in the padcirc.csh. For example, if adcprep.csh says 512 cores and there are PE0511 folders (starts at PE0000) and padcirc.csh has 10 writer cores, then padcirc.csh should say 522 cores. 

Once you are satisfied with the files submit the adcprep.csh first and then check the prepped files. Then run the padcirc.csh

**Step 3: Initial Run with Storm - Coarse Resolution Simulation**

For the first part of the simulation with storm winds, we need all of the files that are used above (tidal spinup-not these exact files) plus a few additional files. The additional files contain the storm wind information. For hindcasts, we would typically use the OWI winds from ADCIRC as they are most accurate. During a real time forecasting mode, we would use the parameterized wind field from ADCIRC. These come in the form of a fort.22* file. The specific number specifies the spatial scale of the winds, with the fort.22 being the winds from ADCIRC. Additionally, we have some other files for the submission scripts and an additional submission script due to this being a hot start from the tidal spin up performed above. Not all submission will need these files, it depends on how you have your submission scripts (adcprep.csh and padcirc.csh) written. Another important thing to note, as of 2/25/21 adcirpolate needs to have binary/localized hotstart files. 

Files needed:

fort.13  
fort.14  
fort.15  
fort.22  
fort.221  
fort.222   
fort.223  
fort.224  
fort.67/68 (from tidal spinup run)  
HSOFS_Matthew_Offsetsurface (if using the HSOFS mesh)  
adcprep  
adcprep.csh  
padcirc  
padcirc.csh  
adchot.csh  
in.prep1  
in.prep2  
in.prephot  

Same for above, check all the files and make sure they are set-up properly. For this first run, assuming it is the initial ‘switch’ before moving to another mesh using adcirpolate, there will be several things that need to be set-up depending on the parameters of the simulation. 

fort.15  

Some of the parameters to note:
IHOT 67 - since there is hotstart file for this simulation, put either 67 or 68, check files
RNDAY 15+X where X is the length of the first part of the hotstart/storm run
DT 1.0 for HSOFS & OW, 0.5 for FEMA-HR
NWS -12 for including OWI winds for Matthew
WTIMINC 900 300 , for every 15 min/900 seconds reading of wind file, 300 ignore. 
NOUTE 3 or 5 for netCDF files 1 for binary* need 1 for adcirpolate
TOUTSE, TOUTFE set for start and end time, here would be 15, starting after 15 of spin up then 15+X, where X is length of current run
NSPOOL/V/E/GV/GE - 3600 means every hour water elevations written
NHSINC - 43200 means hotstart written once every half a day, important if X is not a full day

***Additional point - make sure that the fort.15 is writing binary/ascii output files. This is necessary for adcirpolate to work properly.***

fort.22

In the fort.22 file the second line specifies at which point the simulation starts reading from the wind files. Depending on the data in which your simulation starts, for example, with Hurricane Matthew, we start the storm simulation on October 2nd, 2016 at midnight. In the wind files, the dates start 24 hours before this so we need to shift when the reading begins. In the wind files, time steps are in 15 minute increments. So for 24 hours, that would be 96 - 15 minute increments. The negative sign indicates that we are going forward in time, as opposed to back in time. The first and third lines for this simulation should be 2 (for two sets of OWI wind files) and 1.0 (multiplier if needed) respectively. 

For submitting, there are three submission scripts (these are for Henry2):

adcprep.csh
adchot.csh
padcirc.csh

They involve three other files:

in.prep1
in.prep2
in.prephot

in.prep1

The in.prep1 file has three lines. The first line is the number of cores, the second line is number of wind file sets and the third line is fort.14 for the nodal attributes file. 

in.prep2

The in.prep1 file has two lines. The first line is the number of cores, and the second line is again for the number of sets of wind files.

in.prephot

The in.prephot file has three lines. The first line is the number of cores, the second line is set of wind files and the third line is 67/68 for the hotstart file.

Examples of the submission scripts are below. Main things that would need changing from simulation to simulation would be the number of cores (-n), wall clock time (-W) and the job name (-J) to help you remember the run details.

adcprep.csh

	#! /bin/csh
	#BSUB -J adcprep
	#BSUB -o adcprep.%J
	#BSUB -e adcprep.%J
	#BSUB -W 15
	#BSUB -n 1
	#BSUB -q queuename

	./adcprep <in.prep1
	./adcprep <in.prep2

adchot.csh

	#! /bin/csh
	#BSUB -J adcprep
	#BSUB -o adcprep.%J
	#BSUB -e adcprep.%J
	#BSUB -W 15
	#BSUB -n 1
	#BSUB -q queue name

	./adcprep <in.prephot

padcirc.csh

	#! /bin/csh
	#BSUB -J padcirc
	#BSUB -o padcirc.%J
	#BSUB -e padcirc.%J
	#BSUB -W 2880
	#BSUB -n 320
	#BSUB -q queuename

	mpirun ./padcirc -W 10

Once all of the files are set and ready. Submit adcprep.csh and check that all the PE directories are correctly generated. After adcprep is run, we submit adchot.csh to put the hotstart file into the PE directories. This is only for ascii type hotstart files, if we have a netCDF type file, we do not need to perform this step. When that has run successfully, we then submit padcirc.csh

**Step 4: Adcirpolate**

Set up the files in a specific way, there needs to be a directory named coarse and a directory named fine. Within the coarse directory, place all of the step three run files. It is easiest to just rename the directory they are already in and then create an empty fine directory. 

When the above run (part three) is done there will be a hotstart file in the PE0000 directory, check this directory for the latest written file, should be a 67 or 68 file. This needs to be a binary file for adcirpolate to work. Take this fort.6* file and move it out of the PE directory and into the coarse directory. Then, make sure to delete whichever hot start file is already in the coarse directory and replace it with this one and name it fort.67. Once this is complete, run the adchot.csh in this directory to get this hotstart file in the PE directories. 

Next we set up the fine directory. Within this directory we need the information for the next run which is typically a higher resolution simulation. Files will include those necessary for above parts of the ADCIRC run (submission scripts, executables, etc) as well as the new mesh files. For this example the list includes: 

fort.13
fort.14
fort.15
Highres_mattthew_Offsetsurface
fort.22* (same as above, but need to change the fort.22 second line and add the time spent in the coarse simulation)

Make sure all of these files are set up, look in the fort.15 and make sure the time has been set to include the entire time of the simulation on this mesh, i.e. 15+X+Y (15 is the spin up time, X is time from first (coarse) simulation with Matthew winds, and then Y is the time for this part of the simulation). When all is set, submit adcprep.csh while in the fine directory to create the PE directories. Once this is complete, we can perform adcirpolate. Go up one directory above the fine and coarse directories and pull in the adcirpolate executable as well as the adcirpolate submission script adcirpolate.csh (might be named hot_reader.csh)and then submit adcirpolate.csh (this will need to be run in parallel so there are likely 500+ cores needed) 

Cores

The number of cores should match what is in the adcprep.csh submission script, not the padcirc.csh submission script. So it does not include the writer cores just the number of PE directories. 

hot_reader.csh

	#! /bin/csh
	#BSUB -J adcirpolate
	#BSUB -o adcirpolate.%J
	#BSUB -e adcirpolate.%J
	#BSUB -W 2880
	#BSUB -n 522
	#BSUB -q queue

	mpirun ./adcirpolate

**Step 5: Second Run with Storm - High-Resolution Simulation**

When the above adcirpolate is complete, the fine simulation is ready to be performed. Go to the fine directory for this simulation and delete all of the PE directories. Then submit the submissions scripts in the following order:

adcprep.csh
adchot.csh
padcirc.csh

Submission scripts are very similar to the ones previously used for the coarse and switching steps. Double check to make sure the number of cores and the job name are correct. 
