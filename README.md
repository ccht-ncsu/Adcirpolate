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

Make sure to set up the netCDF libraries for your specific machine. For example for Henry2 and STAMPEDE2: 

     ifeq ($(MACHINENAME),henry2)
        NETCDFHOME     :=/usr/local/apps/netcdf-centos7/4.6.1-intel2017/
        HDF5HOME       :=/usr/local/apps/hdf-centos7/hdf5/v1.10.2/
        FLIBS          := ${FLIBS} -I${NETCDFHOME}/include -L${NETCDFHOME}/lib -lnetcdff -lnetcdf -L${HDF5HOME}/lib -lhdf5_hl -lhdf5
     endif
     
     ifeq ($(MACHINENAME),stampede2)
        NETCDFHOME :=${TACC_NETCDF_DIR}
        FLIBS      := $(FLIBS) -L$(NETCDFHOME)/lib -lnetcdff -lnetcdf
     endif

Then we compile. An example for Henry2:

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

	setenv ESMF_CONFIG_FILE /location/of/esmf/installl_debug/lib/libO/Linux.intel.64.intelmpi.default/esmf.mk

**Step 2: Gather Run Input Files - Tidal Spinup**

For a tidal spinup file we just need tidal constituents for the time span of the combined simulations. Tidal spinup runs are typically 15 days in length. Files needed:

Mesh file - fort.13  
Nodal attributes file - fort.14  
Main control file - fort.15  

For Hurricane Matthew in particular (and possibly others) there is an HSOFS offset file for the vertical datum. 

HSOFS_Matthew_Offsetsurface  

Then the executable files and the associated submission script files. 

adcprep  
adcprep.csh/.sub or similar  
padcirc  
padcirc.csh/.sub  

Once all of these files have been collected, check them over and make sure they are set up specifically for your storm and time span and that you have the correct mesh. 

fort.15  

Some of the parameters to note:  
IHOT 0 - since no hotstart  
RNDAY 15.0 - typical for a tidal spinup  
NOUTE 3 - for netCDF files, 5 for netCDF4  
NSPOOL/V/E/GV/GE 3600 - output is written every hour. 
NHSINC 86400 - means hotstart written once a day, every 24 hours  

Cores

Also, mind the number of cores within both submissions scripts as well as the writer cores in padcirc submission script. The number of cores in the padcirc submission script should match the number of generated subdirectories or PE files from adcprep and number of writer cores specified in the padcirc submission script. For example, if adcprep uses 512 cores it will create 512 PE directories (PE0511 folders, starts at PE0000) and say you are using 10 writer cores, then the padcirc submission script should ask for 522 cores. 

Once you are satisfied with the files, submit adcprep first and then check the prepped files. Then run padcirc.

**Step 3: Initial Run with Storm - Coarse Resolution Simulation**

For the first part of the simulation with storm winds, we need all of the files that are used above (tidal spinup-not these exact files) plus a few additional files. The additional files contain the storm wind information. For hindcasts, we would typically use the OWI winds from ADCIRC as they are most accurate. During real time forecasting mode, we typically use the parameterized wind field from ADCIRC. Wind files come in the form of fort.22* files. Depending on the wind model being used, different fort.22 files will be needed. Please refer to adcirc.org for more information.  

Files needed for a standard ADCIRC simulation with wind forcings:

fort.13  
fort.14  
fort.15  
fort.22* (wind files) 
fort.67/68 (from tidal spinup run)  
Offsetsurface (if necessary)  
adcprep  
adcprep.csh  
padcirc  
padcirc.csh  
adchot.csh  
in.prep1  
in.prep2  
in.prephot  

Before submitting, check all the files and make sure they are set-up properly for your storm and machine. For this first run, assuming it is the initial ‘switch’ before moving to another mesh using Adcirpolate, there will be several things that need to be set-up in a specific manner. 

fort.15  

Some parameters to note:  
IHOT 67 - since there is hotstart file for this simulation, put either 67 or 68 (check files for specific number)   
RNDAY 15+X - where X is the length of the first part of the hotstart/storm run and 15 is the spinup run time   
WTIMINC 900 300 - for every 15 min/900 seconds reading of wind file, 300 ignore   
NOUTE 3 or 5 for netCDF files 1 for binary* need 1 for adcirpolate  
TOUTSE, TOUTFE - set for start and end time, here would be 15, starting after 15 of spin up then 15+X, where X is length of current run  
NSPOOL/V/E/GV/GE - 3600 output written every hour  
NHSINC - 43200 means hotstart written once every half a day, important if X is not a full day  

***Additional point - make sure that the fort.15 is writing binary/ascii output files. This is necessary for adcirpolate to work properly.***

fort.22  

Read if using OWI winds. In the fort.22 file the second line specifies at which point the simulation starts reading from the OWI wind and pressure files. Depending on the data in which your simulation starts. For example, with Hurricane Matthew, we start the storm simulation on October 2nd, 2016 at midnight. In the wind files, the data is provided beginning 24 hours before this so we need to shift when to start reading the data. In the OWI wind files, time steps are in 15 minute increments. So for 24 hours, that would be 96 - 15 minute increments. The negative sign indicates that we are going forward in time, as opposed to back in time. The first and third lines for this simulation should be 2 (for two sets of OWI wind files) and 1.0 (multiplier if needed) respectively. 

Once all of the files are set and ready. Submit adcprep and check that all the PE directories are correctly generated. After adcprep is run, we submit adchot to put the hotstart file into the PE directories. This is only for ascii type hotstart files, if we have a netCDF type file, we do not need to perform this step. When that has run successfully, we then submit padcirc.

**Step 4: Adcirpolate**

Set up the files in a specific way, there needs to be a directory named coarse and a directory named fine. Within the coarse directory, place all of the step three run files. It is easiest to just rename the directory they are already in and then create an empty fine directory. 

When the above run (part three) is done there will be a hotstart file in the PE0000 directory, check this directory for the latest written file, should be a fort.67 or fort.68 file. This needs to be a binary file for adcirpolate to work. Take this fort.6* file and move it out of the PE directory and into the main coarse directory. Then, make sure to delete whichever hot start file is already in the coarse directory and replace it with this one and name it fort.67. Once this is complete, submit adchot to get this new hotstart file intp the PE directories. 

Next we set up the fine directory. Within this directory we need the information for the next run which is typically a higher resolution simulation. Files will include those necessary for a typical ADCIRC run (submission scripts, executables, etc.) as well as the new mesh files. For this example the list includes: 

fort.13  
fort.14 (new higher resolution mesh)  
fort.15  
Offsetsurface  
fort.22* (same as above, but need to change the fort.22 second line to read winds starting with the fine simulation)  

Make sure all of these files are set up then submit adcprep while in the fine directory to create the PE directories. Once this is complete, we can perform adcirpolate. Go up one directory above the fine and coarse directories and pull in the adcirpolate executable as well as the adcirpolate submission script and submit adcirpolate. The number of cores specified in the adcirpolate submission script needs to match the number of cores in both the fine simulation. Meaning, if 512 PE directories were generated then the adcirpolate submission script needs to ask for 512 cores.

**Step 5: Second Run with Storm - High-Resolution Simulation**

When the above adcirpolate is complete, the fine simulation is ready to be performed. Go to the fine directory for this simulation and delete all of the PE directories. Then submit the submissions scripts in the following order:

adcprep  
adchot  
padcirc  

Submission scripts are very similar to the ones previously used for the coarse and switching steps. Double check to make sure the number of cores and the job name are correct. 
