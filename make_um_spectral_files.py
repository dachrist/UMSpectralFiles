from calc_column_density import get_column_density
import convert_atmo_to_socrates_input_spectrum as convert
import os
import tempfile
from pathlib import Path
import shutil

__all__ = ["make_um_spectral_files"]

def make_um_spectral_files(verbose=False,workdir=None,specs=['H2O','CO2',"CO","CH4","NH3","H2","He","Na","K","Li","Rb","Cs","HCN"],stellar_spectrum=None,stellar_spectrum_format=None,bands=32,
                        estf_dir = "/data/jgoyal/kabs/",grav=None,file_prefix=None,file_suffix="",atmo_profile=None,skel_file=None,run=True):
    """
    Creates spectral files in the current working directory.

    verbose - If True, it will provide more information about what it is doing.
    workdir - If specified, it will do alll of its temporary work in this directory.  If left as None, it will create a temporary work directory
    specs   - An array of strings specifying the species involved
    stellar_spectrum - The file containing the stellar spectrum.  It can either be in 'ATMO' format or 'UM' format.
    stellar_spectrum_format - See the previous comment.
    bands   - The number of bands
    estf_dir - The directory with all of the k-data
    grav    - The planet gravity, in m/s^2
    file_prefix - The prefix for the spectral file names (required).
    file_suffix - The suffix for the spectral file names (option).
    atmo_profile - The NetCDF containing the ATMO atmospheric profile.
    skel_file - The skeleton file.  Currently this file doesn't create it.
    run - If true, it runs Socrates automatically.  Useful if you just want to see the scripts.

    """

    # Check that we have access to Socrates
    if (shutil.which("tidy_90") is None):
        raise Exception("Can't find tidy_90.  Socrates not in path?")
    if (shutil.which("prep_spec") is None):
        raise Exception("Can't find prep_spec.  Socrates not in path?")
    
    # Generic files and directories

    # Temperature resolution for the thermal source function.  Doesn't need to be changed.
    Tmin = 70
    Tmax = 4000
    Tsteps = 1000
    
    # Verify inputs

    if (workdir is None):
        # No specified work directory, 
        tempdir = tempfile.TemporaryDirectory()
        workdir = tempdir.name + "/"
    else:                        
        # Check that work directory exists
        if (not os.path.isdir(workdir)):
            if (verbose):
                print("Specified work directory doesn't exist.  Creating it.")
            try:
                os.mkdir(workdir)
            except:
                raise Exception("Can't create the specifield work directory.")
    # Gravity
    if (grav is None):
        raise Exception("The surface gravity needs to be specified in m/s^2.")
    if (grav < 0.):
        raise Exception("The gravity parameter is negative.")
    if (grav > 50.):
        raise Warning("The specified gravity is rather high.  Are you sure it is in m/s^2")

    # Output file prefix/suffix
    if (file_prefix is None):
        raise Exception("file_prefix needs to be set.")

    # ATMO Profile
    if (atmo_profile is None):
        raise Exception("atmo_profile needs to be specified.")
    if (not os.path.isfile(atmo_profile)):
        raise Exception("Specified atmo_profile doesn't exist.")

    # Skeleton file
    if (skel_file is None):
        raise Exception("A skeleton file needs to be specified.")
    if (not os.path.isfile(skel_file)):
        raise Exception("Specified skeleton file doesn't exist.")

    # Spectrum
    if (not os.path.isfile(stellar_spectrum)):
        raise Exception("Specified stellar_spectrum doesn't exist.")
    if (stellar_spectrum_format == "UM"):
        # Nothing needs to be done beyond copying over the spectrum
        if (verbose):
            print("Copying over stellar spectrum.")
        try:
            os.popen("cp {} {}/stellar_spectrum.txt".format(stellar_spectrum,workdir))
        except:
            raise Exception("Error copying stellar spectrum to work directory.")
    elif (stellar_spectrum_format == "ATMO"):
        # Since this is only a temporary conversion that will be deleted, the star name and the source aren't that important.
        convert.write_um_spec(fname=stellar_spectrum, fname_out='{}/stellar_spectrum.txt'.format(workdir), star_name='Star', source='Some source')
    else:
        raise Exception("Unsupported stellar spectrum format. UM and ATMO are the only options.")
    

    # Copy over the skeleton of the spectral file.
    if (verbose):
        print("Copying over the skeleton spectral file to the work directory.")
    try:
        os.popen("cp {} {}/skel_file".format(skel_file,workdir))
    except:
        raise Exception("Error copying skeleton file to work directory.")
    
    # Column densities:
    if (verbose):
        print("Computing column densities from the ATMO profile.")
    # Compute the column densities from the ATMO profile
    cols_in = get_column_density(fname=atmo_profile,gravity=grav,species=specs)

    # Convert column densities to uppercase
    cols  = {k.upper():v for k,v in cols_in.items()}

    # Check that column densities exist for each species
    for x in specs:
        if x.upper() not in cols:
            raise Exception("Error: can't find key for {} in column density dictionary.".format(x.upper()))
        

    # Generate the SOCRATES scripts for each file.
    
    def generate_scripts(bandtype):

        # output file name
        if (file_suffix != ""):
            fname_out = "{}_{}_{}_{}".format(file_prefix,bandtype,bands,file_suffix)
        else:
            fname_out = "{}_{}_{}".format(file_prefix,bandtype,bands)

    
        # STEP 1: Generate the skeleton

        # Check that the script is on the same page about absorption indices
        
        # STEP 2: Add k-terms and thermal
        try:
            f = open(workdir + "/step2_{}".format(bandtype),"w")
        except:
            raise Exception("Error opening script file in work directory.")
        
        f.write("#!/bin/bash\n")
        f.write("# 1: Add k-terms for each gas\n")
        f.write("# 2: Add thermal source function: tabulated between 70 and 4000 K with 1000 points\n")
        f.write("touch spec_file_step2_{}\n".format(bandtype))
        f.write("rm spec_file_step2_{}\n".format(bandtype))
        f.write("prep_spec << EOF\n")
        f.write("skel_file\n")
        f.write("n\n")
        f.write("spec_file_step2_{}\n".format(bandtype))
    
        # Move ESTF data to temporary dir, replacing absorber index, and write script
    
        abs_idx = 1

        for spec in specs:
            if (spec == "He" or spec == "H2"):
                prefix = "h2-{}".format(spec.lower())
            else:
                prefix = spec.lower()
            estf_file = estf_dir + '{}_{}_t5e-3_uw1116'.format(prefix,bands)
            estf_temp = "estf_{}".format(spec.lower())
            fi = open(estf_file,"r")
            fo = open((workdir + "/"+estf_temp).strip(),"w")

            lines = fi.readlines()

            # Replace the indices so that they are appropriate
            for line in lines:
                if (": for absorber index" in line):
                    s = line.replace("in band", " ")
                    s = s.replace(": for absorber index"," ")
                    band,idx = s.split()
                    a,b = line.split(":")
                    c = b.replace("{}".format(idx),"{}".format(abs_idx))
                    out =a + ":" + c
                    fo.write(out)
                else:
                    fo.write(line)

            fi.close()
            fo.close()

            f.write("5\n")
            if (abs_idx > 1):
                f.write("y\n")
            f.write("{}\n".format(estf_temp))
            abs_idx = abs_idx + 1

        # Write the thermal stuff
        f.write("6\n")
        f.write("n\n")
        f.write("T\n")
        f.write("{} {}\n".format(Tmin,Tmax))
        f.write("{}\n".format(Tsteps))

        # If shortwave, add in the spectrum
        if (bandtype == "sw"):
            f.write("2\n")
            f.write("N\n")
            f.write("stellar_spectrum.txt\n")
            f.write("Y\n")
            f.write("3\n")
            f.write("H\n")
            f.write("/home/bd257/Software/socrates/trunk/data/gases/refract_H2\n")

        # Write the file and exit
        f.write("-1\n")
        f.write("EOF\n")
        f.close()

        # STEP 3: Tidy
        try:
            f = open(workdir+"/step3_{}".format(bandtype),"w")
        except:
            raise Exception("Error writing step 3 script to work directory.")
        
        # Write column densities
        f.write("#!/bin/bash\n")
        for x in specs:
            f.write("{}=\"{}\"\n".format(x.upper(),cols[x.upper()]))
        f.write("\n")
        f.write("tidy_90 << EOF\n")
        f.write("spec_file_step2_{}\n".format(bandtype))
        f.write("n\n")
        f.write("{}\n".format(fname_out))
        f.write("6\n")
        for x in specs:
            f.write("${}\n".format(x.upper()))
        f.write("1\n")
        for x in specs:
            f.write("${}\n".format(x.upper()))
        f.write("0.99990\n")
        f.write("-1\n")
        f.write("EOF\n")

        os.system("ls " + workdir)
        f.close()

        
    # Create the scripts
    if (verbose):
        print("Generating scripts.")
    generate_scripts("lw")
    generate_scripts("sw")

    # Run the scripts
    if (run):

        if (verbose):
            print("Running scripts.")
        old_cwd = os.getcwd()
        os.chdir(workdir)
        os.system(". ./step2_lw")
        os.system(". ./step3_lw")
        os.system(". ./step2_sw")
        os.system(". ./step3_sw")

        # Move the results out of the work directory
        if (verbose):
            print("Moving final spectral files from working directory to CWD.")
        os.system("mv {}* {}".format(file_prefix,old_cwd))
        os.chdir(old_cwd)

    return
