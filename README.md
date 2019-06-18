# Fluorify
Perform Fluorine scanning analysis on trajectories.

# Installation

- Download from git

  git clone https://github.com/adw62/Fluorify

- Install via pip

  pip install path/to/directory
  
- Install yank and multiprocess to get all dependencies
 
  conda install -c omnia yank
  
  conda install -c conda-forge multiprocess
  
# Options

[--job_type=STRING] Scanning and full FEP jobs that can be performed,

    default: 'F'
    options:
            F  - Fluorine
            Cl - Chlorine
            N  - Nitrogen
        
[--output_folder=STRING] Directory for output,

    default: './mol_name_job_type'
            
[--ligand_name=STRING] String ligand is adressed by in input mol2 file,

    default: 'MOL'
    
[--mol_name=STRING] Name of input mol2 file containg ligand,

    default: 'ligand'
    
[--complex_name=STRING] Name of input pdb file containg ligand,

    default: 'complex'
                     
[--solvent_name=STRING] Name of input pdb file containg ligand,

    default: 'solvent'
                      
[--yaml_path=STRING] Path to yaml file containg options for yank exsperiment builder,

    default: './setup.yaml'
      
[--c_atom_list=LIST] List of indices of carbon atoms in ligand.mol2 to be replaced with N,

    default: None
    
[--h_atom_list=LIST] List of indices of hydrogen atoms in ligand.mol2 to be replaced with F or Cl

    default: None
    
[--auto_select=STRING] Automatic selction of indicies for mutation based on input,

     default: None
     options:
             1  - C.1 carbons or their associated H
             2  - C.2 carbons or their associated H
             3  - C.2 carbons or their associated H
             ar - C.ar carbons or their associated H
                       
[--num_frames=INT] Number of frames of trajectory to collect for objective, frames spaced by 5ps,

    default: 10,000

[--net_charge=INT] Net charge of ligand to be passed to antechamber for paramterisation, net_charge should also be set in setup.yaml,

    default: 0                   
            
[--gaff_ver=INT] Gaff version to use in paramterisation,

    default: 2           
    options: 1, 2
    
[--equi=INT] Number of steps of equilibriation, each step is 2fs,

    default: 250000

[--num_fep=INT] Number of the best mutants to test with full FEP,

    default: 1             

[--charge_only=BOOL] Boolean to determine if only charge parameters should be changed,

    default: False                    

[--vdw_only=BOOL] Boolean to determine if only van der Waals parameters should be changed,

    default: False
    
[--num_gpu=INT] Number of GPU for the node where the calculation is run,

    note: This software is not configured to use MPI and should only be run on one node, however this node may have multiple GPUs
    default: 1

# Example usage

Performs a fluorine scan of a ligand replacing hydrogens 13, 14 and 15. Uses 25ps of sampling equilibriated for 10fs to calculate the single step pertubation ddGs. Concludes with full FEP calculations for top three fluorinated mutants.

    Fluorify --job_type='F' --num_frames='5' --equi='5' --h_atom_list='13, 14, 15' --yaml_path='./setup.yaml' --num_fep='3'

