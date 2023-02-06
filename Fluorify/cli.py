#!/usr/local/bin/env python

import os
from simtk import unit
import shutil

from .fluorify import Fluorify, SysBuilder
from docopt import docopt

# =============================================================================================
# COMMAND-LINE INTERFACE
# =============================================================================================

usage = """
Fluorify - Calculate change in binding free energy for fluorinated analogues of molecules.

Usage:
  Fluorify [--output_folder=STRING] [--mol_name=STRING] [--ligand_name=STRING] [--complex_name=STRING] [--solvent_name=STRING]
            [--yaml_path=STRING] [--setup_path=STRING] [--c_atom_list=STRING] [--h_atom_list=STRING] [--num_frames=INT] [--net_charge=INT]
            [--gaff_ver=INT] [--equi=INT] [--num_fep=INT] [--auto_select=STRING] [--param=STRING]
            [--num_gpu=INT] [--exclude_dualtopo=BOOL] [--job_type=STRING]...
    
Options:
    --job_type=STRING
        Scanning and full FEP jobs that can be performed,
        default: 'F'
        options:
                F  - Fluorine
                Cl - Chlorine
                N  - Nitrogen
                
    --mol_name=STRING
        Name of input mol2 file containing ligand,
        default: 'ligand'
    
    --output_folder=STRING
        Directory for output,
        default: './mol_name_job_type'
     
    --ligand_name=STRING
        String ligand is addressed by in input mol2 file,
        default: 'MOL'
      
    --complex_name=STRING
        Name of input pdb file containing ligand,
        default: 'complex'
     
    --solvent_name=STRING 
        Name of input pdb file containing ligand,   
        default: 'solvent'
        
    --net_charge=INT
        Net charge of ligand to be passed to antechamber for parameterization, net_charge should also be set in setup.yaml,
        default: 0 
        
    --param=STRING 
        String for what parameters to modify in nonbonded force of alchemical system
        default: all
        options: 'charge', 'sigma', 'vdw', 'all'
    
    --gaff_ver=INT
        Gaff version to use in paramterisation,
        default: 2           
        options: 1, 2
        
    --yaml_path=STRING
        Path to yaml file containing options for yank experiment builder,
        default: './setup.yaml'
     
    --setup_path=STRING
        Dummy option to use Fluorify system builder of YANK
        default: 'None'
        
    --c_atom_list=STRING
        List of indices of carbon atoms in ligand.mol2 to be replaced with N,
        default: None
     
    --h_atom_list=STRING
        List of indices of hydrogen atoms in ligand.mol2 to be replaced with F or Cl,
        default: None
        
    --auto_select=STRING
        Automatic selection of indices for mutation based on input,
         default: None
         options:
                 1  - C.1 carbons or their associated H
                 2  - C.2 carbons or their associated H
                 3  - C.2 carbons or their associated H
                 ar - C.ar carbons or their associated H
     
    --num_frames=INT
         Number of frames of trajectory to collect for objective, frames spaced by 5ps,
         default: 10,000
    
    --equi=INT
        Number of steps of equilibration, each step is 2fs,
        default: 250000
     
    --num_fep=INT
        Number of the best mutants to test with full FEP,
        default: 1   
        Note: Changes in OpenMM 7.7 may make full FEP calculations slow. https://github.com/openmm/openmm/issues/252
        
    --num_gpu=INT 
        Number of GPU for the node where the calculation is run,
        default: 1
        Note: This software is not configured to use MPI and should only be run on one node, however this node may have multiple GPUs
     
    --exclude_dualtopo=BOOL
        Excludes any atoms in dual topology from seeing each other.
        default: 1
        Note:  Fluorines are added to a dual topology as typically the hydrogens they are mutated from are constrained. It is not possible to alchemically interpolate the hydrogen constraint into a C-F harmonic bond.
    
"""

def run_automatic_pipeline(yaml_file_path, complex_name, solvent_name):
    """Run YANK's automatic pipeline."""
    print('Using YANK system builder')
    from yank.experiment import ExperimentBuilder
    exp_builder = ExperimentBuilder(yaml_file_path)

    # Modify the output directory of the setup to be consistent
    # with the hardcoded paths in Fluorify and FSim. The searched
    # path is 'input/complex_name/complex_name.pdb' so we also
    # need to modify the name of the system.
    exp_builder.output_dir = '.'
    exp_builder.setup_dir = 'input'

    # Run the automatic pipeline.
    exp_builder.setup_experiments()
    assert len(exp_builder._db.systems) == 1, 'Setting up multiple systems is not currently supported'
    system_name = next(iter(list(exp_builder._db.systems.keys())))

    # Copy YANK setup files to match the Fluorify folder structure.
    for phase_name, user_phase_name in zip(['complex', 'solvent'], [complex_name, solvent_name]):
        # Create Fluorify directory structure.
        fluorify_phase_dir = os.path.join('input', user_phase_name)
        os.makedirs(fluorify_phase_dir, exist_ok=True)
        for extension in ['.prmtop', '.pdb']:
            yank_file_path = os.path.join(exp_builder.setup_dir, 'systems', system_name, phase_name + extension)
            fluorify_file_path = os.path.join(fluorify_phase_dir, user_phase_name + extension)
            shutil.copyfile(yank_file_path, fluorify_file_path)

def main(argv=None):
    args = docopt(usage, argv=argv, options_first=True)

    msg = 'No {0} specified using default {1}'

    #Name assigned to the complex leg of the simulation
    if args['--complex_name']:
        complex_name = args['--complex_name']
    else:
        complex_name = 'complex'
        print(msg.format('complex name', complex_name))

    #Name applied to the solvent leg of the simulation
    if args['--solvent_name']:
        solvent_name = args['--solvent_name']
    else:
        solvent_name = 'solvent'
        print(msg.format('solvent name', solvent_name))

    #Net charge of the ligand in the simulation
    if args['--net_charge']:
        net_charge = int(args['--net_charge'])
    else:
        net_charge = None
        print(msg.format('net charge', net_charge))

    #What ver. of the Gaff force feild to use
    if args['--gaff_ver']:
        gaff_ver = int(args['--gaff_ver'])
        if gaff_ver != 1 and gaff_ver != 2:
            raise ValueError('Can only use gaff ver. 1 or 2')
    else:
        gaff_ver = 2
        print(msg.format('gaff version', gaff_ver))

    # Dir path to a Yank setup file
    if args['--yaml_path']:
        # Use yank system builder
        run_automatic_pipeline(args['--yaml_path'], complex_name, solvent_name)
        # All these variables passed are dummies we are using yank to prep system.
        systems = SysBuilder('./input/', './receptor.pdb', './ligand.mol2', 'amber14/protein.ff14SB.xml',
                             'amber14/spce.xml', './gaff.xml', 1.0 * unit.nanometers, 0.15 * unit.molar,
                             using_yank=True)

    # Dummy arg to ask for Fluorify system builder with fixed args.
    elif args['--setup_path']:
        #READ OPTIONS
        if net_charge is None:
            net_charge = 0
        systems = SysBuilder('./input/', './receptor.pdb', './ligand.mol2', 'amber14/protein.ff14SB.xml',
                             'amber14/spce.xml', './gaff.xml', 1.0 * unit.nanometers, 0.15 * unit.molar, ligand_charge=net_charge)
    else:
        raise ValueError('No set up script provided. Set setup_path or yaml_path')

    #Name for mol2 file with ligand
    if args['--mol_name']:
        mol_name = args['--mol_name']
    else:
        mol_name = 'ligand'
        print(msg.format('mol file', mol_name + '.mol2'))

    #resname of ligand
    if args['--ligand_name']:
        ligand_name = args['--ligand_name']
    else:
        ligand_name = 'MOL'
        print(msg.format('ligand residue name', ligand_name))

    #Number of snapshots to collect
    if args['--num_frames']:
        num_frames = int(args['--num_frames'])
    else:
        num_frames = 10000
        print(msg.format('number of frames', num_frames))

    #Number of steps to equilibriate
    if args['--equi']:
        equi = int(args['--equi'])
    else:
        equi = 250000
        print(msg.format('Number of equilibration steps', equi))

    #What params to alchmically transform
    if args['--param']:
        param = str(args['--param'])
        accepted_param = ['charge', 'sigma', 'vdw', 'all']
        if param not in accepted_param:
            raise ValueError('param selected not in accepted params: {}'.format(accepted_param))
    else:
        param = ['all']

    for x in param:
        print('Mutating {} ligand parameters...'.format(x))

    if args['--exclude_dualtopo']:
        exclude_dualtopo = int(args['--exclude_dualtopo'])
    else:
        exclude_dualtopo = True
        print('Excluding dual topology from seeing itself')

    # Charge optimisation no longer supported see https://github.com/adw62/Ligand_Charge_Optimiser
    opt = False

    print('Scanning ligand...')
    if args['--c_atom_list']:
        c_atom_list = []
        pairs = args['--c_atom_list']
        pairs = pairs.replace(" ", "")
        c_name = pairs.replace(",", "")
        pairs = pairs.split('and')
        for pair in pairs:
            tmp = []
            pair = pair.split(',')
            for atom in pair:
                tmp.append(atom)
            c_atom_list.append(tmp)
    else:
        c_atom_list = None

    if args['--h_atom_list']:
        h_atom_list = []
        pairs = args['--h_atom_list']
        pairs = pairs.replace(" ", "")
        h_name = pairs.replace(",", "")
        pairs = pairs.split('and')
        for pair in pairs:
            tmp = []
            pair = pair.split(',')
            for atom in pair:
                tmp.append(atom)
            h_atom_list.append(tmp)
    else:
        h_atom_list = None

    #vestiage of sulphur mutations not compatable with SSP
    o_atom_list = None

    if args['--auto_select']:
        auto_select = args['--auto_select']
        auto = ['1', '2', '3', 'ar']
        if auto_select not in auto:
            raise ValueError('Allowed automatic selections {}'.format(auto))
        if c_atom_list is not None or h_atom_list is not None or o_atom_list is not None:
            raise ValueError('Automatic target atom selection will conflict with populated atom lists')
    else:
        if c_atom_list is None and h_atom_list is None and o_atom_list is None:
            raise ValueError('No target atoms specified')
        else:
            auto_select = None

    if args['--job_type']:
        job_type = args['--job_type'][0]
        allowed_jobs = ['F', 'Cl', 'N', 'NxF', 'NxCl', 'VDW']
        if job_type not in allowed_jobs:
            raise ValueError('Allowed elements {}'.format(allowed_jobs))
    else:
        job_type = 'F'
        print(msg.format('job_type', job_type))

    if args['--output_folder']:
        output_folder = args['--output_folder']
    else:
        id = ''
        if auto_select is not None:
            id += '_auto_' + auto_select
        if h_atom_list is not None:
            id += '_H' + h_name
        if c_atom_list is not None:
            id += '_C' + c_name
        if o_atom_list is not None:
            id += '_O' + o_name
        output_folder = './' + mol_name + '_' + job_type + id + '/'

    if args['--num_gpu']:
        num_gpu = int(args['--num_gpu'])
    else:
        num_gpu = 1
        print(msg.format('number of GPUs per node', num_gpu))

    if args['--num_fep']:
        num_fep = args['--num_fep']
    else:
        num_fep = 1
        print(msg.format('number of FEP calculations', num_fep))

    Fluorify(output_folder, mol_name, ligand_name, net_charge, complex_name, solvent_name, job_type, auto_select,
             c_atom_list, h_atom_list, num_frames, param, gaff_ver, num_gpu, num_fep, equi, exclude_dualtopo,
             opt, o_atom_list, systems)

