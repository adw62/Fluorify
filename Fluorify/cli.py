#!/usr/local/bin/env python

import os
import shutil
import logging

from .fluorify import Fluorify
from docopt import docopt

logger = logging.getLogger(__name__)

# =============================================================================================
# COMMAND-LINE INTERFACE
# =============================================================================================

usage = """
FLUORIFY
Usage:
  Fluorify [--output_folder=STRING] [--mol_name=STRING] [--ligand_name=STRING] [--complex_name=STRING] 
           [--solvent_name=STRING] [--yaml_path=STRING] [--c_atom_list=LIST] [--h_atom_list=LIST] [--num_frames=INT]
           [--net_charge=INT] [--gaff_ver=INT] [--auto_select=STRING] [--charge_only=BOOL] [--optimize=BOOL] [--num_gpu=INT] [--job_type=STRING]...
"""


def run_automatic_pipeline(yaml_file_path, complex_name, solvent_name):
    """Run YANK's automatic pipeline."""
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

    if args['--mol_name']:
        mol_name = args['--mol_name']
    else:
        mol_name = 'ligand'
        logger.debug(msg.format('mol file', mol_name + '.mol2'))

    if args['--ligand_name']:
        ligand_name = args['--ligand_name']
    else:
        ligand_name = 'MOL'
        logger.debug(msg.format('ligand residue name', ligand_name))

    if args['--complex_name']:
        complex_name = args['--complex_name']
    else:
        complex_name = 'complex'
        logger.debug(msg.format('complex name', complex_name))

    if args['--solvent_name']:
        solvent_name = args['--solvent_name']
    else:
        solvent_name = 'solvent'
        logger.debug(msg.format('solvent name', solvent_name))

    if args['--num_frames']:
        num_frames = int(args['--num_frames'])
    else:
        num_frames = 10000
        logger.debug(msg.format('number of frames', num_frames))

    if args['--net_charge']:
        net_charge = int(args['--net_charge'])
    else:
        net_charge = None
        logger.debug(msg.format('net charge', net_charge))

    if args['--gaff_ver']:
        gaff_ver = int(args['--gaff_ver'])
        if gaff_ver != 1 and gaff_ver != 2:
            raise ValueError('Can only use gaff ver. 1 or 2')
    else:
        gaff_ver = 2
        logger.debug(msg.format('gaff version', net_charge))

    if args['--charge_only']:
        charge_only = int(args['--charge_only'])
    else:
        charge_only = False
    if charge_only == True:
        logger.debug('Mutating ligand charges only...')
    else:
        logger.debug('Mutating all ligand parameters...')

    if args['--optimize']:
        opt = int(args['--optimize'])
    else:
        opt = False
    if opt == True:
        logger.debug('Optimizing ligand parameters...')
        c_atom_list = None
        h_atom_list = None
        auto_select = None
        job_type = 'optimize'
    else:
        logger.debug('Scanning ligand...')
        if args['--c_atom_list']:
            c_atom_list = []
            pairs = args['--c_atom_list']
            pairs = pairs.replace(" ", "")
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
            pairs = pairs.split('and')
            for pair in pairs:
                tmp = []
                pair = pair.split(',')
                for atom in pair:
                    tmp.append(atom)
                h_atom_list.append(tmp)
        else:
            h_atom_list = None

        if args['--auto_select']:
            auto_select = args['--auto_select']
            auto = ['1', '2', '3', 'ar']
            if auto_select not in auto:
                raise ValueError('Allowed automatic selections {}'.format(auto))
            if c_atom_list is not None or h_atom_list is not None:
                raise ValueError('Automatic target atom selection will conflict with populated atom lists')
        else:
            if c_atom_list is None and h_atom_list is None:
                raise ValueError('No target atoms specified')
            else:
                auto_select = None

        if args['--job_type']:
            job_type = args['--job_type'][0]
            allowed_jobs = ['F', 'Cl', 'N', 'NxF', 'NxCl', 'OH']
            if job_type not in allowed_jobs:
                raise ValueError('Allowed elements {}'.format(allowed_jobs))
        else:
            job_type = 'F'
            logger.debug(msg.format('job_type', job_type))

    if args['--output_folder']:
        output_folder = args['--output_folder']
    else:
        output_folder = './' + mol_name + '_' + job_type + '/'
        logger.debug(msg.format('output folder', output_folder))

    # Run the setup pipeline.
    if args['--yaml_path']:
        run_automatic_pipeline(args['--yaml_path'], complex_name, solvent_name)

    if args['--num_gpu']:
        num_gpu = int(args['--num_gpu'])
    else:
        num_gpu = 1

    Fluorify(output_folder, mol_name, ligand_name, net_charge, complex_name, solvent_name,
         job_type, auto_select, c_atom_list, h_atom_list, num_frames, charge_only, gaff_ver, opt, num_gpu)

