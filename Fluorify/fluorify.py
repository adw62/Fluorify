#!/usr/bin/env python

from .energy import FSim
from .mol2 import Mol2, MutatedLigand
from .mutants import Mutants

import os
import time
import itertools
import mdtraj as md
import numpy as np
import shutil
from simtk import unit

import openmoltools as moltools
import parmed as pmd
import simtk.openmm as mm


#CONSTANTS
e = unit.elementary_charges

class Fluorify(object):
    def __init__(self, output_folder, mol_name, ligand_name, net_charge, complex_name, solvent_name, job_type,
                 auto_select, c_atom_list, h_atom_list, num_frames, param, gaff_ver, num_gpu,
                 num_fep, equi, exclude_dualtopo, opt, o_atom_list, systems):
        """

        :param output_folder: STRING, path to output folder.
        :param mol_name: STRING, name of ligand mol2 file.
        :param ligand_name: STRING, name of ligand resname.
        :param net_charge: INT, value for ligand netcharge.
        :param complex_name: STRING, name for complex system.
        :param solvent_name: STRING, name for solvent system.
        :param job_type: STRING,
        :param auto_select:
        :param c_atom_list:
        :param h_atom_list:
        :param num_frames:
        :param param:
        :param gaff_ver:
        :param num_gpu:
        :param num_fep:
        :param equi:
        :param exclude_dualtopo:
        :param opt:
        :param o_atom_list:
        :param systems:
        """

        self.output_folder = output_folder
        self.net_charge = net_charge
        self.job_type = job_type
        self.num_frames = num_frames
        self.gaff_ver = gaff_ver
        self.num_fep = int(num_fep)

        # Prepare directories/files and read in ligand from mol2 file
        mol_file = mol_name + '.mol2'
        input_folder = './input/'

        complex_sim_dir = input_folder + complex_name + '/'
        solvent_sim_dir = input_folder + solvent_name + '/'

        if os.path.isdir(self.output_folder):
            print('Output folder {} already exists. '
                  'Will attempt to skip ligand parametrisation, proceed with caution...'.format(self.output_folder))
        else:
            try:
                os.makedirs(self.output_folder)
            except:
                print('Could not create output folder {}'.format(self.output_folder))
        shutil.copy2(input_folder+mol_file, self.output_folder)
        self.mol = Mol2()
        try:
            Mol2.get_data(self.mol, input_folder, mol_file)
        except:
            raise ValueError('Could not load molecule {}'.format(self.input_folder + mol_file))

        # Check ligand atom order is consistent across input topologies.
        input_files = [input_folder + mol_file, complex_sim_dir + complex_name + '.pdb',
                       solvent_sim_dir + solvent_name + '.pdb']
        self.mol2_ligand_atoms, complex_ligand_atoms, solvent_ligand_atoms = get_atom_list(input_files, ligand_name)
        if self.mol2_ligand_atoms != complex_ligand_atoms:
            raise ValueError('Names and or name casing and or atom order of ligand not matched across input files.'
                             'Parameters will not be applied where expected')
        if complex_ligand_atoms != solvent_ligand_atoms:
            raise ValueError('Names and or name casing and or atom order of ligand not matched across input files.'
                             'Parameters will not be applied where expected')

        input_files = input_files[1:3]
        self.complex_offset, self.solvent_offset = get_ligand_offset(input_files, self.mol2_ligand_atoms, ligand_name)
        print('Parametrize wild type ligand...')
        wt_ligand = MutatedLigand(file_path=self.output_folder, mol_name=mol_name,
                                  net_charge=self.net_charge, gaff=self.gaff_ver)

        print('Loading complex and solvent systems...')

        #COMPLEX
        self.complex_sys = []
        self.complex_sys.append(FSim(ligand_name=ligand_name, sim_name=complex_name, input_folder=input_folder,
                                     param=param, num_gpu=num_gpu, offset=self.complex_offset, opt=opt,
                                     exclude_dualtopo=exclude_dualtopo, system=systems.complex))
        self.complex_sys.append([complex_sim_dir + complex_name + '.dcd'])
        self.complex_sys.append(complex_sim_dir + complex_name + '.pdb')

        if not os.path.isfile(self.complex_sys[1][0]):
            self.complex_sys[1] = [complex_sim_dir + complex_name + '_gpu' + str(x) + '.dcd' for x in range(num_gpu)]
            for name in self.complex_sys[1]:
                if not os.path.isfile(name):
                    self.complex_sys[1] = self.complex_sys[0].run_parallel_dynamics(complex_sim_dir, complex_name,
                                                                                    self.num_frames, equi, None)
        #SOLVENT
        self.solvent_sys = []
        self.solvent_sys.append(FSim(ligand_name=ligand_name, sim_name=solvent_name, input_folder=input_folder,
                                     param=param, num_gpu=num_gpu, offset=self.solvent_offset, opt=opt,
                                     exclude_dualtopo=exclude_dualtopo, system=systems.solvent))
        self.solvent_sys.append([solvent_sim_dir + solvent_name + '.dcd'])
        self.solvent_sys.append(solvent_sim_dir + solvent_name + '.pdb')
        if not os.path.isfile(self.solvent_sys[1][0]):
            self.solvent_sys[1] = [solvent_sim_dir + solvent_name + '_gpu' + str(x) + '.dcd' for x in range(num_gpu)]
            for name in self.solvent_sys[1]:
                if not os.path.isfile(name):
                    self.solvent_sys[1] = self.solvent_sys[0].run_parallel_dynamics(solvent_sim_dir, solvent_name,
                                                                                    self.num_frames, equi, None)

        Fluorify.scanning(self, wt_ligand, auto_select, c_atom_list, h_atom_list, o_atom_list)

    def scanning(self, wt_ligand, auto_select, c_atom_list, h_atom_list, o_atom_list):
        """

        :param wt_ligand:
        :param auto_select:
        :param c_atom_list:
        :param h_atom_list:
        :param o_atom_list:
        :return:
        """

        #Generate mutant systems with selected pertibations
        mutated_systems, mutations = Fluorify.element_perturbation(self, auto_select, c_atom_list, h_atom_list, o_atom_list)

        """
        Write Mol2 files with substitutions of selected atoms.
        Run antechamber and tleap on mol2 files to get prmtop.
        Create OpenMM systems of ligands from prmtop files.
        Extract ligand parameters from OpenMM systems.
        """
        print('Parametrize mutant ligands...')
        t0 = time.time()

        mutated_ligands = []
        for index, sys in enumerate(mutated_systems):
            mol_name = 'molecule'+str(index)
            Mol2.write_mol2(sys, self.output_folder, mol_name)
            mutated_ligands.append(MutatedLigand(file_path=self.output_folder, mol_name=mol_name,
                                                 net_charge=self.net_charge, gaff=self.gaff_ver))

        wt_parameters = wt_ligand.get_parameters()
        mutant_parameters = []
        for i, ligand in enumerate(mutated_ligands):
            mute = mutations[i]['subtract']
            mutant_parameters.append(ligand.get_parameters(mute))

        #last entry of mutant is wildtype
        mutant_parameters.append(wt_parameters)
        mutations.append({'add': [], 'subtract': [], 'replace': [None], 'replace_insitu': [None]})

        mutant_params = Mutants(mutant_parameters, mutations, self.complex_sys[0], self.solvent_sys[0])
        del mutant_parameters

        t1 = time.time()
        print('Took {} seconds'.format(t1 - t0))

        """
        Apply ligand charges to OpenMM complex and solvent systems.
        Calculate potential energy of simulation with mutant charges.
        Calculate free energy change from wild type to mutant.
        """
        print('Calculating free energies...')
        t0 = time.time()


        print('Computing complex potential energies...')
        complex_free_energy = FSim.treat_phase(self.complex_sys[0], mutant_params.complex_params, self.complex_sys[1],
                                               self.complex_sys[2], self.num_frames)
        print('Computing solvent potential energies...')
        solvent_free_energy = FSim.treat_phase(self.solvent_sys[0], mutant_params.solvent_params, self.solvent_sys[1],
                                               self.solvent_sys[2], self.num_frames)

        #RESULT
        best_mutants = []
        for i, energy in enumerate(complex_free_energy):
            atom_names = []
            replace = mutations[i]['replace']
            replace.extend(mutations[i]['replace_insitu'])
            for atom in replace:
                atom_index = int(atom)-1
                atom_names.append(self.mol2_ligand_atoms[atom_index])
            binding_free_energy = energy - solvent_free_energy[i]
            best_mutants.append([binding_free_energy, atom_names, i])
            '''
            print('dGs for molecule{}.mol2 with'
                  ' {} substituted for {} complex dG = {}, solvent dG = {}'.format(str(i), atom_names, self.job_type,
                                                                                   energy, solvent_free_energy[i]))
            '''
            print('ddG for molecule{}.mol2 with'
                  ' {} substituted for {} = {}'.format(str(i), atom_names, self.job_type, binding_free_energy))
        best_mutants = sorted(best_mutants)
        t1 = time.time()
        print('Took {} seconds'.format(t1 - t0))

        x_best = min(self.num_fep, len(best_mutants))

        print('Calculating FEP for {} best mutants...'.format(x_best))
        t0 = time.time()
        for x in range(x_best):
            complex_dg, complex_error = self.complex_sys[0].run_parallel_fep(mutant_params, 0, best_mutants[x][2],
                                                                             20000, 50, 12)
            solvent_dg, solvent_error = self.solvent_sys[0].run_parallel_fep(mutant_params, 1, best_mutants[x][2],
                                                                             20000, 50, 12)
            ddg_fep = complex_dg - solvent_dg
            ddg_error = (complex_error**2+solvent_error**2)**0.5
            print('Mutant {}:'.format(best_mutants[x][1]))
            print('ddG Fluorine Scanning = {}'.format(best_mutants[x][0]))
            print('ddG FEP = {} +- {}'.format(ddg_fep, ddg_error))
        t1 = time.time()
        print('Took {} seconds'.format(t1 - t0))

    def element_perturbation(self, auto_select, c_atom_list, h_atom_list, o_atom_list):
        """
        Takes a job_type with an auto selection or atom_lists and turns this into the correct mutated_systems
        and a list of mutations to keep track of the mutation applied to each system.

        if single mutation:
            add fluorine, pyridine or hydroxyl
        else double mutation:
            add fluorine and pyridine
        """
        job_type = self.job_type.split('x')
        c_atoms = atom_selection(c_atom_list)
        h_atoms = atom_selection(h_atom_list)
        o_atoms = atom_selection(o_atom_list)
        if len(job_type) == 1:
            if job_type == ['N']:
                job_type = 'N.ar'
                if h_atoms is not None:
                    raise ValueError('hydrogen or oxygen atom list provided but pyridine scanning job requested')
                mutated_systems, mutations = add_nitrogens(self.mol, job_type, auto_select, c_atoms)
            elif job_type == ['F'] or job_type == ['Cl']:
                job_type = job_type[0]
                if c_atoms is not None or o_atom_list is not None:
                    raise ValueError('carbon or oxygen atom list provided but fluorine scanning job requested')
                mutated_systems, mutations = add_fluorines(self.mol, job_type, auto_select, h_atoms)
            elif job_type == ['OH']:
                job_type = 'O'
                if c_atoms is not None or o_atom_list is not None:
                    raise ValueError('carbon atom list provided but hydroxyl scanning job requested')
                mutated_systems, mutations = add_hydroxyl(self.mol, job_type, auto_select, h_atoms)
            elif job_type == ['S']:
                job_type = job_type[0]
                if c_atom_list is not None or h_atom_list is not None:
                    raise ValueError('carbon or hydrogen atom list provided but sulphur scanning job requested')
                mutated_systems, mutations = add_sulphurs(self.mol, job_type, auto_select, o_atoms)
            elif job_type == ['VDW']:
                job_type = 'H'
                if c_atoms is not None or o_atom_list is not None:
                    raise ValueError('carbon or oxygen atom list provided but VDW scanning job requested')
                mutated_systems, mutations = add_fluorines(self.mol, job_type, auto_select, h_atoms)
            else:
                raise ValueError('No recognised job type provided')

        else:
            if h_atoms is None or c_atoms is None:
                raise ValueError('Mixed substitution requested must provide carbon and hydrogen atom lists')
            job_type[0] = 'N.ar'
            mutated_systems = []
            mutations = []
            f_mutated_systems, f_mutations = add_fluorines(self.mol, job_type[1], auto_select, h_atoms)
            for i, mol in enumerate(f_mutated_systems):
                p_mutated_systems, p_mutations = add_nitrogens(mol, job_type[0],
                                                               auto_select, c_atoms, modified_atom_type=job_type[1])
                #compound fluorination and pyridination
                for j in range(len(p_mutations)):
                    p_mutations[j]['add'].extend(f_mutations[i]['add'])
                    p_mutations[j]['subtract'].extend(f_mutations[i]['subtract'])
                    p_mutations[j]['replace'].extend(f_mutations[i]['replace'])
                    p_mutations[j]['replace_insitu'].extend(f_mutations[i]['replace_insitu'])
                mutated_systems.extend(p_mutated_systems)
                mutations.extend(p_mutations)
        return mutated_systems, mutations


class SysBuilder(object):
    def __init__(self, cwd, protein, ligand, protein_FF, water_FF, gaff_FF, boxPadding, ionicStrength, ligand_charge=0,
                 gaff=2, using_yank=False):

        if using_yank is True:
            #if using yank then this class is a dummy and contains nothing
            self.solvent = None
            self.complex = None
            return

        print('Using Fluorify system builder...')

        supported_water_modles = ['tip3p', 'tip4p', 'tip4pew', 'spce']
        for model in supported_water_modles:
            if model in water_FF:
                water_model = model

        # make file locations
        self.protein_file = os.path.join(cwd, protein)
        self.ligand_file = os.path.join(cwd, ligand)

        prmtop_filename, inpcrd_filename, ligand_xml = SysBuilder.paramaterise_ligand(self, ligand_charge, gaff)
        position_file = mm.app.AmberInpcrdFile(inpcrd_filename)
        topology_file = mm.app.AmberPrmtopFile(prmtop_filename)

        FF = mm.app.ForceField(protein_FF, water_FF, gaff_FF, ligand_xml)
        pdb = mm.app.PDBFile(self.protein_file)

        # Create COMPLEX System
        # add protein
        complex = mm.app.Modeller(pdb.topology, pdb.positions)
        # add ligand
        complex.add(topology_file.topology, position_file.positions)
        # solvate system
        complex.addSolvent(FF, model=water_model, ionicStrength=ionicStrength, neutralize=True, padding=boxPadding)
        # create openmm system object of complex system
        self.complex = FF.createSystem(complex.topology, nonbondedMethod=mm.app.PME,
                                       nonbondedCutoff=1.0 * unit.nanometers, constraints=mm.app.HBonds,
                                       rigidWater=True, ewaldErrorTolerance=0.0005)

        # Create Solvent system
        # add ligand
        solvent = mm.app.Modeller(topology_file.topology, position_file.positions)
        # solvate system
        solvent.addSolvent(FF, model=water_model, ionicStrength=ionicStrength, neutralize=True, padding=boxPadding)
        # create openmm system object of complex system
        self.solvent = FF.createSystem(solvent.topology, nonbondedMethod=mm.app.PME,
                                       nonbondedCutoff=1.0 * unit.nanometers, constraints=mm.app.HBonds,
                                       rigidWater=True, ewaldErrorTolerance=0.0005)

        # Make dirs
        paths = [os.path.join(cwd, 'complex'), os.path.join(cwd, 'solvent')]
        for path in paths:
            try:
                os.mkdir(path)
            except OSError:
                print("Creation of the directory {} failed".format(path))
            else:
                print("Successfully created the directory {} ".format(path))

        outfile = open(os.path.join(paths[0], 'complex.pdb'), 'w')
        mm.app.PDBFile.writeFile(complex.topology, complex.positions, outfile)
        outfile.close()

        outfile = open(os.path.join(paths[1], 'solvent.pdb'), 'w')
        mm.app.PDBFile.writeFile(solvent.topology, solvent.positions, outfile)
        outfile.close()

    def paramaterise_ligand(self, ligand_charge, gaff):
        if gaff == 1:
            gaff_version = 'gaff'
            leaprc = 'leaprc.gaff'
        elif gaff == 2:
            gaff_version = 'gaff2'
            leaprc = 'leaprc.gaff2'

        ligand_name = 'ligand'
        ligand_xml = "{}.xml".format(ligand_name)
        #Could use mol2.run_ante
        new_ligand_file, ligand_frcmod = moltools.amber.run_antechamber(ligand_name, self.ligand_file,
                                                                        charge_method="bcc",
                                                                        net_charge=ligand_charge,
                                                                        gaff_version=gaff_version)
        prmtop_filename, inpcrd_filename = moltools.amber.run_tleap(ligand_name, new_ligand_file, ligand_frcmod,
                                                                    leaprc=leaprc)

        # Convert frcmod to XML
        ff = pmd.openmm.OpenMMParameterSet.from_parameterset(pmd.amber.AmberParameterSet(ligand_frcmod))
        mol2 = pmd.load_file(new_ligand_file)
        ff.residues[mol2.name] = mol2
        ff.write(ligand_xml)

        return prmtop_filename, inpcrd_filename, ligand_xml


def atom_selection(atom_list):

    if atom_list is None:
        return None

    #generate all permutations
    atom_list = list(itertools.product(*atom_list))
    #remove permutations with duplicates
    tmp = []
    for x in atom_list:
        if len(x) == len(set(x)):
            tmp.append(x)
    atom_list = tmp
    atom_list = [list(sorted(q)) for q in atom_list]
    #remove equivilant permutations
    tmp = []
    for x in atom_list:
        if x not in tmp:
            tmp.append(x)
    atom_list = tmp
    return atom_list


def get_coordinate_system(mutant, neighbours, h_atom):
    """
    get an origin and axis based on atom loacations
    """
    axis = []
    positions = []
    positions.append(Mol2.get_atom_position(mutant, int(h_atom)))
    neighbours = [x for x in neighbours if int(x) is not int(h_atom)]
    if len(neighbours) != 2:
        raise ValueError('cant hydroxyl')
    for atom in neighbours:
        positions.append(Mol2.get_atom_position(mutant, int(atom)))
    axis.append(positions[0])                   #r1
    axis.append(positions[1] - positions[0])    #r12
    axis.append(positions[2] - positions[0])    #r13
    axis.append(np.cross(axis[1], axis[2]))     #rcross
    return axis


def add_hydroxyl(mol, new_element, auto_select, atom_list):
    """
    Not fully implemented
    """
    raise ValueError('Hydoxyl not fully implemented')
    mutations = []
    if atom_list is None:
        ValueError('No auto for hydroxyl')
    else:
        hydrogens = Mol2.get_atom_by_string(mol, 'H')
        for pair in atom_list:
            tmp = [x for x in pair if x not in hydrogens]
            if len(tmp) > 0:
                raise ValueError('Atoms {} are not recognised as hydrogens and therefore can not be hydroxylated'.format(tmp))
        bonded_h = atom_list
    mutated_systems = Mol2.mutate(mol, bonded_h, new_element)

    for pair, mutant in zip(bonded_h, mutated_systems):
        for atom in pair:
            h_neigh = Mol2.get_bonded_neighbours(mutant, atom)
            c_neigh = Mol2.get_bonded_neighbours(mutant, h_neigh[0])
            axis = get_coordinate_system(mutant, c_neigh, atom)
            weights = [[-0.2, -0.2, -0.2], [-0.2, -0.2, -0.2], [0.2, 0.2, 0.2]] #defines the position of hydrogen relative to ring
            xyz = weights[0]*axis[1] + weights[1]*axis[2] + weights[2]*axis[3]
            #axis[0] is origin
            position = [axis[0][0]+xyz[0], axis[0][1]+xyz[1], axis[0][2]+xyz[2]]
            mutant.add_atom(int(atom), 'H', position)

    #Atom indexes should be 0 indexed to interface with openmm
    #Build list of dictionaries each dict describes mutation applied corresponding system
    #TODO
    for i, mutant in enumerate(mutated_systems):
        mutations.append({'add': [], 'subtract': [], 'replace': [], 'replace_insitu': []})
        for atom in bonded_h[i]:
            mutations[i]['add'].append([int(atom), c_neigh, weights])
    return mutated_systems, mutations


def add_fluorines(mol, new_element, auto_select, atom_list):
    mutations = []
    if atom_list is None:
        carbon_type = 'C.' + auto_select
        carbons = Mol2.get_atom_by_string(mol, carbon_type)
        carbons_neighbours = []
        for atom in carbons:
            carbons_neighbours.extend(Mol2.get_bonded_neighbours(mol, atom))
        hydrogens = Mol2.get_atom_by_string(mol, 'H')
        bonded_h = [[x] for x in hydrogens if x in carbons_neighbours]
    else:
        hydrogens = Mol2.get_atom_by_string(mol, 'H')
        for pair in atom_list:
            tmp = [x for x in pair if x not in hydrogens]
            if len(tmp) > 0:
                raise ValueError('Atoms {} are not recognised as hydrogens and therefore can not be fluorinated'.format(tmp))
        bonded_h = atom_list
    mutated_systems = Mol2.mutate(mol, bonded_h, new_element)

    #Build list of dictionaries each dict describes mutation applied corresponding system
    for i, mutant in enumerate(mutated_systems):
        mutations.append({'add': [], 'subtract': [], 'replace': [], 'replace_insitu': []})
        for atom in bonded_h[i]:
            mutations[i]['replace'].append(int(atom))
    return mutated_systems, mutations

def add_sulphurs(mol, new_element, auto_select, atom_list):
    mutations = []
    if atom_list is None:
        oxygen_type = 'O.' + auto_select
        atom_list = Mol2.get_atom_by_string(mol, oxygen_type)
        atom_list = [[x] for x in atom_list]
    mutated_systems = Mol2.mutate(mol, atom_list, new_element)
    #should check this is a double bonded sulphur with only carbon neighbours

    #Build list of dictionaries each dict describes mutation applied corresponding system
    for i, mutant in enumerate(mutated_systems):
        mutations.append({'add': [], 'subtract': [], 'replace': [], 'replace_insitu': []})
        for atom in atom_list[i]:
            mutations[i]['replace_insitu'].append(int(atom))
    return mutated_systems, mutations


def add_nitrogens(mol, new_element, auto_select, atom_list, modified_atom_type=None):
    """
    Look for carbons with one hydrogen neighbour.
    Reduce list of carbons to those with one hydrogen neighbour.
    Make a note of which hydrogen is the neighbour.
    Swap the carbon for a nitrogen and label hydrogen to be muted
    or use user provided list of carbons
    """
    mutations = []
    carbons, bonded_h = get_single_neighbour_carbons(mol, modified_atom_type, auto_select)
    if atom_list is None:
        carbons = [[x] for x in carbons]
        # need to index from 0 to interface with openmm
        hydrogens_to_remove = [[int(x) - 1] for x in bonded_h]
    else:
        hydrogens_to_remove = []
        for pair in atom_list:
            tmp = [x for x in pair if x not in carbons]
            if len(tmp) > 0:
                raise ValueError('Atoms {} are not recognised as carbons with'
                                 ' one hydrogen neighbour and therefore can not be pyridinated'.format(tmp))
        carbons = atom_list
        for pair in carbons:
            h_tmp = []
            for atom in pair:
                neighbours = Mol2.get_bonded_neighbours(mol, atom)
                for neighbour in neighbours:
                    #This bonded_h includes atom which have been pyridinated.
                    if neighbour in bonded_h:
                        h_tmp.append(neighbour)
            h_tmp.sort(reverse=True)
            # need to index from 0 to interface with openmm.
            h_tmp = [int(x) - 1 for x in h_tmp]
            hydrogens_to_remove.append(h_tmp)

    #Mutate carbons to nitrogen and remove hydrogens.
    mutated_systems = Mol2.mutate(mol, carbons, new_element)
    for i, mutant in enumerate(mutated_systems):
        for atom in hydrogens_to_remove[i]:
            #mol2 indexed from 1 so +1
            mutant.remove_atom(atom+1)

    # Build list of dictionaries, each dict describes mutation applied corresponding system
    for i, mutant in enumerate(mutated_systems):
        mutations.append({'add': [], 'subtract': [], 'replace': [], 'replace_insitu': []})
        for atom in carbons[i]:
            mutations[i]['replace_insitu'].append(int(atom))
        for atom in hydrogens_to_remove[i]:
            mutations[i]['subtract'].append(int(atom))
    return mutated_systems, mutations


def get_single_neighbour_carbons(mol, modified_atom_type, carbon_type=None):
    if carbon_type is not None:
        carbon_type = 'C.' + carbon_type
        carbons = Mol2.get_atom_by_string(mol, carbon_type)
    else:
        carbons = Mol2.get_atom_by_string(mol, 'C.', wild_card=True)

    if modified_atom_type is not None:
        modified_atoms = Mol2.get_atom_by_string(mol, modified_atom_type)
    else:
        modified_atoms = []

    hydrogens = Mol2.get_atom_by_string(mol, 'H')
    bonded_h = []
    c_tmp = []
    for atom in carbons:
        h_neigh = []
        neighbours = Mol2.get_bonded_neighbours(mol, atom)
        for neighbour in neighbours:
            if neighbour in hydrogens:
                h_neigh.append(neighbour)
            if neighbour in modified_atoms:
                h_neigh.append(neighbour)
        if len(h_neigh) == 1:
            bonded_h.extend(h_neigh)
            c_tmp.append(atom)
    carbons = c_tmp
    return carbons, bonded_h


def get_ligand_offset(input_files, mol2_ligand_atoms, ligand_name):
    """
    almost certain bond will be inconsistent between inputs need to map bonds between inputs
    :return:
    """
    offset = []
    for file in input_files:
        snapshot = md.load(file)
        offset.append(snapshot.topology.select('resname {} and name {}'.format(ligand_name, mol2_ligand_atoms[0])))
    return tuple(offset)


def get_atom_list(input_files, resname):
    atom_lists = []
    for file in input_files:
        atoms = []
        traj = md.load(file)
        top = traj.topology
        mol = top.select('resname '+resname)
        for idx in mol:
            atoms.append(str(top.atom(idx)).split('-')[1])
        atom_lists.append(atoms)
    return tuple(atom_lists)

