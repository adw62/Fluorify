#!/usr/bin/env python

from .energy import FSim
from .mol2 import Mol2, MutatedLigand
from .optimize import Optimize

import os
import time
import itertools
import mdtraj as md
import numpy as np
import shutil
from simtk import unit
import logging
import copy

logger = logging.getLogger(__name__)

#CONSTANTS
e = unit.elementary_charges

class Mutants(object):
    def __init__(self, params, atom_order, virt_atoms, bond_order, exception_order, virt_exceptions, mutations, offset, virtual_offset):
        """
        A class for applying mutant parameters to openmm system
        :param params: List of lists for mutant parameters
        :param atom_order: list of lists for order of atoms in each sys
        :param exception_order: list of lists for order of exceptions in each sys
        :param bond_order: list of lists for order of bonds in each sys
        :param sys_offset: number of addition atoms between begining of ligand in each sys
        :param virtual_offset: number of total atoms in each sys, such that this will be the
                                first index of any vitrual atoms added to represent dual topology
        """
        self.offset = offset
        self.virtual_offset = virtual_offset

        self.nonbonded_params = Mutants.build_nonbonded(self, params, mutations,
                                                                               virt_atoms, atom_order)
        self.exception_params = Mutants.build_exceptions(self, params, mutations,
                                                                                virt_exceptions, exception_order)
        self.bond_params = Mutants.build_bonds(self, params, bond_order)

    def build_nonbonded(self, params, mutations, virt_atoms, atom_order):
        #Reorder params from ligand only system to match solvent and complex systems
        nonbonded_params = [[], []]
        for i, (sys_atom_order, sys_offset) in enumerate(zip(atom_order, self.offset)):
            sys_nonbonded_params = copy.deepcopy([x[0] for x in params])
            for j, mutant_params in enumerate(sys_nonbonded_params):
                map = {x['id']: x for x in mutant_params}
                sys_nonbonded_params[j] = [map[int(atom - sys_offset)] for atom in sys_atom_order]
                nonbonded_params[i].extend(sys_nonbonded_params)

        # Build nonbonded ghosts which handle parmas for dual topology
        ghosts = [[0.0 * e, 0.26 * unit.nanometer, 0.0 * unit.kilojoules_per_mole] for i in
                     range(len(virt_atoms[0]))]
        ghosts = [copy.deepcopy(ghosts) for i in range(len(nonbonded_params[0]))]

        # transfer params from original topology to ghost topology
        nonbonded_ghosts = [[], []]
        for i, (sys, sys_virt_atoms) in enumerate(zip(nonbonded_params, virt_atoms)):
            for mutant_params, mutant, ghost in zip(sys, mutations, ghosts):
                atom_idxs = mutant['replace']
                if atom_idxs[0] is not None:
                    for atom in atom_idxs:
                        atom = int(atom-1)
                        transfer_params = copy.deepcopy(mutant_params[atom])
                        if int(transfer_params['id']) != int(atom):
                            raise ValueError('Mismatch between indexes whilst building mutant parameters')
                        transfer_params = transfer_params['data']
                        mutant_params[atom] = [0.0*e, 0.26*unit.nanometer, 0.0*unit.kilojoules_per_mole]
                        transfer_index = sys_virt_atoms.index(atom)
                        ghost[transfer_index] = transfer_params
            nonbonded_ghosts[i].extend(ghosts)

        print(nonbonded_ghosts[0][0])

        return [nonbonded_params, nonbonded_ghosts]

    def build_exceptions(self, params, mutations, virt_exceptions, exception_order):
        """

        :param params:
        :param mutations:
        :param virt_exceptions: ghost from complex and solsvent systems
        :param exception_order: order of vannila exceptions for complex and solvent systems
        :return:
        """
        #reorder
        exception_params = [[], []]
        for i, (sys_excep_order, sys_offset) in enumerate(zip(exception_order, self.offset)):
            sys_exception_params = copy.deepcopy([x[3] for x in params])
            for j, mutant_parmas in enumerate(sys_exception_params):
                map = {x['id']: x for x in mutant_parmas}
                sys_exception_params[j] = [map[frozenset(int(x-sys_offset) for x in atom)] for atom in sys_excep_order]
                exception_params[i].extend(sys_exception_params)

        """
        #TEST
        for x,y in zip(sys_exception_params[1][0], exception_order[1]):
            if x['id'] != frozenset(int(atom-0) for atom in y):
                raise ValueError(x['id'], frozenset(int(atom-0) for atom in y))
        """

        #Build ghost exceptions
        exception_ghosts = copy.deepcopy(exception_params)
        for i, (sys_exception_params, sys_virt_order, sys_offset) in enumerate(zip(exception_ghosts, virt_exceptions, self.offset)):
            for j, mutant_parmas in enumerate(sys_exception_params):
                map = {x['id']: x for x in mutant_parmas}
                exception_ghosts[i][j] = [map[frozenset(int(x-sys_offset) for x in atom)] for atom in sys_virt_order]

        #CORRECT PLZ
        zero = 0.0

        # zero flourines in original topology aka exception params
        for i, sys_exception_params in enumerate(exception_params):
            for j, (mutant_parmas, mutant) in enumerate(zip(sys_exception_params, mutations)):
                atom_idxs = mutant['replace']
                if None not in atom_idxs:
                    for atom in atom_idxs:
                        atom = int(atom-1)
                        for k, excep1 in enumerate(mutant_parmas):
                            if atom in excep1['id']:
                                exception_params[i][j][k] = zero

        # zero everything but fluorines in dual topology aka exception_ghosts
        for i, sys_exception_ghosts in enumerate(exception_ghosts):
            for j, (mutant_parmas, mutant) in enumerate(zip(sys_exception_ghosts, mutations)):
                atom_idxs = mutant['replace']
                if None not in atom_idxs:
                    for atom in atom_idxs:
                        atom = int(atom-1)
                        for k, excep1 in enumerate(mutant_parmas):
                            if atom not in excep1['id']:
                                exception_ghosts[i][j][k] = zero

        #remove ids

        return [exception_params, exception_ghosts]

    def build_bonds(self, bond_order):
        pass

    def get_complex_mutants(self):
        return [[x, y, z] for x, y, z in zip(self.non_bonded_params[0], self.exception_params[0], self.bond_params[0])]

    def get_solvent_mutants(self):
        return [[x, y, z] for x, y, z in zip(self.non_bonded_params[1], self.exception_params[1], self.bond_params[1])]


def trim_mutant_exceptions(ligand_exception_order, exceptions):
    ans = []
    for exception in exceptions:
        tmp = []
        for excep1, excep2 in zip(ligand_exception_order, exception):
            for excep3 in h_exceptions:
                if excep1 == excep3:
                    tmp.append(excep2)
        ans.append(tmp)
    return ans