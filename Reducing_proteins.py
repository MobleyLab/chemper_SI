import os
from openeye import oechem
import pickle
import copy
import glob
import itertools
from chemper.chemper_utils import check_smirks_to_reference, get_typed_molecules, create_tuples_for_clusters
from making_proteins import everything_from_fastas, print_order_type_data, at_least_one_passed, ParameterSystem, by_biggest_size, by_biggest_smirks
import cmiles
import json
from parmed.modeller import ResidueTemplate
from chemper.smirksify import Reducer, print_smirks


def convert_json_and_oeb(json_file, mol_dir='./mol_files/'):
    with open(json_file, 'r') as inputf:
        d = json.load(inputf)

    mol_dir = os.path.abspath(mol_dir)
    mol_files = [os.path.join(mol_dir, m) for m in d['mol_files']]
    mols = list()
    for mol_file in mol_files:
        mol = oechem.OEMol()
        ifs = oechem.oemolistream(mol_file)
        while oechem.OEReadMolecule(ifs,mol):
            mols.append(oechem.OEMol(mol))

    return mols, d['smirks_lists'], d['clusters']

# # Final Dictionary
#
# The final results I want are:
#
# initial SMIRKS, reduced SMIRKS (None if failed), four sorting, all fragments so the dictionary will have:
#
# **Fragement**
# * ordering
#     - 'initial smirks'
#     - 'final smirks'
#

final_dict = dict()
iterations = 1000

file_keys = [
    ('big', ['big_smirks', 'biggest_size']),
    ('small', ['small_smirks', 'small_size'])
]
for fn_label, cluster_orders in file_keys:
    fns = glob.glob('./all_aminos_together/allIn1_%s_99sbildn_*_1mols.json' % fn_label)
    print('='*80)
    print(' '*20,fn_label)
    print('='*80)
    for f in fns:
        frag = f.split('_')[-2]
        if frag == 'torsion':
            prefix = f.split('_')[-3]
            frag = '%s_%s' % (prefix, frag)
        print('-'*80)
        print(' '*30,frag)
        print('-'*80)
        mols, dsmirks, dclusters = convert_json_and_oeb(f)

        if frag not in final_dict:
            final_dict[frag] = dict()

        for order in cluster_orders:
            if order in final_dict[frag]:
                if 'output' in final_dict[frag][order] or 'output_10k' in final_dict[frag][order]:
                    print('Already in dictionary ', frag, order)
                    continue

            final_dict[frag][order] = dict()
            d = dsmirks[order][frag]
            type_list = [(l, s) for l,s in d['type_list']]
            final_dict[frag][order]['initial'] = type_list

            print('ORIGINAL', order)
            print_smirks(type_list)

            if not d['checked']:
                final_dict[frag][order]['output_10k'] = None
                continue

            red = Reducer(type_list, mols, verbose=False)
            final_dict[frag][order]['output_1k'] = red.run(iterations)
            print('REDUCED 1k ', order)
            print_smirks(final_dict[frag][order]['output_1k'])

            final_dict[frag][order]['output_5k'] = red.run(4000)
            print('REDUCED 5k ', order)
            print_smirks(final_dict[frag][order]['output_5k'])



pickle.dump(final_dict, open('./all_aminos_together/reduced_smirks_dict_5k.p', 'wb'))


